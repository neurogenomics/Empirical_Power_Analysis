#!/usr/bin/env bash
# =============================================================================
# pipeline/run.sh -- one-command reproduction of the Empirical Power Analysis
# paper on a PBS Pro cluster (Imperial CX3) or, for testing, on a single host.
# -----------------------------------------------------------------------------
# Usage:
#   pipeline/run.sh WORKDIR [options]
#
#   WORKDIR must contain (see pipeline/README.md, "Preparing a work directory"):
#     raw/         PsychAD-HBCC.h5ad PsychAD-MSSM.h5ad PsychAD-RADC.h5ad media-1.csv
#     reference/   the offline reference bundle (pipeline/reference/README.md)
#     epa.sif      the Apptainer image (or set EPA_SIF in WORKDIR/epa.env)
#
# Options:
#   --dry-run       write the job scripts and print what would be submitted; submit nothing
#   --status        print the state of every stage and exit
#   --local         run every stage sequentially on this host instead of submitting
#                   (testing only; needs apptainer/singularity and a lot of RAM)
#   --only S1,S2    restrict this invocation to the named stages (dependencies
#                   must already be complete)
#   --allow-dirty   snapshot uncommitted changes too (the run is then not tied
#                   to a clean commit; recorded in code/COMMIT)
#   -h, --help      show this help
#
# Behaviour:
#   * The first invocation freezes the code (git archive of HEAD, minus Datalink/)
#     into WORKDIR/code/code.tar. Every job runs that snapshot, so editing the
#     repository while the pipeline runs changes nothing.
#   * Each stage is a PBS job; dependencies are expressed with
#     `-W depend=afterok:<id>`. Every setting is written into the generated job
#     script itself (WORKDIR/jobs/<stage>.pbs); nothing is passed via `qsub -v`.
#   * Re-running the same command resumes: stages whose results are complete are
#     skipped, stages still queued/running are reused as dependencies, and only
#     missing or failed work is resubmitted (for the grid array: only the
#     missing chunks). This check happens here, on the login node, so no large
#     allocation is ever requested just to discover there is nothing to do.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$PIPELINE_DIR/.." && pwd)"

usage() { sed -n '2,/^# =\{70,\}$/p' "$0" | sed '1d;$d' | sed 's/^# \{0,1\}//'; }

# ---- arguments ----------------------------------------------------------------
WORK="" MODE="submit" ONLY="" ALLOW_DIRTY=false
while [ "$#" -gt 0 ]; do
  case "$1" in
    --dry-run)     MODE="dry-run"; shift ;;
    --status)      MODE="status"; shift ;;
    --local)       MODE="local"; shift ;;
    --only)        ONLY="$2"; shift 2 ;;
    --allow-dirty) ALLOW_DIRTY=true; shift ;;
    -h|--help)     usage; exit 0 ;;
    -*)            echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
    *)             [ -z "$WORK" ] || { echo "only one WORKDIR allowed" >&2; exit 2; }; WORK="$1"; shift ;;
  esac
done
[ -n "$WORK" ] || { usage >&2; exit 2; }
[ -d "$WORK" ] || { echo "WORKDIR does not exist: $WORK" >&2; exit 2; }
WORK="$(cd -- "$WORK" && pwd)"
case "$WORK" in *[[:space:]]*) echo "WORKDIR path must not contain spaces: $WORK" >&2; exit 2 ;; esac

die()  { echo "run.sh: FATAL: $*" >&2; exit 1; }
info() { echo "run.sh: $*"; }

# ---- configuration --------------------------------------------------------------
# Defaults come from pipeline/epa.env.defaults; WORKDIR/epa.env (created on first
# use from the defaults) holds this run's settings. Edit it before the first run.
if [ ! -f "$WORK/epa.env" ]; then
  cp "$PIPELINE_DIR/epa.env.defaults" "$WORK/epa.env"
  info "created $WORK/epa.env from the defaults (edit it to change settings)"
fi
# shellcheck source=/dev/null
source "$PIPELINE_DIR/epa.env.defaults"
# shellcheck source=/dev/null
source "$WORK/epa.env"
: "${EPA_SIF:=$WORK/epa.sif}"
EPA_RAW_DIR="$WORK/raw"
EPA_REF_DIR="$WORK/reference"

mkdir -p "$WORK"/{code,jobs,logs,results,state}

# ---- stage table ---------------------------------------------------------------
# pipeline/stages.tsv: name  deps(comma|-)  ncpus  mem  walltime  kind(single|array:grid|array:N)
STAGES=() declare -A S_DEPS S_NCPUS S_MEM S_WALL S_KIND
while IFS=$'\t' read -r name deps ncpus mem wall kind; do
  case "$name" in ''|'#'*) continue ;; esac
  STAGES+=("$name")
  S_DEPS[$name]="$deps"; S_NCPUS[$name]="$ncpus"; S_MEM[$name]="$mem"
  S_WALL[$name]="$wall"; S_KIND[$name]="$kind"
  # Per-stage resource override: EPA_RES_<stage>="ncpus mem walltime"
  ovr="EPA_RES_${name}"
  if [ -n "${!ovr:-}" ]; then
    read -r S_NCPUS[$name] S_MEM[$name] S_WALL[$name] <<< "${!ovr}"
  fi
done < "$PIPELINE_DIR/stages.tsv"
for s in "${STAGES[@]}"; do
  [ "${S_DEPS[$s]}" = "-" ] && continue
  IFS=',' read -r -a _d <<< "${S_DEPS[$s]}"
  for d in "${_d[@]}"; do [ -n "${S_KIND[$d]:-}" ] || die "stage $s depends on unknown stage $d"; done
done

# ---- stage state -----------------------------------------------------------------
# kind = single | array:grid | array:<N>. An array stage's task list is
# "<task-id> <output-tar-name>" per line; the grid has one task per (mode, chunk).
expected_tasks() {
  local s="$1" k="${S_KIND[$1]#array:}" m c
  if [ "$k" = grid ]; then
    for m in reads cells; do
      for c in $(seq 1 "$EPA_GRID_CHUNKS"); do echo "$m:$c epa_${m}_chunk_${c}"; done
    done
  else
    for c in $(seq 1 "$k"); do echo "$c task_${c}"; done
  fi
}
todo_tasks() {  # array tasks whose result tar is not there yet
  local id out
  expected_tasks "$1" | while read -r id out; do
    [ -f "$WORK/results/$1/$out.tar" ] || echo "$id $out"
  done
}
is_array() { case "${S_KIND[$1]}" in array:*) return 0 ;; *) return 1 ;; esac; }
stage_complete() {
  if is_array "$1"; then [ -z "$(todo_tasks "$1")" ]; else [ -f "$WORK/results/$1/.complete" ]; fi
}
# A job is "live" if the scheduler still has it queued, held or running.
# Unknown to the scheduler (history expired) counts as not live; any other qstat
# failure stops run.sh rather than risk submitting a duplicate of a running job.
job_live() {
  local id="$1" out st
  [ "$MODE" = "submit" ] || return 1
  if ! out="$(qstat -f -x "$id" 2>&1)"; then
    case "$out" in *"Unknown Job"*|*"unknown job"*) return 1 ;; esac
    die "could not query job $id with qstat (output: $out); not resubmitting anything"
  fi
  st="$(awk -F' = ' '/job_state/{print $2; exit}' <<< "$out")"
  case "$st" in Q|R|H|W|B|T|E|S) return 0 ;; *) return 1 ;; esac
}
stage_state() {
  local s="$1"
  if stage_complete "$s"; then echo complete; return; fi
  if [ -f "$WORK/state/$s.jobid" ] && job_live "$(cat "$WORK/state/$s.jobid")"; then echo live; return; fi
  if [ -f "$WORK/state/$s.jobid" ]; then echo "failed-or-killed"; return; fi
  echo pending
}

# Task lists of every array stage (read by jobs that consume array outputs).
for s in "${STAGES[@]}"; do
  if is_array "$s"; then expected_tasks "$s" > "$WORK/jobs/$s.expected"; fi
done

if [ "$MODE" = "status" ]; then
  printf '%-18s %-18s %s\n' STAGE STATE DETAIL
  for s in "${STAGES[@]}"; do
    detail=""
    if is_array "$s"; then
      n_all="$(expected_tasks "$s" | grep -c .)"
      n_todo="$(todo_tasks "$s" | grep -c . || true)"
      detail="$((n_all - n_todo))/$n_all tasks done"
    fi
    [ -f "$WORK/state/$s.jobid" ] && detail="$detail job $(cat "$WORK/state/$s.jobid")"
    printf '%-18s %-18s %s\n' "$s" "$(stage_state "$s")" "$detail"
  done
  exit 0
fi

# ---- preflight -------------------------------------------------------------------
info "work dir   : $WORK"
info "repository : $REPO_ROOT"
missing=0
while read -r f; do
  case "$f" in ''|'#'*) continue ;; esac
  [ -f "$EPA_RAW_DIR/$f" ] || { echo "  missing raw input: raw/$f" >&2; missing=1; }
done < <(awk '{print $2}' "$PIPELINE_DIR/inputs.sha256")
while read -r _ f; do
  case "$f" in ''|'#'*) continue ;; esac
  [ -f "$EPA_REF_DIR/$f" ] || { echo "  missing reference file: reference/$f" >&2; missing=1; }
done < "$PIPELINE_DIR/reference/MANIFEST.sha256"
[ -f "$EPA_SIF" ] || { echo "  missing container image: $EPA_SIF" >&2; missing=1; }
[ "$missing" -eq 0 ] || die "inputs incomplete (see above; pipeline/README.md explains how to prepare WORKDIR)"
# The reference bundle is small: verify it fully on every invocation. The raw
# inputs (~57 GB) are verified inside the first job (stage shard), not here.
( cd "$EPA_REF_DIR" && sha256sum --quiet -c "$PIPELINE_DIR/reference/MANIFEST.sha256" ) \
  || die "reference/ does not match pipeline/reference/MANIFEST.sha256 (rebuild it with pipeline/reference/fetch_reference.sh)"
if [ "$MODE" = "submit" ]; then command -v qsub >/dev/null || die "qsub not found (use --local or --dry-run off-cluster)"; fi
if [ "$MODE" = "local" ] && [ -z "${EPA_CONTAINER_BIN:-}" ]; then command -v apptainer >/dev/null || command -v singularity >/dev/null || die "apptainer/singularity not found"; fi

# ---- code snapshot (once per work directory) -------------------------------------
if [ ! -f "$WORK/code/code.tar" ]; then
  dirty="$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no -- . ':!Datalink' | head -1)"
  commit="$(git -C "$REPO_ROOT" rev-parse HEAD)"
  if [ -n "$dirty" ]; then
    $ALLOW_DIRTY || die "repository has uncommitted changes; commit them or pass --allow-dirty"
    info "snapshotting working tree with uncommitted changes (--allow-dirty)"
    stash="$(git -C "$REPO_ROOT" stash create)"
    git -C "$REPO_ROOT" archive --format=tar -o "$WORK/code/code.tar" "$stash" -- . ':!Datalink'
    echo "$commit+dirty" > "$WORK/code/COMMIT"
  else
    git -C "$REPO_ROOT" archive --format=tar -o "$WORK/code/code.tar" HEAD -- . ':!Datalink'
    echo "$commit" > "$WORK/code/COMMIT"
  fi
  sha256sum "$EPA_SIF" | cut -d' ' -f1 > "$WORK/code/SIF_SHA256"
  if [ -n "${EPA_SIF_SHA256:-}" ] && [ "$EPA_SIF_SHA256" != "$(cat "$WORK/code/SIF_SHA256")" ]; then
    rm -f "$WORK/code/code.tar"
    die "container image checksum does not match EPA_SIF_SHA256 in epa.env"
  fi
  cp "$WORK/epa.env" "$WORK/code/epa.env.at-snapshot"
  info "code snapshot: $(cat "$WORK/code/COMMIT")  (image sha256 $(cut -c1-12 "$WORK/code/SIF_SHA256"))"
else
  info "code snapshot: $(cat "$WORK/code/COMMIT")  (frozen at first run; delete $WORK/code/ to re-snapshot)"
  if [ "$(git -C "$REPO_ROOT" rev-parse HEAD)" != "$(cut -d+ -f1 "$WORK/code/COMMIT")" ]; then
    info "NOTE: the repository HEAD has moved since the snapshot; jobs still run the snapshot"
  fi
fi

# Resolved settings every job sources (absolute paths, nothing inherited from the
# login shell).
{
  echo "# generated by pipeline/run.sh $(date -u +%FT%TZ); do not edit, edit epa.env"
  echo "EPA_WORK='$WORK'"
  echo "EPA_SIF='$EPA_SIF'"
  echo "EPA_RAW_DIR='$EPA_RAW_DIR'"
  echo "EPA_REF_DIR='$EPA_REF_DIR'"
  compgen -v EPA_ | grep -vE '^EPA_(WORK|SIF|RAW_DIR|REF_DIR)$' | while read -r v; do
    printf "%s='%s'\n" "$v" "${!v}"
  done
} > "$WORK/epa.env.resolved"

# ---- job script generation -------------------------------------------------------
write_job() {  # stage [ntasks]
  local s="$1" ntasks="${2:-0}" f="$WORK/jobs/$1.pbs"
  mkdir -p "$WORK/logs/$s"
  {
    echo "#!/usr/bin/env bash"
    echo "#PBS -N epa-$s"
    echo "#PBS -l select=1:ncpus=${S_NCPUS[$s]}:mem=${S_MEM[$s]}"
    echo "#PBS -l walltime=${S_WALL[$s]}"
    echo "#PBS -j oe"
    echo "#PBS -o $WORK/logs/$s/"
    [ -n "${EPA_PBS_ACCOUNT:-}" ] && echo "#PBS -P $EPA_PBS_ACCOUNT"
    if [ "$ntasks" -gt 1 ]; then echo "#PBS -J 1-${ntasks}%${EPA_GRID_CONCURRENCY}"; fi
    cat <<EOF
# Generated by pipeline/run.sh for stage '$s'. All settings are baked in here;
# nothing is passed through qsub -v.
set -euo pipefail
export EPA_WORK='$WORK' EPA_STAGE='$s'
export EPA_JOB_KIND='${S_KIND[$s]}' EPA_JOB_NTASKS='$ntasks'
source "\$EPA_WORK/epa.env.resolved"
export EPA_LOCAL_NCPUS='${S_NCPUS[$s]}'
BOOT="\${TMPDIR:?TMPDIR not set}/epa-boot-$s-\${PBS_ARRAY_INDEX:-0}-\$\$"
mkdir -p "\$BOOT"
tar -xf "\$EPA_WORK/code/code.tar" -C "\$BOOT" pipeline/lib pipeline/stages
source "\$BOOT/pipeline/lib/common.sh"
source "\$BOOT/pipeline/stages/$s.sh"
rm -rf "\$BOOT"
if [ "\$EPA_JOB_KIND" != single ]; then
  # Array task: line PBS_ARRAY_INDEX of the to-do list written at submission
  # ("<task-id> <output-name>"); a single remaining task is submitted as a plain job.
  read -r EPA_TASK EPA_TASK_OUT < <(sed -n "\${PBS_ARRAY_INDEX:-1}p" "\$EPA_WORK/jobs/$s.todo")
  [ -n "\${EPA_TASK:-}" ] || { echo "no task for index \${PBS_ARRAY_INDEX:-1}" >&2; exit 2; }
  export EPA_TASK EPA_TASK_OUT
fi
stage_main
EOF
  } > "$f"
  chmod +x "$f"
}

# ---- submission ------------------------------------------------------------------
declare -A JOBID      # stage -> id of a live or newly submitted job this round
FAKE_ID=1000
submit() {  # stage ntasks depstring
  local s="$1" ntasks="$2" dep="$3" id
  local args=()
  [ -n "$dep" ] && args+=(-W "depend=afterok:$dep")
  if [ "$MODE" = "dry-run" ]; then
    id="DRY$FAKE_ID[]"; [ "$ntasks" -gt 1 ] || id="DRY$FAKE_ID"; FAKE_ID=$((FAKE_ID + 1))
    echo "  qsub ${args[*]} $WORK/jobs/$s.pbs   -> $id"
  else
    id="$(qsub "${args[@]}" "$WORK/jobs/$s.pbs")"
    echo "$id" > "$WORK/state/$s.jobid"
    echo "  submitted $s -> $id${dep:+ (after $dep)}"
  fi
  JOBID[$s]="$id"
}

run_local() {  # stage ntasks
  local s="$1" ntasks="$2" i
  local log="$WORK/logs/$s/local.log"
  export TMPDIR_LOCAL_BASE="${EPA_LOCAL_SCRATCH:-$WORK/scratch}"
  mkdir -p "$TMPDIR_LOCAL_BASE"
  info "local: running $s (log $log)"
  if [ "$ntasks" -gt 1 ]; then
    for i in $(seq 1 "$ntasks"); do
      TMPDIR="$TMPDIR_LOCAL_BASE" PBS_ARRAY_INDEX="$i" NCPUS="${S_NCPUS[$s]}" \
        bash "$WORK/jobs/$s.pbs" >> "$log" 2>&1 || die "local stage $s task $i failed (see $log)"
    done
  else
    TMPDIR="$TMPDIR_LOCAL_BASE" NCPUS="${S_NCPUS[$s]}" \
      bash "$WORK/jobs/$s.pbs" >> "$log" 2>&1 || die "local stage $s failed (see $log)"
  fi
}

selected() {
  [ -z "$ONLY" ] && return 0
  case ",$ONLY," in *",$1,"*) return 0 ;; *) return 1 ;; esac
}

info "mode       : $MODE"
for s in "${STAGES[@]}"; do
  selected "$s" || continue
  state="$(stage_state "$s")"
  if [ "$state" = complete ]; then echo "  $s: complete, skipping"; continue; fi
  if [ "$state" = live ]; then
    JOBID[$s]="$(cat "$WORK/state/$s.jobid")"; echo "  $s: already queued/running as ${JOBID[$s]}"; continue
  fi

  # Dependencies: complete ones need nothing; live/just-submitted ones become
  # afterok edges; anything else means an upstream stage is not available.
  deps=()
  if [ "${S_DEPS[$s]}" != "-" ]; then
    IFS=',' read -r -a _d <<< "${S_DEPS[$s]}"
    for d in "${_d[@]}"; do
      if stage_complete "$d"; then continue; fi
      if [ -n "${JOBID[$d]:-}" ]; then deps+=("${JOBID[$d]}"); continue; fi
      if [ -n "$ONLY" ]; then die "stage $s needs $d, which is neither complete nor queued"; fi
      die "internal: dependency $d of $s was not scheduled"
    done
  fi
  dep="$(IFS=:; echo "${deps[*]:-}")"

  ntasks=1
  if is_array "$s"; then
    todo_tasks "$s" > "$WORK/jobs/$s.todo"
    ntasks="$(grep -c . "$WORK/jobs/$s.todo" || true)"
    echo "  $s: $ntasks of $(grep -c . "$WORK/jobs/$s.expected") tasks to run"
  fi
  write_job "$s" "$ntasks"
  case "$MODE" in
    local)   [ "${#deps[@]}" -eq 0 ] || die "local mode: unexpected live dependency for $s"
             run_local "$s" "$ntasks" ;;
    *)       submit "$s" "$ntasks" "$dep" ;;
  esac
done

if [ "$MODE" = "submit" ]; then
  echo
  info "submitted. Check progress with:  pipeline/run.sh $WORK --status"
  info "if anything fails, fix it and run the same command again; finished work is kept."
fi
