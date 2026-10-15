#!/usr/bin/env bash
# =============================================================================
# pipeline/tests/test_pipeline.sh -- behavioural tests for the orchestration
# layer (run.sh + lib/common.sh), with qsub / qstat / the container stubbed.
# Runs anywhere in a few seconds; no cluster, no container, no data.
#
#   bash pipeline/tests/test_pipeline.sh
#
# Covers: fresh submission (dependency edges, array sizes, generated job
# scripts with no qsub -v), dry run, resume after partial completion (only the
# missing grid chunks are resubmitted, complete stages are skipped and dropped
# from dependency lists), reuse of queued/running jobs, refusal on missing
# inputs or a dirty repository, the in-job helpers (env file with commas,
# stage-in checksum verification, atomic stage-out, array completeness), and a
# full local-mode job with a fake container.
# =============================================================================
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SRC_PIPELINE="$(cd -- "$HERE/.." && pwd)"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
PASS=0 FAIL=0
ok()   { PASS=$((PASS + 1)); echo "  ok   $*"; }
bad()  { FAIL=$((FAIL + 1)); echo "  FAIL $*"; }
check() { local d="$1"; shift; if "$@"; then ok "$d"; else bad "$d"; fi; }

# ---- a throwaway repository containing this pipeline -----------------------------
REPO="$T/repo"
mkdir -p "$REPO/R" "$REPO/Datalink"
cp -r "$SRC_PIPELINE" "$REPO/pipeline"
rm -rf "$REPO/pipeline/tests"
echo 'cat("hello\n")' > "$REPO/R/hello.R"
# fake reference bundle + raw inputs with matching manifests
REFSRC="$T/refsrc"; mkdir -p "$REFSRC/biomart"
for f in biomart/a.csv b.txt; do echo "ref $f" > "$REFSRC/$f"; done
( cd "$REFSRC" && sha256sum biomart/a.csv b.txt ) > "$REPO/pipeline/reference/MANIFEST.sha256"
RAWSRC="$T/rawsrc"; mkdir -p "$RAWSRC"
for f in PsychAD-HBCC.h5ad media-1.csv; do echo "raw $f" > "$RAWSRC/$f"; done
( cd "$RAWSRC" && sha256sum PsychAD-HBCC.h5ad media-1.csv ) > "$REPO/pipeline/inputs.sha256"
( cd "$REPO" && git init -q && git add -A && git -c user.email=t@t -c user.name=t commit -qm init )

new_work() {
  local w="$T/work-$1"
  mkdir -p "$w"
  cp -r "$RAWSRC" "$w/raw"; cp -r "$REFSRC" "$w/reference"
  echo "fake image" > "$w/epa.sif"
  cat > "$w/epa.env" <<EOF
EPA_GRID_CHUNKS=3
EPA_AD_GRID_CHUNKS=2
EPA_GRID_CONCURRENCY=2
EPA_N_RUNS=2
EPA_DE_METHODS="deseq2,edger,limma"
EOF
  echo "$w"
}

# ---- stub scheduler ----------------------------------------------------------------
BIN="$T/bin"; mkdir -p "$BIN"
cat > "$BIN/qsub" <<'EOF'
#!/usr/bin/env bash
# records "<id> <args...>" and prints a PBS-style id (arrays get the [] suffix)
n=$(( $(wc -l < "$QSUB_LOG" 2>/dev/null || echo 0) + 1000 ))
script="${@: -1}"
if grep -q '^#PBS -J' "$script"; then id="${n}[].pbs"; else id="${n}.pbs"; fi
echo "$id $*" >> "$QSUB_LOG"
echo "$id"
EOF
cat > "$BIN/qstat" <<'EOF'
#!/usr/bin/env bash
# job_state comes from $QSTAT_STATE/<id> (default: finished)
id="${@: -1}"
if [ -n "${QSTAT_FAIL:-}" ]; then echo "qstat: cannot connect to server" >&2; exit 1; fi
if [ -n "${QSTAT_UNKNOWN:-}" ]; then echo "qstat: Unknown Job Id $id" >&2; exit 153; fi
st="F"; [ -f "$QSTAT_STATE/$id" ] && st="$(cat "$QSTAT_STATE/$id")"
echo "Job Id: $id"; echo "    job_state = $st"
EOF
chmod +x "$BIN"/*
export PATH="$BIN:$PATH" QSTAT_STATE="$T/qstat"; mkdir -p "$QSTAT_STATE"

RUN="$REPO/pipeline/run.sh"
submitted() { awk '{print $NF}' "$QSUB_LOG" | xargs -n1 basename | sed 's/\.pbs$//'; }
dep_of() {  # stage -> its depend string ("" if none)
  grep "/jobs/$1.pbs\$" "$QSUB_LOG" | grep -o 'depend=afterok:[^ ]*' | sed 's/depend=afterok://' || true
}
id_of() { grep "/jobs/$1.pbs\$" "$QSUB_LOG" | awk '{print $1}'; }

echo "== fresh submission"
W="$(new_work fresh)"; export QSUB_LOG="$T/qsub-fresh.log"; : > "$QSUB_LOG"
"$RUN" "$W" > "$T/out-fresh.txt" 2>&1 || { cat "$T/out-fresh.txt"; exit 1; }
n_stages="$(grep -cvE '^(#|$)' "$REPO/pipeline/stages.tsv")"
check "every stage submitted once ($n_stages)" test "$(submitted | wc -l)" -eq "$n_stages"
check "shard has no dependency" test -z "$(dep_of shard)"
check "grid depends on shard and gt_manifest" test "$(dep_of grid)" = "$(id_of shard):$(id_of gt_manifest)"
check "pcf_aggregate depends on the grid array id" grep -qF "$(id_of grid)" <<< "$(dep_of pcf_aggregate)"
check "grid array has 2 x chunks tasks" grep -q '^#PBS -J 1-6%2$' "$W/jobs/grid.pbs"
check "xgb_search array of 108" grep -q '^#PBS -J 1-108%' "$W/jobs/xgb_search.pbs"
check "ad_grid array has AD chunk tasks" grep -q '^#PBS -J 1-2%2$' "$W/jobs/ad_grid.pbs"
check "ad_grid task list is cells-only" test "$(cat "$W/jobs/ad_grid.todo")" = "$(printf 'cells:1 epa_cells_chunk_1\ncells:2 epa_cells_chunk_2')"
ok_dep=true
for s in ad_cartesian ad_de; do grep -qF "$(id_of "$s")" <<< "$(dep_of finalise)" || ok_dep=false; done
check "finalise waits for the AD stages" "$ok_dep"
check "commit time recorded for SOURCE_DATE_EPOCH" grep -Eq '^[0-9]+$' "$W/code/COMMIT_TIME"
check "no job uses qsub -v" bash -c "! grep -q -- ' -v ' '$QSUB_LOG'"
check "job scripts carry their resources" grep -q '^#PBS -l select=1:ncpus=24:mem=200gb$' "$W/jobs/shard.pbs"
check "code snapshot excludes Datalink" bash -c "! tar -tf '$W/code/code.tar' | grep -q '^Datalink'"
check "code snapshot has the pipeline" bash -c "tar -tf '$W/code/code.tar' | grep -q '^pipeline/stages/grid.sh$'"
check "commas survive into resolved env" grep -q "^EPA_DE_METHODS='deseq2,edger,limma'$" "$W/epa.env.resolved"
check "grid to-do lists all tasks" test "$(wc -l < "$W/jobs/grid.todo")" -eq 6
check "job ids recorded" test -f "$W/state/knit.jobid"

echo "== re-run while everything is still queued: nothing resubmitted"
for s in $(submitted); do echo Q > "$QSTAT_STATE/$(id_of "$s")"; done
before="$(wc -l < "$QSUB_LOG")"
"$RUN" "$W" > "$T/out-requeue.txt" 2>&1
check "no new submissions" test "$(wc -l < "$QSUB_LOG")" -eq "$before"
check "reports queued jobs" grep -q "already queued/running" "$T/out-requeue.txt"

echo "== qstat failure stops run.sh; expired job history allows resubmission"
before="$(wc -l < "$QSUB_LOG")"
if QSTAT_FAIL=1 "$RUN" "$W" > "$T/out-qfail.txt" 2>&1; then bad "qstat failure ignored"; else check "qstat failure refused" grep -q "could not query job" "$T/out-qfail.txt"; fi
check "nothing submitted on qstat failure" test "$(wc -l < "$QSUB_LOG")" -eq "$before"
QSTAT_UNKNOWN=1 "$RUN" "$W" > "$T/out-qunk.txt" 2>&1
check "unknown jobs resubmitted" test "$(wc -l < "$QSUB_LOG")" -eq "$((before + n_stages))"

echo "== resume: upstream done, 2 grid chunks missing, grid job gone"
W2="$(new_work resume)"; export QSUB_LOG="$T/qsub-resume.log"; : > "$QSUB_LOG"
for s in shard gt_manifest eda fgt deg_loss; do mkdir -p "$W2/results/$s"; touch "$W2/results/$s/.complete"; done
mkdir -p "$W2/results/grid"
for t in epa_reads_chunk_1 epa_reads_chunk_2 epa_cells_chunk_1 epa_cells_chunk_3; do touch "$W2/results/grid/$t.tar"; done
"$RUN" "$W2" > "$T/out-resume.txt" 2>&1 || { cat "$T/out-resume.txt"; exit 1; }
check "complete stages skipped" bash -c "! grep -q '/jobs/shard.pbs\$' '$QSUB_LOG'"
check "grid resubmitted with only missing chunks" test "$(cat "$W2/jobs/grid.todo")" = "$(printf 'reads:3 epa_reads_chunk_3\ncells:2 epa_cells_chunk_2')"
check "grid now has no dependency (upstream complete)" test -z "$(dep_of grid)"
check "2-task array" grep -q '^#PBS -J 1-2%2$' "$W2/jobs/grid.pbs"

echo "== single missing chunk is submitted as a plain job"
W3="$(new_work single)"; export QSUB_LOG="$T/qsub-single.log"; : > "$QSUB_LOG"
for s in shard gt_manifest; do mkdir -p "$W3/results/$s"; touch "$W3/results/$s/.complete"; done
mkdir -p "$W3/results/grid"
for m in reads cells; do for c in 1 2 3; do touch "$W3/results/grid/epa_${m}_chunk_${c}.tar"; done; done
rm "$W3/results/grid/epa_cells_chunk_2.tar"
"$RUN" "$W3" --only grid > "$T/out-single.txt" 2>&1
check "no -J for one task" bash -c "! grep -q '^#PBS -J' '$W3/jobs/grid.pbs'"
check "--only submits just that stage" test "$(submitted)" = "grid"

echo "== --only refuses when a dependency is neither complete nor queued"
W4="$(new_work only)"; export QSUB_LOG="$T/qsub-only.log"; : > "$QSUB_LOG"
if "$RUN" "$W4" --only grid > "$T/out-only.txt" 2>&1; then bad "--only should fail"; else ok "--only fails on missing upstream"; fi

echo "== dry run submits nothing"
W5="$(new_work dry)"; export QSUB_LOG="$T/qsub-dry.log"; : > "$QSUB_LOG"
"$RUN" "$W5" --dry-run > "$T/out-dry.txt" 2>&1
check "qsub never called" test ! -s "$QSUB_LOG"
check "dry run prints the chain" grep -q "depend=afterok:DRY" "$T/out-dry.txt"

echo "== preflight refusals"
W6="$(new_work missing)"; rm "$W6/raw/media-1.csv"
if "$RUN" "$W6" > "$T/out-missing.txt" 2>&1; then bad "missing raw input accepted"; else check "missing raw input refused" grep -q "missing raw input" "$T/out-missing.txt"; fi
W7="$(new_work badref)"; echo tampered >> "$W7/reference/b.txt"
if "$RUN" "$W7" > "$T/out-badref.txt" 2>&1; then bad "tampered reference accepted"; else ok "tampered reference refused"; fi
W8="$(new_work dirty)"; echo "# edit" >> "$REPO/R/hello.R"
if "$RUN" "$W8" > "$T/out-dirty.txt" 2>&1; then bad "dirty repo accepted"; else check "dirty repo refused" grep -q "uncommitted changes" "$T/out-dirty.txt"; fi
( cd "$REPO" && git checkout -q R/hello.R )

echo "== status"
"$RUN" "$W2" --status > "$T/out-status.txt"
check "status shows grid progress" grep -Eq '^grid +(live|pending|failed-or-killed) +4/6 tasks done' "$T/out-status.txt"

echo "== in-job helpers (lib/common.sh) under set -euo pipefail"
H="$T/helpers"; mkdir -p "$H/work/results/up" "$H/work/code" "$H/tmp" "$H/work/jobs"
( mkdir -p "$H/src/Datalink/A" && echo 1 > "$H/src/Datalink/A/x.txt" && cd "$H/src/Datalink" && tar -cf "$H/work/results/up/data.tar" A )
( cd "$H/work/results/up" && sha256sum data.tar > data.tar.sha256 && touch .complete )
( cd "$REPO" && git archive --format=tar -o "$H/work/code/code.tar" HEAD ) && echo testcommit > "$H/work/code/COMMIT"
cat > "$H/fake-apptainer" <<'EOF'
#!/usr/bin/env bash
# fake container: record the env file, then run "Rscript <file>" as a no-op that
# creates Datalink/out/<basename>.done inside the bound workspace
ws=""; envf=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    exec|--cleanenv|--no-home) shift ;;
    --bind) case "$2" in *:/workspace) ws="${2%%:/workspace}" ;; esac; shift 2 ;;
    --env-file) envf="$2"; shift 2 ;;
    --pwd) shift 2 ;;
    *.sif) shift; break ;;
    *) shift ;;
  esac
done
cp "$envf" "$ws/.last_env"
mkdir -p "$ws/Datalink/out"; touch "$ws/Datalink/out/$(basename "${2:-x}").done"
EOF
chmod +x "$H/fake-apptainer"
(
  set -euo pipefail
  export EPA_WORK="$H/work" EPA_STAGE=demo TMPDIR="$H/tmp" EPA_CONTAINER_BIN="$H/fake-apptainer"
  export EPA_SIF="$H/x.sif" EPA_RAW_DIR="$H" EPA_REF_DIR="$H"
  # shellcheck source=/dev/null
  source "$REPO/pipeline/lib/common.sh"
  epa_init_ws
  epa_env DE_METHODS "deseq2,edger,limma"
  epa_env DE_METHODS "deseq2,edger"
  epa_in up
  epa_R R/hello.R
  epa_require A/x.txt out/hello.R.done
  epa_out result A out
  epa_complete
  grep -q "^DE_METHODS='deseq2,edger'$" "$WS/.last_env"
  [ "$(grep -c '^DE_METHODS=' "$WS/.last_env")" -eq 1 ]
) > "$T/out-helpers.txt" 2>&1 && ok "stage-in, container env (commas, overrides), stage-out, complete" \
  || { bad "helper lifecycle"; cat "$T/out-helpers.txt"; }
check "stage-out tar has both paths" bash -c "tar -tf '$H/work/results/demo/result.tar' | grep -q '^A/x.txt$' && tar -tf '$H/work/results/demo/result.tar' | grep -q '^out/hello.R.done$'"
check "stage-out checksum valid" bash -c "cd '$H/work/results/demo' && sha256sum -c --quiet result.tar.sha256"
check "no partial files left" bash -c "! ls -a '$H/work/results/demo' | grep -q partial"
check "scratch cleaned up" bash -c "! ls '$H/tmp' | grep -q epa-demo"
echo corrupt >> "$H/work/results/up/data.tar"
( set -euo pipefail; export EPA_WORK="$H/work" EPA_STAGE=demo2 TMPDIR="$H/tmp" EPA_SIF=x EPA_RAW_DIR="$H" EPA_REF_DIR="$H"
  source "$REPO/pipeline/lib/common.sh"; epa_init_ws; epa_in up ) > "$T/out-corrupt.txt" 2>&1 \
  && bad "corrupt tar accepted" || check "corrupt upstream tar refused" grep -q "checksum mismatch" "$T/out-corrupt.txt"
( set -euo pipefail; export EPA_WORK="$H/work" EPA_STAGE=demo3 TMPDIR="$H/tmp" EPA_SIF=x EPA_RAW_DIR="$H" EPA_REF_DIR="$H"
  source "$REPO/pipeline/lib/common.sh"; epa_init_ws; epa_out nothing does/not/exist ) > "$T/out-noout.txt" 2>&1 \
  && bad "missing output accepted" || ok "missing declared output fails the job"
printf 'reads:1 epa_reads_chunk_1\nreads:2 epa_reads_chunk_2\n' > "$H/work/jobs/arr.expected"
mkdir -p "$H/work/results/arr"; ( cd "$H/src/Datalink" && tar -cf "$H/work/results/arr/epa_reads_chunk_1.tar" A ); \
  ( cd "$H/work/results/arr" && sha256sum epa_reads_chunk_1.tar > epa_reads_chunk_1.tar.sha256 )
( set -euo pipefail; export EPA_WORK="$H/work" EPA_STAGE=demo4 TMPDIR="$H/tmp" EPA_SIF=x EPA_RAW_DIR="$H" EPA_REF_DIR="$H"
  source "$REPO/pipeline/lib/common.sh"; epa_init_ws; epa_in_array arr ) > "$T/out-arr.txt" 2>&1 \
  && bad "incomplete array accepted" || check "incomplete array refused" grep -q "incomplete" "$T/out-arr.txt"

echo "== local mode runs a real generated job script end to end"
W9="$(new_work local)"; export QSUB_LOG="$T/qsub-local.log"; : > "$QSUB_LOG"
mkdir -p "$REPO/pipeline/stages"
cat > "$REPO/pipeline/stages/zz_demo.sh" <<'EOF'
stage_main() {
  epa_init_ws
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_R R/hello.R
  epa_out demo out
  epa_complete
}
EOF
cp "$REPO/pipeline/stages.tsv" "$T/stages.tsv.bak"
printf '# test\nzz_demo\t-\t1\t1gb\t00:10:00\tsingle\n' > "$REPO/pipeline/stages.tsv"
( cd "$REPO" && git add -A && git -c user.email=t@t -c user.name=t commit -qm demo )
echo "EPA_CONTAINER_BIN=$H/fake-apptainer" >> "$W9/epa.env"
EPA_CONTAINER_BIN="$H/fake-apptainer" PATH="$H:$PATH" "$RUN" "$W9" --local > "$T/out-local.txt" 2>&1 \
  || { cat "$T/out-local.txt"; cat "$W9/logs/zz_demo/local.log"; }
check "local job completed" test -f "$W9/results/zz_demo/.complete"
check "local job output tar" bash -c "tar -tf '$W9/results/zz_demo/demo.tar' | grep -q 'out/hello.R.done'"
"$RUN" "$W9" --local > "$T/out-local2.txt" 2>&1
check "local re-run skips completed stage" grep -q "zz_demo: complete, skipping" "$T/out-local2.txt"

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
