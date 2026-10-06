# shellcheck shell=bash
# =============================================================================
# pipeline/lib/common.sh -- helpers sourced by every stage job.
# -----------------------------------------------------------------------------
# A stage job runs on ONE compute node and follows the same lifecycle every
# time (see pipeline/README.md, "Anatomy of a stage"):
#
#   1. epa_init_ws        unpack the frozen code snapshot into node-local
#                         scratch ($TMPDIR) and create an empty Datalink/
#   2. epa_in / epa_ref   unpack upstream result tars + reference files into
#                         that local Datalink/ (never onto the shared filesystem)
#   3. epa_R / epa_py     run the analysis inside the container, bound to the
#                         local workspace
#   4. epa_out            tar the declared outputs back to results/<stage>/,
#                         written atomically (.partial -> mv) with a .sha256
#   5. epa_complete       mark the stage complete (results/<stage>/.complete)
#
# The shared work directory therefore only ever sees a handful of large tars,
# never thousands of small files.
#
# Conventions:
#   * Every function fails loudly (non-zero exit) on any problem; stage scripts
#     run under `set -euo pipefail`, so a failed step fails the PBS job and PBS
#     cancels the stages that depend on it (afterok).
#   * Paths given to epa_out / epa_require are relative to the workspace
#     Datalink/ directory, exactly as the R code sees them (/workspace/Datalink).
# =============================================================================

# ---- logging -----------------------------------------------------------------
epa_log()  { printf '[%s] [%s%s] %s\n' "$(date '+%F %T')" "$EPA_STAGE" "${EPA_TASK:+#$EPA_TASK}" "$*"; }
epa_die()  { epa_log "FATAL: $*" >&2; exit 1; }
epa_step() { epa_log "---- $* ----"; }

# ---- workspace ------------------------------------------------------------------
# Creates $WS (the container's /workspace) on node-local scratch and unpacks the
# frozen code snapshot into it.
epa_init_ws() {
  : "${EPA_WORK:?}" "${EPA_STAGE:?}"
  local scratch="${TMPDIR:-}"
  [ -n "$scratch" ] || epa_die "TMPDIR is not set (expected node-local scratch inside a PBS job)"
  WS="$scratch/epa-${EPA_STAGE}${EPA_TASK:+-$EPA_TASK}"
  rm -rf "$WS"
  mkdir -p "$WS"
  tar -xf "$EPA_WORK/code/code.tar" -C "$WS"
  mkdir -p "$WS/Datalink/raw_data" "$WS/.home" "$WS/.cache"
  EPA_CENV="$WS/.epa_container_env"
  : > "$EPA_CENV"
  epa_env POWER_ANALYSIS_ROOT /workspace
  epa_env HOME /workspace/.home
  epa_env XDG_CACHE_HOME /workspace/.cache
  epa_env RENV_CONFIG_AUTOLOADER_ENABLED FALSE
  epa_env POWER_ANALYSIS_FORCE_RECALC false
  epa_env TZ UTC
  epa_env LC_ALL C.UTF-8
  epa_env MPLBACKEND Agg
  # Dates written into outputs (e.g. the exported model bundle) come from the
  # code commit, not the wall clock, so reruns give byte-identical files.
  [ -f "$EPA_WORK/code/COMMIT_TIME" ] && epa_env SOURCE_DATE_EPOCH "$(cat "$EPA_WORK/code/COMMIT_TIME")"
  # Pin BLAS/OpenMP threading so floating-point reductions do not depend on
  # how many cores a node happens to expose.
  epa_env OMP_NUM_THREADS 1
  epa_env OPENBLAS_NUM_THREADS 1
  epa_env MKL_NUM_THREADS 1
  # Number of cores this job may use (PBS exports NCPUS; local mode sets it).
  EPA_NCPUS="${NCPUS:-${EPA_LOCAL_NCPUS:-1}}"
  epa_log "workspace $WS  (ncpus=$EPA_NCPUS, host=$(hostname))"
  trap 'epa_cleanup_ws' EXIT
}

epa_cleanup_ws() {
  if [ "${EPA_KEEP_SCRATCH:-false}" = "true" ]; then
    epa_log "EPA_KEEP_SCRATCH=true: leaving $WS in place"
  else
    rm -rf "${WS:?}" 2>/dev/null || true
  fi
}

# Add KEY=VALUE to the container environment. Values may contain commas
# (e.g. DE_METHODS=deseq2,edger,limma): they are passed via --env-file, never via
# --env or qsub -v, both of which split on commas.
# Setting a key again replaces its previous value.
epa_env() {
  local key="$1" val="$2"
  case "$val" in *"'"*) epa_die "epa_env: single quote not allowed in $key" ;; esac
  grep -v "^${key}=" "$EPA_CENV" > "$EPA_CENV.tmp" || true
  printf "%s='%s'\n" "$key" "$val" >> "$EPA_CENV.tmp"
  mv -f "$EPA_CENV.tmp" "$EPA_CENV"
}

# ---- container execution -------------------------------------------------------
_epa_container_bin() {
  if [ -n "${EPA_CONTAINER_BIN:-}" ]; then echo "$EPA_CONTAINER_BIN"; return; fi
  command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null \
    || epa_die "neither apptainer nor singularity is on PATH"
}

# epa_exec [--cwd DIR] CMD ARGS...  -- run CMD in the container, cwd /workspace[/DIR]
epa_exec() {
  local cwd="/workspace"
  if [ "${1:-}" = "--cwd" ]; then cwd="/workspace/$2"; shift 2; fi
  local bin; bin="$(_epa_container_bin)"
  local binds=(--bind "$WS:/workspace")
  binds+=(--bind "$EPA_RAW_DIR:/workspace/Datalink/raw_data:ro")
  epa_log "exec: $*  (cwd $cwd)"
  "$bin" exec --cleanenv --no-home "${binds[@]}" \
    --env-file "$EPA_CENV" --pwd "$cwd" \
    "$EPA_SIF" "$@"
}

# epa_R [--cwd DIR] script.R [args...]  (script path relative to /workspace or to --cwd)
epa_R()  { if [ "${1:-}" = "--cwd" ]; then epa_exec --cwd "$2" Rscript "${@:3}"; else epa_exec Rscript "$@"; fi; }
# epa_py [--cwd DIR] script.py [args...]  -- the container's pinned venv Python
epa_py() { if [ "${1:-}" = "--cwd" ]; then epa_exec --cwd "$2" /opt/python-env/bin/python "${@:3}"; else epa_exec /opt/python-env/bin/python "$@"; fi; }

# ---- staging in ---------------------------------------------------------------
# epa_in STAGE [NAME...]  -- unpack results/STAGE/NAME.tar (default: every tar of
# that stage) into $WS/Datalink. Each tar's checksum is verified first.
epa_in() {
  local stage="$1"; shift
  local dir="$EPA_WORK/results/$stage"
  [ -f "$dir/.complete" ] || [ "${EPA_IN_ALLOW_INCOMPLETE:-false}" = "true" ] \
    || epa_die "upstream stage '$stage' is not complete ($dir/.complete missing)"
  local tars=()
  if [ "$#" -gt 0 ]; then
    local n; for n in "$@"; do tars+=("$dir/$n.tar"); done
  else
    local t; for t in "$dir"/*.tar; do [ -e "$t" ] && tars+=("$t"); done
  fi
  [ "${#tars[@]}" -gt 0 ] || epa_die "no result tars found for stage '$stage' in $dir"
  local t
  for t in "${tars[@]}"; do
    [ -f "$t" ] || epa_die "missing result tar $t"
    [ "${EPA_VERIFY_TARS:-true}" = true ] && _epa_verify_sha "$t"
    epa_log "stage-in $(basename "$t") ($(du -h "$t" | cut -f1))"
    tar -xf "$t" -C "$WS/Datalink"
  done
}

# epa_in_array STAGE [PATTERN]  -- unpack the task tars of a completed array
# stage (all of them, or those whose name contains PATTERN). Fails unless every
# expected task tar is present.
epa_in_array() {
  local stage="$1" pat="${2:-}" dir="$EPA_WORK/results/$1" id out n=0
  [ -f "$EPA_WORK/jobs/$stage.expected" ] || epa_die "no task list for array stage $stage"
  while read -r id out; do
    [ -f "$dir/$out.tar" ] || epa_die "array stage $stage incomplete: $out.tar missing"
  done < "$EPA_WORK/jobs/$stage.expected"
  while read -r id out; do
    case "$out" in *"$pat"*) ;; *) continue ;; esac
    [ "${EPA_VERIFY_TARS:-true}" = true ] && _epa_verify_sha "$dir/$out.tar"
    tar -xf "$dir/$out.tar" -C "$WS/Datalink"
    n=$((n + 1))
  done < "$EPA_WORK/jobs/$stage.expected"
  [ "$n" -gt 0 ] || epa_die "epa_in_array: no $stage tasks match '$pat'"
  epa_log "stage-in $stage: $n task tars${pat:+ matching '$pat'}"
}

# epa_ref SRC DEST  -- copy a reference-bundle file (path relative to the
# reference/ dir) to a path relative to $WS/Datalink. Small files only.
epa_ref() {
  local src="$EPA_REF_DIR/$1" dest="$WS/Datalink/$2"
  [ -f "$src" ] || epa_die "reference file missing: $src"
  mkdir -p "$(dirname "$dest")"
  cp "$src" "$dest"
}

# epa_ref_dir -- make the whole reference bundle visible at /workspace/reference.
epa_ref_dir() {
  mkdir -p "$WS/reference"
  cp -r "$EPA_REF_DIR/." "$WS/reference/"
}

_epa_verify_sha() {
  local t="$1"
  [ -f "$t.sha256" ] || epa_die "checksum file missing for $t"
  local want got
  want="$(cut -d' ' -f1 < "$t.sha256")"
  got="$(sha256sum "$t" | cut -d' ' -f1)"
  [ "$want" = "$got" ] || epa_die "checksum mismatch for $t (want $want got $got)"
}

# ---- checks ----------------------------------------------------------------------
# epa_require PATH...  -- fail unless each path (relative to Datalink) exists.
epa_require() {
  local p missing=0
  for p in "$@"; do
    if [ ! -e "$WS/Datalink/$p" ]; then epa_log "missing expected output: Datalink/$p" >&2; missing=1; fi
  done
  [ "$missing" -eq 0 ] || epa_die "expected outputs missing (see above)"
}

# ---- staging out ---------------------------------------------------------------
# epa_out NAME PATH...  -- tar PATHs (relative to Datalink, globs allowed) into
# results/<stage>/NAME.tar, atomically, with a sha256 sidecar.
epa_out() {
  local name="$1"; shift
  local dir="$EPA_WORK/results/$EPA_STAGE"
  mkdir -p "$dir"
  local final="$dir/$name.tar"
  local tmp="$dir/.$name.$$.partial.tar"
  local list="$WS/.out_$name.lst"
  : > "$list"
  local p m found
  ( cd "$WS/Datalink"
    for p in "$@"; do
      found=0
      # shellcheck disable=SC2086
      for m in $p; do
        if [ -e "$m" ]; then printf '%s\0' "$m" >> "$list"; found=1; fi
      done
      [ "$found" -eq 1 ] || { echo "epa_out: nothing matches '$p'" >&2; exit 1; }
    done
  ) || epa_die "epa_out $name: declared output missing"
  ( cd "$WS/Datalink" && tar --null -T "$list" --sort=name --owner=0 --group=0 --numeric-owner -cf "$tmp" )
  sha256sum "$tmp" | sed "s#$tmp#$name.tar#" > "$tmp.sha256"
  mv -f "$tmp.sha256" "$final.sha256"
  mv -f "$tmp" "$final"
  epa_log "stage-out $name.tar ($(du -h "$final" | cut -f1))"
}

# epa_complete -- mark the stage complete. Only call after every epa_out.
epa_complete() {
  local dir="$EPA_WORK/results/$EPA_STAGE"
  mkdir -p "$dir"
  {
    echo "stage=$EPA_STAGE"
    echo "finished=$(date -u +%FT%TZ)"
    echo "host=$(hostname)"
    echo "code_commit=$(cat "$EPA_WORK/code/COMMIT")"
    echo "sif_sha256=$(cat "$EPA_WORK/code/SIF_SHA256" 2>/dev/null || echo unknown)"
    ( cd "$dir" && cat ./*.sha256 2>/dev/null ) || true
  } > "$dir/.complete.tmp"
  mv -f "$dir/.complete.tmp" "$dir/.complete"
  epa_log "stage complete"
}
