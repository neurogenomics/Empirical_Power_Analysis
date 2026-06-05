#!/usr/bin/env bash
# =============================================================================
# docker-all.sh — build -> smoke -> run, end to end.
#
# Thin wrapper that chains the three standalone scripts in order, stopping at
# the first failure. Each phase can be skipped, so the common flows are simple:
#
#   ./docker-all.sh                          # fresh: build + smoke + run (pcf)
#   ./docker-all.sh --skip-build             # already built: smoke + run
#   ./docker-all.sh --skip-build --skip-smoke   # just run
#   ./docker-all.sh --skip-run               # build + smoke only (verify image)
#   ./docker-all.sh --pipeline all --max-cores 24
#
# All build/smoke/run options are forwarded to the relevant phase:
#   --image NAME:TAG     (all phases)
#   --datalink PATH      (smoke, run)
#   --pipeline pcf|all   (run)
#   --max-cores N        (run)
#   --no-overwrite|--force, --mode cells|reads|both, --extra-arg "STR" (run)
#   --dry-run            print each phase's command, run nothing
#
# Unattended:
#   nohup ./docker-all.sh --pipeline all --max-cores 24 > nohup.docker-all.out 2>&1 & disown
#
# Exit code: 0 if all run phases succeed; otherwise the failing phase's code.
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

IMAGE="empirical-power-analysis:overnight"
DATALINK="${SCRIPT_DIR}/Datalink"
PIPELINE="pcf"
MAX_CORES="auto"
OVERWRITE="--no-overwrite"
MODE="both"
DRY_RUN=0
SKIP_BUILD=0
SKIP_SMOKE=0
SKIP_RUN=0
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-build)   SKIP_BUILD=1; shift ;;
        --skip-smoke)   SKIP_SMOKE=1; shift ;;
        --skip-run)     SKIP_RUN=1; shift ;;
        --image)        IMAGE="$2"; shift 2 ;;
        --datalink)     DATALINK="$2"; shift 2 ;;
        --pipeline)     PIPELINE="$2"; shift 2 ;;
        --max-cores)    MAX_CORES="$2"; shift 2 ;;
        --no-overwrite) OVERWRITE="--no-overwrite"; shift ;;
        --force)        OVERWRITE="--force"; shift ;;
        --mode)         MODE="$2"; shift 2 ;;
        --extra-arg)    EXTRA_ARGS+=("$2"); shift 2 ;;
        --dry-run)      DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown flag: $1 (try --help)" >&2; exit 2 ;;
    esac
done

phase() { printf '\n=== %s ===\n' "$1"; }

# ---- Phase 1: build ---------------------------------------------------------
if [[ "${SKIP_BUILD}" -eq 0 ]]; then
    phase "build"
    CMD=( "${SCRIPT_DIR}/docker-build.sh" --image "${IMAGE}" )
    if [[ "${DRY_RUN}" -eq 1 ]]; then echo "[dry-run] ${CMD[*]}"; else "${CMD[@]}" || { rc=$?; echo "docker-all: build failed (exit $rc)" >&2; exit "$rc"; }; fi
else
    echo "Skipping build (--skip-build)."
fi

# ---- Phase 2: smoke ---------------------------------------------------------
if [[ "${SKIP_SMOKE}" -eq 0 ]]; then
    phase "smoke"
    CMD=( "${SCRIPT_DIR}/docker-smoke.sh" --image "${IMAGE}" --datalink "${DATALINK}" )
    if [[ "${DRY_RUN}" -eq 1 ]]; then echo "[dry-run] ${CMD[*]}"; else "${CMD[@]}" || { rc=$?; echo "docker-all: smoke failed (exit $rc)" >&2; exit "$rc"; }; fi
else
    echo "Skipping smoke (--skip-smoke)."
fi

# ---- Phase 3: run -----------------------------------------------------------
if [[ "${SKIP_RUN}" -eq 0 ]]; then
    phase "run"
    CMD=( "${SCRIPT_DIR}/docker-run.sh"
          --image "${IMAGE}" --datalink "${DATALINK}"
          --pipeline "${PIPELINE}" --max-cores "${MAX_CORES}"
          "${OVERWRITE}" --mode "${MODE}" )
    for a in ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}; do CMD+=( --extra-arg "$a" ); done
    if [[ "${DRY_RUN}" -eq 1 ]]; then CMD+=( --dry-run ); fi
    "${CMD[@]}" || { rc=$?; echo "docker-all: run finished with exit $rc" >&2; exit "$rc"; }
else
    echo "Skipping run (--skip-run)."
fi

echo
echo "docker-all: done."
