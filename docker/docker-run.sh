#!/usr/bin/env bash
# =============================================================================
# docker-run.sh — run the pipeline inside the built image (Phase 3).
#
# The standalone run step that pairs with docker-build.sh + docker-smoke.sh.
# Runs the image with the same runtime contract as the smoke test (repo bound
# on /workspace, Datalink on /workspace/Datalink, dropped to the host UID:GID
# via gosu), dispatching the requested pipeline. Foreground `docker run --rm`,
# so PID-1 tini forwards SIGINT/SIGTERM into the container and the exit code is
# the pipeline's own — no detached-container bookkeeping needed.
#
# Writes a provenance manifest, translates the container exit code into a
# meaningful status, and tees the run to a log file. No heartbeat / healthcheck
# / auto-restart (deliberately omitted).
#
# Usage (run from the repo root):
#   ./docker/docker-run.sh [--pipeline pcf|all] [--max-cores N] [--no-overwrite|--force]
#                   [--mode cells|reads|both] [--image NAME:TAG] [--datalink PATH]
#                   [--extra-arg "STR" ...] [--mini-run] [--dry-run]
#
#   --mini-run  Tiny end-to-end smoke test: forces --pipeline pcf and runs
#               cartesian_downsampling -> depth_equalization -> the two depth-eq
#               variant cartesian stages (DepthEq / ShallowFewer) -> power_curve_fitting
#               on a capped donor universe (<=10 donors/experiment, Astro cell
#               type). Variant shards are written to MINI_-prefixed dataset dirs
#               and all outputs land under Datalink/MINI_*, so production data is
#               never touched. Use to verify the container works before committing
#               to a real multi-week run. Implies --force and skips only sharding
#               (needs the raw atlas) and figure_composition (needs all panels).
#
# Unattended:
#   nohup ./docker/docker-run.sh --pipeline all --max-cores 24 > nohup.docker-run.out 2>&1 & disown
#
# Exit codes:
#   0  success
#   1  pipeline ran but reported partial failure (container exit 2)
#   2  pre-flight failed (no daemon, missing image, bad args)
#   6  container crashed (OOM / signal / other non-zero)
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# This script lives in <repo>/docker/; the repo root (bind-mounted on /workspace,
# and where logs/ + the git metadata live) is the directory above it.
REPO_ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"

IMAGE="empirical-power-analysis:overnight"
DATALINK="${REPO_ROOT}/Datalink"
PIPELINE="pcf"
MAX_CORES="auto"
OVERWRITE="--no-overwrite"
MODE="both"
DRY_RUN=0
MINI_RUN=0
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --image)        IMAGE="$2"; shift 2 ;;
        --datalink)     DATALINK="$2"; shift 2 ;;
        --pipeline)     PIPELINE="$2"; shift 2 ;;
        --max-cores)    MAX_CORES="$2"; shift 2 ;;
        --no-overwrite) OVERWRITE="--no-overwrite"; shift ;;
        --force)        OVERWRITE="--force"; shift ;;
        --mode)         MODE="$2"; shift 2 ;;
        --extra-arg)    EXTRA_ARGS+=("$2"); shift 2 ;;
        --mini-run)     MINI_RUN=1; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown flag: $1 (try --help)" >&2; exit 2 ;;
    esac
done

# ---- mini-run: force a tiny, isolated cartesian -> ... -> pcf smoke test -----
# Forces the power-curve pipeline, regenerates (so the smoke always executes),
# and skips only the stages that can't run from existing Roussos shards:
# sharding (needs the raw atlas) and figure_composition (needs all panels). This
# leaves cartesian_downsampling -> depth_equalization -> the two depth-eq variant
# cartesian stages -> power_curve_fitting, exercising the DepthEq / ShallowFewer
# paths too. MINI_RUN=true is passed into the container, where config.R shrinks
# the grid, caps the donor universe, writes variant shards to MINI_-prefixed
# dataset dirs, and redirects all outputs to Datalink/MINI_*.
MINI_SKIP="sharding,figure_composition"
if [[ "${MINI_RUN}" -eq 1 ]]; then
    PIPELINE="pcf"
    OVERWRITE="--force"
    echo "Mini-run: forcing --pipeline pcf, --force, and MINI_RUN=true in the container."
fi

# ---- validate args ----------------------------------------------------------
case "${PIPELINE}" in pcf|all) ;; *) echo "ERROR: --pipeline must be pcf|all, got '${PIPELINE}'" >&2; exit 2 ;; esac
case "${MODE}"     in cells|reads|both) ;; *) echo "ERROR: --mode must be cells|reads|both, got '${MODE}'" >&2; exit 2 ;; esac

NPROC="$(nproc 2>/dev/null || echo 1)"
if [[ "${MAX_CORES}" == "auto" ]]; then MAX_CORES=$(( NPROC < 24 ? NPROC : 24 )); fi
if ! [[ "${MAX_CORES}" =~ ^[1-9][0-9]*$ ]]; then echo "ERROR: --max-cores must be a positive integer, got '${MAX_CORES}'" >&2; exit 2; fi
if (( MAX_CORES > NPROC )); then echo "WARN: --max-cores ${MAX_CORES} > nproc ${NPROC}; clamping to ${NPROC}." >&2; MAX_CORES="${NPROC}"; fi

# ---- preflight --------------------------------------------------------------
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: cannot reach the Docker daemon (see docker-build.sh for hints)." >&2
    exit 2
fi
if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
    echo "ERROR: image '${IMAGE}' not found. Build it first: ./docker/docker-build.sh" >&2
    exit 2
fi
mkdir -p "${DATALINK}"

# ---- logging + provenance ---------------------------------------------------
TS="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_ROOT="${REPO_ROOT}/logs/docker-run/${TS}"
mkdir -p "${LOG_ROOT}"
ln -sfn "${TS}" "${REPO_ROOT}/logs/docker-run/latest"
RUN_LOG="${LOG_ROOT}/run.log"
MANIFEST="${LOG_ROOT}/manifest.json"

GIT_SHA="$(cd "${REPO_ROOT}" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOST="$(hostname)"
IMAGE_DIGEST="$(docker image inspect --format '{{.Id}}' "${IMAGE}" 2>/dev/null || echo unknown)"
START="" ; FINISH="" ; EXIT_CODE="" ; STATUS="pending"

# Wrapper args forwarded to the in-container pipeline driver.
WRAPPER_ARGS=( "--max-cores" "${MAX_CORES}" "${OVERWRITE}" "--mode" "${MODE}" )
if [[ "${MINI_RUN}" -eq 1 ]]; then WRAPPER_ARGS+=( "--skip" "${MINI_SKIP}" ); fi
if (( ${#EXTRA_ARGS[@]} > 0 )); then WRAPPER_ARGS+=( "${EXTRA_ARGS[@]}" ); fi

# Minimal JSON-array encoder for the (controlled) extra args.
extra_json() {
    if (( ${#EXTRA_ARGS[@]} == 0 )); then printf '[]'; return; fi
    local out="" a
    for a in "${EXTRA_ARGS[@]}"; do out+="\"${a//\"/\\\"}\","; done
    printf '[%s]' "${out%,}"
}

write_manifest() {
    cat > "${MANIFEST}" <<EOF
{
  "schema_version": 1,
  "timestamp_utc":  "${TS}",
  "host":           "${HOST}",
  "git_sha":        "${GIT_SHA}",
  "pipeline":       "${PIPELINE}",
  "image_tag":      "${IMAGE}",
  "image_digest":   "${IMAGE_DIGEST}",
  "datalink_path":  "${DATALINK}",
  "max_cores":      ${MAX_CORES},
  "host_nproc":     ${NPROC},
  "overwrite_flag": "${OVERWRITE}",
  "mode":           "${MODE}",
  "mini_run":       $( [[ "${MINI_RUN}" -eq 1 ]] && echo true || echo false ),
  "extra_args":     $(extra_json),
  "started":        "${START}",
  "finished":       "${FINISH}",
  "exit_code":      "${EXIT_CODE}",
  "exit_status":    "${STATUS}",
  "log_file":       "${RUN_LOG}"
}
EOF
}

# If we're interrupted, record it before exiting (foreground docker run will
# already have forwarded the signal to the container via tini).
on_signal() {
    STATUS="aborted"; FINISH="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; write_manifest
    echo "docker-run.sh: interrupted; manifest at ${MANIFEST}" >&2
    exit 130
}
trap on_signal INT TERM

# ---- summary ----------------------------------------------------------------
echo "Pipeline:   ${PIPELINE}$( [[ "${MINI_RUN}" -eq 1 ]] && echo '  (MINI smoke test)' )"
echo "Image:      ${IMAGE} (${IMAGE_DIGEST})"
echo "Datalink:   ${DATALINK}"
echo "Args:       ${WRAPPER_ARGS[*]}"
echo "Log:        ${RUN_LOG}"
echo "Manifest:   ${MANIFEST}"

DOCKER_CMD=(
    docker run --rm
    -e "PIPELINE=${PIPELINE}"
    -e "POWER_ANALYSIS_ROOT=/workspace"
    -e "HOME=/tmp"
    -e "TARGET_UID=$(id -u)"
    -e "TARGET_GID=$(id -g)"
)
if [[ "${MINI_RUN}" -eq 1 ]]; then DOCKER_CMD+=( -e "MINI_RUN=true" ); fi
DOCKER_CMD+=(
    -v "${REPO_ROOT}:/workspace"
    -v "${DATALINK}:/workspace/Datalink"
    "${IMAGE}"
    "${WRAPPER_ARGS[@]}"
)

if [[ "${DRY_RUN}" -eq 1 ]]; then
    STATUS="dry_run"; write_manifest
    echo
    echo "[dry-run] ${DOCKER_CMD[*]}"
    exit 0
fi

# ---- run --------------------------------------------------------------------
START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; STATUS="running"; write_manifest
echo "Streaming container output to ${RUN_LOG} ..."
"${DOCKER_CMD[@]}" 2>&1 | tee "${RUN_LOG}"
EXIT_CODE="${PIPESTATUS[0]}"
FINISH="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# ---- translate exit code ----------------------------------------------------
# run_all.sh / run_power_curve_fitting.sh exit 2 == "ran, some stages failed",
# distinct from a hard crash (OOM 137, signal 143, other non-zero).
case "${EXIT_CODE}" in
    0)   STATUS="success";              SCRIPT_EXIT=0 ;;
    2)   STATUS="pipeline_partial_fail"; SCRIPT_EXIT=1 ;;
    137) STATUS="container_oom_killed"; SCRIPT_EXIT=6 ;;
    143) STATUS="container_terminated"; SCRIPT_EXIT=6 ;;
    *)   STATUS="container_exit_${EXIT_CODE}"; SCRIPT_EXIT=6 ;;
esac
write_manifest

trap - INT TERM
echo
echo "Final status: ${STATUS} (container exit ${EXIT_CODE})"
echo "Manifest:     ${MANIFEST}"
exit "${SCRIPT_EXIT}"
