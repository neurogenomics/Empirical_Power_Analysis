#!/usr/bin/env bash
# =============================================================================
# docker-overnight.sh
#
# Single-entry, unattended overnight orchestrator: builds the Docker image,
# smoke-tests it, then launches the chosen pipeline inside the container with
# bind-mounted code and Datalink. Designed to be left running for days.
#
# Pipeline phases:
#   1. Build       — DOCKER_BUILDKIT=1 docker build (uses cache mounts inside
#                    the Dockerfile; rebuilds in <5 min once warm).
#   2. Smoke test  — R session + critical package load. Bails before the long
#                    run if anything in the image is broken.
#   3. Pipeline    — docker run -d, --restart=on-failure:3, bind-mounts the
#                    repo into /workspace and Datalink into /workspace/Datalink,
#                    streams `docker logs -f` to a host file so you can tail it.
#
# Usage:
#   nohup ./docker-overnight.sh > nohup.docker-overnight.out 2>&1 &
#   disown
#
# Flags (all optional):
#   --pipeline pcf|all              Default: pcf (power-curve-fitting only)
#   --max-cores N                   Default: 24
#   --no-overwrite | --force        Default: --no-overwrite
#   --mode cells|reads|both         Default: both (cartesian sweep)
#   --image NAME:TAG                Default: empirical-power-analysis:overnight
#   --datalink PATH                 Default: ${SCRIPT_DIR}/Datalink
#   --extra-arg "STRING"            Pass extra args to the wrapper inside the
#                                   container (repeatable; see examples).
#   --skip-build                    Use existing image, don't rebuild.
#   --build-only                    Build + smoke-test, do NOT run the pipeline.
#   --no-restart                    Don't apply --restart=on-failure (useful for
#                                   debugging so a crash doesn't auto-relaunch).
#   --dry-run                       Print what would happen, run nothing.
#   -h | --help
#
# Example with extra args (e.g. --skip on the wrapper):
#   ./docker-overnight.sh --extra-arg "--skip" --extra-arg "full_ground_truth"
#
# Logs land in:
#   logs/docker-overnight/<UTC-timestamp>/
#       01_build.log        full docker build output
#       02_smoke.log        smoke-test stdout/stderr
#       03_run.log          live `docker logs -f` of the running container
#       manifest.json       image tag/digest, host, git SHA, args, start/end
#       container.id        the container ID for manual `docker stop` if needed
# =============================================================================

set -u
set -o pipefail

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PIPELINE="pcf"
MAX_CORES=24
OVERWRITE_FLAG="--no-overwrite"
MODE="both"
IMAGE_TAG="empirical-power-analysis:overnight"
DATALINK_PATH="${SCRIPT_DIR}/Datalink"
SKIP_BUILD=0
BUILD_ONLY=0
NO_RESTART=0
DRY_RUN=0
EXTRA_ARGS=()

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
print_help() { sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pipeline)        PIPELINE="$2"; shift 2 ;;
        --max-cores)       MAX_CORES="$2"; shift 2 ;;
        --no-overwrite)    OVERWRITE_FLAG="--no-overwrite"; shift ;;
        --force)           OVERWRITE_FLAG="--force"; shift ;;
        --mode)            MODE="$2"; shift 2 ;;
        --image)           IMAGE_TAG="$2"; shift 2 ;;
        --datalink)        DATALINK_PATH="$2"; shift 2 ;;
        --extra-arg)       EXTRA_ARGS+=("$2"); shift 2 ;;
        --skip-build)      SKIP_BUILD=1; shift ;;
        --build-only)      BUILD_ONLY=1; shift ;;
        --no-restart)      NO_RESTART=1; shift ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         print_help ;;
        *)
            echo "[docker-overnight] Unknown flag: $1" >&2
            echo "Run with --help for usage." >&2
            exit 64
            ;;
    esac
done

case "${PIPELINE}" in
    pcf|all) ;;
    *) echo "[docker-overnight] --pipeline must be 'pcf' or 'all', got '${PIPELINE}'" >&2; exit 64 ;;
esac
case "${MODE}" in
    cells|reads|both) ;;
    *) echo "[docker-overnight] --mode must be cells|reads|both, got '${MODE}'" >&2; exit 64 ;;
esac

# -----------------------------------------------------------------------------
# Logging dirs
# -----------------------------------------------------------------------------
TS="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_ROOT="${SCRIPT_DIR}/logs/docker-overnight/${TS}"
mkdir -p "${LOG_ROOT}"
LATEST_LINK="${SCRIPT_DIR}/logs/docker-overnight/latest"
ln -sfn "${TS}" "${LATEST_LINK}"

BUILD_LOG="${LOG_ROOT}/01_build.log"
SMOKE_LOG="${LOG_ROOT}/02_smoke.log"
RUN_LOG="${LOG_ROOT}/03_run.log"
MANIFEST="${LOG_ROOT}/manifest.json"
CID_FILE="${LOG_ROOT}/container.id"

# -----------------------------------------------------------------------------
# Banner
# -----------------------------------------------------------------------------
banner() { printf '\n=== %s ===\n' "$1" | tee -a "${LOG_ROOT}/00_orchestrator.log"; }
log()    { printf '[%s] %s\n' "$(date -u +%H:%M:%SZ)" "$1" | tee -a "${LOG_ROOT}/00_orchestrator.log"; }

GIT_SHA="$(cd "${SCRIPT_DIR}" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOSTNAME_LOCAL="$(hostname)"

banner "docker-overnight.sh — overnight Docker orchestrator"
log "timestamp:     ${TS}"
log "host:          ${HOSTNAME_LOCAL}"
log "git SHA:       ${GIT_SHA}"
log "pipeline:      ${PIPELINE}"
log "max-cores:     ${MAX_CORES}"
log "overwrite:     ${OVERWRITE_FLAG}"
log "mode:          ${MODE}"
log "image tag:     ${IMAGE_TAG}"
log "datalink:      ${DATALINK_PATH}"
log "extra args:    ${EXTRA_ARGS[*]:-<none>}"
log "log dir:       ${LOG_ROOT}"
log "skip-build:    ${SKIP_BUILD}"
log "build-only:    ${BUILD_ONLY}"
log "no-restart:    ${NO_RESTART}"
log "dry-run:       ${DRY_RUN}"

if [[ ! -d "${DATALINK_PATH}" ]]; then
    log "WARN: Datalink path '${DATALINK_PATH}' does not exist. Creating an empty one."
    mkdir -p "${DATALINK_PATH}"
fi

# -----------------------------------------------------------------------------
# Manifest stub (filled in as we go)
# -----------------------------------------------------------------------------
write_manifest() {
    cat > "${MANIFEST}" <<EOF
{
  "timestamp_utc":   "${TS}",
  "host":            "${HOSTNAME_LOCAL}",
  "git_sha":         "${GIT_SHA}",
  "pipeline":        "${PIPELINE}",
  "image_tag":       "${IMAGE_TAG}",
  "image_digest":    "${IMAGE_DIGEST:-unknown}",
  "datalink_path":   "${DATALINK_PATH}",
  "max_cores":       ${MAX_CORES},
  "overwrite_flag":  "${OVERWRITE_FLAG}",
  "mode":            "${MODE}",
  "extra_args":      $(printf '%s\n' "${EXTRA_ARGS[@]:-}" | jq -R . 2>/dev/null | jq -s . 2>/dev/null || echo '[]'),
  "container_id":    "${CONTAINER_ID:-pending}",
  "phase":           "${PHASE:-init}",
  "build_started":   "${BUILD_START:-}",
  "build_finished":  "${BUILD_END:-}",
  "smoke_started":   "${SMOKE_START:-}",
  "smoke_finished":  "${SMOKE_END:-}",
  "run_started":     "${RUN_START:-}",
  "run_finished":    "${RUN_END:-}",
  "exit_status":     "${EXIT_STATUS:-pending}"
}
EOF
}
PHASE="init"; write_manifest

# -----------------------------------------------------------------------------
# Pre-flight: docker available?
# -----------------------------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
    log "ERROR: docker CLI not found. Install Docker or run on a host that has it."
    EXIT_STATUS="docker_missing"; PHASE="failed"; write_manifest
    exit 2
fi
DOCKER_VERSION="$(docker --version 2>&1 || echo unknown)"
log "docker:        ${DOCKER_VERSION}"

# -----------------------------------------------------------------------------
# Phase 1 — Build
# -----------------------------------------------------------------------------
if [[ "${SKIP_BUILD}" -eq 0 ]]; then
    banner "Phase 1 — docker build"
    BUILD_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    PHASE="building"; write_manifest

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log "[dry-run] DOCKER_BUILDKIT=1 docker build -t ${IMAGE_TAG} ${SCRIPT_DIR}"
    else
        log "Build log: ${BUILD_LOG}"
        if ! DOCKER_BUILDKIT=1 docker build \
                --progress=plain \
                -t "${IMAGE_TAG}" \
                "${SCRIPT_DIR}" \
                > "${BUILD_LOG}" 2>&1; then
            log "ERROR: docker build failed. Tail of build log:"
            tail -n 40 "${BUILD_LOG}" | sed 's/^/    /'
            EXIT_STATUS="build_failed"; PHASE="failed"; write_manifest
            exit 3
        fi
    fi
    BUILD_END="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    log "Build complete."
else
    log "Skipping build (--skip-build); using existing image ${IMAGE_TAG}"
fi

# Capture image digest (id) for manifest
if [[ "${DRY_RUN}" -eq 0 ]]; then
    IMAGE_DIGEST="$(docker image inspect --format='{{.Id}}' "${IMAGE_TAG}" 2>/dev/null || echo unknown)"
else
    IMAGE_DIGEST="dry-run"
fi
write_manifest

# -----------------------------------------------------------------------------
# Phase 2 — Smoke test
# -----------------------------------------------------------------------------
banner "Phase 2 — smoke test"
SMOKE_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
PHASE="smoke_testing"; write_manifest

SMOKE_R='
critical <- c("DESeq2","apeglm","glmGamPoi","SingleCellExperiment",
              "qs","reticulate","zellkonverter","future.apply",
              "ggplot2","dplyr","tidyr","Matrix")
miss <- character()
for (p in critical) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-20s %s %s\n", p,
              if (ok) "OK" else "FAIL",
              if (ok) as.character(packageVersion(p)) else ""))
  if (!ok) miss <- c(miss, p)
}
if (length(miss) > 0) stop("Missing critical packages: ", paste(miss, collapse=", "))
cat("\nSMOKE TEST PASSED\n")
'

if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "[dry-run] docker run --rm ${IMAGE_TAG} R -e '<smoke checks>'"
else
    log "Smoke log: ${SMOKE_LOG}"
    # --workdir / so renv activates from the project root (.Rprofile lives at /,
    # not at /workspace). Without this, R uses the system library only and
    # cannot see the renv-installed packages.
    if ! docker run --rm --workdir / --entrypoint R "${IMAGE_TAG}" \
            -e "${SMOKE_R}" > "${SMOKE_LOG}" 2>&1; then
        log "ERROR: smoke test failed. Tail of smoke log:"
        tail -n 40 "${SMOKE_LOG}" | sed 's/^/    /'
        EXIT_STATUS="smoke_failed"; PHASE="failed"; write_manifest
        exit 4
    fi
    log "Smoke test PASSED."
fi
SMOKE_END="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
write_manifest

if [[ "${BUILD_ONLY}" -eq 1 ]]; then
    log "--build-only set; not running pipeline."
    EXIT_STATUS="build_only_ok"; PHASE="build_only_complete"; write_manifest
    exit 0
fi

# -----------------------------------------------------------------------------
# Phase 3 — Pipeline run
# -----------------------------------------------------------------------------
banner "Phase 3 — pipeline run inside container"
RUN_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
PHASE="running"; write_manifest

CONTAINER_NAME="empirical-power-${TS}"

# Compose wrapper args (these are what the entrypoint forwards to
# run_power_curve_fitting.sh / run_all.sh).
WRAPPER_ARGS=( "--max-cores" "${MAX_CORES}" "${OVERWRITE_FLAG}" "--mode" "${MODE}" )
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    WRAPPER_ARGS+=( "${EXTRA_ARGS[@]}" )
fi

RESTART_FLAG="--restart=on-failure:3"
if [[ "${NO_RESTART}" -eq 1 ]]; then
    RESTART_FLAG="--restart=no"
fi

# Run as the host user so files written into the bind-mounted Datalink end up
# owned by the host user, not root. This matters for syncthing/git on the host.
HOST_UID="$(id -u)"
HOST_GID="$(id -g)"

DOCKER_RUN_CMD=(
    docker run -d
    --name "${CONTAINER_NAME}"
    "${RESTART_FLAG}"
    -e "PIPELINE=${PIPELINE}"
    -e "POWER_ANALYSIS_ROOT=/workspace"
    # Host UID has no /etc/passwd entry inside the image, so glibc falls back
    # HOME=/ which is read-only. renv tries to mkdir //.cache/R/renv and dies.
    # Pin HOME to a writable tmpfs path so renv has a cache dir.
    -e "HOME=/tmp"
    -u "${HOST_UID}:${HOST_GID}"
    -v "${SCRIPT_DIR}:/workspace"
    -v "${DATALINK_PATH}:/workspace/Datalink"
    --entrypoint "/workspace/docker/entrypoint.sh"
    "${IMAGE_TAG}"
    "${WRAPPER_ARGS[@]}"
)

log "Container: ${CONTAINER_NAME}"
log "Cmd: ${DOCKER_RUN_CMD[*]}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "[dry-run] (would launch container; would tail logs to ${RUN_LOG})"
    EXIT_STATUS="dry_run_ok"; PHASE="dry_run_complete"; write_manifest
    exit 0
fi

if ! CONTAINER_ID="$("${DOCKER_RUN_CMD[@]}" 2>>"${LOG_ROOT}/00_orchestrator.log")"; then
    log "ERROR: docker run failed."
    EXIT_STATUS="run_launch_failed"; PHASE="failed"; write_manifest
    exit 5
fi
echo "${CONTAINER_ID}" > "${CID_FILE}"
log "Container ID: ${CONTAINER_ID}"
write_manifest

# Stream container logs to file in foreground. This call blocks until the
# container exits (or the docker daemon goes away). Because we used `docker
# run -d`, the container is detached from this script — if you SIGTERM this
# orchestrator, the container keeps running. Use `docker stop ${CONTAINER_NAME}`
# to actually stop the work.
log "Streaming container logs to ${RUN_LOG}"
log "(this orchestrator now blocks until the container exits — safe to nohup)"
docker logs -f "${CONTAINER_ID}" >"${RUN_LOG}" 2>&1 || true

# Wait for the container to exit and grab its status.
EXIT_CODE="$(docker inspect --format='{{.State.ExitCode}}' "${CONTAINER_ID}" 2>/dev/null || echo unknown)"
RUN_END="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
log "Container exited with code ${EXIT_CODE}"

if [[ "${EXIT_CODE}" == "0" ]]; then
    EXIT_STATUS="success"
    PHASE="complete"
else
    EXIT_STATUS="container_exit_${EXIT_CODE}"
    PHASE="failed"
fi
write_manifest

banner "docker-overnight.sh complete"
log "Final status: ${EXIT_STATUS}"
log "Manifest: ${MANIFEST}"
log "Run log:  ${RUN_LOG}"

exit $([[ "${EXIT_STATUS}" == "success" ]] && echo 0 || echo 1)
