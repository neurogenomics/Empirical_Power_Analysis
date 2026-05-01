#!/usr/bin/env bash
# =============================================================================
# docker-overnight.sh
#
# Single-entry, unattended overnight orchestrator.
#
# Phases:
#   1. Pre-flight  — Docker present? BuildKit? disk? cores?
#   2. Build       — DOCKER_BUILDKIT=1 docker build (cache mounts inside the
#                    Dockerfile keep rebuilds <5 min once warm).
#   3. Smoke test  — runs the SAME image, with the SAME bind-mount config and
#                    user as the real run, with PIPELINE=smoke. If anything is
#                    broken with the renv library / Python env / passwd fixup,
#                    we fail here, before the multi-day pipeline starts.
#   4. Pipeline    — docker run -d, bind-mounts host repo on /workspace and
#                    Datalink on /workspace/Datalink, streams logs via
#                    `docker logs -f --timestamps` to a host file.
#
# Usage:
#   nohup ./docker-overnight.sh > nohup.docker-overnight.out 2>&1 &
#   disown
#
# Flags (all optional):
#   --pipeline pcf|all         Default: pcf
#   --max-cores N              Default: min(24, nproc). Capped to nproc.
#   --no-overwrite | --force   Default: --no-overwrite
#   --mode cells|reads|both    Default: both
#   --image NAME:TAG           Default: empirical-power-analysis:overnight
#   --datalink PATH            Default: ${SCRIPT_DIR}/Datalink
#   --extra-arg "STRING"       Pass extra args to the wrapper (repeatable)
#   --skip-build               Use existing image, don't rebuild
#   --build-only               Build + smoke-test, do NOT run the pipeline
#   --auto-restart             Add --restart=on-failure:3 (off by default —
#                              auto-restart masks errors and was historically
#                              the source of confusing repeated runs)
#   --min-free-gb N            Refuse to start if Datalink filesystem has
#                              less than N GiB free. Default 50.
#   --dry-run                  Print what would happen, run nothing
#   -h | --help
#
# Logs land in:
#   logs/docker-overnight/<UTC-timestamp>/
#       00_orchestrator.log    timestamped orchestrator events
#       01_build.log           docker build stdout/stderr
#       02_smoke.log           smoke-test container stdout/stderr
#       03_run.log             docker logs --timestamps -f ${container}
#       manifest.json          full provenance: image digest, git SHA,
#                              host, args, phase timings, exit status
#       container.id           container ID for manual `docker stop`
#       orchestrator.pid       this script's PID (kill the orchestrator
#                              without killing the container)
#
# Exit codes:
#    0  success
#    1  pipeline reported failures (run_all.sh exit 2 — partial failure)
#    2  pre-flight check failed (no docker, low disk, bad args)
#    3  build failed
#    4  smoke test failed
#    5  docker run could not start the container
#    6  container died with a non-zero/non-2 exit code (real crash)
#  130  orchestrator received SIGINT
#  143  orchestrator received SIGTERM
# =============================================================================

set -u
set -o pipefail

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PIPELINE="pcf"
MAX_CORES="auto"          # resolved to min(24, nproc) below
OVERWRITE_FLAG="--no-overwrite"
MODE="both"
IMAGE_TAG="empirical-power-analysis:overnight"
DATALINK_PATH="${SCRIPT_DIR}/Datalink"
SKIP_BUILD=0
BUILD_ONLY=0
AUTO_RESTART=0
MIN_FREE_GB=50
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
        --auto-restart)    AUTO_RESTART=1; shift ;;
        --min-free-gb)     MIN_FREE_GB="$2"; shift 2 ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         print_help ;;
        *)
            echo "[docker-overnight] Unknown flag: $1" >&2
            echo "Run with --help for usage." >&2
            exit 2
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Validation of arg values (keep all checks here so the rest of the script is
# free to assume the values are sane).
# -----------------------------------------------------------------------------
case "${PIPELINE}" in
    pcf|all) ;;
    *) echo "[docker-overnight] --pipeline must be 'pcf' or 'all', got '${PIPELINE}'" >&2; exit 2 ;;
esac
case "${MODE}" in
    cells|reads|both) ;;
    *) echo "[docker-overnight] --mode must be cells|reads|both, got '${MODE}'" >&2; exit 2 ;;
esac
case "${OVERWRITE_FLAG}" in
    --no-overwrite|--force) ;;
    *) echo "[docker-overnight] internal: invalid OVERWRITE_FLAG ${OVERWRITE_FLAG}" >&2; exit 2 ;;
esac

NPROC="$(nproc 2>/dev/null || echo 1)"
if [[ "${MAX_CORES}" == "auto" ]]; then
    MAX_CORES=$(( NPROC < 24 ? NPROC : 24 ))
fi
if ! [[ "${MAX_CORES}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[docker-overnight] --max-cores must be a positive integer, got '${MAX_CORES}'" >&2
    exit 2
fi
if (( MAX_CORES > NPROC )); then
    echo "[docker-overnight] WARN: --max-cores ${MAX_CORES} > host nproc ${NPROC}; clamping to ${NPROC}." >&2
    MAX_CORES="${NPROC}"
fi

if ! [[ "${MIN_FREE_GB}" =~ ^[0-9]+$ ]]; then
    echo "[docker-overnight] --min-free-gb must be a non-negative integer, got '${MIN_FREE_GB}'" >&2
    exit 2
fi

# -----------------------------------------------------------------------------
# Logging dirs
# -----------------------------------------------------------------------------
TS="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_ROOT="${SCRIPT_DIR}/logs/docker-overnight/${TS}"
mkdir -p "${LOG_ROOT}"
LATEST_LINK="${SCRIPT_DIR}/logs/docker-overnight/latest"
ln -sfn "${TS}" "${LATEST_LINK}"

ORCH_LOG="${LOG_ROOT}/00_orchestrator.log"
BUILD_LOG="${LOG_ROOT}/01_build.log"
SMOKE_LOG="${LOG_ROOT}/02_smoke.log"
RUN_LOG="${LOG_ROOT}/03_run.log"
MANIFEST="${LOG_ROOT}/manifest.json"
CID_FILE="${LOG_ROOT}/container.id"
PID_FILE="${LOG_ROOT}/orchestrator.pid"
echo $$ > "${PID_FILE}"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
banner() { printf '\n=== %s ===\n' "$1" | tee -a "${ORCH_LOG}"; }
log()    { printf '[%s] %s\n' "$(date -u +%H:%M:%SZ)" "$1" | tee -a "${ORCH_LOG}"; }
err()    { printf '[%s] ERROR: %s\n' "$(date -u +%H:%M:%SZ)" "$1" | tee -a "${ORCH_LOG}" >&2; }

GIT_SHA="$(cd "${SCRIPT_DIR}" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOSTNAME_LOCAL="$(hostname)"

CONTAINER_ID=""
PHASE="init"
BUILD_START=""; BUILD_END=""
SMOKE_START=""; SMOKE_END=""
RUN_START=""; RUN_END=""
EXIT_STATUS="pending"
EXIT_CODE_RAW=""
IMAGE_DIGEST=""

# -----------------------------------------------------------------------------
# Manifest writer. Idempotent — call as often as needed; always reflects the
# current state of the orchestrator. jq is optional (used only to format the
# extra_args array); we fall back to a string-join if jq is missing on host.
# -----------------------------------------------------------------------------
write_manifest() {
    local extra_json
    if command -v jq >/dev/null 2>&1 && (( ${#EXTRA_ARGS[@]} > 0 )); then
        extra_json="$(printf '%s\n' "${EXTRA_ARGS[@]}" | jq -R . | jq -s -c .)"
    elif (( ${#EXTRA_ARGS[@]} > 0 )); then
        # Minimal manual JSON encoding; we control these values so this is safe.
        local joined=""
        for a in "${EXTRA_ARGS[@]}"; do
            joined+="\"${a//\"/\\\"}\","
        done
        extra_json="[${joined%,}]"
    else
        extra_json="[]"
    fi
    cat > "${MANIFEST}" <<EOF
{
  "schema_version":  2,
  "timestamp_utc":   "${TS}",
  "host":            "${HOSTNAME_LOCAL}",
  "git_sha":         "${GIT_SHA}",
  "pipeline":        "${PIPELINE}",
  "image_tag":       "${IMAGE_TAG}",
  "image_digest":    "${IMAGE_DIGEST:-unknown}",
  "datalink_path":   "${DATALINK_PATH}",
  "max_cores":       ${MAX_CORES},
  "host_nproc":      ${NPROC},
  "overwrite_flag":  "${OVERWRITE_FLAG}",
  "mode":            "${MODE}",
  "auto_restart":    ${AUTO_RESTART},
  "min_free_gb":     ${MIN_FREE_GB},
  "extra_args":      ${extra_json},
  "container_id":    "${CONTAINER_ID:-pending}",
  "phase":           "${PHASE}",
  "build_started":   "${BUILD_START}",
  "build_finished":  "${BUILD_END}",
  "smoke_started":   "${SMOKE_START}",
  "smoke_finished":  "${SMOKE_END}",
  "run_started":     "${RUN_START}",
  "run_finished":    "${RUN_END}",
  "exit_status":     "${EXIT_STATUS}",
  "exit_code":       "${EXIT_CODE_RAW:-}"
}
EOF
}

# -----------------------------------------------------------------------------
# Signal trap: forward to the running container so `kill orchestrator` actually
# stops the work, instead of leaving an orphan container running.
# -----------------------------------------------------------------------------
shutdown() {
    local sig="$1"
    err "Orchestrator received ${sig} — stopping container if running."
    if [[ -n "${CONTAINER_ID}" ]]; then
        docker stop --timeout 30 "${CONTAINER_ID}" 2>>"${ORCH_LOG}" || true
    fi
    EXIT_STATUS="orchestrator_${sig}"
    PHASE="aborted"
    write_manifest
    case "${sig}" in
        INT)  exit 130 ;;
        TERM) exit 143 ;;
        *)    exit 1   ;;
    esac
}
trap 'shutdown INT'  INT
trap 'shutdown TERM' TERM

# -----------------------------------------------------------------------------
# Banner / preamble
# -----------------------------------------------------------------------------
banner "docker-overnight.sh — overnight Docker orchestrator"
log "timestamp:     ${TS}"
log "host:          ${HOSTNAME_LOCAL} (${NPROC} cores)"
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
log "auto-restart:  ${AUTO_RESTART}"
log "min-free-gb:   ${MIN_FREE_GB}"
log "dry-run:       ${DRY_RUN}"
log "PID:           $$  (saved to ${PID_FILE})"

write_manifest

# -----------------------------------------------------------------------------
# Phase 0 — Pre-flight checks
# -----------------------------------------------------------------------------
banner "Phase 0 — pre-flight"

if ! command -v docker >/dev/null 2>&1; then
    err "docker CLI not found. Install Docker or run on a host that has it."
    EXIT_STATUS="docker_missing"; PHASE="failed"; write_manifest
    exit 2
fi
DOCKER_VERSION="$(docker --version 2>&1 || echo unknown)"
log "docker:        ${DOCKER_VERSION}"

if ! docker info >/dev/null 2>&1; then
    err "docker daemon not reachable (permissions? socket missing?)."
    EXIT_STATUS="docker_daemon_unreachable"; PHASE="failed"; write_manifest
    exit 2
fi

# BuildKit is required by the Dockerfile (cache mounts, syntax directive).
# Detect early — easier to fix here than to fish through `docker build` errors.
if ! docker buildx version >/dev/null 2>&1; then
    log "WARN: docker buildx not present. We will set DOCKER_BUILDKIT=1 anyway,"
    log "WARN: but if the build fails with 'Unknown directive', install buildx."
fi

if [[ ! -d "${DATALINK_PATH}" ]]; then
    log "WARN: Datalink path '${DATALINK_PATH}' does not exist. Creating it."
    mkdir -p "${DATALINK_PATH}"
fi

# Disk-space pre-flight: refuse to launch if the Datalink filesystem is
# critically low. A multi-day cartesian sweep can write 100+ GiB; running out
# halfway leaves a poison-pill of partial outputs that confuses resume.
DF_AVAIL_KB="$(df -P "${DATALINK_PATH}" | awk 'NR==2 { print $4 }')"
DF_AVAIL_GB=$(( DF_AVAIL_KB / 1024 / 1024 ))
log "disk:          ${DF_AVAIL_GB} GiB free at ${DATALINK_PATH}"
if (( DF_AVAIL_GB < MIN_FREE_GB )); then
    err "Free space ${DF_AVAIL_GB} GiB < required ${MIN_FREE_GB} GiB. Aborting."
    err "Free up disk or pass --min-free-gb ${DF_AVAIL_GB} to proceed anyway."
    EXIT_STATUS="low_disk"; PHASE="failed"; write_manifest
    exit 2
fi

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
            err "docker build failed. Tail of build log:"
            tail -n 60 "${BUILD_LOG}" | sed 's/^/    /' | tee -a "${ORCH_LOG}"
            EXIT_STATUS="build_failed"; PHASE="failed"; write_manifest
            exit 3
        fi
    fi
    BUILD_END="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    log "Build complete."
else
    log "Skipping build (--skip-build); using existing image ${IMAGE_TAG}"
fi

# Capture image digest for the manifest.
if [[ "${DRY_RUN}" -eq 0 ]]; then
    IMAGE_DIGEST="$(docker image inspect --format='{{.Id}}' "${IMAGE_TAG}" 2>/dev/null || echo unknown)"
else
    IMAGE_DIGEST="dry-run"
fi
write_manifest

# -----------------------------------------------------------------------------
# Phase 2 — Smoke test
#
# CRITICAL: this MUST run with the same bind-mount + user config as the real
# pipeline run. The previous orchestrator's smoke test ran from / with no bind
# mount — which is a different environment from the real run, so it could pass
# while the real run was broken. This version dispatches PIPELINE=smoke
# through the real entrypoint with the real bind-mounts.
# -----------------------------------------------------------------------------
banner "Phase 2 — smoke test (real-runtime config)"
SMOKE_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
PHASE="smoke_testing"; write_manifest

HOST_UID="$(id -u)"
HOST_GID="$(id -g)"

SMOKE_RUN_CMD=(
    docker run --rm
    -e "PIPELINE=smoke"
    -e "POWER_ANALYSIS_ROOT=/workspace"
    -e "HOME=/tmp"
    # The container starts as root (entrypoint passwd-fixup needs it) and
    # drops to TARGET_UID:TARGET_GID via gosu before running anything.
    -e "TARGET_UID=${HOST_UID}"
    -e "TARGET_GID=${HOST_GID}"
    -v "${SCRIPT_DIR}:/workspace"
    -v "${DATALINK_PATH}:/workspace/Datalink"
    "${IMAGE_TAG}"
)

if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "[dry-run] ${SMOKE_RUN_CMD[*]}"
else
    log "Smoke log: ${SMOKE_LOG}"
    log "Cmd: ${SMOKE_RUN_CMD[*]}"
    if ! "${SMOKE_RUN_CMD[@]}" > "${SMOKE_LOG}" 2>&1; then
        err "Smoke test failed. Tail of smoke log:"
        tail -n 60 "${SMOKE_LOG}" | sed 's/^/    /' | tee -a "${ORCH_LOG}"
        EXIT_STATUS="smoke_failed"; PHASE="failed"; write_manifest
        exit 4
    fi
    if ! grep -q "SMOKE TEST PASSED" "${SMOKE_LOG}"; then
        err "Smoke test exited 0 but did not print SMOKE TEST PASSED. Suspicious."
        tail -n 60 "${SMOKE_LOG}" | sed 's/^/    /' | tee -a "${ORCH_LOG}"
        EXIT_STATUS="smoke_inconclusive"; PHASE="failed"; write_manifest
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

# Wrapper args forwarded to run_power_curve_fitting.sh / run_all.sh.
WRAPPER_ARGS=( "--max-cores" "${MAX_CORES}" "${OVERWRITE_FLAG}" "--mode" "${MODE}" )
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    WRAPPER_ARGS+=( "${EXTRA_ARGS[@]}" )
fi

# Auto-restart is OFF by default. The previous default of `--restart=on-failure:3`
# masked real failures (a stage failure exits 2; container restarts; the
# restart skips the stage via --no-overwrite; the container exits 0 — looks
# successful even though the stage was never repaired). Opt-in with
# --auto-restart if you really want it.
RESTART_FLAG="--restart=no"
if [[ "${AUTO_RESTART}" -eq 1 ]]; then
    RESTART_FLAG="--restart=on-failure:3"
    log "WARN: --auto-restart ON. Container will be restarted up to 3 times on non-zero exit."
fi

DOCKER_RUN_CMD=(
    docker run -d
    --name "${CONTAINER_NAME}"
    "${RESTART_FLAG}"
    -e "PIPELINE=${PIPELINE}"
    -e "POWER_ANALYSIS_ROOT=/workspace"
    -e "HOME=/tmp"
    # Container starts as root so the entrypoint can write to /etc/passwd, then
    # drops to TARGET_UID:TARGET_GID via gosu. Files written into the
    # bind-mounted Datalink end up owned by the host user, as before.
    -e "TARGET_UID=${HOST_UID}"
    -e "TARGET_GID=${HOST_GID}"
    -v "${SCRIPT_DIR}:/workspace"
    -v "${DATALINK_PATH}:/workspace/Datalink"
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

if ! CONTAINER_ID="$("${DOCKER_RUN_CMD[@]}" 2>>"${ORCH_LOG}")"; then
    err "docker run failed."
    EXIT_STATUS="run_launch_failed"; PHASE="failed"; write_manifest
    exit 5
fi
echo "${CONTAINER_ID}" > "${CID_FILE}"
log "Container ID: ${CONTAINER_ID}"
write_manifest

log "Streaming container logs (--timestamps) to ${RUN_LOG}"
log "(orchestrator now blocks until the container exits — safe to nohup)"
log "(SIGTERM/SIGINT to this orchestrator will stop the container)"

# `docker logs --timestamps -f` prefixes every line with an RFC-3339 timestamp,
# which is invaluable for tracing slow stages. The `|| true` ensures we don't
# bail if the docker daemon goes away mid-stream — we'll still inspect the
# container's exit code below.
docker logs --timestamps -f "${CONTAINER_ID}" >"${RUN_LOG}" 2>&1 || true

# Wait for the container to actually exit (logs -f may return early).
EXIT_CODE_RAW="$(docker inspect --format='{{.State.ExitCode}}' "${CONTAINER_ID}" 2>/dev/null || echo unknown)"
RUN_END="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
log "Container exited with code ${EXIT_CODE_RAW}"

# Translate the raw exit code into a meaningful status. Exit 2 from
# run_all.sh / run_power_curve_fitting.sh means "pipeline ran, some stages
# failed" — surface that distinctly from a hard crash (signal kill, OOM).
case "${EXIT_CODE_RAW}" in
    0)
        EXIT_STATUS="success";              PHASE="complete"; ORCH_EXIT=0 ;;
    2)
        EXIT_STATUS="pipeline_partial_fail"; PHASE="partial"; ORCH_EXIT=1 ;;
    137)
        EXIT_STATUS="container_oom_killed"; PHASE="failed";  ORCH_EXIT=6 ;;
    143)
        EXIT_STATUS="container_terminated"; PHASE="failed";  ORCH_EXIT=6 ;;
    *)
        EXIT_STATUS="container_exit_${EXIT_CODE_RAW}"; PHASE="failed"; ORCH_EXIT=6 ;;
esac
write_manifest

banner "docker-overnight.sh complete"
log "Final status: ${EXIT_STATUS}"
log "Manifest:     ${MANIFEST}"
log "Run log:      ${RUN_LOG}"
log "Container:    ${CONTAINER_NAME} (id ${CONTAINER_ID})"

# Clear traps so we don't double-handle on natural exit.
trap - INT TERM
exit "${ORCH_EXIT}"
