#!/usr/bin/env bash
# =============================================================================
# docker/entrypoint.sh
#
# Container entrypoint. The container starts as root (required so we can
# fix up /etc/passwd / chown /tmp / write the heartbeat file regardless of
# the host UID), then drops to the TARGET_UID:TARGET_GID requested by the
# caller (docker-run.sh / docker-smoke.sh) via `gosu` for the actual run.
#
# Privilege model:
#   * tini is PID 1 (set in the Dockerfile ENTRYPOINT). It forwards SIGTERM
#     and reaps zombies.
#   * This script runs as root.
#   * If TARGET_UID / TARGET_GID are set, we add a /etc/passwd entry for
#     them and `exec gosu` to drop to that user before starting the pipeline.
#   * If they are not set (e.g. the container was started with `docker run -u`
#     directly, bypassing this contract), we proceed as the existing user
#     and emit a warning.
#
# Pipelines (PIPELINE env var, default = pcf):
#   pcf      -> ./run_power_curve_fitting.sh "$@"
#   all      -> ./run_all.sh                "$@"
#   shell    -> drop into bash (debug)
#   smoke    -> embedded R smoke test (no pipeline)
# =============================================================================

set -euo pipefail

log()  { printf '[entrypoint %(%H:%M:%SZ)T] %s\n' -1 "$*"; }
warn() { printf '[entrypoint %(%H:%M:%SZ)T] WARN: %s\n' -1 "$*" >&2; }
fail() { printf '[entrypoint %(%H:%M:%SZ)T] FATAL: %s\n' -1 "$*" >&2; exit "${2:-1}"; }

# -----------------------------------------------------------------------------
# Step 1 — Privilege handoff. If we're root and TARGET_UID/GID are requested,
# create a /etc/passwd entry for them, then re-exec ourselves as that user via
# gosu. The second invocation goes straight to step 2 below.
# -----------------------------------------------------------------------------
if [[ "$(id -u)" -eq 0 ]]; then
    : "${TARGET_UID:=}"
    : "${TARGET_GID:=}"
    if [[ -n "${TARGET_UID}" && -n "${TARGET_GID}" ]]; then
        if ! getent passwd "${TARGET_UID}" >/dev/null 2>&1; then
            log "Adding /etc/passwd entry for UID=${TARGET_UID} GID=${TARGET_GID}"
            echo "epa:x:${TARGET_UID}:${TARGET_GID}:Empirical Power Analysis runtime user:${HOME:-/tmp}:/bin/bash" >> /etc/passwd
        fi
        if ! getent group "${TARGET_GID}" >/dev/null 2>&1; then
            echo "epa:x:${TARGET_GID}:" >> /etc/group
        fi

        # Make HOME writable by the target user. /tmp is the default and is
        # already 1777, but be defensive in case HOME has been overridden.
        : "${HOME:=/tmp}"
        if [[ ! -w "${HOME}" ]]; then
            mkdir -p "${HOME}"
            chown "${TARGET_UID}:${TARGET_GID}" "${HOME}" || true
        fi

        # Carry the existing env across the privilege drop. gosu does this by
        # default, but we exec under env to be explicit about which vars matter.
        log "Dropping privileges via gosu to ${TARGET_UID}:${TARGET_GID}"
        exec gosu "${TARGET_UID}:${TARGET_GID}" "$0" "$@"
    else
        warn "Running as root (TARGET_UID/TARGET_GID not set). Files written"
        warn "to bind-mounted Datalink will be owned by root on the host."
    fi
fi

# -----------------------------------------------------------------------------
# Step 2 — Sanity checks (we're now the runtime user).
# -----------------------------------------------------------------------------
HOME="${HOME:-/tmp}"
if [[ ! -d "${HOME}" ]]; then
    mkdir -p "${HOME}"
fi
if [[ ! -w "${HOME}" ]]; then
    fail "HOME='${HOME}' is not writable. Set HOME to a writable path via 'docker run -e HOME=...'." 65
fi
mkdir -p "${HOME}/.cache/R/renv" "${HOME}/.local/share"

WORKSPACE="${POWER_ANALYSIS_ROOT:-/workspace}"
if [[ ! -d "${WORKSPACE}/R" ]]; then
    fail "POWER_ANALYSIS_ROOT='${WORKSPACE}' is missing the R/ directory.
Did you forget to bind-mount the repo? Expected:
    docker run -v <repo>:${WORKSPACE} ..." 66
fi
cd "${WORKSPACE}"

if [[ ! -d "${RENV_PATHS_LIBRARY:-/opt/renv-library}" ]]; then
    fail "renv library missing at ${RENV_PATHS_LIBRARY:-/opt/renv-library}.
This image is broken — rebuild from a clean state." 67
fi

# Disable renv's autoloader inside the container. The bind-mounted host repo
# carries a .Rprofile that sources renv/activate.R; if it ran, renv would
# isolate .libPaths() to the (host-platform) project library and shadow the
# image's pinned library. We instead rely on the R_LIBS_SITE baked into the
# image's Renviron.site (see docker/install_packages.R), which makes the pinned
# packages visible deterministically without renv re-activating at runtime.
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

# -----------------------------------------------------------------------------
# Step 3 — Heartbeat ticker for HEALTHCHECK.
# -----------------------------------------------------------------------------
HEARTBEAT_FILE="${CBERT_HEARTBEAT_FILE:-/tmp/empirical-power-analysis.heartbeat}"
touch "${HEARTBEAT_FILE}"

heartbeat_loop() {
    while true; do
        sleep 60
        touch "${HEARTBEAT_FILE}" 2>/dev/null || exit 0
    done
}
heartbeat_loop &
HEARTBEAT_PID=$!

# -----------------------------------------------------------------------------
# Step 4 — Signal forwarding.
# -----------------------------------------------------------------------------
CHILD_PID=""
shutdown() {
    local sig="$1"
    log "Received ${sig}; shutting down children."
    if [[ -n "${CHILD_PID}" ]] && kill -0 "${CHILD_PID}" 2>/dev/null; then
        kill -s "${sig}" "${CHILD_PID}" 2>/dev/null || true
    fi
    if kill -0 "${HEARTBEAT_PID}" 2>/dev/null; then
        kill "${HEARTBEAT_PID}" 2>/dev/null || true
    fi
}
trap 'shutdown TERM' TERM
trap 'shutdown INT'  INT

# -----------------------------------------------------------------------------
# Step 5 — Banner.
# -----------------------------------------------------------------------------
log "=== empirical-power-analysis container starting ==="
log "user:      $(id 2>/dev/null || echo unknown)"
log "host:      $(hostname)"
log "workspace: ${WORKSPACE}"
log "renv lib:  ${RENV_PATHS_LIBRARY:-/opt/renv-library}"
log "py env:    ${RETICULATE_PYTHON:-unset}"
log "pipeline:  ${PIPELINE:-pcf}"
log "args:      $*"
log "heartbeat: ${HEARTBEAT_FILE}"

# -----------------------------------------------------------------------------
# Step 6 — Dispatch.
# -----------------------------------------------------------------------------
PIPELINE="${PIPELINE:-pcf}"
case "${PIPELINE}" in
    pcf)
        if [[ ! -x ./run_power_curve_fitting.sh ]]; then
            fail "./run_power_curve_fitting.sh missing or not executable in ${WORKSPACE}." 68
        fi
        ./run_power_curve_fitting.sh "$@" &
        ;;
    all)
        if [[ ! -x ./run_all.sh ]]; then
            fail "./run_all.sh missing or not executable in ${WORKSPACE}." 68
        fi
        ./run_all.sh "$@" &
        ;;
    shell)
        # exec replaces this script; heartbeat ticker becomes orphan but is
        # killed when the container exits. No `wait` afterwards.
        kill "${HEARTBEAT_PID}" 2>/dev/null || true
        exec /bin/bash "$@"
        ;;
    smoke)
        Rscript -e '
critical <- c("DESeq2","apeglm","glmGamPoi","SingleCellExperiment",
              "qs","reticulate","zellkonverter","future.apply",
              "ggplot2","dplyr","tidyr","Matrix","Rcpp")
bad <- character()
for (p in critical) {
  ok <- tryCatch({ suppressMessages(library(p, character.only = TRUE)); TRUE },
                 error = function(e) { cat(sprintf("FAIL %s: %s\n", p, conditionMessage(e))); FALSE })
  if (ok) cat(sprintf("OK   %-22s %s\n", p, packageVersion(p)))
  else bad <- c(bad, p)
}
cat("\nlibPaths:\n"); for (l in .libPaths()) cat(" ", l, "\n")
cat("\nReticulate Python:\n")
tryCatch({
  reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
  reticulate::py_run_string("import anndata, sys; print(f\"  anndata {anndata.__version__} on {sys.version.split()[0]}\")")
}, error = function(e) cat("  RETICULATE FAILED:", conditionMessage(e), "\n"))
if (length(bad) > 0) stop("Smoke-test failed: ", paste(bad, collapse=", "))
cat("\nSMOKE TEST PASSED\n")
' &
        ;;
    *)
        fail "Unknown PIPELINE='${PIPELINE}'. Expected: pcf | all | shell | smoke" 64
        ;;
esac

CHILD_PID=$!
set +e
wait "${CHILD_PID}"
CHILD_EXIT=$?
set -e
kill "${HEARTBEAT_PID}" 2>/dev/null || true

log "Pipeline child (${PIPELINE}) exited with code ${CHILD_EXIT}"
exit "${CHILD_EXIT}"
