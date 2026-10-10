#!/usr/bin/env bash
# =============================================================================
# docker/entrypoint.sh -- container entrypoint for `docker run`.
#
# The reproduction pipeline does not use this: pipeline/run.sh runs the
# Apptainer image with `apptainer exec`, which bypasses the entrypoint. This
# script only serves `docker run` (docker/docker-smoke.sh and interactive use).
#
# tini is PID 1 (Dockerfile ENTRYPOINT). The container starts as root; if
# TARGET_UID / TARGET_GID are set, we add a passwd entry and drop to that user
# with gosu, so files written to a bind-mounted directory belong to the caller.
#
# Commands:
#   (none) | smoke   run docker/smoke_test.R (the default)
#   shell            interactive bash
#   anything else    executed as given, e.g. `docker run IMAGE Rscript -e ...`
# =============================================================================
set -euo pipefail

log()  { printf '[entrypoint] %s\n' "$*"; }
warn() { printf '[entrypoint] WARN: %s\n' "$*" >&2; }

if [[ "$(id -u)" -eq 0 && -n "${TARGET_UID:-}" && -n "${TARGET_GID:-}" ]]; then
    getent group "${TARGET_GID}" >/dev/null 2>&1 || echo "epa:x:${TARGET_GID}:" >> /etc/group
    getent passwd "${TARGET_UID}" >/dev/null 2>&1 \
        || echo "epa:x:${TARGET_UID}:${TARGET_GID}::/tmp:/bin/bash" >> /etc/passwd
    export HOME=/tmp
    log "dropping privileges to ${TARGET_UID}:${TARGET_GID}"
    exec gosu "${TARGET_UID}:${TARGET_GID}" "$0" "$@"
elif [[ "$(id -u)" -eq 0 ]]; then
    warn "running as root (TARGET_UID/TARGET_GID not set); files written to bind mounts will be root-owned"
fi

export HOME="${HOME:-/tmp}"
mkdir -p "${HOME}"
if [[ -d "${POWER_ANALYSIS_ROOT:-/workspace}" ]]; then
    cd "${POWER_ANALYSIS_ROOT:-/workspace}"
fi

case "${1:-smoke}" in
    smoke) exec Rscript /usr/local/share/epa/smoke_test.R ;;
    shell) shift || true; exec /bin/bash "$@" ;;
    *)     exec "$@" ;;
esac
