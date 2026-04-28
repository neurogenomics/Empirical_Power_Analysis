#!/usr/bin/env bash
# =============================================================================
# docker/entrypoint.sh
#
# Container entrypoint. Sets POWER_ANALYSIS_ROOT to /workspace (the bind-mounted
# repo) and dispatches to the requested pipeline wrapper. Anything passed as
# arguments goes straight through to the wrapper.
#
# Selection: the PIPELINE env var ("pcf" | "all" | "shell") chooses the wrapper.
# Default is "pcf" (matches docker-overnight.sh's default).
# =============================================================================

set -euo pipefail

export POWER_ANALYSIS_ROOT="${POWER_ANALYSIS_ROOT:-/workspace}"
cd "${POWER_ANALYSIS_ROOT}"

PIPELINE="${PIPELINE:-pcf}"

case "${PIPELINE}" in
    pcf)
        exec ./run_power_curve_fitting.sh "$@"
        ;;
    all)
        exec ./run_all.sh "$@"
        ;;
    shell)
        # Drop into interactive bash for debugging; ignore any args.
        exec /bin/bash
        ;;
    *)
        echo "[entrypoint] Unknown PIPELINE='${PIPELINE}'. Expected: pcf | all | shell" >&2
        exit 64
        ;;
esac
