#!/usr/bin/env bash
# =============================================================================
# docker-smoke.sh — verify the built image's environment is correct.
#
# Runs the image with PIPELINE=smoke, using the SAME runtime contract as a real
# pipeline run (repo bind-mounted on /workspace, Datalink on
# /workspace/Datalink, dropped to the host UID:GID via gosu). The container's
# embedded smoke test loads every critical R package plus reticulate/anndata
# and prints "SMOKE TEST PASSED" on success.
#
# No raw data is required; this only checks the environment, not the pipeline.
#
# Usage:
#   ./docker-smoke.sh [--image NAME:TAG] [--datalink PATH]
#
# Flags:
#   --image NAME:TAG   Image to test (default: empirical-power-analysis:overnight)
#   --datalink PATH    Host Datalink dir to mount (default: <repo>/Datalink)
#   -h | --help        Show this help
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IMAGE="empirical-power-analysis:overnight"
DATALINK="${SCRIPT_DIR}/Datalink"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --image)    IMAGE="$2"; shift 2 ;;
        --datalink) DATALINK="$2"; shift 2 ;;
        -h|--help) sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown flag: $1 (try --help)" >&2; exit 2 ;;
    esac
done

# ---- Preflight --------------------------------------------------------------
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: cannot reach the Docker daemon (see docker-build.sh for hints)." >&2
    exit 2
fi
if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
    echo "ERROR: image '${IMAGE}' not found. Build it first: ./docker-build.sh" >&2
    exit 2
fi
mkdir -p "${DATALINK}"

# ---- Run smoke test ---------------------------------------------------------
echo "Smoke-testing '${IMAGE}' (host UID:GID $(id -u):$(id -g))"
LOG="$(mktemp)"
trap 'rm -f "${LOG}"' EXIT

set +e
docker run --rm \
    -e PIPELINE=smoke \
    -e POWER_ANALYSIS_ROOT=/workspace \
    -e HOME=/tmp \
    -e "TARGET_UID=$(id -u)" \
    -e "TARGET_GID=$(id -g)" \
    -v "${SCRIPT_DIR}:/workspace" \
    -v "${DATALINK}:/workspace/Datalink" \
    "${IMAGE}" 2>&1 | tee "${LOG}"
RUN_RC="${PIPESTATUS[0]}"
set -e

echo
if [[ "${RUN_RC}" -ne 0 ]]; then
    echo "SMOKE TEST FAILED: container exited ${RUN_RC}." >&2
    exit 4
fi
if ! grep -q "SMOKE TEST PASSED" "${LOG}"; then
    echo "SMOKE TEST FAILED: container exited 0 but never printed 'SMOKE TEST PASSED'." >&2
    exit 4
fi
echo "SMOKE TEST PASSED — image environment is good."
