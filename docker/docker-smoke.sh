#!/usr/bin/env bash
# =============================================================================
# docker-smoke.sh -- rerun the image's smoke test (docker/smoke_test.R).
#
# Checks, with no data and no repository mount, that every R package matches
# renv.lock exactly, every Python package matches python-requirements.txt, the
# packages the pipeline uses load, and the plotting fonts exist. The same test
# already runs at the end of docker-build.sh and inside build-sif.sh; this
# script is for checking an image you pulled or built earlier.
#
# Usage (run from the repo root):
#   ./docker/docker-smoke.sh [--image NAME:TAG]
#
# Flags:
#   --image NAME:TAG   Image to test (default: empirical-power-analysis:v1.0.0)
#   -h | --help        Show this help
# =============================================================================
set -euo pipefail

IMAGE="empirical-power-analysis:v1.0.0"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --image) IMAGE="$2"; shift 2 ;;
        -h|--help) sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown flag: $1 (try --help)" >&2; exit 2 ;;
    esac
done

if ! docker info >/dev/null 2>&1; then
    echo "ERROR: cannot reach the Docker daemon (see docker-build.sh for hints)." >&2
    exit 2
fi
if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
    echo "ERROR: image '${IMAGE}' not found. Build it (./docker/docker-build.sh) or pull it." >&2
    exit 2
fi

echo "Smoke-testing ${IMAGE} ($(docker image inspect -f '{{.Id}}' "${IMAGE}"))"
docker run --rm -e "TARGET_UID=$(id -u)" -e "TARGET_GID=$(id -g)" "${IMAGE}" smoke
