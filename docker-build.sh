#!/usr/bin/env bash
# =============================================================================
# docker-build.sh — build the Empirical Power Analysis image.
#
# The simple counterpart to docker-overnight.sh: it does ONE thing — build the
# image — with clear preflight checks and readable errors. No manifest, no
# heartbeat, no pipeline run. Pair it with docker-smoke.sh to verify the build.
#
# Usage:
#   ./docker-build.sh [--image NAME:TAG]
#
# Flags:
#   --image NAME:TAG   Image tag to build (default: empirical-power-analysis:overnight)
#   -h | --help        Show this help
#
# The Dockerfile is BuildKit-optional, so this builds on any host: if the
# buildx plugin is present we use BuildKit; otherwise we fall back to the
# legacy builder. Either way the image is identical.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IMAGE="empirical-power-analysis:overnight"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --image) IMAGE="$2"; shift 2 ;;
        -h|--help) sed -n '2,/^# =\{60,\}/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown flag: $1 (try --help)" >&2; exit 2 ;;
    esac
done

# ---- Preflight --------------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker CLI not found on PATH. Install Docker first." >&2
    exit 2
fi
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: cannot reach the Docker daemon." >&2
    echo "  - Is it running?           (e.g. 'sudo service docker start')" >&2
    echo "  - Are you in the 'docker' group?  'sudo usermod -aG docker \$USER'" >&2
    echo "    then start a new login shell (or 'newgrp docker')." >&2
    exit 2
fi
# Pick the builder: BuildKit if buildx is available, else the legacy builder.
# The Dockerfile builds correctly under both. --progress=plain is a BuildKit-only
# flag, so only pass it when BuildKit is in use (the legacy builder rejects it).
BUILD_FLAGS=()
if docker buildx version >/dev/null 2>&1; then
    BUILDKIT=1
    BUILD_FLAGS+=(--progress=plain)
    echo "Builder: BuildKit (buildx present)"
else
    BUILDKIT=0
    echo "Builder: legacy (buildx not installed; this is fine, just slower)"
fi

# ---- Build ------------------------------------------------------------------
echo "Building image '${IMAGE}' from ${SCRIPT_DIR}"
echo "(this can take 15-30 min on a cold cache; R/Bioconductor packages compile)"
DOCKER_BUILDKIT="${BUILDKIT}" docker build "${BUILD_FLAGS[@]}" -t "${IMAGE}" "${SCRIPT_DIR}"

echo
echo "Build complete: ${IMAGE}"
docker image inspect --format '  id:      {{.Id}}' "${IMAGE}" 2>/dev/null || true
docker image inspect --format '  size:    {{.Size}} bytes' "${IMAGE}" 2>/dev/null || true
echo "Next: ./docker-smoke.sh --image ${IMAGE}"
