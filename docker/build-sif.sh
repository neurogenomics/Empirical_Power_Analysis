#!/usr/bin/env bash
# =============================================================================
# docker/build-sif.sh -- turn the Docker image into the Apptainer image (.sif)
# that pipeline/run.sh runs on the cluster.
# -----------------------------------------------------------------------------
# Usage (repo root):
#   ./docker/build-sif.sh [--image NAME:TAG | --from-registry REF] [--out FILE]
#
#   --image NAME:TAG     local Docker image to convert
#                        (default: empirical-power-analysis:v1.0.0, the
#                        docker-build.sh default)
#   --from-registry REF  pull a published image instead, ideally by digest, e.g.
#                        ghcr.io/neurogenomics/empirical_power_analysis@sha256:<digest>
#                        (needs only apptainer, no Docker)
#   --out FILE           output path (default: empirical-power-analysis.sif)
#
# The image contains software only, never data. Writes FILE and FILE.sha256,
# then runs docker/smoke_test.R inside the new .sif (with --no-home, as the
# pipeline does) and deletes nothing: a failed smoke test exits non-zero.
# Put the checksum in WORKDIR/epa.env as EPA_SIF_SHA256 to pin a run to it.
#
# Where to build: anywhere with apptainer (or singularity) and, for --image,
# Docker. On a cluster without Docker, either build on a workstation and copy
# the .sif, or use --from-registry.
# =============================================================================
set -euo pipefail
IMAGE="empirical-power-analysis:v1.0.0" REF="" OUT="empirical-power-analysis.sif"
while [ "$#" -gt 0 ]; do
  case "$1" in
    --image)         IMAGE="$2"; shift 2 ;;
    --from-registry) REF="$2"; shift 2 ;;
    --out)           OUT="$2"; shift 2 ;;
    -h|--help)       sed -n '2,/^# =\{70,\}$/p' "$0" | sed '1d;$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
done
APP="$(command -v apptainer || command -v singularity || true)"
[ -n "$APP" ] || { echo "ERROR: apptainer/singularity not found" >&2; exit 2; }
[ -e "$OUT" ] && { echo "ERROR: $OUT exists; remove it first" >&2; exit 2; }

if [ -n "$REF" ]; then
  echo "building $OUT from registry image $REF"
  "$APP" build "$OUT" "docker://$REF"
else
  command -v docker >/dev/null || { echo "ERROR: docker not found (use --from-registry)" >&2; exit 2; }
  docker image inspect "$IMAGE" >/dev/null || { echo "ERROR: no local image $IMAGE (run docker/docker-build.sh)" >&2; exit 2; }
  echo "building $OUT from local Docker image $IMAGE ($(docker image inspect -f '{{.Id}}' "$IMAGE"))"
  # Export to an archive first so apptainer never needs access to the Docker daemon socket.
  tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
  docker save -o "$tmp/image.tar" "$IMAGE"
  "$APP" build "$OUT" "docker-archive://$tmp/image.tar"
fi
sha256sum "$OUT" | tee "$OUT.sha256"
echo "smoke test:"
( cd /tmp && "$APP" exec --cleanenv --no-home "$OUT" Rscript /usr/local/share/epa/smoke_test.R )
echo "done: $OUT"
