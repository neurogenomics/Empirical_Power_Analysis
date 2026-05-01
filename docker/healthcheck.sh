#!/usr/bin/env bash
# =============================================================================
# docker/healthcheck.sh
#
# Reports container health based on the heartbeat file's mtime. A healthy
# container has touched the heartbeat within HEARTBEAT_STALE_MIN minutes
# (default 30). Otherwise the container is "unhealthy" and `docker ps` will
# flag it; orchestration layers can react without having to parse logs.
#
# This is intentionally lightweight (no R, no jq) so the healthcheck adds
# negligible CPU even on a container running at 100%.
# =============================================================================

set -euo pipefail

HEARTBEAT_FILE="${CBERT_HEARTBEAT_FILE:-/tmp/empirical-power-analysis.heartbeat}"
HEARTBEAT_STALE_MIN="${HEARTBEAT_STALE_MIN:-30}"

if [[ ! -f "${HEARTBEAT_FILE}" ]]; then
    echo "no heartbeat file at ${HEARTBEAT_FILE}" >&2
    exit 1
fi

# stat -c %Y prints mtime in seconds since epoch.
mtime="$(stat -c %Y "${HEARTBEAT_FILE}")"
now="$(date +%s)"
age_sec=$(( now - mtime ))
limit_sec=$(( HEARTBEAT_STALE_MIN * 60 ))

if (( age_sec > limit_sec )); then
    echo "stale heartbeat: ${age_sec}s old (limit ${limit_sec}s)" >&2
    exit 1
fi

echo "ok: heartbeat ${age_sec}s old"
exit 0
