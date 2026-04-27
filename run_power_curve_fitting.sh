#!/usr/bin/env bash
# =============================================================================
# run_power_curve_fitting.sh
#
# Power-curve-only driver. Runs the minimum DAG required to (re)generate the
# fitted power equations from raw shards, skipping pipelines that don't feed
# into power_curve_fitting:
#
#       sharding      -----------\
#           |                     \
#           v                      v
#  cartesian_downsampling   depth_equalization
#  (cells + reads)                |
#           \                    /
#            \                  /
#             v                v
#           power_curve_fitting
#
# Skipped: eda_roussos, bulk_deg_processing, de_diagnostics,
#          figure_data_generation
#
# Implemented as a thin wrapper around run_all.sh so flags, logging, lock
# handling, and FORCE_RECALC plumbing all behave identically. Anything you
# pass on the command line is forwarded to run_all.sh; the only thing this
# script enforces is the --skip set above.
#
# Usage (same flags as run_all.sh):
#   ./run_power_curve_fitting.sh
#   ./run_power_curve_fitting.sh --max-cores 24 --no-overwrite
#   ./run_power_curve_fitting.sh --force                     # regenerate all
#   ./run_power_curve_fitting.sh --mode cells                # cartesian only cells
#   ./run_power_curve_fitting.sh --from cartesian_downsampling
#   ./run_power_curve_fitting.sh --dry-run
#
# Recommended unattended invocation:
#   nohup ./run_power_curve_fitting.sh --max-cores 24 --no-overwrite \
#       > nohup.run_pcf.out 2>&1 &
#   disown
# =============================================================================

set -u
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Stages that are NOT on the power-curve-fitting critical path
SKIP="full_ground_truth,eda_roussos,bulk_deg_processing,de_diagnostics,figure_data_generation"

# Forward every CLI arg straight through; if the user passes their own --skip
# we merge the two skip lists so we never accidentally enable a skipped stage.
USER_SKIP=""
PASSTHRU=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip) USER_SKIP="$2"; shift 2 ;;
        *)      PASSTHRU+=("$1"); shift ;;
    esac
done

if [[ -n "${USER_SKIP}" ]]; then
    SKIP="${SKIP},${USER_SKIP}"
fi

echo "==> run_power_curve_fitting.sh dispatching to run_all.sh"
echo "    skip stages:       ${SKIP}"
echo "    sharding datasets: roussos  (PCF only consumes Roussos shards)"
echo "    extra args:        ${PASSTHRU[*]:-<none>}"

exec "${SCRIPT_DIR}/run_all.sh" \
    --skip "${SKIP}" \
    --sharding-datasets roussos \
    "${PASSTHRU[@]}"
