#!/usr/bin/env bash
# ==============================================================================
# Fire-and-forget wrapper: XGBoost realistic-grid reweighting sweep.
# ==============================================================================
# Submits Stages B-realistic (search array), C-realistic (select+SHAP refit) and
# D-realistic (SHAP figures) as ONE dependency-chained batch on Imperial CX3.
# Each stage only starts if the previous one completed successfully (afterok /
# afterokarray). No polling, no babysitting -- fire it, walk away.
#
# Precondition (checked here):
#   Stage A must already have written HPC_XGB_Input/xgb-aggregated-<tag>.tar,
#   AND that tar must contain aggregated_fitting_data_realistic.qs (produced by
#   the updated 01e_Filter_Realistic.R step in cx3-xgb-aggregate.pbs).
#   If the tar is missing OR does not include the realistic .qs, this wrapper
#   refuses to submit and tells you to re-run Stage A first.
#
# Usage (from Empirical_Power_Analysis root on CX3):
#   ./hpc/submit-xgb-realistic-chain.sh
#
# Optional environment overrides (same names as the individual PBS scripts):
#   EPA_SIF                  path to the apptainer .sif (default: $PWD/epa-cx3-2026-07-19.sif)
#   EPA_DATALINK             RDS Datalink dir           (default: $PWD/Datalink)
#   EPA_DE_METHODS           comma list                 (default: deseq2,edger,limma)
#   EPA_XGB_TAG              tag for the sweep          (default: <methods>_runs01-10)
#   EPA_XGB_SHAP_MAX_ROWS    SHAP subsample rows        (default: 200000)
# ==============================================================================
set -euo pipefail

# ---- resolve repo root --------------------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
EPA_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
cd "$EPA_ROOT"

# ---- defaults -----------------------------------------------------------------
: "${EPA_SIF:=$EPA_ROOT/epa-cx3-2026-07-19.sif}"
: "${EPA_DATALINK:=$EPA_ROOT/Datalink}"
: "${EPA_DE_METHODS:=deseq2,edger,limma}"
: "${EPA_XGB_SHAP_MAX_ROWS:=200000}"
_safe_methods="${EPA_DE_METHODS//,/-}"
: "${EPA_XGB_TAG:=${_safe_methods}_runs01-10}"

export EPA_SIF EPA_DATALINK EPA_DE_METHODS EPA_XGB_TAG EPA_XGB_SHAP_MAX_ROWS

AGG_TAR="$EPA_DATALINK/HPC_XGB_Input/xgb-aggregated-${EPA_XGB_TAG}.tar"

echo "================ submit-xgb-realistic-chain ================"
echo "repo root     : $EPA_ROOT"
echo "sif           : $EPA_SIF"
echo "datalink      : $EPA_DATALINK"
echo "de methods    : $EPA_DE_METHODS"
echo "xgb tag       : $EPA_XGB_TAG"
echo "aggregated tar: $AGG_TAR"
echo "============================================================"

# ---- sanity checks ------------------------------------------------------------
[ -f "$EPA_SIF" ] || { echo "FATAL: sif missing: $EPA_SIF" >&2; exit 2; }
[ -d "$EPA_DATALINK" ] || { echo "FATAL: datalink dir missing: $EPA_DATALINK" >&2; exit 2; }

if [ ! -f "$AGG_TAR" ]; then
  cat >&2 <<EOF
FATAL: aggregated tar missing:
  $AGG_TAR
Run Stage A first:
  qsub -v EPA_SIF=$EPA_SIF,EPA_DATALINK=$EPA_DATALINK hpc/cx3-xgb-aggregate.pbs
Then re-run this wrapper.
EOF
  exit 3
fi

if ! tar tf "$AGG_TAR" 2>/dev/null | grep -q 'aggregated_fitting_data_realistic\.qs$'; then
  cat >&2 <<EOF
FATAL: aggregated tar does NOT contain aggregated_fitting_data_realistic.qs:
  $AGG_TAR
This tar predates the 01e_Filter_Realistic.R Stage A extension. Regenerate:
  1. Delete (or rename) the tar so Stage A doesn't short-circuit:
       mv "$AGG_TAR" "$AGG_TAR.pre-01e.bak"
  2. Re-submit Stage A:
       qsub -v EPA_SIF=$EPA_SIF,EPA_DATALINK=$EPA_DATALINK hpc/cx3-xgb-aggregate.pbs
  3. Re-run this wrapper.
EOF
  exit 3
fi

echo "OK: aggregated tar present and includes realistic .qs. Submitting chain."

# PBS Pro's -v uses commas as VAR separators; values containing commas
# (EPA_DE_METHODS=deseq2,edger,limma) must be single-quoted inside the -v string.
# Build the shared -v payload once with the methods value pre-quoted.
VARS_COMMON="EPA_SIF=$EPA_SIF,EPA_DATALINK=$EPA_DATALINK,EPA_DE_METHODS='$EPA_DE_METHODS',EPA_XGB_TAG=$EPA_XGB_TAG"

# ---- Stage B-realistic: search array (36 tasks) -------------------------------
JOB_B=$(qsub -v "$VARS_COMMON" \
             "$EPA_ROOT/hpc/cx3-xgb-realistic-search.pbs")
echo "Stage B (realistic search) : $JOB_B"

# Imperial CX3's PBS Pro accepts `afterok:<full-id>` where <full-id> is qsub's
# raw returned id including the `[].server` suffix (e.g. `3516535[].pbs-7`).
# The whole -W value MUST be quoted so the shell doesn't glob the `[]`.
# `afterokarray` was rejected on this cluster ("illegal -W value") -- plain
# `afterok` on the array parent is the working form. Confirmed 2026-07-30.
JOB_B_DEP="$JOB_B"

# ---- Stage C-realistic: select + 5-model SHAP refit ---------------------------
JOB_C=$(qsub -W "depend=afterok:${JOB_B_DEP}" \
             -v "$VARS_COMMON" \
             "$EPA_ROOT/hpc/cx3-xgb-realistic-select.pbs")
echo "Stage C (realistic select) : $JOB_C  (depends on $JOB_B_DEP)"
JOB_C_DEP="$JOB_C"

# ---- Stage D-realistic: SHAP figures -----------------------------------------
JOB_D=$(qsub -W "depend=afterok:${JOB_C_DEP}" \
             -v "$VARS_COMMON,EPA_XGB_SHAP_MAX_ROWS=$EPA_XGB_SHAP_MAX_ROWS" \
             "$EPA_ROOT/hpc/cx3-xgb-realistic-shap.pbs")
echo "Stage D (realistic shap)   : $JOB_D  (depends on $JOB_C_DEP)"

echo "============================================================"
echo "Chain submitted. Monitor with:"
echo "  qstat -tanw1 -u \$USER"
echo "Outputs will land in:"
echo "  Stage B: $EPA_DATALINK/HPC_XGB_Search/incoming/$EPA_XGB_TAG-realistic/config_*.tar"
echo "  Stage C: $EPA_DATALINK/HPC_XGB_Model/xgb-model-${EPA_XGB_TAG}-realistic.tar"
echo "  Stage D: $EPA_DATALINK/HPC_XGB_Figures/shap-figures-${EPA_XGB_TAG}-realistic.tar"
echo "============================================================"
