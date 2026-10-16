# shellcheck shell=bash
# Stage ad_gt_manifest: as gt_manifest, for the AD vs control contrast and
# cells mode only. The ground truth for run r is a 50% donor split stratified
# by AD status x sex, with the same seed (1234 + r) as the sex grid; each
# subsample then draws each group as close to 50:50 male/female as it can.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  epa_init_ws
  epa_in ad_shard shards_Roussos_ADvsCtrl
  epa_cartesian_env ad
  epa_env DOWNSAMPLE_MODE cells
  epa_step "step 1: ground truth"
  epa_R R/cartesian_downsampling/run_pipeline.R 1 "$EPA_NCPUS"
  epa_step "step 1.5: job manifest ($EPA_AD_GRID_CHUNKS chunks)"
  epa_R R/cartesian_downsampling/run_pipeline.R 1.5 "$EPA_NCPUS"

  local runs; runs="$(printf 'runs%02d-%02d' 1 "$EPA_N_RUNS")"
  local nch; nch="$(printf '%03d' "$EPA_AD_GRID_CHUNKS")"
  epa_require "Cartesian_Downsampling/_hpc_job_manifests/manifest_Roussos_ADvsCtrl_cells_${runs}_${nch}_chunks.csv"
  epa_out control \
    "Cartesian_Downsampling/_hpc_job_manifests" \
    "Cartesian_Downsampling/*/Roussos_ADvsCtrl/*/*/Run_*/Ground_Truth*.qs" \
    "Cartesian_Downsampling/*/Roussos_ADvsCtrl/*/*/Run_*/Validation_Donors.qs" \
    "Cartesian_Downsampling/*/Roussos_ADvsCtrl/*/*/.ground_truth_complete"
  epa_complete
}
