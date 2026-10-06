# shellcheck shell=bash
# Stage ad_grid (array): the AD vs control grid, cells mode, one task per chunk.
# Same step-2 code, seeds and design grid as the sex grid; design points that
# need more donors than the AD validation half holds are dropped automatically.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  local chunk="${EPA_TASK##*:}"
  epa_init_ws
  epa_in ad_shard shards_Roussos_ADvsCtrl
  epa_in ad_gt_manifest control
  epa_cartesian_env ad
  epa_env DOWNSAMPLE_MODE cells
  epa_env JOB_CHUNK "$chunk"
  epa_step "step 2: AD cells chunk $chunk of $EPA_AD_GRID_CHUNKS"
  epa_R R/cartesian_downsampling/run_pipeline.R 2 "$EPA_NCPUS"
  epa_out "$EPA_TASK_OUT" "Cartesian_Downsampling/*/Roussos_ADvsCtrl/cells"
}
