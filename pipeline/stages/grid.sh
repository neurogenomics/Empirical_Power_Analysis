# shellcheck shell=bash
# Stage grid (array): one task per (mode, chunk). Each task runs cartesian step 2
# for its share of the donors x cells|reads x iteration design, for every
# robustness run, cell type and DE method. A job's random seed depends only on
# (run, donors, mode, depth, iteration), so results do not depend on how the
# grid is chunked. Output: one tar per task with that task's per-iteration DE
# results (plus the mode's ground-truth copies, identical across tasks).
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  local mode="${EPA_TASK%%:*}" chunk="${EPA_TASK##*:}"
  epa_init_ws
  epa_in shard shards_Roussos
  epa_in gt_manifest control
  epa_cartesian_env
  epa_env DOWNSAMPLE_MODE "$mode"
  epa_env JOB_CHUNK "$chunk"
  epa_step "step 2: $mode chunk $chunk of $EPA_GRID_CHUNKS"
  epa_R R/cartesian_downsampling/run_pipeline.R 2 "$EPA_NCPUS"
  epa_out "$EPA_TASK_OUT" "Cartesian_Downsampling/*/Roussos/$mode"
}
