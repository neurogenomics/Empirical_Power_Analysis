# shellcheck shell=bash
# Shared body of stages cartesian_cells / cartesian_reads: aggregate the grid
# for one mode and draw its figures (cartesian steps 3, 3.5, 4, 5, 6).
#   3    per-method aggregation over all robustness runs -> Plots/Aggregated_Summary_Data.qs
#   3.5  trial-weighted pooling across methods           -> all_methods/Plots/...
#   4    discrete recall/precision grids                 -> <method>/Plots/*.png
#   5    continuous (semi-log) surfaces                  -> <method>/Continuous_Viz_Log/*.png
#   6    tiny/small-effect expression panels             -> <method>/Plots/Figure_TS/*.png
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

cartesian_post() {
  local mode="$1"
  epa_init_ws
  epa_in shard shards_Roussos
  epa_in_array grid "epa_${mode}_"
  epa_cartesian_env
  epa_env DOWNSAMPLE_MODE "$mode"
  epa_env CARTESIAN_INCLUDE_ALL_METHODS true
  epa_step "steps 3,3.5,4,5,6 ($mode)"
  epa_R R/cartesian_downsampling/run_pipeline.R 3,3.5,4,5,6 "$EPA_NCPUS"
  local ct
  for ct in Astro Micro Oligo; do
    epa_require \
      "Cartesian_Downsampling/$ct/Roussos/$mode/all_methods/Plots/Aggregated_Summary_Data.qs" \
      "Cartesian_Downsampling/$ct/Roussos/$mode/all_methods/Plots/Mean_Recall_Total.png" \
      "Cartesian_Downsampling/$ct/Roussos/$mode/all_methods/Continuous_Viz_Log/Mean_Recall_Total_Continuous_SemiLog.png" \
      "Cartesian_Downsampling/$ct/Roussos/$mode/all_methods/Plots/Figure_TS/Figure_T_Recall_${mode}.png"
  done
  epa_out "cartesian_${mode}" \
    "Cartesian_Downsampling/*/Roussos/$mode/*/Plots" \
    "Cartesian_Downsampling/*/Roussos/$mode/*/Continuous_Viz" \
    "Cartesian_Downsampling/*/Roussos/$mode/*/Continuous_Viz_Log"
  epa_complete
}
