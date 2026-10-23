# shellcheck shell=bash
# Stage ad_cartesian: aggregate the AD vs control grid and draw its figures
# (cartesian steps 3, 3.5, 4, 5, 6, as for the sex grid), then the
# method-pooled recall and precision panels (by adjusted-p quartile and by
# mean expression) and the AD manuscript figure (tiny-bin recall and
# precision, recall by adjusted-p quartile, tiny-bin recall by expression) for
# each cell type.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  epa_init_ws
  epa_in ad_shard shards_Roussos_ADvsCtrl
  epa_in_array ad_grid
  epa_cartesian_env ad
  epa_env DOWNSAMPLE_MODE cells
  epa_env CARTESIAN_INCLUDE_ALL_METHODS true
  epa_step "steps 3,3.5,4,5,6 (AD vs control, cells)"
  epa_R R/cartesian_downsampling/run_pipeline.R 3,3.5,4,5,6 "$EPA_NCPUS"
  local ct c=Cartesian_Downsampling d=AD_vs_Control/Composed_Figures
  for ct in Astro Micro Oligo; do
    epa_step "recall / precision panels ($ct)"
    epa_R R/figure_composition/compose_recall_panels.R --cell "$ct" --mode cells \
      --dataset Roussos_ADvsCtrl --composed-dir "$d"
    epa_R R/figure_composition/compose_precision_panels.R --cell "$ct" --mode cells \
      --dataset Roussos_ADvsCtrl --composed-dir "$d"
    epa_R R/figure_composition/compose_ad_panels.R --cell "$ct" \
      --dataset Roussos_ADvsCtrl --composed-dir "$d"
    local lc; lc="$(echo "$ct" | tr '[:upper:]' '[:lower:]')"
    epa_require \
      "$c/$ct/Roussos_ADvsCtrl/cells/all_methods/Plots/Mean_Recall_Total.png" \
      "$c/$ct/Roussos_ADvsCtrl/cells/all_methods/Plots/Mean_Precision_Total.png" \
      "$c/$ct/Roussos_ADvsCtrl/cells/all_methods/Continuous_Viz_Log/Mean_Recall_Total_Continuous_SemiLog.png" \
      "$c/$ct/Roussos_ADvsCtrl/cells/all_methods/Continuous_Viz_Log/Mean_Precision_Total_Continuous_SemiLog.png" \
      "$d/figure_recall_panels_${lc}_cells.png" "$d/figure_precision_panels_${lc}_cells.png" \
      "$d/figure_ad_panels_${lc}_cells.png"
  done
  epa_out ad_cartesian \
    "$c/*/Roussos_ADvsCtrl/cells/*/Plots" \
    "$c/*/Roussos_ADvsCtrl/cells/*/Continuous_Viz" \
    "$c/*/Roussos_ADvsCtrl/cells/*/Continuous_Viz_Log" \
    "$d"
  epa_complete
}
