# shellcheck shell=bash
# Stage deg_loss: astrocyte DEG loss under depth / cell-count equalisation
# (Supplementary Table 3).
#   R/depth_equalization/run_pipeline.R  multinomial downsampling of every astrocyte to
#                                        microglia-equivalent counts per cell (DepthEq) and
#                                        depth + cell-count matching (ShallowFewer)
#   deg_loss_analysis.R                  cell-count-only variant (CellsOnly), then DE for
#                                        all variants x 10 runs x 3 methods
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  epa_step "depth equalisation shards"
  epa_R R/depth_equalization/run_pipeline.R
  epa_step "DEG loss analysis"
  epa_R R/depth_equalization/deg_loss_analysis.R "$EPA_NCPUS"
  epa_require Depth_Eq_DEG_Loss/deg_loss_table.csv
  epa_out deg_loss "Depth_Eq_DEG_Loss/*.csv"
  epa_complete
}
