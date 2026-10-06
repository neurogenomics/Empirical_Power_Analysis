# shellcheck shell=bash
# Stage deg_loss: DEG loss under depth / cell-count equalisation towards
# microglia, for astrocytes and oligodendrocytes (Supplementary Table 3).
# Per cell type:
#   R/depth_equalization/run_pipeline.R  multinomial downsampling of every cell to
#                                        microglia-equivalent counts per cell (DepthEq)
#                                        and depth + cell-count matching (ShallowFewer)
#   deg_loss_analysis.R                  cell-count-only variant (CellsOnly), then DE for
#                                        all variants x 10 runs x 3 methods
# then combine_deg_loss_tables.R stacks the two cell types into one table.
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  local ct
  for ct in Astro Oligo; do
    epa_env DEPTH_EQ_CELL_TYPE "$ct"
    epa_step "$ct: depth equalisation shards"
    epa_R R/depth_equalization/run_pipeline.R
    epa_step "$ct: DEG loss analysis"
    epa_R R/depth_equalization/deg_loss_analysis.R "$EPA_NCPUS"
    epa_require "Depth_Eq_DEG_Loss/$ct/deg_loss_table.csv"
    # The equalised shards are only needed by this cell type's analysis.
    rm -rf "$WS/Datalink/Shards/$ct/Roussos_DepthEq" "$WS/Datalink/Shards/$ct/Roussos_ShallowFewer" \
           "$WS/Datalink/Shards/$ct/Roussos_CellsOnly" "$WS/Datalink/Depth_Eq_DEG_Loss/$ct/GT_Cache"
  done
  epa_step "combine cell types"
  epa_R R/depth_equalization/combine_deg_loss_tables.R
  epa_require Depth_Eq_DEG_Loss/deg_loss_table.csv
  epa_out deg_loss "Depth_Eq_DEG_Loss/*.csv" Depth_Eq_DEG_Loss/Astro Depth_Eq_DEG_Loss/Oligo
  epa_complete
}
