# shellcheck shell=bash
# =============================================================================
# Stage paper_values -- every data-derived number quoted in the manuscript.
#   R/figure_data_generation step 9 -> Paper_Source_Data/paper_values.csv
#   (key, value, unit, definition, source), read by val() at knit time.
# =============================================================================
stage_main() {
  epa_init_ws
  epa_in eda
  epa_in tables
  epa_in deg_loss
  epa_in pcf_aggregate fitting_data
  epa_in xgb_select model
  epa_in xgb_post
  epa_in gene_class
  epa_step "figure data step 9: paper values"
  epa_R R/figure_data_generation/run_pipeline.R 9
  epa_require Paper_Source_Data/paper_values.csv
  epa_out paper_values Paper_Source_Data/paper_values.csv
  epa_complete
}
