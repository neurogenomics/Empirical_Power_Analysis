# shellcheck shell=bash
# Stage tables: source tables that combine the full-cohort DE with the grid.
#   figure_data_generation 6   autosomal vs sex-chromosome truth genes (Supplementary
#                              Table 2; truth-set sizes quoted in the Results)
#   figure_data_generation 8   tidy recall/precision grid (numbers quoted in the Results)
stage_main() {
  epa_init_ws
  epa_in fgt
  epa_in pcf_aggregate gt_dds
  epa_in cartesian_cells
  epa_in cartesian_reads
  epa_env BIOMART_OFFLINE true
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  epa_step "figure data steps 6 and 8"
  epa_R R/figure_data_generation/run_pipeline.R 6,8
  epa_require _scratch_autosomal_vs_sexchr_summary.csv \
              Paper_Source_Data/recall_precision_grid.csv
  epa_out tables "_scratch_autosomal_vs_sexchr_*.csv" Paper_Source_Data/recall_precision_grid.csv
  epa_complete
}
