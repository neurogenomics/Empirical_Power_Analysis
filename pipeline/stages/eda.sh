# shellcheck shell=bash
# Stage eda: analyses that need only the shards.
#   R/eda_roussos step 1          PCA gallery, unfiltered and filtered cohorts (Supp. Fig. PCA)
#   Validate_Sex_Labels.R All     XIST / Y-gene check of sex labels (Supplementary Table 4)
#   figure_data_generation 1      cohort statistics (dataset_statistics.csv: donor counts
#                                 quoted in the text, Supplementary Data 1 sheet 1)
stage_main() {
  epa_init_ws
  epa_in shard
  epa_step "PCA gallery"
  epa_R R/eda_roussos/run_pipeline.R 1
  epa_step "sex label validation"
  epa_R R/eda_roussos/Validate_Sex_Labels.R All
  epa_step "dataset statistics"
  epa_R R/figure_data_generation/run_pipeline.R 1
  local ct ds
  for ct in Astro Micro Oligo; do
    epa_require "Exploratory_Analysis/$ct/Roussos_Full/PCA/PCA_FullCohort_sex_scaled.png" \
                "Exploratory_Analysis/$ct/Roussos/PCA/PCA_FilteredCohort_sex_scaled.png"
  done
  epa_require Sex_Validation/sex_exclusion_combined.csv Paper_Source_Data/dataset_statistics.csv
  epa_out eda Exploratory_Analysis Sex_Validation Paper_Source_Data/dataset_statistics.csv
  epa_complete
}
