# shellcheck shell=bash
# Stage ad_pcf: turn the AD vs control grid into a table with exactly the
# power model's training layout (R/power_curve_fitting 00, 00b, 01 per DE
# method, 01b, 01e with EPA_CONTRAST=ad_vs_control), for use as an external
# test set. Gene dispersions come from DESeq2 fits on each run's AD truth half.
# Outputs land under Datalink/AD_vs_Control/Power_Curve_Fitting/.
stage_main() {
  epa_init_ws
  epa_in ad_shard shards_Roussos_ADvsCtrl
  epa_in_array ad_grid
  epa_env EPA_CONTRAST ad_vs_control
  epa_env SHARD_DATASET Roussos_ADvsCtrl
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  local p=AD_vs_Control/Power_Curve_Fitting

  epa_step "00: truth-half DESeq2 checkpoints"
  epa_env CARTESIAN_DE_METHOD deseq2
  epa_R R/power_curve_fitting/00_Generate_GT_DDS_Checkpoints.R
  epa_step "00b: shard read-count cache"
  epa_R R/power_curve_fitting/00b_Cache_Shard_Read_Counts.R

  local fd="$WS/Datalink/$p/Fitting_Data" m f
  IFS=',' read -r -a methods <<< "$EPA_DE_METHODS"
  for m in "${methods[@]}"; do
    epa_step "01: test rows for $m"
    epa_env CARTESIAN_DE_METHOD "$m"
    epa_R R/power_curve_fitting/01_Build_Fitting_Dataset.R
    while IFS= read -r -d '' f; do
      mv -f "$f" "${f%.qs}__${m}.qs"
    done < <(find "$fd" -type f -name 'unified_fitting_data_per_iteration.qs' -print0)
  done
  epa_step "01b: aggregate; 01e: n_cells >= 100"
  epa_env CARTESIAN_DE_METHOD deseq2
  epa_R R/power_curve_fitting/01b_Aggregate_Fitting_Data.R
  epa_R R/power_curve_fitting/01e_Filter_Design.R
  epa_require "$p/Fitting_Data/aggregated_fitting_data.qs"
  epa_out ad_fitting_data "$p/Fitting_Data/aggregated_fitting_data.qs" "$p/Fitting_Data/filter_*.csv"
  epa_complete
}
