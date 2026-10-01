# shellcheck shell=bash
# Stage pcf_aggregate: turn the grid into the power-model training table
# (R/power_curve_fitting 00, 00b, 01 per DE method, 01b, 01e).
#   00   full-truth-half DESeq2 fits per run (gene dispersion / baseMean covariates)
#   00b  per-donor read-count cache
#   01   per-iteration detection outcomes, once per DE method
#   01b  aggregate to detection rate per (gene, design point)  -> aggregated_fitting_data_all.qs
#   01e  realistic sub-grid (n_cells >= 100)                    -> aggregated_fitting_data_realistic.qs
# Outputs: fitting_data.tar (the two training tables) and gt_dds.tar (the
# DESeq2 truth-half checkpoints and the DESeq2 truth-gene set, which the
# supplementary tables use).
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_in_array grid
  epa_env SHARD_DATASET Roussos
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"

  epa_step "00: truth-half DESeq2 checkpoints"
  epa_env CARTESIAN_DE_METHOD deseq2
  epa_R R/power_curve_fitting/00_Generate_GT_DDS_Checkpoints.R
  local ct r
  for ct in Astro Micro Oligo; do
    for r in $(seq 1 "$EPA_N_RUNS"); do
      epa_require "Power_Curve_Fitting/GT_DDS_Checkpoints/$ct/GT_DDS_Run_$r.qs"
    done
  done

  epa_step "00b: shard read-count cache"
  epa_R R/power_curve_fitting/00b_Cache_Shard_Read_Counts.R

  # 01 writes unified_fitting_data_per_iteration.qs; rename per method so 01b
  # sees all of them. It also rewrites Fitting_Data/shared_truth_genes.qs on
  # every call, so keep the DESeq2 copy (the truth set the supplementary truth
  # characteristics are reported for).
  local fd="$WS/Datalink/Power_Curve_Fitting/Fitting_Data" m f
  IFS=',' read -r -a methods <<< "$EPA_DE_METHODS"
  for m in "${methods[@]}"; do
    epa_step "01: fitting data for $m"
    epa_env CARTESIAN_DE_METHOD "$m"
    epa_R R/power_curve_fitting/01_Build_Fitting_Dataset.R
    while IFS= read -r -d '' f; do
      mv -f "$f" "${f%.qs}__${m}.qs"
    done < <(find "$fd" -type f -name 'unified_fitting_data_per_iteration.qs' -print0)
    if [ "$m" = deseq2 ]; then cp -f "$fd/shared_truth_genes.qs" "$WS/.shared_truth_genes_deseq2.qs"; fi
  done
  cp -f "$WS/.shared_truth_genes_deseq2.qs" "$fd/shared_truth_genes.qs"

  epa_step "01b: aggregate; 01e: realistic sub-grid"
  epa_env CARTESIAN_DE_METHOD deseq2
  epa_R R/power_curve_fitting/01b_Aggregate_Fitting_Data.R
  epa_R R/power_curve_fitting/01e_Filter_Realistic.R

  epa_out fitting_data \
    Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_all.qs \
    Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_realistic.qs \
    "Power_Curve_Fitting/Fitting_Data/filter_*.csv"
  epa_out gt_dds \
    Power_Curve_Fitting/GT_DDS_Checkpoints \
    Power_Curve_Fitting/Fitting_Data/shared_truth_genes.qs
  epa_complete
}
