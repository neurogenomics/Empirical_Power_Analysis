# shellcheck shell=bash
# Stage fgt: full-cohort ("full ground truth") sex DE per cell type and the gene
# annotation audit.
#   R/full_ground_truth steps 1, 4  DESeq2 on all PsychAD donors; gene universe by
#                                   biotype / chromosome (Supplementary Data 1 sheet 2)
#   R/de_diagnostics step 4         DEG chromosome distribution (Supp. Fig.)
# Gene annotation comes from the shipped biomaRt caches (BIOMART_OFFLINE=true).
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_ref biomart/gene_audit_biomart_cache.csv Full_Ground_Truth/gene_audit_biomart_cache.csv
  epa_ref biomart/gene_biotype_cache.csv DE_Diagnostics/gene_biotype_cache.csv
  epa_env BIOMART_OFFLINE true
  local half=$(( EPA_NCPUS / 2 > 0 ? EPA_NCPUS / 2 : 1 ))
  epa_step "full ground truth DE (step 1) and gene audit (step 4)"
  epa_R R/full_ground_truth/run_pipeline.R 1,4 "$half" "$half"
  epa_step "DEG chromosome grid (de_diagnostics step 4)"
  epa_R R/de_diagnostics/run_pipeline.R 4
  local ct
  for ct in Astro Micro Oligo; do
    epa_require "Full_Ground_Truth/$ct/Roussos_DE_results.csv" "Full_Ground_Truth/$ct/Roussos_dds_checkpoint.qs"
  done
  epa_require Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv \
              DE_Diagnostics/Chromosome_Grid/DEG_Chr_Grid_Combined.png
  epa_out fgt Full_Ground_Truth DE_Diagnostics Run_Output_Logs
  epa_complete
}
