# shellcheck shell=bash
# Stage ad_de: full-cohort DE for both contrasts with the grid's DE workers,
# then the effect-size / expression comparison of AD vs control DEGs with
# autosomal sex-biased DEGs (R/ad_contrast 01 and 02).
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_in ad_shard shards_Roussos_ADvsCtrl
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_step "full-cohort DE: sex (all filtered donors)"
  epa_env EPA_CONTRAST sex
  epa_env SHARD_DATASET Roussos
  epa_R R/ad_contrast/01_Full_Cohort_DE.R
  epa_step "full-cohort DE: AD vs control"
  epa_env EPA_CONTRAST ad_vs_control
  epa_env SHARD_DATASET Roussos_ADvsCtrl
  epa_R R/ad_contrast/01_Full_Cohort_DE.R
  epa_step "effect distributions"
  epa_R R/ad_contrast/02_Compare_Effect_Distributions.R
  local d=AD_vs_Control
  epa_require "$d/Effect_Distributions/effect_distributions_pooled.png" \
              "$d/Effect_Distributions/effect_distribution_summary.csv"
  epa_out ad_de "$d/Full_Cohort_DE/*.csv" "$d/Effect_Distributions"
  epa_complete
}
