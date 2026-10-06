# shellcheck shell=bash
# Stage ad_shard: donor shards for the AD vs control contrast. Keeps donors whose
# only recorded diagnosis is AD and donors with no recorded diagnosis (every
# diagnosis evaluated); see R/sharding/Build_AD_Control_Shards.R.
# Outputs: shards_Roussos_ADvsCtrl.tar, ad_donors.tar (group counts).
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_step "select AD-only and no-diagnosis donors"
  epa_R R/sharding/Build_AD_Control_Shards.R "$EPA_NCPUS"
  local ct
  for ct in Astro Micro Oligo; do
    epa_require "Shards/$ct/Roussos_ADvsCtrl/00_MASTER_MANIFEST.qs"
  done
  epa_out shards_Roussos_ADvsCtrl "Shards/*/Roussos_ADvsCtrl"
  epa_out ad_donors AD_vs_Control/ad_control_donor_counts.csv
  epa_complete
}
