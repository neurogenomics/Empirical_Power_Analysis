# shellcheck shell=bash
# Container environment shared by every stage that runs R/cartesian_downsampling.
#   epa_cartesian_env          the paper's sex contrast on Shards/<CT>/Roussos
#   epa_cartesian_env ad       AD vs control on Shards/<CT>/Roussos_ADvsCtrl
#                              (R/shared/contrast.R); its own chunk count
epa_cartesian_env() {
  if [ "${1:-sex}" = ad ]; then
    epa_env EPA_CONTRAST ad_vs_control
    epa_env SHARD_DATASET Roussos_ADvsCtrl
    epa_env JOB_N_CHUNKS "$EPA_AD_GRID_CHUNKS"
  else
    epa_env EPA_CONTRAST sex
    epa_env SHARD_DATASET Roussos
    epa_env JOB_N_CHUNKS "$EPA_GRID_CHUNKS"
  fi
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  epa_env RUN_START 1
  epa_env RUN_END "$EPA_N_RUNS"
}
