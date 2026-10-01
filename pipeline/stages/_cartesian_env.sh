# shellcheck shell=bash
# Container environment shared by every stage that runs R/cartesian_downsampling.
epa_cartesian_env() {
  epa_env SHARD_DATASET Roussos
  epa_env DE_METHODS "$EPA_DE_METHODS"
  epa_env N_ROBUSTNESS_RUNS "$EPA_N_RUNS"
  epa_env RUN_START 1
  epa_env RUN_END "$EPA_N_RUNS"
  epa_env JOB_N_CHUNKS "$EPA_GRID_CHUNKS"
}
