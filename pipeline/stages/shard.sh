# shellcheck shell=bash
# Stage shard: raw PsychAD .h5ad files -> one SingleCellExperiment per donor.
#   R/sharding/Run_Shard_Roussos.R runs, for Astro, Micro and Oligo:
#     Run_Shard_Roussos_Full.R  all donors (Roussos_Full, used for the unfiltered PCA)
#     filter                    donors aged >= 40 with >= 30 cells (Roussos)
#     Annotate_Roussos_PMI.R    post-mortem interval from media-1.csv (a DE covariate)
# Outputs: shards_Roussos.tar (analysis cohort), shards_Roussos_Full.tar.
stage_main() {
  epa_init_ws
  if [ "${EPA_VERIFY_INPUTS:-true}" = true ]; then
    epa_step "verify raw inputs against pipeline/inputs.sha256"
    ( cd "$EPA_RAW_DIR" && sha256sum -c "$WS/pipeline/inputs.sha256" ) \
      || epa_die "raw inputs do not match pipeline/inputs.sha256"
  fi
  epa_step "shard PsychAD (Astro, Micro, Oligo)"
  epa_R R/sharding/Run_Shard_Roussos.R Astro,Micro,Oligo
  local ct
  for ct in Astro Micro Oligo; do
    epa_require "Shards/$ct/Roussos/00_MASTER_MANIFEST.qs" "Shards/$ct/Roussos/.filter_complete" \
                "Shards/$ct/Roussos/.pmi_annotated" "Shards/$ct/Roussos_Full/00_MASTER_MANIFEST.qs"
  done
  epa_out shards_Roussos "Shards/*/Roussos"
  epa_out shards_Roussos_Full "Shards/*/Roussos_Full"
  epa_complete
}
