# shellcheck shell=bash
# Stage gt_manifest: ground truth for every robustness run plus the grid job
# manifests (cartesian steps 1 and 1.5). The ground truth for run r is a
# sex-stratified 50% donor split (seed 1234 + r) analysed with all DE methods.
# Output: control.tar (manifests + Ground_Truth*.qs + Validation_Donors.qs for
# both modes), which every grid task unpacks.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_cartesian_env

  epa_step "step 1: ground truth (written under cells/, content is mode-independent)"
  epa_env DOWNSAMPLE_MODE cells
  epa_R R/cartesian_downsampling/run_pipeline.R 1 "$EPA_NCPUS"

  epa_step "copy ground truth cells/ -> reads/"
  local src dst
  while IFS= read -r -d '' src; do
    dst="${src/\/cells\//\/reads\/}"
    mkdir -p "$(dirname "$dst")"
    cp -f "$src" "$dst"
  done < <(find "$WS/Datalink/Cartesian_Downsampling" -path '*/cells/*' \
             \( -name 'Ground_Truth*.qs' -o -name 'Validation_Donors.qs' \
                -o -name '.ground_truth_complete' \) -print0)

  epa_step "step 1.5: job manifests ($EPA_GRID_CHUNKS chunks per mode)"
  epa_env DOWNSAMPLE_MODE cells
  epa_R R/cartesian_downsampling/run_pipeline.R 1.5 "$EPA_NCPUS"
  epa_env DOWNSAMPLE_MODE reads
  epa_R R/cartesian_downsampling/run_pipeline.R 1.5 "$EPA_NCPUS"

  local runs; runs="$(printf 'runs%02d-%02d' 1 "$EPA_N_RUNS")"
  local nch; nch="$(printf '%03d' "$EPA_GRID_CHUNKS")"
  epa_require \
    "Cartesian_Downsampling/_hpc_job_manifests/manifest_Roussos_cells_${runs}_${nch}_chunks.csv" \
    "Cartesian_Downsampling/_hpc_job_manifests/manifest_Roussos_reads_${runs}_${nch}_chunks.csv"
  epa_out control \
    "Cartesian_Downsampling/_hpc_job_manifests" \
    "Cartesian_Downsampling/*/Roussos/*/*/Run_*/Ground_Truth*.qs" \
    "Cartesian_Downsampling/*/Roussos/*/*/Run_*/Validation_Donors.qs" \
    "Cartesian_Downsampling/*/Roussos/*/*/.ground_truth_complete"
  epa_complete
}
