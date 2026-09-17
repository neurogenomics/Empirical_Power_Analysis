# ==============================================================================
# SCRIPT: Build HPC Job Manifest
# PURPOSE: Create cost-balanced PBS array chunks for the downsampling grid
# INPUT: Validation donor splits from Step 1
# OUTPUT: CSV manifest consumed by Step 2 when JOB_CHUNK is set
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
})

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",
    "config.R",
    "../config.R"
  )

  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }

  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 1.5: Build HPC Job Manifest")
log_info("Manifest path: %s", JOB_MANIFEST_PATH)
log_info("Chunks: %d", JOB_N_CHUNKS)

build_job_manifest(
  n_chunks = JOB_N_CHUNKS,
  manifest_path = JOB_MANIFEST_PATH
)

log_section("STEP 1.5 COMPLETE")
log_info("HPC job manifest ready: %s", JOB_MANIFEST_PATH)
