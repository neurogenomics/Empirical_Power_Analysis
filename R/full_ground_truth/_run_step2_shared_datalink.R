# Wrapper: run step 2 (Cross-Dataset Diagnostics) using the shared Datalink
# where all DE results actually live.
# Usage: Rscript _run_step2_shared_datalink.R

# Set env vars BEFORE config.R is sourced by the step script
Sys.setenv(OUTPUT_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/Full_Ground_Truth")
Sys.setenv(CARTESIAN_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/Cartesian_Downsampling")

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  here::here()
}

cat("PROJECT_ROOT:", PROJECT_ROOT, "\n")
cat("OUTPUT_BASE will be:", Sys.getenv("OUTPUT_BASE_OVERRIDE"), "\n")
cat("CARTESIAN_BASE will be:", Sys.getenv("CARTESIAN_BASE_OVERRIDE"), "\n")

source(file.path(PROJECT_ROOT, "R", "full_ground_truth", "02_Cross_Dataset_Diagnostics.R"))
