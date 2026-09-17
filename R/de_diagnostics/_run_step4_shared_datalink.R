# Wrapper: run step 4 (Chromosome Grid) using the shared Datalink
# Usage: Rscript _run_step4_shared_datalink.R

Sys.setenv(FULL_GT_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/Full_Ground_Truth")
Sys.setenv(CARTESIAN_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/Cartesian_Downsampling")
Sys.setenv(SHARD_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/Shards")
Sys.setenv(DE_DIAG_OUTPUT_BASE_OVERRIDE = "/mnt/data/shared/poweranalysis/DE_Diagnostics")

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  here::here()
}

# Force recalculation since we're regenerating with new labels
FORCE_RECALC <- TRUE

cat("Running chromosome grid with shared Datalink paths\n")
source(file.path(PROJECT_ROOT, "R", "de_diagnostics", "04_DEG_Chromosome_Grid.R"))
