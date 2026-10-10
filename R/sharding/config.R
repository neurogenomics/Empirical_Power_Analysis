# ==============================================================================
# CONFIG: Sharding Pipeline — Shared Configuration
# ==============================================================================
# Source this file at the top of every Run_Shard_*.R script.
#
# Resolves PROJECT_ROOT (same convention as all other pipelines):
#   1. $POWER_ANALYSIS_ROOT environment variable, if set.
#   2. Otherwise via here::here() from the nearest .here / .Rproj / DESCRIPTION.
# ==============================================================================

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("The `here` package is required. Install with: install.packages(\"here\")")
  }
  here::here()
}

SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "sharding")

project_path <- function(...) file.path(PROJECT_ROOT, ...)

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================

# Input: raw single-cell files (not distributed; see Datalink/raw_data_README.md)
RAW_DATA_DIR <- project_path("Datalink", "raw_data")

# Output: donor-level SingleCellExperiment shards
SHARD_BASE   <- project_path("Datalink", "Shards")

# ==============================================================================
# SHARED FILTERS (applied to every dataset unless overridden per script)
# ==============================================================================

MIN_CELLS_PER_DONOR <- 30L
MIN_AGE_YEARS       <- 40L

# Cell-count feasibility levels used by the manifest
CELL_STEPS <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA)
