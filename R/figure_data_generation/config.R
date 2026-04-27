# ==============================================================================
# CONFIG: Figure Data Generation Pipeline v1.0.0
# ==============================================================================
# Generates all CSV source data tables needed by paper_final_draft.rmd.
# Replaces hardcoded data frames with computed values from actual pipeline
# outputs and shard files.
#
# Output: Datalink/Paper_Source_Data/*.csv
# Paper reads these directly via read.csv() instead of hardcoding values.
# ==============================================================================

# ==============================================================================
# 1. PROJECT ROOT
# ==============================================================================

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("The `here` package is required. Install with: install.packages(\"here\")")
  }
  here::here()
}
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "figure_data_generation")

project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. INPUT PATHS
# ==============================================================================

# Shard files (for dataset statistics)
SHARD_BASE <- project_path("Datalink", "Shards")

# Full ground truth results (for truth gene characteristics + FDP)
GT_BASE <- project_path("Datalink", "Full_Ground_Truth")

# Power curve fitting exports (for model coefficients)
PCF_BASE <- project_path("Datalink", "Power_Curve_Fitting")

# DE Diagnostics (for shared gene diagnostics)
DE_DIAG_BASE <- project_path("Datalink", "DE_Diagnostics")

# Glial characteristics (existing summary stats for Roussos/DeJager)
GLIAL_BASE <- project_path("Datalink", "Glial_Characteristics")

# ==============================================================================
# 3. OUTPUT PATH
# ==============================================================================

OUTPUT_BASE <- project_path("Datalink", "Paper_Source_Data")

# ==============================================================================
# 4. DATASET CONFIGURATION (all 5 datasets)
# ==============================================================================

DATASETS <- list(
  Roussos = list(
    label = "Roussos (PsychAD)",
    brain_region = "DLPFC",
    shard_dir_suffix = "Roussos",
    donor_id_col = "donor_id",
    sex_col = "sex",
    sex_encoding = c(male = "male", female = "female")
  ),
  DeJager = list(
    label = "DeJager (ROSMAP)",
    brain_region = "DLPFC",
    shard_dir_suffix = "DeJager",
    donor_id_col = "donor_id",
    sex_col = "msex",
    sex_encoding = c(male = "1", female = "0")
  ),
  Allen = list(
    label = "Allen (SEA-AD)",
    brain_region = "MTG",
    shard_dir_suffix = "Allen",
    donor_id_col = "donor_id",
    sex_col = "sex",
    sex_encoding = c(male = "male", female = "female")
  ),
  Johansen = list(
    label = "Johansen",
    brain_region = "Cerebral Cortex",
    shard_dir_suffix = "Johansen",
    donor_id_col = "donor_id",
    sex_col = "sex",
    sex_encoding = c(male = "male", female = "female")
  ),
  Rexach = list(
    label = "Rexach",
    brain_region = "Insular Cortex",
    shard_dir_suffix = "Rexach",
    donor_id_col = "donor_id",
    sex_col = "sex",
    sex_encoding = c(male = "male", female = "female")
  )
)

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ==============================================================================
# 5. PARALLELIZATION
# ==============================================================================

N_CORES <- 12

# ==============================================================================
# 6. CONTROL FLAGS
# ==============================================================================

# Force recalculation (ignore cached outputs).
# Override via environment variable POWER_ANALYSIS_FORCE_RECALC (used by run_all.sh).
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}

# ==============================================================================
# 7. PIPELINE VERSION
# ==============================================================================

CONFIG_VERSION <- "1.0.0"

message(sprintf("[CONFIG] Figure Data Generation Pipeline v%s", CONFIG_VERSION))
message(sprintf("[CONFIG] Output: %s", OUTPUT_BASE))
