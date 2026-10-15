# ==============================================================================
# CONFIG: Power Curve Fitting Pipeline v4.2.0 - Shared Configuration
# ==============================================================================
# This file contains all shared settings for the power curve fitting pipeline.
# Source this file at the start of each script.
#
# v4.2.0: Extended truth set to include Tiny LFC bin (0.2-1.0) by lowering
#          GT_LFC_THRESHOLD from 1.0 to 0.2. This captures small effect sizes
#          more representative of real-world sex-biased expression.
# v4.0.0: Switched from beta regression to logistic regression (glm, family=binomial)
#          because the response is binary per-iteration detection (0/1).
# v3.1.0: Unified model merging cells + reads modes via actual per-iteration
#          n_cells and reads_per_cell parsed from sampling metadata files.
# ==============================================================================

# ==============================================================================
# 1. PROJECT ROOT (HARDCODED FOR RELIABILITY)
# ==============================================================================

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("The `here` package is required. Install with: install.packages(\"here\")")
  }
  here::here()
}
# Script directory for this pipeline
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "power_curve_fitting")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. PATH CONFIGURATION
# ==============================================================================

# Input: Cartesian downsampling results
CARTESIAN_BASE <- project_path("Datalink", "Cartesian_Downsampling")

# Cartesian now stores DE results under {mode}/{method}/. The fitted power
# equation remains DESeq2-based unless this is overridden explicitly.
CARTESIAN_DE_METHOD <- if (nzchar(Sys.getenv("CARTESIAN_DE_METHOD"))) {
  Sys.getenv("CARTESIAN_DE_METHOD")
} else {
  "deseq2"
}
stopifnot(CARTESIAN_DE_METHOD %in% c("deseq2", "edger", "limma"))

# Input: Donor shards (for Step 0 GT DDS regeneration + Step 00b cache)
SHARD_BASE <- project_path("Datalink", "Shards")

# Input: BioMart gene annotation cache
BIOMART_CACHE <- project_path("Datalink", "Full_Ground_Truth", "gene_audit_biomart_cache.csv")

# Output: Power curve fitting results (Tiny LFC bin 0.2-1.0)
OUTPUT_BASE <- project_path("Datalink", "Power_Curve_Fitting")

# Dataset used in cartesian downsampling (legacy — kept for backward compatibility)
SHARD_DATASET <- "Roussos"

# ==============================================================================
# 2b. MULTI-DATASET CONFIGURATION
# ==============================================================================
# Each entry defines a dataset to process through the pipeline.
# - dataset:       Name of the shard subdirectory (e.g., "Roussos", "Roussos_DepthEq")
# - cell_types:    Which cell types to process for this dataset
# - model_names:   Named vector mapping cell_type -> model name in outputs.
#                  If NULL, cell type name is used as-is (e.g., "Astro" -> "Astro").
# - fit_combined:  Whether to fit combined/additive/interaction models across CTs.
#                  Requires 2+ cell types; set FALSE for single-CT datasets.

DATASET_CONFIGS <- list(
  list(
    dataset = "Roussos",
    cell_types = c("Astro", "Micro", "Oligo"),
    model_names = NULL,
    fit_combined = TRUE,
    special_mode = NULL
  ),
  list(
    dataset = "Roussos_DepthEq",
    cell_types = c("Astro"),
    model_names = c(Astro = "Astro_Shallow"),
    fit_combined = FALSE,
    special_mode = NULL
  ),
  list(
    dataset = "Roussos_ShallowFewer",
    cell_types = c("Astro"),
    model_names = c(Astro = "Astro_Shallow_Fewer"),
    fit_combined = FALSE,
    special_mode = "gt_only",  # GT-only mode: new GT, existing Val/iterations
    
    # For GT-only mode: GT comes from this dataset, but iteration results
    # come from the cartesian_source dataset
    cartesian_source = "Roussos_DepthEq",
    cartesian_target = "Roussos_ShallowFewer"
  )
)

# ==============================================================================
# 3. GRID CONFIGURATION (must match cartesian_downsampling)
# ==============================================================================

# Donor steps (number of donors to sample)
DONOR_STEPS <- c(5, 10, 15, 25, 50, 75, 100, 150, 200, 300, 500, 600)

# Cell steps (cells per donor) - used in cells mode
CELL_STEPS <- c(5, 10, 20, 50, 75, 100, 250, 500, 750, 1000, 1500, 2000, 3000, NA)

# Read steps (total reads per donor) - used in reads mode
READ_STEPS <- c(25000, 50000, 100000, 250000, 500000,
                1000000, 2500000, 5000000, 7500000,
                10000000, 12500000, 15000000, NA)

# Operating modes
MODES <- c("cells", "reads")

# ==============================================================================
# 4. ROBUSTNESS AND ITERATION SETTINGS (must match cartesian_downsampling)
# ==============================================================================

N_ROBUSTNESS_RUNS <- 10
N_ITERATIONS <- 20

# ==============================================================================
# 5. POWER CURVE FITTING SETTINGS
# ==============================================================================

# Ground truth significance thresholds
GT_PADJ_THRESHOLD <- 0.05
GT_LFC_THRESHOLD <- 0.2    # |LFC| > 0.2 defines the truth set (includes Tiny bin 0.2-1.0)

# Detection threshold (for subsample results)
DETECT_PADJ_THRESHOLD <- 0.05

# Minimum gene count filter (must match cartesian_downsampling)
MIN_GENE_COUNT <- 10

# [v4.0.0] Model family: logistic regression (glm, family=binomial)
# Replaced beta regression in v4.0.0 because the response is binary (0/1)
# per-iteration detection, not a continuous proportion. Logistic regression
# uses the same logit link and produces identical coefficient interpretation.
# The Smithson-Verkuilen squeeze (BETA_SQUEEZE) is no longer needed.

# Minimum number of valid observations per gene per grid point
# (across all iterations and runs) to include in fitting dataset
MIN_VALID_OBS <- 10

# ==============================================================================
# 6. v3.1.0 MODEL CONFIGURATION
# ==============================================================================
# v3.1.0 merges cells + reads modes by parsing sampling metadata to extract
# actual n_cells and mean reads_per_cell per iteration. The formula is:
#
#   detection ~ log_donors + log_n_cells + log_reads_per_cell +
#               log_abs_gt_lfc + log_gt_basemean + log_gt_dispersion
#
# 8 models total: (3 per-CT + 1 combined) x 2 gene sets x 1 aggregation level
#
# Gene sets:
#   - "full": all truth genes per cell type (different genes per CT)
#   - "shared": intersection of truth genes across all CTs (same genes)
#
# Aggregation level:
#   - "per_iteration": no aggregation across iterations
#                      → 1 row per gene x grid_point x run x iteration

GENE_SETS <- c("full", "shared")
AGGREGATION_LEVELS <- c("per_iteration")

# LRT significance threshold
LRT_ALPHA <- 0.05

# Confidence level for coefficient CI overlap check
CI_LEVEL <- 0.95

# ==============================================================================
# 6b. MEAN READS PER CELL (from EDA: Roussos dataset)
# ==============================================================================
# Source: Datalink/Exploratory_Analysis/Read_Depth_Summary_All_CellTypes_Roussos.csv
# These are used in Step 04 for cells<->reads cross-validation reference values.
# In v3.1.0, actual per-iteration reads_per_cell replaces these constants in models.

MEAN_READS_PER_CELL <- list(
  Astro = 9456,
  Micro = 4426,
  Oligo = 5742
)

# ==============================================================================
# 7. CELL TYPE CONFIGURATION
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ==============================================================================
# 8. METADATA COLUMN MAPPINGS (must match cartesian_downsampling)
# ==============================================================================

DONOR_ID_COL <- "donor_id"
SEX_COL <- "sex"
AGE_COL <- "development_stage"
MALE_VAL <- "male"
FEMALE_VAL <- "female"
CATEGORICAL_VARS <- c("Source", "self_reported_ethnicity", "disease")
# pmi (post-mortem interval, minutes) annotated onto Roussos shards by
# R/sharding/Annotate_Roussos_PMI.R. Scaled to pmi_scaled in run_deseq_for_dds().
NUMERIC_VARS <- c("pmi")

# Contrast override (sex by default; see R/shared/contrast.R).
source(file.path(PROJECT_ROOT, "R", "shared", "contrast.R"))
invisible(lapply(DATASET_CONFIGS, function(d) check_contrast_dataset(d$dataset)))

# ==============================================================================
# 9. DESEQ2 SETTINGS (for Step 0 GT DDS regeneration)
# ==============================================================================

FDR_THRESHOLD <- 0.05
SEED_BASE <- 1234   # Must match cartesian_downsampling

# ==============================================================================
# 10. PARALLEL PROCESSING
# ==============================================================================

if (!exists("N_CORES")) {
  N_CORES <- 24
}

N_CORES_HEAVY <- min(6, N_CORES)

setup_parallel <- function(n_cores = N_CORES, heavy = FALSE) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel is missing from the pinned R library. Restore renv.lock or rebuild the Docker image.")
  }
  cores_to_use <- if (heavy) min(n_cores, N_CORES_HEAVY) else n_cores
  BiocParallel::MulticoreParam(
    workers = cores_to_use,
    progressbar = TRUE,
    stop.on.error = FALSE
  )
}

suppressPackageStartupMessages({
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BPPARAM <- setup_parallel(N_CORES, heavy = FALSE)
    BPPARAM_HEAVY <- setup_parallel(N_CORES, heavy = TRUE)
  }
})

# ==============================================================================
# 11. VISUALISATION SETTINGS
# ==============================================================================

# Representative gene covariate values for prediction plots
# (will be computed from data in Step 03, these are fallback defaults)
DEFAULT_REPRESENTATIVE_LFC <- 2.0       # |LFC| for "typical truth gene"
DEFAULT_REPRESENTATIVE_BASEMEAN <- 100   # baseMean for "typical expressed gene"
DEFAULT_REPRESENTATIVE_DISPERSION <- 0.5 # dispersion for "typical gene"

# Power target levels for contour lines
POWER_TARGETS <- c(0.5, 0.8, 0.9, 0.95)

# Plot dimensions
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# Colour palette for cell types (including sensitivity analysis variants)
CELL_TYPE_COLOURS <- c(
  Astro = "#440154",              # viridis 1/5 (purple)
  Micro = "#3B528B",              # viridis 2/5 (blue)
  Oligo = "#21908C",              # viridis 3/5 (teal)
  Astro_Shallow = "#5DC863",      # viridis 4/5 (green) -- depth-equalized
  Astro_Shallow_Fewer = "#FDE725" # viridis 5/5 (yellow) -- depth + count matched
)

# ==============================================================================
# 12. LOOKUP TABLE SETTINGS (for Step 05 export)
# ==============================================================================

# Grid for lookup table generation (v3.1.0: uses n_cells + reads_per_cell)
LOOKUP_DONOR_STEPS <- c(10, 25, 50, 100, 200, 500)
LOOKUP_CELL_STEPS <- c(50, 100, 250, 500, 1000)
LOOKUP_RPC_STEPS <- c(2000, 4000, 6000, 8000, 10000)
LOOKUP_LFC_VALUES <- c(1.5, 2.0, 3.0, 4.0, 5.0)
LOOKUP_BASEMEAN_VALUES <- c(10, 50, 100, 500, 1000)
LOOKUP_DISPERSION_VALUES <- c(0.1, 0.5, 1.0)

# ==============================================================================
# 12b. MINI RUN OVERRIDES (smoke test — activate with env var MINI_RUN=true)
# ==============================================================================
# Mirrors the cartesian_downsampling MINI_RUN switch so the power-curve-fitting
# stage of a mini smoke run reads the mini cartesian outputs, re-derives ground
# truth from the same capped donor universe, and writes to an isolated MINI_-
# prefixed directory. See R/cartesian_downsampling/config.R for the rationale.
# All grids must match the cartesian MINI grid for the two stages to line up.
MINI_RUN <- tolower(Sys.getenv("MINI_RUN", "false")) %in% c("true", "1", "yes")

MINI_CELL_TYPE <- if (nzchar(Sys.getenv("MINI_CELL_TYPE"))) {
  Sys.getenv("MINI_CELL_TYPE")
} else {
  "Astro"
}

MINI_MAX_DONORS <- if (nzchar(Sys.getenv("MINI_MAX_DONORS"))) {
  as.integer(Sys.getenv("MINI_MAX_DONORS"))
} else {
  20L
}

if (MINI_RUN) {
  # Read mini cartesian outputs; write isolated PCF outputs.
  CARTESIAN_BASE    <- project_path("Datalink", "MINI_Cartesian_Downsampling")
  OUTPUT_BASE       <- project_path("Datalink", "MINI_Power_Curve_Fitting")
  # Single cell type (Astro) across three datasets: the native run plus the two
  # depth-equalisation variants, which live under MINI_-prefixed shard/dataset
  # dirs (see depth_equalization/config.R) so production data is never touched.
  # No combined model is fitted (combined needs 2+ cell types).
  CELL_TYPES        <- MINI_CELL_TYPE
  DATASET_CONFIGS   <- list(
    list(
      dataset      = "Roussos",
      cell_types   = MINI_CELL_TYPE,
      model_names  = NULL,
      fit_combined = FALSE,
      special_mode = NULL
    ),
    list(
      dataset      = "MINI_Roussos_DepthEq",
      cell_types   = MINI_CELL_TYPE,
      model_names  = c(Astro = "Astro_Shallow"),
      fit_combined = FALSE,
      special_mode = NULL
    ),
    list(
      dataset          = "MINI_Roussos_ShallowFewer",
      cell_types       = MINI_CELL_TYPE,
      model_names      = c(Astro = "Astro_Shallow_Fewer"),
      fit_combined     = FALSE,
      special_mode     = "gt_only",
      cartesian_source = "MINI_Roussos_DepthEq",
      cartesian_target = "MINI_Roussos_ShallowFewer"
    )
  )
  DONOR_STEPS       <- c(5, 10)
  CELL_STEPS        <- c(50, 100, NA)
  READ_STEPS        <- c(100000, 250000, NA)
  N_ROBUSTNESS_RUNS <- 1
  N_ITERATIONS      <- 2
  GENE_SETS         <- c("full")  # "shared" == "full" with a single cell type
  message(sprintf("[MINI_RUN] power_curve_fitting: cell type '%s', %d datasets, isolated outputs -> %s",
                  MINI_CELL_TYPE, length(DATASET_CONFIGS), OUTPUT_BASE))
}

# ==============================================================================
# 13. PATH HELPERS
# ==============================================================================

#' Get paths for a cell type (mode-agnostic in v3.1.0)
#'
#' Legacy function — uses SHARD_DATASET global. For multi-dataset support,
#' prefer get_paths_for_dataset() which accepts an explicit dataset name.
#'
#' @param cell_type Cell type name
#' @param mode "cells" or "reads" (for cartesian input paths only)
#' @return Named list of paths
get_paths <- function(cell_type, mode = "cells") {
  list(
    # Input: cartesian downsampling results
    cartesian_dir = file.path(CARTESIAN_BASE, cell_type, SHARD_DATASET, mode, CARTESIAN_DE_METHOD),
    shard_dir = file.path(SHARD_BASE, cell_type, SHARD_DATASET),
    manifest_path = file.path(SHARD_BASE, cell_type, SHARD_DATASET, "00_MASTER_MANIFEST.qs"),

    # Output: power curve fitting results (mode-agnostic in v3.1.0)
    gt_dds_dir = file.path(OUTPUT_BASE, "GT_DDS_Checkpoints", cell_type),
    shard_cache_dir = file.path(OUTPUT_BASE, "Shard_Cache", cell_type),
    fitting_dir = file.path(OUTPUT_BASE, "Fitting_Data", cell_type),
    model_dir = file.path(OUTPUT_BASE, "Models"),
    plots_dir = file.path(OUTPUT_BASE, "Plots"),
    cross_val_dir = file.path(OUTPUT_BASE, "Plots", "Cross_Validation"),
    export_dir = file.path(OUTPUT_BASE, "Exports"),
    summary_dir = file.path(OUTPUT_BASE, "Summary")
  )
}

#' Get the model name for a cell type within a dataset config
#'
#' @param cell_type Cell type name (e.g., "Astro")
#' @param ds_config A single DATASET_CONFIGS entry (list)
#' @return Model name string (e.g., "Astro" or "Astro_Shallow")
get_dataset_model_name <- function(cell_type, ds_config) {
  if (!is.null(ds_config$model_names) && cell_type %in% names(ds_config$model_names)) {
    return(ds_config$model_names[[cell_type]])
  }
  return(cell_type)
}

#' Get paths for a cell type within a specific dataset
#'
#' Multi-dataset-aware version of get_paths(). Uses the dataset name for
#' input paths (shards, cartesian results) and the model name for output
#' paths (shard cache, GT DDS, fitting data, models).
#'
#' @param cell_type Cell type name (e.g., "Astro")
#' @param ds_config A single DATASET_CONFIGS entry (list)
#' @param mode "cells" or "reads" (for cartesian input paths only)
#' @return Named list of paths
get_paths_for_dataset <- function(cell_type, ds_config, mode = "cells") {
  dataset <- ds_config$dataset
  model_name <- get_dataset_model_name(cell_type, ds_config)

  list(
    # Input: cartesian downsampling results (keyed by dataset name)
    cartesian_dir = file.path(CARTESIAN_BASE, cell_type, dataset, mode, CARTESIAN_DE_METHOD),
    shard_dir = file.path(SHARD_BASE, cell_type, dataset),
    manifest_path = file.path(SHARD_BASE, cell_type, dataset, "00_MASTER_MANIFEST.qs"),

    # Output: power curve fitting results (keyed by model name)
    gt_dds_dir = file.path(OUTPUT_BASE, "GT_DDS_Checkpoints", model_name),
    shard_cache_dir = file.path(OUTPUT_BASE, "Shard_Cache", model_name),
    fitting_dir = file.path(OUTPUT_BASE, "Fitting_Data", model_name),
    model_dir = file.path(OUTPUT_BASE, "Models"),
    plots_dir = file.path(OUTPUT_BASE, "Plots"),
    cross_val_dir = file.path(OUTPUT_BASE, "Plots", "Cross_Validation"),
    export_dir = file.path(OUTPUT_BASE, "Exports"),
    summary_dir = file.path(OUTPUT_BASE, "Summary"),

    # Extra metadata for downstream use
    dataset = dataset,
    model_name = model_name,
    cell_type = cell_type
  )
}

#' Get model output directory for a v3.1.0 model variant
#'
#' @param gene_set "full" or "shared"
#' @param agg_level "per_iteration"
#' @return Path to model directory
get_model_dir <- function(gene_set = "full", agg_level = "per_iteration") {
  file.path(OUTPUT_BASE, "Models", paste0(gene_set, "_", agg_level))
}

#' Get the file pattern regex for a mode's metadata files
#'
#' @param mode "cells" or "reads"
#' @return List with pattern and parse_regex
get_metadata_file_pattern <- function(mode) {
  if (mode == "reads") {
    list(
      pattern = "^N\\d+_R.*_I\\d+_sampling_metadata\\.qs$",
      parse_regex = "N(\\d+)_R([0-9\\.e\\+]+|NA)_I(\\d+)_sampling_metadata"
    )
  } else {
    list(
      pattern = "^N\\d+_C.*_I\\d+_sampling_metadata\\.qs$",
      parse_regex = "N(\\d+)_C(\\d+|NA)_I(\\d+)_sampling_metadata"
    )
  }
}

#' Get the file pattern regex for a mode's DE result files
#'
#' @param mode "cells" or "reads"
#' @return List with pattern and parse_regex
get_file_pattern <- function(mode) {
  if (mode == "reads") {
    list(
      pattern = "^N\\d+_R.*_I\\d+\\.qs$",
      parse_regex = "N(\\d+)_R([0-9\\.e\\+]+|NA)_I(\\d+)"
    )
  } else {
    list(
      pattern = "^N\\d+_C.*_I\\d+\\.qs$",
      parse_regex = "N(\\d+)_C(\\d+|NA)_I(\\d+)"
    )
  }
}

#' Get the Y-axis column name for a mode
#'
#' @param mode "cells" or "reads"
#' @return Column name
get_y_col_name <- function(mode) {
  if (mode == "reads") "Target_Reads" else "Cells"
}

#' Get the Y-axis steps for a mode
#'
#' @param mode "cells" or "reads"
#' @return Vector of steps
get_y_steps <- function(mode) {
  if (mode == "reads") READ_STEPS else CELL_STEPS
}

# ==============================================================================
# 14. CONTROL FLAGS
# ==============================================================================

# Force recalculation (ignore cached fits and outputs).
# Override via environment variable POWER_ANALYSIS_FORCE_RECALC (used by run_all.sh).
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}

# ==============================================================================
# 15. LOGGING HELPERS (fallback if shared library not loaded)
# ==============================================================================

if (!exists("log_info")) {
  log_info <- function(msg, ...) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s", timestamp, sprintf(msg, ...)))
  }
}

if (!exists("log_section")) {
  log_section <- function(title) {
    cat("\n")
    log_info(paste(rep("=", 60), collapse = ""))
    log_info(title)
    log_info(paste(rep("=", 60), collapse = ""))
  }
}

if (!exists("log_error")) {
  log_error <- function(msg, ...) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] [ERROR] %s", timestamp, sprintf(msg, ...)))
  }
}

if (!exists("log_warning")) {
  log_warning <- function(msg, ...) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] [WARNING] %s", timestamp, sprintf(msg, ...)))
  }
}

# ==============================================================================
# 16. VERSION INFO
# ==============================================================================

PIPELINE_VERSION <- "4.2.0"
PIPELINE_DATE <- "2026-02-26"

log_info("Power Curve Fitting Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Cell types: %s | Cores: %d", paste(CELL_TYPES, collapse = ", "), N_CORES)
log_info("GT thresholds: padj < %.2f, |LFC| > %.1f | Detection: padj < %.2f",
         GT_PADJ_THRESHOLD, GT_LFC_THRESHOLD, DETECT_PADJ_THRESHOLD)
log_info("Gene sets: %s | Aggregation: %s",
         paste(GENE_SETS, collapse = ", "), paste(AGGREGATION_LEVELS, collapse = ", "))
log_info("Datasets: %d configured", length(DATASET_CONFIGS))
for (.ds_cfg in DATASET_CONFIGS) {
  .ds_model_names <- sapply(.ds_cfg$cell_types, function(ct) get_dataset_model_name(ct, .ds_cfg))
  log_info("  %s: %s (combined=%s)",
           .ds_cfg$dataset,
           paste(sprintf("%s->%s", .ds_cfg$cell_types, .ds_model_names), collapse = ", "),
           .ds_cfg$fit_combined)
}
rm(.ds_cfg, .ds_model_names)
