# ==============================================================================
# CONFIG: Glial Characteristics Cross-Dataset Comparison Pipeline
# ==============================================================================
# This file contains all shared settings for comparing cell-level and donor-level
# characteristics between Roussos and DeJager datasets (both DLPFC) across
# the three glial cell types.
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
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "glial_characteristics")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. PATH CONFIGURATION
# ==============================================================================

# Input: Donor shards location
SHARD_BASE <- project_path("Datalink", "Shards")

# Output: Analysis results location
OUTPUT_BASE <- project_path("Datalink", "Glial_Characteristics")

# Output subdirectories
OUTPUT_METRICS <- file.path(OUTPUT_BASE, "Metrics")
OUTPUT_PLOTS <- file.path(OUTPUT_BASE, "Plots")
OUTPUT_TESTS <- file.path(OUTPUT_BASE, "Tests")
OUTPUT_SUMMARY <- file.path(OUTPUT_BASE, "Summary")

# Existing EDA outputs (for reusing Roussos data)
EDA_BASE <- project_path("Datalink", "Exploratory_Analysis")

# ==============================================================================
# 3. DATASET CONFIGURATION
# ==============================================================================

# Each dataset has specific metadata column names and encodings
DATASETS <- list(
  Roussos = list(
    shard_dir_suffix = "Roussos",
    donor_id_col = "donor_id",
    sex_col = "sex",
    sex_encoding = c(male = "male", female = "female"),
    label = "Roussos (DLPFC)",
    short_label = "Roussos",
    brain_region = "DLPFC"
  ),
  DeJager = list(
    shard_dir_suffix = "DeJager",
    donor_id_col = "donor_id",  # In manifest; shards use individualID internally
    sex_col = "msex",
    sex_encoding = c(male = "1", female = "0"),
    label = "De Jager (DLPFC)",
    short_label = "DeJager",
    brain_region = "DLPFC"
  )
)

# Dataset display order (for consistent plotting)
DATASET_ORDER <- c("Roussos", "DeJager")

# Dataset colours for plotting
DATASET_COLOURS <- c(
  Roussos = "#4E79A7",   # Steel blue
  DeJager = "#E15759"    # Salmon
)

# ==============================================================================
# 4. CELL TYPES
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# Reference cell type for computing ratios (denominator)
REFERENCE_CELL_TYPE <- "Micro"

# Cell type colours (from shared/prettify_labels.R)
# These will be loaded from prettify_labels.R

# ==============================================================================
# 5. METRICS TO COMPUTE
# ==============================================================================

# Metrics to extract and compare
METRICS <- c(
  "umis_per_cell",       # Library size: colSums(counts)
  "cells_per_donor",     # Cell abundance: ncol(sce) or manifest cell_count
  "genes_per_cell"       # Complexity: colSums(counts > 0)
)

# Metric labels for plots
METRIC_LABELS <- c(
  umis_per_cell = "UMIs per Cell",
  cells_per_donor = "Cells per Donor",
  genes_per_cell = "Genes Detected per Cell"
)

# Metric short labels for tables
METRIC_SHORT <- c(
  umis_per_cell = "UMIs/cell",
  cells_per_donor = "Cells/donor",
  genes_per_cell = "Genes/cell"
)

# ==============================================================================
# 6. BOOTSTRAP PARAMETERS
# ==============================================================================

# Number of bootstrap resamples for ratio CIs
N_BOOTSTRAP <- 10000

# Bootstrap seed for reproducibility
BOOTSTRAP_SEED <- 42

# Confidence level for bootstrap CIs
CI_LEVEL <- 0.95

# ==============================================================================
# 7. PARALLELIZATION
# ==============================================================================

# Number of cores for parallel processing
N_CORES <- 12

# ==============================================================================
# 8. SAMPLING PARAMETERS
# ==============================================================================
# For memory-efficient plotting, subsample cells when creating distributions

# Maximum cells per cell type for histogram/violin plots
MAX_CELLS_PER_CT <- 100000

# Random seed for subsampling
SUBSAMPLE_SEED <- 42

# ==============================================================================
# 9. PLOTTING PARAMETERS
# ==============================================================================

# Default plot dimensions (inches)
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 6

# Number of bins for histograms
HISTOGRAM_BINS <- 100

# Base font size for publication-quality plots
BASE_FONT_SIZE <- 16

# DPI for saved plots
PLOT_DPI <- 300

# ==============================================================================
# 10. PIPELINE STEP SELECTION
# ==============================================================================

# Which steps to run (1=Extract Metrics, 2=Comparison Plots, 3=Bootstrap Tests)
STEPS_TO_RUN <- 1:3

# ==============================================================================
# 11. FORCE RECALCULATION
# ==============================================================================

# If FALSE, skip steps whose output files already exist
# Set to TRUE to regenerate all outputs unconditionally
FORCE_RECALC <- TRUE

# ==============================================================================
# 12. LOGGING
# ==============================================================================

# Log file location
LOG_DIR <- project_path("Datalink", "Run_Output_Logs")

# Pipeline name for logging
PIPELINE_NAME <- "glial_characteristics"

# ==============================================================================
# 13. VERSION
# ==============================================================================

CONFIG_VERSION <- "1.0.0"
LAST_UPDATED <- "2026-02-27"

message(sprintf("[CONFIG] Glial Characteristics Pipeline v%s (Updated: %s)", 
                CONFIG_VERSION, LAST_UPDATED))
