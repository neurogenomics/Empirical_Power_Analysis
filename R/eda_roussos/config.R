# ==============================================================================
# CONFIG: EDA Roussos Pipeline - Shared Configuration
# ==============================================================================
# This file contains all shared settings for the exploratory data analysis
# pipeline for Roussos dataset. Source this file at the start of each script.
# ==============================================================================

# ==============================================================================
# 1. PROJECT ROOT (HARDCODED FOR RELIABILITY)
# ==============================================================================
# This ensures scripts work when sourced in RStudio, not just via Rscript.
# DO NOT use commandArgs() or sys.frame() for path detection.

PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("The `here` package is required. Install with: install.packages(\"here\")")
  }
  here::here()
}
# Script directory for this pipeline
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "eda_roussos")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. PATH CONFIGURATION
# ==============================================================================

# Input: Donor shards location
SHARD_BASE <- project_path("Datalink", "Shards")

# Datasets to analyze (both unfiltered and filtered Roussos shards)
# NOTE: "Roussos" = filtered (age >= 40, cells >= 30) — used by power analysis pipelines
#       "Roussos_Full" = unfiltered — complete EDA picture
SHARD_DATASETS <- c("Roussos_Full", "Roussos")

# Backward compatibility: default single dataset for any code still using SHARD_DATASET
SHARD_DATASET <- SHARD_DATASETS[1]

# Human-readable labels for plot titles, keyed by dataset name
COHORT_LABELS <- c(
  "Roussos_Full" = "Full Cohort",
  "Roussos"      = "Filtered Cohort"
)

# Output: Analysis results location
OUTPUT_BASE <- project_path("Datalink", "Exploratory_Analysis")

# Full Ground Truth: DE results and pseudobulk checkpoints
FULL_GT_BASE <- project_path("Datalink", "Full_Ground_Truth")

# Helper function to get output paths for a specific cell type and dataset
get_output_paths <- function(cell_type, dataset = SHARD_DATASET) {
  base_dir <- file.path(OUTPUT_BASE, cell_type, dataset)
  
  list(
    base = base_dir,
    pca = file.path(base_dir, "PCA"),
    histograms = file.path(base_dir, "Histograms"),
    qc = file.path(base_dir, "QC"),
    tables = file.path(base_dir, "Tables")
  )
}

# Paper skeleton location
PAPER_SKELETON <- project_path("docs", "paper_skeleton.Rmd")

# ==============================================================================
# 3. CELL TYPES TO ANALYZE
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ==============================================================================
# 4. METADATA COLUMN NAMES
# ==============================================================================

# Core columns (required)
DONOR_ID_COL <- "donor_id"
SEX_COL <- "sex"

# Optional columns (will check existence before using)
AGE_COL <- "development_stage"
SOURCE_COL <- "Source"
DISEASE_COL <- "disease"
ETHNICITY_COL <- "self_reported_ethnicity"
# NOTE: Cell subtype deliberately excluded from PCA covariates. Subtypes are
# defined by transcriptomic signature (collinear with the expression data PCA
# operates on), exist on a continuous spectrum rather than as discrete categories
# (Trapnell 2015; Tanay & Regev 2017), and their inclusion as covariates would
# create a circular dependency (double-dipping). Additionally, the pseudobulk
# aggregation assigns each donor the subtype of its first cell — an arbitrary
# and uninformative summary of a mixed-subtype donor.

# All metadata columns to check for PCA
# NOTE: "Age_Years" is a numeric column created by parse_age_stage_to_years()
#       in 01_PCA_Diagnostics.R. Using this instead of AGE_COL ("development_stage")
#       ensures PCA plots color age on a continuous viridis scale, not as a
#       categorical variable with 40+ discrete levels.
PCA_METADATA_COLS <- c(
  SEX_COL,
  "Age_Years",
  SOURCE_COL,
  DISEASE_COL,
  ETHNICITY_COL,
  "Log2_Cells"  # Will be calculated from n_cells
)

# ==============================================================================
# 5. PCA PARAMETERS
# ==============================================================================

# Number of top variable genes to use for PCA
N_PCA_GENES <- 10000

# Minimum gene filter (genes with rowSums < this are removed)
MIN_GENE_COUNTS <- 10

# ==============================================================================
# 6. HISTOGRAM PARAMETERS
# ==============================================================================

# Number of bins for histograms
HISTOGRAM_BINS <- 100

# ==============================================================================
# 7. GENE EXPRESSION PARAMETERS
# ==============================================================================

# Threshold for zero-inflated genes (proportion of zeros)
ZERO_INFLATED_THRESHOLD <- 0.90

# Number of top expressed genes to report
TOP_GENES_N <- 50

# ==============================================================================
# 8. QC PARAMETERS
# ==============================================================================

# Mitochondrial gene patterns (regex)
MT_PATTERN <- "^MT-"

# Ribosomal gene patterns (regex)
RIBO_PATTERN <- "^RP[SL]"

# ==============================================================================
# 9. PARALLELIZATION
# ==============================================================================

# Number of cores for parallel processing
N_CORES <- 12

# ==============================================================================
# 10. PIPELINE STEP SELECTION
# ==============================================================================

# Which steps to run by default (override via command-line args)
# Steps: 1=PCA, 3=Read Depth, 6=DE Filtering Summary
STEPS_TO_RUN <- c(1, 3, 6)

# ==============================================================================
# 11. FORCE RECALCULATION
# ==============================================================================

# If FALSE, skip steps whose output files already exist for a given dataset.
# Set to TRUE to regenerate all outputs unconditionally.
FORCE_RECALC <- TRUE

# ==============================================================================
# 12. DE FILTERING SUMMARY PARAMETERS (Step 06)
# ==============================================================================

# Dataset name for DE results (must match Full Ground Truth pipeline)
DE_DATASET <- "Roussos"

# LFC magnitude bins for filtering summary (|LFC| thresholds)
DE_FILTER_LFC_BINS <- c(1, 2, 3, Inf)
DE_FILTER_LFC_LABELS <- c("1-2", "2-3", "3+")

# Filter category labels and colours
DE_FILTER_CATEGORIES <- c(
  "Low-Count Pre-Filter",
  "Cook's Distance Outlier",
  "Independent Filtering",
  "Fully Tested"
)

DE_FILTER_COLOURS <- c(
  "Low-Count Pre-Filter"       = "#999999",  # grey
  "Cook's Distance Outlier"    = "#E41A1C",  # red
  "Independent Filtering"      = "#FF7F00",  # orange
  "Fully Tested"               = "#4DAF4A"   # green
)

# ==============================================================================
# 13. PLOTTING PARAMETERS
# ==============================================================================

# Default plot dimensions (inches)
PLOT_WIDTH <- 7
PLOT_HEIGHT <- 6

# Color palettes
COLOR_PALETTE_CATEGORICAL <- "Set1"  # For discrete variables
COLOR_PALETTE_CONTINUOUS <- "viridis"  # For continuous variables

# ==============================================================================
# 14. LOGGING
# ==============================================================================

# Log file location
LOG_DIR <- project_path("Datalink", "Run_Output_Logs")

# Pipeline name for logging
PIPELINE_NAME <- "eda_roussos"

# ==============================================================================
# 15. VERSION
# ==============================================================================

CONFIG_VERSION <- "1.4.0"
LAST_UPDATED <- "2026-02-20"

message(sprintf("[CONFIG] EDA Roussos Pipeline v%s (Updated: %s)", 
                CONFIG_VERSION, LAST_UPDATED))
