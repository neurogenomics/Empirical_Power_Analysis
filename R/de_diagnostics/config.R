# ==============================================================================
# CONFIG: DE Diagnostics Pipeline - Shared Configuration
# ==============================================================================
# This file contains all shared settings for the DE diagnostics pipeline.
# Six steps:
#   Step 1: Density plots of DE result metrics across cell types
#   Step 2: GT-vs-Validation PCA diagnostics for cartesian downsampling splits
#   Step 3: Split DEG recovery — compare split GT vs full GT DE results
#   Step 4: DEG chromosome grid — counts by chromosome type × LFC magnitude
#   Step 5: Shared gene diagnostics — paired cross-CT comparisons for hypothesis
#   Step 6: Subtype expression divergence — per-gene split-violin figures (Micro/Oligo/Astro)
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
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "de_diagnostics")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. PATH CONFIGURATION
# ==============================================================================

# Input: Donor shards location
SHARD_BASE <- project_path("Datalink", "Shards")

# Input: Full ground truth DE results
FULL_GT_BASE <- project_path("Datalink", "Full_Ground_Truth")

# Input: Cartesian downsampling results (for GT/Val split metadata)
CARTESIAN_BASE <- project_path("Datalink", "Cartesian_Downsampling")

# Output: Diagnostic results
OUTPUT_BASE <- project_path("Datalink", "DE_Diagnostics")

# ==============================================================================
# 3. CELL TYPES
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ==============================================================================
# 4. STEP 1: DENSITY PLOT SETTINGS
# ==============================================================================

# Dataset to use for DE results (from full ground truth pipeline)
DE_DATASET <- "Roussos"

# Significance filter for density plots
PADJ_CUTOFF <- 0.05
LFC_MIN_HIST <- 2.0

# Metrics to plot (column name, file suffix, axis label, x-axis limits)
# lim = NULL means auto-scale from data (with DENSITY_PADDING on each side).
# lim = c(lower, upper) overrides with fixed bounds.
# lim_lower = value sets a fixed lower bound while auto-scaling the upper.
DENSITY_METRICS <- list(
  list(col = "Abs_LFC",           suffix = "LFC",        lab = "Abs(LFC)",           lim = NULL, lim_lower = 2),
  list(col = "lfcSE",             suffix = "SE",         lab = "Standard Error",     lim = NULL),
  list(col = "Log10_SNR",         suffix = "LogSNR",     lab = "Log10(SNR)",         lim = NULL),
  list(col = "Log10_BaseMean",    suffix = "BaseMean",   lab = "Log10(Mean Counts)", lim = NULL),
  list(col = "Log10_Dispersion",  suffix = "Dispersion", lab = "Log10(Dispersion)",  lim = NULL)
)

# Padding fraction for auto-scaled density plot x-axes (5% on each side)
DENSITY_PADDING <- 0.05

# BioMart settings
BIOMART_CACHE_FILE <- "gene_biotype_cache.csv"
BIOMART_DATASET <- "hsapiens_gene_ensembl"
BIOMART_HOST <- "https://useast.ensembl.org"
BIOMART_HOST_FALLBACK <- "https://www.ensembl.org"

# ==============================================================================
# 5. STEP 2: PCA SETTINGS
# ==============================================================================

# Shard dataset for PCA (filtered Roussos used by cartesian downsampling)
PCA_SHARD_DATASET <- "Roussos"

# Number of top variable genes for PCA
N_PCA_GENES <- 10000

# Minimum gene count filter for VST
MIN_GENE_COUNTS <- 10

# Number of robustness runs to process
N_ROBUSTNESS_RUNS <- 6

# Cartesian downsampling modes to process
CARTESIAN_MODES <- c("cells", "reads")

# PCA covariates to plot (each generates a separate plot per run)
# "Split" is always first — it gets density contours
PCA_COVARIATES <- c("Split", "sex", "Age_Years", "Source")

# Seed base for cartesian downsampling (must match cartesian config)
CARTESIAN_SEED_BASE <- 1234

# ==============================================================================
# 6. STEP 3: SPLIT DEG RECOVERY SETTINGS
# ==============================================================================

# Significance thresholds for DEG recovery comparison
RECOVERY_FDR_THRESHOLD <- 0.05
RECOVERY_LFC_THRESHOLD <- 1.0

# LFC magnitude bins for stratified recovery (matches cartesian downsampling)
RECOVERY_LFC_BINS <- list(
  "Small (1-2)"  = c(1.0, 2.0),
  "Medium (2-3)" = c(2.0, 3.0),
  "Large (3+)"   = c(3.0, Inf)
)

# Chromosome categories
CHR_AUTOSOMAL <- as.character(1:22)
CHR_SEX <- c("X", "Y")

# Recovery bar chart colours
RECOVERY_BAR_COLOURS <- c(
  "Full_GT" = "#2166AC",
  "Split_GT" = "#B2182B",
  "Overlap" = "#4DAF4A"
)

# UpSet plot settings
UPSET_NSETS <- 7       # Full GT + 6 splits
UPSET_ORDER_BY <- "freq"

# LFC scatter settings
SCATTER_ALPHA <- 0.4
SCATTER_POINT_SIZE <- 1.5
SCATTER_COLOURS <- c(
  "Both Significant" = "#4DAF4A",
  "Full GT Only"     = "#2166AC",
  "Split Only"       = "#B2182B",
  "Neither"          = "grey70"
)

# ==============================================================================
# 7. STEP 4: DEG CHROMOSOME GRID SETTINGS
# ==============================================================================

# Grid plot dimensions
CHR_GRID_WIDTH <- 8
CHR_GRID_HEIGHT <- 4
CHR_GRID_COMBINED_WIDTH <- 14

# Grid cell text size
CHR_GRID_TEXT_SIZE <- 6

# Heatmap fill gradient
CHR_GRID_FILL_LOW <- "#F7FBFF"
CHR_GRID_FILL_HIGH <- "#08519C"

# ==============================================================================
# 8. STEP 5: SHARED GENE DIAGNOSTICS SETTINGS
# ==============================================================================

# Path to shared truth genes (23 genes, intersection across all 3 CTs)
# Produced by power_curve_fitting Step 01
SHARED_GENES_PATH <- project_path("Datalink", "Power_Curve_Fitting",
                                   "Fitting_Data", "shared_truth_genes.qs")

# Significance thresholds for full truth set summary (matches density plots)
SHARED_GENE_PADJ_CUTOFF <- 0.05
SHARED_GENE_LFC_MIN <- 2.0

# Metrics for paired Wilcoxon tests
# Each: column name in the comparison table, display label, log-scale flag
SHARED_GENE_METRICS <- list(
  list(col = "dispersion",  lab = "Dispersion",  log_scale = TRUE),
  list(col = "abs_lfc",     lab = "|LFC|",       log_scale = FALSE),
  list(col = "lfcSE",       lab = "lfcSE",       log_scale = TRUE),
  list(col = "snr",         lab = "SNR",         log_scale = TRUE),
  list(col = "baseMean",    lab = "baseMean",    log_scale = TRUE)
)

# Metrics for Pearson correlation (subset of above)
SHARED_GENE_CORR_METRICS <- c("dispersion", "baseMean")

# Cell type pairs for pairwise comparisons
SHARED_GENE_CT_PAIRS <- list(
  c("Astro", "Micro"),
  c("Astro", "Oligo"),
  c("Micro", "Oligo")
)

# Scatter plot dimensions
SHARED_SCATTER_WIDTH <- 6
SHARED_SCATTER_HEIGHT <- 6
SHARED_SCATTER_DPI <- 150
SHARED_SCATTER_POINT_SIZE <- 2.5
SHARED_SCATTER_ALPHA <- 0.8

# ==============================================================================
# 8.5. STEP 6: SUBTYPE EXPRESSION DIVERGENCE SETTINGS
# ==============================================================================

# Number of random donors to sample per cell type
SUBTYPE_DIV_N_DONORS <- 500

# Seed for reproducible donor sampling
SUBTYPE_DIV_SEED <- 42

# Subtype column name in shard colData
SUBTYPE_DIV_SUBTYPE_COL <- "subtype"

# Cell types to analyse and their subtype definitions
# Each entry: list of subtypes (character vector) and short display labels (named vector)
# Micro has no subtypes — use NULL subtypes to indicate "pool all cells"
SUBTYPE_DIV_CT_CONFIG <- list(
  Oligo = list(
    subtypes = c("Oligo_OPALIN", "Oligo_RBFOX1"),
    labels   = c("Oligo_OPALIN" = "OPALIN", "Oligo_RBFOX1" = "RBFOX1"),
    test     = "wilcoxon"
  ),
  Astro = list(
    subtypes = c("Astro_ADAMTSL3", "Astro_GRIA1", "Astro_PLSCR1", "Astro_WIF1"),
    labels   = c("Astro_ADAMTSL3" = "ADAMTSL3", "Astro_GRIA1" = "GRIA1",
                  "Astro_PLSCR1" = "PLSCR1", "Astro_WIF1" = "WIF1"),
    test     = "kruskal"
  ),
  Micro = list(
    subtypes = NULL,
    labels   = c("All" = "All Micro"),
    test     = "none"
  )
)

# Cell type display order (for combined figure column ordering)
SUBTYPE_DIV_CT_ORDER <- c("Oligo", "Astro", "Micro")

# --- Per-gene split-violin figure layout ---

# X-axis positions: Micro | Oligo | Astro1 | Astro2
# Each position maps to a cell type and a set of subtypes for the split violin
SUBTYPE_DIV_XPOS <- list(
  list(
    label    = "Micro\n(All)",
    ct       = "Micro",
    subtypes = NULL,          # no split — single pooled violin
    type     = "pooled"
  ),
  list(
    label    = "Oligo\n(OPALIN | RBFOX1)",
    ct       = "Oligo",
    subtypes = c("Oligo_OPALIN", "Oligo_RBFOX1"),
    type     = "split"
  ),
  list(
    label    = "Astro\n(ADAMTSL3 | GRIA1)",
    ct       = "Astro",
    subtypes = c("Astro_ADAMTSL3", "Astro_GRIA1"),
    type     = "split"
  ),
  list(
    label    = "Astro\n(PLSCR1 | WIF1)",
    ct       = "Astro",
    subtypes = c("Astro_PLSCR1", "Astro_WIF1"),
    type     = "split"
  )
)

# Per-gene figure dimensions (inches)
SUBTYPE_DIV_FIG_WIDTH <- 8
SUBTYPE_DIV_FIG_HEIGHT <- 5
SUBTYPE_DIV_DPI <- 150

# Donor-level mean overlay points
SUBTYPE_DIV_POINT_SIZE <- 0.8
SUBTYPE_DIV_POINT_ALPHA <- 0.3

# Violin aesthetics
SUBTYPE_DIV_VIOLIN_ALPHA <- 0.7
SUBTYPE_DIV_VIOLIN_SCALE <- "width"

# Annotation text sizes
SUBTYPE_DIV_ANNOT_SIZE <- 2.5
SUBTYPE_DIV_PVAL_SIZE <- 3.0

# Parallel processing cores for shard loading
SUBTYPE_DIV_N_CORES <- 40

# Gene symbol cache path (for human-readable filenames)
SUBTYPE_DIV_SYMBOL_CACHE <- project_path("Datalink", "DE_Diagnostics",
                                          "gene_symbol_cache.csv")

# ==============================================================================
# 9. METADATA COLUMN NAMES
# ==============================================================================

DONOR_ID_COL <- "donor_id"
SEX_COL <- "sex"
AGE_COL <- "development_stage"
SOURCE_COL <- "Source"

# ==============================================================================
# 10. PLOT SETTINGS
# ==============================================================================

# Plot dimensions
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 6
PLOT_WIDTH_STRATIFIED <- 10
PLOT_HEIGHT_STRATIFIED <- 8
PLOT_DPI <- 150

# PCA plot dimensions
PCA_PLOT_WIDTH <- 7
PCA_PLOT_HEIGHT <- 7

# Colour palettes
DENSITY_PALETTE <- "Set1"
SPLIT_COLOURS <- c("Ground_Truth" = "#E69F00", "Validation" = "#56B4E9")
SOURCE_PALETTE <- "Set1"

# ==============================================================================
# 11. PARALLEL PROCESSING
# ==============================================================================

if (!exists("N_CORES")) {
  N_CORES <- 12
}

# ==============================================================================
# 12. CONTROL FLAGS
# ==============================================================================

FORCE_RECALC <- TRUE

# ==============================================================================
# 13. HELPER FUNCTIONS: Paths
# ==============================================================================

#' Get output paths for density plots
get_density_output_dir <- function() {
  file.path(OUTPUT_BASE, "Density_Plots", "Cross_CellType")
}

#' Get output paths for PCA plots
#' @param mode "cells" or "reads"
#' @param cell_type Cell type name
get_pca_output_dir <- function(mode, cell_type) {
  file.path(OUTPUT_BASE, "PCA", mode, cell_type)
}

#' Get path to full ground truth DE results CSV
#' @param cell_type Cell type name
#' @param dataset Dataset name
get_de_results_path <- function(cell_type, dataset = DE_DATASET) {
  file.path(FULL_GT_BASE, cell_type, paste0(dataset, "_DE_results.csv"))
}

#' Get path to full ground truth dds checkpoint
#' @param cell_type Cell type name
#' @param dataset Dataset name
get_dds_checkpoint_path <- function(cell_type, dataset = DE_DATASET) {
  file.path(FULL_GT_BASE, cell_type, paste0(dataset, "_dds_checkpoint.qs"))
}

#' Get path to cartesian downsampling GT metadata
#' @param cell_type Cell type name
#' @param mode "cells" or "reads"
#' @param run_id Run number
get_gt_metadata_path <- function(cell_type, mode, run_id) {
  file.path(CARTESIAN_BASE, cell_type, PCA_SHARD_DATASET, mode,
            paste0("Run_", run_id), "Ground_Truth_sampling_metadata.qs")
}

#' Get path to shard directory
#' @param cell_type Cell type name
get_shard_dir <- function(cell_type) {
  file.path(SHARD_BASE, cell_type, PCA_SHARD_DATASET)
}

#' Get path to cartesian downsampling split GT DE results
#' @param cell_type Cell type name
#' @param mode "cells" or "reads"
#' @param run_id Run number
get_cartesian_gt_path <- function(cell_type, mode, run_id) {
  file.path(CARTESIAN_BASE, cell_type, PCA_SHARD_DATASET, mode,
            paste0("Run_", run_id), "Ground_Truth.qs")
}

#' Get output directory for split recovery plots
#' @param cell_type Cell type name (optional — omit for top-level)
get_split_recovery_output_dir <- function(cell_type = NULL) {
  if (is.null(cell_type)) {
    file.path(OUTPUT_BASE, "Split_Recovery")
  } else {
    file.path(OUTPUT_BASE, "Split_Recovery", cell_type)
  }
}

#' Get output directory for chromosome grid plots
get_chr_grid_output_dir <- function() {
  file.path(OUTPUT_BASE, "Chromosome_Grid")
}

#' Get output directory for shared gene diagnostics
get_shared_gene_output_dir <- function() {
  file.path(OUTPUT_BASE, "Shared_Gene_Diagnostics")
}

#' Get output directory for subtype expression divergence
get_subtype_div_output_dir <- function() {
  file.path(OUTPUT_BASE, "Subtype_Expression_Divergence")
}

# ==============================================================================
# 14. LOGGING HELPERS (fallback if shared library not loaded)
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

if (!exists("log_warn")) {
  log_warn <- function(msg, ...) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] [WARNING] %s", timestamp, sprintf(msg, ...)))
  }
}

# ==============================================================================
# 15. VERSION INFO
# ==============================================================================

PIPELINE_VERSION <- "1.5.0"
PIPELINE_DATE <- "2026-02-23"
PIPELINE_NAME <- "de_diagnostics"

log_info("DE Diagnostics Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Cell Types: %s", paste(CELL_TYPES, collapse = ", "))
log_info("DE Dataset: %s | PCA Modes: %s | PCA Runs: %d",
         DE_DATASET, paste(CARTESIAN_MODES, collapse = "+"), N_ROBUSTNESS_RUNS)
