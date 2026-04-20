# ==============================================================================
# BULK DEG PROCESSING PIPELINE - CONFIGURATION
# ==============================================================================
# Version: 1.0.0
# Last updated: 2026-02-10
#
# This file contains all paths and parameters for the bulk DEG processing
# pipeline. It uses hardcoded PROJECT_ROOT for reliability across execution
# contexts (RStudio Source, Background Jobs, command-line Rscript).
# ==============================================================================

# ==============================================================================
# PROJECT ROOT (HARDCODED FOR RELIABILITY)
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
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "bulk_deg_processing")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# INPUT PATHS
# ==============================================================================

# Raw bulk DEG data file
INPUT_FILE <- project_path("Datalink", "Bulk_Data", "Bulk_Significant_DEGs.txt")

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================

# Output directory for pipeline results
OUTPUT_DIR <- project_path("Datalink", "Bulk_Analysis")

# Intermediate files (cached)
BULK_ALL_FILE <- file.path(OUTPUT_DIR, "Bulk_DE_genes_All.qs")
FILTERED_FILE <- file.path(OUTPUT_DIR, "Filtered_Consistent_Genes.csv")

# Final outputs
ANNOTATED_FILE <- file.path(OUTPUT_DIR, "Filtered_Consistent_Genes_with_Descriptions.csv")
SEX_CHR_FILE <- file.path(OUTPUT_DIR, "Sex_Chr_Proportions_Intervals.csv")

# Recovery by expression output directory
RECOVERY_EXPR_DIR <- file.path(OUTPUT_DIR, "recovery_by_expression")

# Toy test output directory
TOY_OUTPUT_DIR <- file.path(OUTPUT_DIR, "toy")

# ==============================================================================
# FILTERING PARAMETERS
# ==============================================================================

# Tissue consistency threshold
# Genes must be significant in at least this proportion of tissues
TISSUE_CONSISTENCY_THRESHOLD <- 0.90

# Significance threshold (lfsr)
P_VALUE_THRESHOLD <- 0.05

# Log fold-change threshold (absolute value)
# Set to 0 to include all significant genes regardless of effect size
LFC_THRESHOLD <- 0

# ==============================================================================
# SEX CHROMOSOME ANALYSIS PARAMETERS
# ==============================================================================

# LFC intervals for sex chromosome analysis
# Format: list("label" = c(min, max))
# Use Inf for open-ended intervals
LFC_INTERVALS <- list(
  "0 to 0.5" = c(0,0.5),
  "> 0.5"     = c(0.5, Inf)
)

# ==============================================================================
# BIOMART PARAMETERS
# ==============================================================================

# BioMart dataset
BIOMART_DATASET <- "hsapiens_gene_ensembl"

# Attributes to fetch
BIOMART_ATTRIBUTES <- c("ensembl_gene_id", "external_gene_name", "description")

# ==============================================================================
# RUNTIME PARAMETERS
# ==============================================================================

# Number of cores for parallel operations (currently not used, but reserved)
N_CORES <- 2

# Force recalculation (ignore cached files)
FORCE_RECALC <- TRUE

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

# Validate that required input files exist
validate_inputs <- function() {
  if (!file.exists(INPUT_FILE)) {
    stop(sprintf("Input file not found: %s", INPUT_FILE))
  }
  return(TRUE)
}

# Create output directories if they don't exist
create_output_dirs <- function() {
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
    message(sprintf("Created output directory: %s", OUTPUT_DIR))
  }
  if (!dir.exists(TOY_OUTPUT_DIR)) {
    dir.create(TOY_OUTPUT_DIR, recursive = TRUE)
    message(sprintf("Created toy output directory: %s", TOY_OUTPUT_DIR))
  }
  return(TRUE)
}

# ==============================================================================
# CONFIGURATION SUMMARY
# ==============================================================================

print_config <- function() {
  cat("\n")
  cat("==============================================================================\n")
  cat("BULK DEG PROCESSING PIPELINE - CONFIGURATION\n")
  cat("==============================================================================\n")
  cat(sprintf("Project root: %s\n", PROJECT_ROOT))
  cat(sprintf("Script directory: %s\n", SCRIPT_DIR))
  cat(sprintf("Input file: %s\n", INPUT_FILE))
  cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
  cat("\n")
  cat("FILTERING PARAMETERS:\n")
  cat(sprintf("  Tissue consistency: %.0f%%\n", TISSUE_CONSISTENCY_THRESHOLD * 100))
  cat(sprintf("  P-value threshold: %.2f\n", P_VALUE_THRESHOLD))
  cat(sprintf("  LFC threshold: %.2f\n", LFC_THRESHOLD))
  cat("\n")
  cat("SEX CHROMOSOME INTERVALS:\n")
  for (label in names(LFC_INTERVALS)) {
    interval <- LFC_INTERVALS[[label]]
    cat(sprintf("  %s: [%.1f, %.1f)\n", label, interval[1], interval[2]))
  }
  cat("==============================================================================\n")
  cat("\n")
}

# ==============================================================================
# VISUALIZATION PARAMETERS (STEP 05)
# ==============================================================================

# Which single-cell datasets to use for concordance plots
# Options: "Roussos", "Allen", "DeJager", "Johansen", "Rexach"
# Can specify multiple: c("Roussos", "Allen")
SC_DATASETS_FOR_CONCORDANCE <- c("Roussos")

# Path to Full Ground Truth results (SC DE results)
FULL_GT_BASE <- project_path("Datalink", "Full_Ground_Truth")

# Gene labelling in plots
VIS_LABEL_GENES <- TRUE           # Label genes in concordance/volcano plots
VIS_CALC_SUBGROUP_LINES <- FALSE  # Calculate separate trend lines for autosomes/sex chr

# Histogram axis limits (for LFC and SE distributions)
VIS_HIST_XLIM <- c(2, 8)      # X-axis limits for |LFC| > 2 histogram
VIS_HIST_YLIM <- c(0, 5)      # Y-axis limits for |LFC| > 2 histogram
VIS_SE_XLIM <- c(0, 0.45)     # X-axis limits for SE distribution histogram
VIS_SE_YLIM <- c(0, 5)        # Y-axis limits for SE distribution histogram

# Histogram bin widths
VIS_LFC_BIN_WIDTH <- 0.5      # Bin width for LFC histograms
VIS_SE_BIN_WIDTH <- 0.01      # Bin width for SE histograms

# BioMart connection (for chromosome and biotype annotation)
# GTEx v8 used GENCODE v26 = Ensembl 88 (March 2017), but Ensembl 88 is no longer
# available in the archive. Ensembl 109 (Feb 2023) is the best available compromise:
# resolves 77/79 bulk genes (vs 74/79 with latest Ensembl). The 2 unresolvable genes
# (ENSG00000228906, ENSG00000205664) are retired IDs that fall back to "Autosome".
ENSEMBL_VERSION <- 109

# ==============================================================================
# RECOVERY CLASSIFICATION (STEPS 05 & 06)
# ==============================================================================

# 4-category recovery scheme for bulk→SC concordance
# Recovered:    FDR < 0.05 in SC AND matching LFC sign with bulk
# Insignificant: Gene present in SC with a valid padj, but padj >= 0.05 or wrong LFC direction
# Low-Count Filtered: Gene present in SC but padj = NA (DESeq2 independent filtering or outlier flagging)
# Absent:       Gene not present in SC dataset at all
RECOVERY_COLORS <- c(
  "Recovered"    = "#00BA38",
  "Insignificant" = "#E7B800",
  "Low-Count Filtered"  = "#F8766D",
  "Absent"       = "#999999"
)

# Ordered factor levels for consistent plotting (bottom to top in stacked bars)
RECOVERY_LEVELS <- c("Recovered", "Insignificant", "Low-Count Filtered", "Absent")

# ==============================================================================
# INITIALIZATION
# ==============================================================================

# Run validation and setup when config is sourced
if (!interactive()) {
  # Only print config when running non-interactively (e.g., via Rscript)
  # Suppress during interactive sessions to avoid clutter
  validate_inputs()
  create_output_dirs()
}
