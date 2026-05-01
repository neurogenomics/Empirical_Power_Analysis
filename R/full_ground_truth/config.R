# ==============================================================================
# CONFIG: Full Ground Truth Pipeline - Shared Configuration
# ==============================================================================
# This file contains all shared settings for the full ground truth DE analysis
# pipeline. Source this file at the start of each script.
#
# Pipeline runs DESeq2 on all 7 datasets x 3 cell types using ALL donors,
# then performs cross-dataset diagnostics (UpSet, correlation, concordance).
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
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "full_ground_truth")

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
OUTPUT_BASE <- project_path("Datalink", "Full_Ground_Truth")

# ==============================================================================
# 3. CELL TYPES
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ==============================================================================
# 4. DATASET CONFIGURATION
# ==============================================================================
# Each dataset has different column names, covariates, sex encodings, and
# age parsing requirements. This is the single source of truth.

# --- Age Parsing Functions ---

# Complex string parser (e.g. "20 year-old stage")
parse_verbose_age <- function(x) {
  is_unknown_string <- function(s_in) {
    s <- tolower(trimws(as.character(s_in)))
    grepl("^(unknown|not available|na|n/a|nan)$", s)
  }

  s <- tolower(trimws(as.character(x)))
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))

  is_year    <- grepl("year[ -]?old.*stage", s)
  is_month   <- grepl("month[ -]?old.*stage", s)
  is_week    <- grepl("week[ -]?old.*stage", s)
  is_day     <- grepl("day[ -]?old.*stage", s)
  is_unknown <- is_unknown_string(s)

  y <- rep(NA_real_, length(s))
  y[is_year]    <- n[is_year]
  y[is_month]   <- n[is_month] / 12
  y[is_week]    <- n[is_week]  / 52.1429
  y[is_day]     <- n[is_day]   / 365.25
  y[is_unknown] <- NA_real_

  # Fallback: if no specific unit found but it is a number, assume years
  if (all(is.na(y)) && !any(is_unknown) && any(!is.na(n))) y <- n
  return(y)
}

# Simple numeric parser (e.g. "80", "80.5")
parse_simple_numeric <- function(x) {
  as.numeric(as.character(x))
}

# --- Dataset Definitions ---

DATASET_CONFIG <- list(

  "Roussos" = list(
    cols = list(
      donor = "donor_id",
      sex   = "sex",
      age   = "development_stage"
    ),
    covariates = c("Source", "self_reported_ethnicity", "disease"),
    # pmi (post-mortem interval, minutes) annotated onto shards from the
    # PsychAD supplementary table (media-1.csv) by R/sharding/Annotate_Roussos_PMI.R.
    # Standardised in run_deseq_for_dataset() to enter the design as pmi_scaled.
    numeric_covariates = c("pmi"),
    sex_values = list(male = "male", female = "female"),
    age_parser = parse_verbose_age
  ),

  "Allen" = list(
    cols = list(
      donor = "donor_id",
      sex   = "sex",
      age   = "age"
    ),
    covariates = c("ethnicity", "disease"),
    numeric_covariates = c("PMI"),
    sex_values = list(male = "male", female = "female"),
    age_parser = parse_simple_numeric
  ),

  "Rexach" = list(
    cols = list(
      donor = "donor_id",
      sex   = "sex",
      age   = "development_stage"
    ),
    covariates = c("self_reported_ethnicity", "disease"),
    numeric_covariates = c(),
    sex_values = list(male = "male", female = "female"),
    age_parser = parse_verbose_age
  ),

  # FIX: batch_vendor_name removed — 14 levels across 22 Astro donors caused
  # rank-deficient model matrix (confirmed DESeq2 error, log 2026-02-10).
  # library_prep_set had 1 level across all cell types (auto-filtered by pipeline).
  # medical_conditions retained (3–4 levels, no rank issues in any cell type).
  "Johansen" = list(
    cols = list(
      donor = "donor_name",
      sex   = "sex",
      age   = "age_numeric"
    ),
    covariates = c("medical_conditions"),
    numeric_covariates = c(),
    sex_values = list(male = "M", female = "F"),
    age_parser = parse_simple_numeric
  ),

  # DeJager metadata (age_death, msex, pmi, etc.) is now baked into shards
  # via Run_Shard_DeJager.R (joined from DeJagerMeta.csv at sharding time)
  "DeJager" = list(
    cols = list(
      donor = "individualID",
      sex   = "msex",
      age   = "age_death"
    ),
    covariates = c("study"),
    numeric_covariates = c("pmi"),
    sex_values = list(male = "1", female = "0"),
    age_parser = parse_simple_numeric
  )
)

# Which datasets to process (all by default)
RUN_DATASETS <- names(DATASET_CONFIG)

# ==============================================================================
# 5. DESEQ2 SETTINGS
# ==============================================================================

# Minimum gene count filter (genes with fewer total counts are excluded)
MIN_GENE_COUNT <- 10

# FDR threshold for significance
FDR_THRESHOLD <- 0.05

# LFC threshold for significance
LFC_THRESHOLD <- 1

# ==============================================================================
# 6. CROSS-DATASET DIAGNOSTIC SETTINGS
# ==============================================================================

# Reference dataset for correlation and concordance analyses
REFERENCE_DATASET <- "Roussos"

# Number of top genes for LFC correlation analysis
N_TOP_GENES <- 100

# Target gene for sanity check (XIST - X-chromosome inactivation)
TARGET_GENE_ID   <- "ENSG00000229807"
TARGET_GENE_NAME <- "XIST"

# ==============================================================================
# 6b. EXPANDED LFC AGREEMENT MATRIX SETTINGS
# ==============================================================================

# Cartesian downsampling output base (for reading 50:50 split metadata)
CARTESIAN_BASE <- file.path(PROJECT_ROOT, "Datalink", "Cartesian_Downsampling")

# Cartesian downsampling dataset and mode (must match cartesian config)
CARTESIAN_DATASET <- "Roussos"
CARTESIAN_MODE    <- "cells"

# Upper bound: number of cartesian runs to use for split-half comparison
N_UPPER_BOUND_RUNS <- 3

# Lower bound: number of null (MI-minimised label-flip) runs
N_NULL_RUNS <- 3

# Seed offset for null runs (null run i uses SEED_BASE + NULL_SEED_OFFSET + i)
NULL_SEED_OFFSET <- 101

# Cartesian seed base (must match cartesian config for donor split reproducibility)
CARTESIAN_SEED_BASE <- 1234

# ==============================================================================
# 7. PARALLEL PROCESSING
# ==============================================================================

# Cores for pseudobulk aggregation
N_CORES_AGG <- 4

# Cores for DESeq2
N_CORES_DE <- 4

# Setup BiocParallel backend
setup_parallel <- function(n_cores) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel required. Install with: BiocManager::install('BiocParallel')")
  }
  BiocParallel::MulticoreParam(
    workers = n_cores,
    progressbar = TRUE,
    stop.on.error = FALSE
  )
}

# ==============================================================================
# 8. GENE AUDIT MATRIX SETTINGS (Step 04)
# ==============================================================================

# BioMart cache for audit annotations (chromosome, biotype, band)
AUDIT_BIOMART_CACHE_FILE <- "gene_audit_biomart_cache.csv"

# Significance thresholds for the "significant gene universe" audit
AUDIT_FDR_THRESHOLD <- FDR_THRESHOLD   # 0.05
AUDIT_LFC_THRESHOLD <- LFC_THRESHOLD   # 1.0

# BioMart connection settings (reuse from existing pipeline if available)
AUDIT_BIOMART_DATASET <- "hsapiens_gene_ensembl"
AUDIT_BIOMART_HOST <- "https://useast.ensembl.org"
AUDIT_BIOMART_HOST_FALLBACK <- "https://www.ensembl.org"

# PAR region identification via BioMart chromosome band
# PAR1: Xp22.33 / Yp11.32 (tip of short arm)
# PAR2: Xq28 / Yq12 (tip of long arm)
# We match genes on X or Y whose band starts with these prefixes
PAR1_BANDS_X <- "p22.33"
PAR1_BANDS_Y <- "p11.32"
PAR2_BANDS_X <- "q28"
PAR2_BANDS_Y <- "q12"

# Pseudogene biotypes (comprehensive list from Ensembl)
PSEUDOGENE_BIOTYPES <- c(
  "processed_pseudogene", "unprocessed_pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
  "translated_processed_pseudogene", "translated_unprocessed_pseudogene",
  "polymorphic_pseudogene", "pseudogene",
  "IG_pseudogene", "TR_pseudogene",
  "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
  "TR_J_pseudogene", "TR_V_pseudogene",
  "rRNA_pseudogene", "unitary_pseudogene"
)

# Autosomal chromosomes
AUTOSOMAL_CHROMOSOMES <- as.character(1:22)

# Gene audit output file names
AUDIT_FULL_UNIVERSE_FILE <- "Gene_Audit_Full_Universe.csv"
AUDIT_SIGNIFICANT_FILE   <- "Gene_Audit_Significant.csv"
AUDIT_COMBINED_FULL_FILE <- "Gene_Audit_Full_Universe_Combined.csv"
AUDIT_COMBINED_SIG_FILE  <- "Gene_Audit_Significant_Combined.csv"

# ==============================================================================
# 9. REPRODUCIBILITY
# ==============================================================================

# No stochastic operations in this pipeline (uses all donors, no sampling)
# Seeds not required but documented for completeness
SEED_BASE <- 42

# ==============================================================================
# 10. CONTROL FLAGS
# ==============================================================================

# Force recalculation even if output files exist.
# Override via environment variable POWER_ANALYSIS_FORCE_RECALC (used by run_all.sh).
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}

# ==============================================================================
# 11. HELPER FUNCTION: Get Paths for Cell Type
# ==============================================================================

get_paths <- function(cell_type, output_base = OUTPUT_BASE) {
  list(
    # Input
    shard_base = file.path(SHARD_BASE, cell_type),

    # Output: {OUTPUT_BASE}/{CellType}/
    output_dir      = file.path(output_base, cell_type),
    diagnostic_dir  = file.path(output_base, cell_type, "Diagnostic_Results")
  )
}

# Get shard directory for a specific dataset and cell type
get_shard_dir <- function(cell_type, dataset) {
  file.path(SHARD_BASE, cell_type, dataset)
}

# ==============================================================================
# 12. LOGGING SETUP
# ==============================================================================

# Source shared logging infrastructure

# Fallback logging helpers (used if shared logging not available)
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
# 13. SANITY CHECKS
# ==============================================================================

# Source shared sanity check functions
SANITY_CHECKS_SCRIPT <- project_path("R", "shared", "sanity_checks.R")

# ==============================================================================
# 14. VERSION INFO
# ==============================================================================

PIPELINE_VERSION <- "1.1.0"
PIPELINE_DATE    <- "2026-02-19"

log_info("Full Ground Truth Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Datasets: %s", paste(RUN_DATASETS, collapse = ", "))
log_info("Cell Types: %s", paste(CELL_TYPES, collapse = ", "))
log_info("Cores: AGG=%d, DE=%d", N_CORES_AGG, N_CORES_DE)
