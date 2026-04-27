# ==============================================================================
# CONFIG: Cartesian Downsampling Pipeline - Shared Configuration
# ==============================================================================
# This file contains all shared settings for the Cartesian downsampling power
# analysis pipeline. Source this file at the start of each script.
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
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "cartesian_downsampling")

# Helper function to get paths relative to project root
project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# 2. PATH CONFIGURATION
# ==============================================================================

# Input: Donor shards location
SHARD_BASE <- project_path("Datalink", "Shards")

# Dataset to use for downsampling (subdirectory within each cell type's shard folder)
# Options: "Roussos", "Roussos_DepthEq", "Allen", "DeJager", "Johansen", "Rexach"
SHARD_DATASET <- "Roussos"

# Output: Analysis results location
OUTPUT_BASE <- project_path("Datalink", "Cartesian_Downsampling")

# Note: Output structure is now {OUTPUT_BASE}/{CellType}/{Dataset}/
# The SUFFIX variables are no longer used - dataset name is used directly

# Feasibility data for overlay on plots (optional)
FEASIBILITY_DIR <- project_path("Datalink", "Cross-Dataset-Analysis", "Dataset_Sufficiency_Matrices")

# ==============================================================================
# 3. GRID CONFIGURATION
# ==============================================================================

# Downsampling mode: "cells" or "reads"
# - "cells": Downsample by cells per donor (original behavior)
# - "reads": Downsample by total reads (UMI counts) per donor
# >>> DEPTH EQUALIZATION: Run cells first, then change to "reads" and rerun <<<
# Can be overridden by environment variable DOWNSAMPLE_MODE (used by run_all.sh).
DOWNSAMPLE_MODE <- if (nzchar(Sys.getenv("DOWNSAMPLE_MODE"))) {
  Sys.getenv("DOWNSAMPLE_MODE")
} else {
  "cells"
}
stopifnot(DOWNSAMPLE_MODE %in% c("cells", "reads"))

# Donor steps (number of donors to sample)
DONOR_STEPS <- c(5, 10, 15, 25, 50, 75, 100, 150, 200, 300, 500, 600)

# Cell steps (cells per donor to sample) - used when DOWNSAMPLE_MODE == "cells"
# NA represents "Max" (use all available cells)
CELL_STEPS <- c(5, 10, 20, 50, 75, 100, 250, 500, 750, 1000, 1500, 2000, 3000, NA)

# Read steps (total reads per donor to sample) - used when DOWNSAMPLE_MODE == "reads"
# NA represents "Max" (use all available reads)
# Values chosen to roughly correspond to CELL_STEPS assuming ~5000 reads/cell
READ_STEPS <- c(25000, 50000, 100000, 250000, 500000, 
                1000000, 2500000, 5000000, 7500000, 
                10000000, 12500000, 15000000, NA)

# ==============================================================================
# 4. ROBUSTNESS SETTINGS
# ==============================================================================

# Number of independent robustness runs (each with different ground truth split)
N_ROBUSTNESS_RUNS <- 6

# Number of subsample iterations per grid point
N_SUBSAMPLES <- 20

# Ground truth split proportion (fraction of donors used for ground truth)
GT_SPLIT <- 0.5

# ==============================================================================
# 5. DESEQ2 SETTINGS
# ==============================================================================
 
# Minimum gene count filter (genes with fewer total counts are excluded)
MIN_GENE_COUNT <- 10
# DESeq2 dispersion fitType.
# - "parametric" (DESeq2 default): slower but is what the v8 paper used.
# - "glmGamPoi":  5-10x faster; published & first-party DESeq2 option
#                 (Ahlmann-Eltze & Huber 2020, Bioinformatics).
#                 Requires BiocManager::install("glmGamPoi").
# Change to "glmGamPoi" only after running the equivalence benchmark
# (tests/test_fittype_equivalence.R) and updating paper Methods.
DE_FIT_TYPE <- "glmGamPoi"

# FDR threshold for significance
FDR_THRESHOLD <- 0.05

# LFC intervals for binning results
LFC_INTERVALS <- list(
  "Tiny (0.2-1.0)"    = c(0.2, 1.0),
  "Small (1.0-2.0)"   = c(1.0, 2.0),
  "Medium (2.0-3.0)"  = c(2.0, 3.0),
  "Large (3.0+)"      = c(3.0, Inf)
)

# Chromosome categories for analysis
CHR_CATEGORIES <- c("Total", "Autosomal", "SexChr")

# Gene of interest for sanity check (XIST)
GENE_OF_INTEREST <- "ENSG00000229807"

# ==============================================================================
# 6. METADATA COLUMN MAPPINGS
# ==============================================================================

DONOR_ID_COL <- "donor_id"
SEX_COL <- "sex"
AGE_COL <- "development_stage"

# Sex values in the data
MALE_VAL <- "male"
FEMALE_VAL <- "female"

# Additional categorical covariates for DESeq2 model
CATEGORICAL_VARS <- c("Source", "self_reported_ethnicity", "disease")

# ==============================================================================
# 7. PARALLEL PROCESSING
# ==============================================================================

# Default number of cores for parallel processing
# Can be overridden by setting N_CORES before sourcing this file.
# Default lowered from 24 -> 20 (2026-04-21) per Chris Crampton: on the
# neurogenomics shared box, 64-core Multicore runs hit OOM because R's GC
# breaks copy-on-write semantics on donor-shard lapply. 20 runs stable;
# push higher only with Snow backend (see BPPARAM_BACKEND below) or after
# migrating to Nextflow per-grid-cell SLURM dispatch.
if (!exists("N_CORES")) {
  N_CORES <- 20
}

# Parallel backend selector. Two options:
#   "Multicore" (DEFAULT) — fork-based, fast spawn. On Linux, R's GC can
#                           break copy-on-write semantics and inflate RSS
#                           under load; capped at N_CORES ~= 20 on this box.
#   "Snow"                — PSOCK-based, separate R processes with predictable
#                           per-worker RSS. Slower spawn, stable steady-state.
#                           Worth trying if you need > 20 cores without OOM.
# Override by setting BPPARAM_BACKEND before sourcing this file.
if (!exists("BPPARAM_BACKEND")) {
  BPPARAM_BACKEND <- "Multicore"
}

# Maximum cores for memory-heavy operations
N_CORES_HEAVY <- min(6, N_CORES)

# Setup BiocParallel backend. Dispatches on BPPARAM_BACKEND.
setup_parallel <- function(n_cores = N_CORES, heavy = FALSE,
                           backend = if (exists("BPPARAM_BACKEND")) BPPARAM_BACKEND
                                      else "Multicore") {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel required. Install with: BiocManager::install('BiocParallel')")
  }

  cores_to_use <- if (heavy) min(n_cores, N_CORES_HEAVY) else n_cores

  if (backend == "Snow") {
    # PSOCK workers — separate R processes, predictable RSS, slower spawn.
    # Try this if 20-core Multicore still OOMs on large grid points.
    BiocParallel::SnowParam(
      workers = cores_to_use,
      type = "SOCK",
      progressbar = TRUE,
      stop.on.error = FALSE
    )
  } else {
    # Fork workers — fast spawn but memory can inflate under GC on Linux
    BiocParallel::MulticoreParam(
      workers = cores_to_use,
      progressbar = TRUE,
      stop.on.error = FALSE
    )
  }
}

# Create default BPPARAM objects
suppressPackageStartupMessages({
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BPPARAM <- setup_parallel(N_CORES, heavy = FALSE)
    BPPARAM_HEAVY <- setup_parallel(N_CORES, heavy = TRUE)
  }
})

# ==============================================================================
# 8. REPRODUCIBILITY
# ==============================================================================

# Base seed for reproducibility
# Each run uses SEED_BASE + run_id
SEED_BASE <- 1234

# ==============================================================================
# 9. CONTROL FLAGS
# ==============================================================================

# Force recalculation even if output files exist.
# Can be overridden by environment variable POWER_ANALYSIS_FORCE_RECALC
# (used by run_all.sh wrapper). Accepts: true/false/1/0/yes/no (case-insensitive).
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}

# ==============================================================================
# 10. CELL TYPE CONFIGURATION
# ==============================================================================

# Cell types are auto-detected from SHARD_BASE
# This function returns cell types that have the specified dataset
get_cell_types <- function(shard_base = SHARD_BASE, dataset = SHARD_DATASET) {
  subdirs <- list.dirs(shard_base, full.names = FALSE, recursive = FALSE)
  
  # Filter: keep only folders that contain the dataset subdirectory
  cell_types <- Filter(function(x) {
    check_path <- file.path(shard_base, x, dataset)
    return(dir.exists(check_path))
  }, subdirs)
  
  if (length(cell_types) == 0) {
    stop(sprintf("No cell type folders found containing '%s' in %s", dataset, shard_base))
  }
  
  return(cell_types)
}

# ==============================================================================
# 11. HELPER FUNCTIONS: Paths and Mode-Aware Utilities
# ==============================================================================

#' Get paths for a cell type, with mode-specific output directories
#' 
#' @param cell_type Cell type name
#' @param shard_base Base directory for shards
#' @param dataset Dataset name
#' @param output_base Base directory for outputs
#' @param mode Downsampling mode ("cells" or "reads")
#' @return List of paths
get_paths <- function(cell_type, 
                      shard_base = SHARD_BASE, 
                      dataset = SHARD_DATASET,
                      output_base = OUTPUT_BASE,
                      mode = DOWNSAMPLE_MODE) {
  list(
    # Input
    shard_dir = file.path(shard_base, cell_type, dataset),
    manifest_path = file.path(shard_base, cell_type, dataset, "00_MASTER_MANIFEST.qs"),
    
    # Output: {OUTPUT_BASE}/{CellType}/{Dataset}/{mode}/
    # Separate subdirectories for cells and reads modes
    results_dir = file.path(output_base, cell_type, dataset, mode),
    plots_dir = file.path(output_base, cell_type, dataset, mode, "Plots"),
    
    # Feasibility (optional)
    feasibility_path = file.path(FEASIBILITY_DIR, paste0("Feasibility_", cell_type, ".csv"))
  )
}

#' Get the Y-axis steps for the current mode
#' 
#' @param mode Downsampling mode ("cells" or "reads")
#' @return Vector of steps (CELL_STEPS or READ_STEPS)
get_y_steps <- function(mode = DOWNSAMPLE_MODE) {
  if (mode == "reads") {
    return(READ_STEPS)
  } else {
    return(CELL_STEPS)
  }
}

#' Get the Y-axis column name for the current mode
#' 
#' @param mode Downsampling mode ("cells" or "reads")
#' @return Column name ("Cells" or "Target_Reads")
get_y_col_name <- function(mode = DOWNSAMPLE_MODE) {
  if (mode == "reads") {
    return("Target_Reads")
  } else {
    return("Cells")
  }
}

#' Get the file pattern regex for the current mode
#' 
#' @param mode Downsampling mode ("cells" or "reads")
#' @return List with pattern and parse_regex
get_file_pattern <- function(mode = DOWNSAMPLE_MODE) {
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

# ==============================================================================
# 12. LOGGING HELPERS (fallback if shared library not loaded)
# ==============================================================================
# takes precedence when sourced by run_pipeline.R

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
# 13. VERSION INFO
# ==============================================================================

PIPELINE_VERSION <- "2.2.1"
PIPELINE_DATE <- "2026-03-31"

log_info("Cartesian Downsampling Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Mode: %s | Dataset: %s | Cores: %d | Runs: %d | Subsamples: %d", 
         DOWNSAMPLE_MODE, SHARD_DATASET, N_CORES, N_ROBUSTNESS_RUNS, N_SUBSAMPLES)
log_info("Grid: %d donors x %d %s steps = %d points",
         length(DONOR_STEPS), length(get_y_steps()), 
         ifelse(DOWNSAMPLE_MODE == "reads", "read", "cell"),
         length(DONOR_STEPS) * length(get_y_steps()))
