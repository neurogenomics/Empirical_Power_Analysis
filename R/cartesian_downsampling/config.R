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
# Can be overridden by environment variable SHARD_DATASET (used by run_all.sh).
SHARD_DATASET <- if (nzchar(Sys.getenv("SHARD_DATASET"))) {
  Sys.getenv("SHARD_DATASET")
} else {
  "Roussos"
}

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

# Number of independent robustness runs (each with a different ground truth
# split). This is the TOTAL number of runs the project targets. Env-overridable
# so the total can grow (e.g. 10 -> 20) without a code edit. To compute only a
# sub-range incrementally (e.g. runs 6-10 after 1-5 are banked) set RUN_START /
# RUN_END in section 12c below: runs are seeded deterministically
# (SEED_BASE + run_id), so any sub-range is reproducible and disjoint.
N_ROBUSTNESS_RUNS <- if (nzchar(Sys.getenv("N_ROBUSTNESS_RUNS"))) {
  as.integer(Sys.getenv("N_ROBUSTNESS_RUNS"))
} else {
  10L
}
stopifnot(is.finite(N_ROBUSTNESS_RUNS), N_ROBUSTNESS_RUNS >= 1)

# Number of subsample iterations per grid point
N_SUBSAMPLES <- 20

# Ground truth split proportion (fraction of donors used for ground truth)
GT_SPLIT <- 0.5

# ==============================================================================
# 5. DIFFERENTIAL EXPRESSION SETTINGS
# ==============================================================================
 
# Minimum gene count filter (genes with fewer total counts are excluded)
MIN_GENE_COUNT <- 10

# Phase 1 runs all enabled fixed-effect DE methods on the same pseudobulk
# aggregation. dreamlet is reserved for a later mixed-model implementation.
IMPLEMENTED_DE_METHODS <- c("deseq2", "edger", "limma")
PLANNED_DE_METHODS <- c("dreamlet")
ALL_DE_METHODS <- c(IMPLEMENTED_DE_METHODS, PLANNED_DE_METHODS)
DE_METHODS <- if (nzchar(Sys.getenv("DE_METHODS"))) {
  trimws(strsplit(Sys.getenv("DE_METHODS"), ",", fixed = TRUE)[[1]])
} else {
  IMPLEMENTED_DE_METHODS
}
DE_METHODS <- DE_METHODS[nzchar(DE_METHODS)]
stopifnot(length(DE_METHODS) > 0, all(DE_METHODS %in% IMPLEMENTED_DE_METHODS))

# DESeq2 dispersion fitType.
# - "parametric" (DESeq2 default): slower but is what the v8 paper used.
# - "glmGamPoi":  5-10x faster; published & first-party DESeq2 option
#                 (Ahlmann-Eltze & Huber 2020, Bioinformatics).
#                 Must be present in renv.lock and the pinned Docker image.
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

# BaseMean intervals for binning by expression level
BASEMEAN_INTERVALS <- list(
  "Very Low (<10)"    = c(0, 10),
  "Low (10-50)"       = c(10, 50),
  "Medium (50-200)"   = c(50, 200),
  "High (200+)"       = c(200, Inf)
)

# Coarse bins for LFC x BaseMean cartesian product analysis
LFC_COARSE <- list(
  "|LFC| < 2"  = c(0.2, 2.0),
  "|LFC| >= 2" = c(2.0, Inf)
)

# Fine LFC x BaseMean cross-classification for the small-effect "figure of 4"
# (Figure T / Figure S). Each of these two LFC bins is crossed with the four
# BASEMEAN_INTERVALS categories so Step 3 can emit per-(donor, cells/reads) grid
# recall/precision within a fixed small-effect LFC band, stratified by
# expression. Names must exactly match keys of LFC_INTERVALS.
FINECROSS_LFC_BINS <- c("Tiny (0.2-1.0)", "Small (1.0-2.0)")

# Label used for the cross-method-pooled (paper "all-methods") outputs. Not a
# real DE method: it names a results/plots subdirectory that lives alongside the
# per-method ones (deseq2/edger/limma), holding the trial-weighted pool produced
# by 08_Aggregate_CrossMethod.R.
CROSS_METHOD_LABEL <- "all_methods"

BASEMEAN_COARSE <- list(
  "Low Expr (<50)"  = c(0, 50),
  "High Expr (50+)" = c(50, Inf)
)

# Expression filter threshold for the "high-expression only" LFC-split analysis.
# Genes with baseMean <= this value are excluded from recall (GT side) and
# precision (predicted side) to test how much low-expression genes affect metrics.
EXPR_FILTER_THRESHOLD <- 30

# Facets (e.g. a BaseMean expression category within an LFC bin) whose mean
# ground-truth gene count falls below this threshold carry essentially no genes,
# so any metric shown for them is misleading. Such facets are greyed out
# uniformly (all tiles set to NA -> grey) rather than displaying spurious values.
EMPTY_FACET_GENE_THRESHOLD <- 1

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

# Additional categorical covariates for the DE model
CATEGORICAL_VARS <- c("Source", "self_reported_ethnicity", "disease")

# Numeric covariates standardised (mean 0, sd 1) before entering the design.
# pmi (post-mortem interval, minutes) is annotated onto Roussos shards by
# R/sharding/Annotate_Roussos_PMI.R from the PsychAD supplementary table.
NUMERIC_VARS <- c("pmi")

# Contrast: sex by default; EPA_CONTRAST=ad_vs_control re-points SEX_COL /
# MALE_VAL / FEMALE_VAL / CATEGORICAL_VARS at the AD vs control subset.
source(file.path(PROJECT_ROOT, "R", "shared", "contrast.R"))
check_contrast_dataset(SHARD_DATASET)

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
    stop("BiocParallel is missing from the pinned R library. Restore renv.lock or rebuild the Docker image.")
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

# HPC array chunking. JOB_CHUNK is intentionally unset for normal local/server
# runs. When set by PBS, Step 2 loads the manifest generated by Step 1.5 and
# runs only that chunk.
JOB_N_CHUNKS <- if (nzchar(Sys.getenv("JOB_N_CHUNKS"))) {
  as.integer(Sys.getenv("JOB_N_CHUNKS"))
} else {
  40L
}
stopifnot(is.finite(JOB_N_CHUNKS), JOB_N_CHUNKS >= 1)

JOB_CHUNK <- if (nzchar(Sys.getenv("JOB_CHUNK"))) {
  as.integer(Sys.getenv("JOB_CHUNK"))
} else {
  NA_integer_
}
if (!is.na(JOB_CHUNK)) {
  stopifnot(is.finite(JOB_CHUNK), JOB_CHUNK >= 1, JOB_CHUNK <= JOB_N_CHUNKS)
}

JOB_MANIFEST_PATH_FROM_ENV <- nzchar(Sys.getenv("JOB_MANIFEST_PATH"))
JOB_MANIFEST_PATH <- Sys.getenv("JOB_MANIFEST_PATH")

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

  # MINI_RUN: restrict to the single configured cell type for a fast smoke test.
  if (exists("MINI_RUN") && isTRUE(MINI_RUN)) {
    if (MINI_CELL_TYPE %in% cell_types) {
      cell_types <- MINI_CELL_TYPE
    } else {
      log_warning("MINI_RUN: configured MINI_CELL_TYPE '%s' not found in shards; falling back to first available '%s'",
                  MINI_CELL_TYPE, cell_types[1])
      cell_types <- cell_types[1]
    }
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
                      mode = DOWNSAMPLE_MODE,
                      method = DE_METHODS[1]) {
  if (!is.null(method)) {
    stopifnot(length(method) == 1,
              method %in% c(IMPLEMENTED_DE_METHODS, CROSS_METHOD_LABEL))
  }
  results_base_dir <- file.path(output_base, cell_type, dataset, mode)
  results_dir <- if (is.null(method)) results_base_dir else file.path(results_base_dir, method)

  list(
    # Input
    shard_dir = file.path(shard_base, cell_type, dataset),
    manifest_path = file.path(shard_base, cell_type, dataset, "00_MASTER_MANIFEST.qs"),
    
    # Output: {OUTPUT_BASE}/{CellType}/{Dataset}/{mode}/{method}/
    # `results_base_dir` is the shared mode directory containing method folders.
    results_base_dir = results_base_dir,
    results_dir = results_dir,
    plots_dir = file.path(results_dir, "Plots"),
    method = method,
    
    # Feasibility (optional)
    feasibility_path = file.path(FEASIBILITY_DIR, paste0("Feasibility_", cell_type, ".csv"))
  )
}

get_paths_by_method <- function(cell_type,
                                methods = DE_METHODS,
                                shard_base = SHARD_BASE,
                                dataset = SHARD_DATASET,
                                output_base = OUTPUT_BASE,
                                mode = DOWNSAMPLE_MODE) {
  paths <- lapply(methods, function(method) {
    get_paths(
      cell_type = cell_type,
      shard_base = shard_base,
      dataset = dataset,
      output_base = output_base,
      mode = mode,
      method = method
    )
  })
  names(paths) <- methods
  paths
}

#' Methods to iterate over in the visualisation steps (04 / 07).
#'
#' Defaults to the real DE methods only, so local/host-native runs are
#' unchanged. When CARTESIAN_INCLUDE_ALL_METHODS is truthy (set by the HPC
#' Stage A figure phase after 08_Aggregate_CrossMethod.R has written the pooled
#' checkpoint), the cross-method label is appended so 04/07 also render the
#' paper "all-methods" heatmaps from the same plotting code.
get_viz_methods <- function() {
  methods <- DE_METHODS
  if (tolower(Sys.getenv("CARTESIAN_INCLUDE_ALL_METHODS", "false")) %in%
        c("true", "1", "yes")) {
    methods <- c(methods, CROSS_METHOD_LABEL)
  }
  methods
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
# 12b. MINI RUN OVERRIDES (smoke test — activate with env var MINI_RUN=true)
# ==============================================================================
# A tiny end-to-end run for verifying the full pipeline (and the Docker / HPC
# container) works without launching the multi-week production grid. When
# active it: shrinks the donor / cell / read grids, drops robustness and
# subsample counts, restricts to a single cell type (via get_cell_types above),
# and isolates ALL outputs under a MINI_-prefixed directory so production
# results are never touched. Driven entirely by environment variables so the
# committed defaults are unchanged for normal runs.
MINI_RUN <- tolower(Sys.getenv("MINI_RUN", "false")) %in% c("true", "1", "yes")

# PVAL_Q_FULL_PLOTS: if TRUE, step 4 generates Precision + Recall + F1 for all
# p-value quartile plots. If FALSE (default), only Recall is produced.
PVAL_Q_FULL_PLOTS <- tolower(Sys.getenv("PVAL_Q_FULL_PLOTS", "false")) %in% c("true", "1", "yes")

# Cell type the mini run is restricted to (env-overridable). get_cell_types()
# filters to this when MINI_RUN is active. Default is Astro so the smoke test
# also exercises the depth-equalisation variant pipelines (DepthEq /
# ShallowFewer), which are Astro-only.
MINI_CELL_TYPE <- if (nzchar(Sys.getenv("MINI_CELL_TYPE"))) {
  Sys.getenv("MINI_CELL_TYPE")
} else {
  "Astro"
}

# Cap on the donor universe loaded before the ground-truth split (used by
# load_donor_metadata() in utils.R). With GT_SPLIT = 0.5 this yields roughly
# MINI_MAX_DONORS/2 ground-truth donors and half for the downsampling pool, so
# no experiment ever exceeds ~MINI_MAX_DONORS/2 donors. Env-overridable.
MINI_MAX_DONORS <- if (nzchar(Sys.getenv("MINI_MAX_DONORS"))) {
  as.integer(Sys.getenv("MINI_MAX_DONORS"))
} else {
  40L
}

if (MINI_RUN) {
  DONOR_STEPS       <- c(10, 20)
  CELL_STEPS        <- c(50, 100, NA)
  READ_STEPS        <- c(100000, 250000, NA)
  N_ROBUSTNESS_RUNS <- 1
  N_SUBSAMPLES      <- 2
  # Isolate outputs so the smoke test never overwrites production results.
  OUTPUT_BASE       <- project_path("Datalink", "MINI_Cartesian_Downsampling")
  log_warning(paste0("MINI_RUN active: grid shrunk to %d donor x %d cell steps, ",
                     "%d run(s) x %d subsample(s), cell type '%s', max donors %d, outputs -> %s"),
              length(DONOR_STEPS), length(CELL_STEPS), N_ROBUSTNESS_RUNS, N_SUBSAMPLES,
              MINI_CELL_TYPE, MINI_MAX_DONORS, OUTPUT_BASE)
}

# ==============================================================================
# 12c. ROBUSTNESS RUN RANGE (incremental / parallel-grid support)
# ==============================================================================
# RUN_START..RUN_END selects which robustness runs THIS invocation generates and
# computes. Defaults to the full 1..N_ROBUSTNESS_RUNS set, so leaving them unset
# reproduces today's behaviour. Set RUN_START=6 RUN_END=10 to extend an
# already-banked 1..5 grid with ZERO recompute: seeds are SEED_BASE + run_id, so
# runs 1..5 stay byte-identical and runs 6..10 get fresh disjoint seeds.
# Generalises cleanly: 11..20 next time is just RUN_START=11 RUN_END=20.
# Defined after the MINI_RUN block so MINI's N_ROBUSTNESS_RUNS <- 1 is respected.
RUN_START <- if (nzchar(Sys.getenv("RUN_START"))) {
  as.integer(Sys.getenv("RUN_START"))
} else {
  1L
}
RUN_END <- if (nzchar(Sys.getenv("RUN_END"))) {
  as.integer(Sys.getenv("RUN_END"))
} else {
  N_ROBUSTNESS_RUNS
}
stopifnot(is.finite(RUN_START), is.finite(RUN_END),
          RUN_START >= 1, RUN_END >= RUN_START, RUN_END <= N_ROBUSTNESS_RUNS)
RUN_IDS <- RUN_START:RUN_END

JOB_MANIFEST_DIR <- file.path(OUTPUT_BASE, "_hpc_job_manifests")
if (!JOB_MANIFEST_PATH_FROM_ENV) {
  # Manifest name carries the run-range so an incremental (e.g. 6-10) manifest
  # never clobbers the 1-5 manifest, and the gt-manifest + grid jobs derive the
  # same path from the same RUN_START/RUN_END env. Method is intentionally NOT
  # in the name: the manifest lists grid points (donor x cell x iter), which are
  # identical across DE methods, so a dreamlet grid reuses the same manifest.
  JOB_MANIFEST_PATH <- file.path(
    JOB_MANIFEST_DIR,
    sprintf("manifest_%s_%s_runs%02d-%02d_%03d_chunks.csv",
            SHARD_DATASET, DOWNSAMPLE_MODE, RUN_START, RUN_END, JOB_N_CHUNKS)
  )
}

# ==============================================================================
# 13. VERSION INFO
# ==============================================================================

PIPELINE_VERSION <- "2.3.0"
PIPELINE_DATE <- "2026-05-28"

log_info("Cartesian Downsampling Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Mode: %s | Dataset: %s | Cores: %d | Runs: %d (active range %d-%d) | Subsamples: %d",
         DOWNSAMPLE_MODE, SHARD_DATASET, N_CORES, N_ROBUSTNESS_RUNS,
         RUN_START, RUN_END, N_SUBSAMPLES)
log_info("Grid: %d donors x %d %s steps = %d points",
         length(DONOR_STEPS), length(get_y_steps()), 
         ifelse(DOWNSAMPLE_MODE == "reads", "read", "cell"),
         length(DONOR_STEPS) * length(get_y_steps()))
