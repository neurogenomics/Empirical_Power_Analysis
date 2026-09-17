# ==============================================================================
# CONFIG: Depth Equalization Pipeline
# ==============================================================================
# Multinomially downsample all cells in a source shard directory to match a
# target mean reads-per-cell, producing depth-equalized shards for sensitivity
# analysis of sequencing depth confounding in power curve fitting.
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
project_path <- function(...) file.path(PROJECT_ROOT, ...)

# Base shard directory (needed for Step 7)
SHARD_BASE <- project_path("Datalink", "Shards")

# ==============================================================================
# 2. SOURCE AND TARGET CONFIGURATION
# ==============================================================================

# Cell type to equalize (only Astro for this sensitivity analysis)
SOURCE_CELL_TYPE <- "Astro"

# Source shard directory
SOURCE_DATASET <- "Roussos"
SOURCE_SHARD_DIR <- project_path("Datalink", "Shards", SOURCE_CELL_TYPE, SOURCE_DATASET)

# Target output directory
TARGET_DATASET <- "Roussos_DepthEq"
TARGET_SHARD_DIR <- project_path("Datalink", "Shards", SOURCE_CELL_TYPE, TARGET_DATASET)

# ==============================================================================
# 3. DOWNSAMPLING PARAMETERS
# ==============================================================================

# Source mean reads per cell (from EDA: Read_Depth_Summary_All_CellTypes_Roussos.csv)
SOURCE_MEAN_RPC <- 9456

# Target mean reads per cell (match Microglia)
TARGET_MEAN_RPC <- 4426

# Downsampling fraction: every cell's UMI count is multiplied by this
DOWNSAMPLE_FRACTION <- TARGET_MEAN_RPC / SOURCE_MEAN_RPC
# = 0.4679...

# ==============================================================================
# 3b. ASTRO_SHALLOW_FEWER EXTENSION CONFIGURATION
# ==============================================================================
# Extension: Match cell count distribution to microglia on top of depth
# equalization. This tests whether ground truth composition affects power curves
# when the GT is "easier" (fewer cells = less power = fewer genes detected).

FEWER_EXTENSION <- list(
  enabled = TRUE,
  
  # Target cell type to match distribution
  target_cell_type = "Micro",
  target_dataset = "Roussos",
  
  # Source: Start from depth-equalized shards
  source_dataset = "Roussos_DepthEq",
  
  # Matching strategy: "mean_ratio" = scale all donors by Micro_mean/Astro_mean
  scaling_mode = "mean_ratio",
  
  # Minimum cells threshold: donors below this after scaling are dropped
  min_cells_threshold = 30,
  
  # Output naming
  output_dataset = "Roussos_ShallowFewer",
  model_name = "Astro_Shallow_Fewer",
  
  # Seed for cell subsampling (different from depth equalization seed)
  seed_base = 9998
)

# ==============================================================================
# 4. REPRODUCIBILITY
# ==============================================================================

# Base seed for multinomial downsampling
# Each donor shard uses SEED_BASE + donor_index for reproducibility
SEED_BASE <- 9999

# ==============================================================================
# 5. PARALLEL PROCESSING
# ==============================================================================

if (!exists("N_CORES")) {
  N_CORES <- 24
}

# ==============================================================================
# 6. LOGGING
# ==============================================================================

if (!exists("log_info")) {
  log_info <- function(msg, ...) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s", timestamp, sprintf(msg, ...)))
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

if (!exists("log_section")) {
  log_section <- function(title) {
    cat("\n")
    log_info(paste(rep("=", 60), collapse = ""))
    log_info(title)
    log_info(paste(rep("=", 60), collapse = ""))
  }
}

# ==============================================================================
# 6b. MINI RUN OVERRIDES (smoke test — activate with env var MINI_RUN=true)
# ==============================================================================
# Mirrors the cartesian_downsampling / power_curve_fitting MINI switch. When
# active the depth-equalisation pipeline (a) caps the source Astro shard
# universe to a handful of donors (see run_pipeline.R), (b) writes the
# depth-equalised and shallow-fewer shards to MINI_-prefixed dataset dirs so the
# production 1,293-shard sets and their .depth_eq_complete / .shallow_fewer_complete
# sentinels are never touched, and (c) relaxes the shallow-fewer min-cells drop
# threshold so the tiny capped set isn't discarded entirely.
MINI_RUN <- tolower(Sys.getenv("MINI_RUN", "false")) %in% c("true", "1", "yes")

MINI_MAX_DONORS <- if (nzchar(Sys.getenv("MINI_MAX_DONORS"))) {
  as.integer(Sys.getenv("MINI_MAX_DONORS"))
} else {
  20L
}

if (MINI_RUN) {
  TARGET_DATASET   <- "MINI_Roussos_DepthEq"
  TARGET_SHARD_DIR <- project_path("Datalink", "Shards", SOURCE_CELL_TYPE, TARGET_DATASET)
  FEWER_EXTENSION$source_dataset      <- "MINI_Roussos_DepthEq"
  FEWER_EXTENSION$output_dataset      <- "MINI_Roussos_ShallowFewer"
  FEWER_EXTENSION$min_cells_threshold <- 5L
  log_warning(paste0("MINI_RUN active: depth-eq target -> %s, shallow-fewer -> %s, ",
                     "min_cells_threshold -> %d, source shards capped to %d (see run_pipeline.R)"),
              TARGET_DATASET, FEWER_EXTENSION$output_dataset,
              FEWER_EXTENSION$min_cells_threshold, 2L * MINI_MAX_DONORS)
}

# ==============================================================================
# 7. VERSION INFO
# ==============================================================================

# Force recalculation (ignore cached outputs).
# Override via environment variable POWER_ANALYSIS_FORCE_RECALC (used by run_all.sh).
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}

PIPELINE_VERSION <- "1.0.0"
PIPELINE_DATE <- "2026-02-25"

log_info("Depth Equalization Pipeline Config v%s loaded", PIPELINE_VERSION)
log_info("Source: %s/%s (mean RPC: %d)", SOURCE_CELL_TYPE, SOURCE_DATASET, SOURCE_MEAN_RPC)
log_info("Target: %s/%s (mean RPC: %d)", SOURCE_CELL_TYPE, TARGET_DATASET, TARGET_MEAN_RPC)
log_info("Downsampling fraction: %.4f", DOWNSAMPLE_FRACTION)
