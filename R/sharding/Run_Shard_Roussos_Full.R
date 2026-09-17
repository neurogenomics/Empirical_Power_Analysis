# ==============================================================================
# Shared config (resolves PROJECT_ROOT via here::here() or POWER_ANALYSIS_ROOT env var).
# ==============================================================================
if (!exists("PROJECT_ROOT")) {
  .sharding_config <- tryCatch({
    # commandArgs path works under `Rscript R/sharding/<this>.R`
    .args <- commandArgs(trailingOnly = FALSE)
    .file_arg <- grep("--file=", .args, value = TRUE)
    if (length(.file_arg)) {
      file.path(dirname(sub("--file=", "", .file_arg[1])), "config.R")
    } else if (requireNamespace("here", quietly = TRUE)) {
      file.path(here::here(), "R", "sharding", "config.R")
    } else {
      file.path(getwd(), "R", "sharding", "config.R")
    }
  }, error = function(e) file.path(getwd(), "R", "sharding", "config.R"))
  source(.sharding_config)
  rm(.sharding_config)
}

#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Unfiltered Roussos Shard Creator
# PURPOSE: Reads PsychAD H5AD files and creates donor-level .qs shards
#          WITHOUT age or cell-count filtering. Produces the complete
#          (unfiltered) Roussos donor set for use by the EDA pipeline.
# INPUT: PsychAD H5AD files (HBCC, RADC, MSSM) from Datalink/raw_data/
# OUTPUT: Datalink/Shards/{CellType}/Roussos_Full/donor_*.qs + manifest
# USAGE:
#   Rscript R/sharding/Run_Shard_Roussos_Full.R Astro
#   Rscript R/sharding/Run_Shard_Roussos_Full.R Micro
#   Rscript R/sharding/Run_Shard_Roussos_Full.R Oligo
# NOTE: The FILTERED shards at Datalink/Shards/{CellType}/Roussos/ are
#       NOT modified. Those remain the correct input for power analysis
#       pipelines (cartesian_downsampling, simulations, full_ground_truth).
# ==============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(qs)
  library(dplyr)
  library(Matrix)
})

# --- SETUP PYTHON ---
if (!py_module_available("anndata")) {
  stop(
    "Python module 'anndata' is unavailable. ",
    "Use the pinned Docker image or rebuild it from docker/python-requirements.txt; ",
    "do not install Python packages at runtime."
  )
}
ad <- import("anndata")

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# --- Cell type from command line ---
args <- commandArgs(trailingOnly = TRUE)

# Map short names to the cell_type values used in the H5AD metadata
CELL_TYPE_MAP <- list(
  "Astro"  = "astrocyte",
  "Micro"  = "microglial cell",
  "Oligo"  = "oligodendrocyte"
)

if (length(args) < 1 || !args[1] %in% names(CELL_TYPE_MAP)) {
  stop(sprintf(
    "Usage: Rscript Run_Shard_Roussos_Full.R <CellType>\n  Valid cell types: %s",
    paste(names(CELL_TYPE_MAP), collapse = ", ")
  ))
}

CELL_TYPE_SHORT <- args[1]                        # e.g. "Astro"
TARGET_CELL_TYPE <- CELL_TYPE_MAP[[CELL_TYPE_SHORT]]  # e.g. "astrocyte"

message(sprintf("\n[CONFIG] Cell type: %s (%s)", CELL_TYPE_SHORT, TARGET_CELL_TYPE))
message("[CONFIG] Mode: UNFILTERED (no age or cell-count filters)")

# PATHS
INPUT_FILES <- c(
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-HBCC.h5ad"),
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-RADC.h5ad"),
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-MSSM.h5ad")
)

# Output Directory — Roussos_Full (separate from filtered Roussos)
BASE_OUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")
SHARD_DIR <- file.path(BASE_OUT_DIR, CELL_TYPE_SHORT, "Roussos_Full")
PROGRESS_LOG <- file.path(SHARD_DIR, "sharding_progress.rds")

# DATA COLUMNS
DONOR_ID_COL  <- "donor_id"
CELL_TYPE_COL <- "cell_type"

# PROCESSING
CHUNK_SIZE <- 50000

if (!dir.exists(SHARD_DIR)) dir.create(SHARD_DIR, recursive = TRUE)

message(sprintf("[CONFIG] Output: %s", SHARD_DIR))

# ==============================================================================
# 2. HELPER FUNCTIONS (Resume Logic)
# ==============================================================================

get_progress <- function() {
  if (file.exists(PROGRESS_LOG)) return(readRDS(PROGRESS_LOG))
  return(data.frame(file = character(), chunk = integer(), stringsAsFactors = FALSE))
}

update_progress <- function(file_name, chunk_idx) {
  progress <- get_progress()
  progress <- rbind(progress, data.frame(file = file_name, chunk = chunk_idx))
  saveRDS(progress, PROGRESS_LOG)
}

# ==============================================================================
# 3. CORE PROCESSING LOOP
# ==============================================================================

for (f_path in INPUT_FILES) {
  if (!file.exists(f_path)) {
    message(sprintf("!!! [SKIP] File not found: %s", f_path))
    next
  }

  f_name <- basename(f_path)
  message(sprintf("\n>>> Opening File: %s", f_name))

  # Open in backed mode to save RAM
  adata <- ad$read_h5ad(f_path, backed = "r")

  # --------------------------------------------------------------------------
  # STEP A: Pre-Scan Metadata to Identify ALL Donors of Target Cell Type
  # --------------------------------------------------------------------------
  message("   [PRE-SCAN] Loading metadata to determine donors...")

  # Load only metadata into R
  full_meta_df <- as.data.frame(adata$obs)

  # Filter by Cell Type only — NO age or cell-count filtering
  type_matches <- full_meta_df[[CELL_TYPE_COL]] == TARGET_CELL_TYPE
  target_meta  <- full_meta_df[type_matches, ]

  if (nrow(target_meta) == 0) {
    message(sprintf("   [SKIP] No cells of type '%s' found in this file.", TARGET_CELL_TYPE))
    adata$file$close()
    next
  }

  # Get ALL donor IDs (no filtering)
  valid_donor_ids <- unique(target_meta[[DONOR_ID_COL]])

  # Report donor stats for transparency
  donor_stats <- target_meta %>%
    group_by(!!sym(DONOR_ID_COL)) %>%
    summarise(n_cells = n(), .groups = "drop")

  message(sprintf("   [INFO] Total donors with %s cells: %d",
                  TARGET_CELL_TYPE, length(valid_donor_ids)))
  message(sprintf("   [INFO] Cell count range: %d - %d (median: %.0f)",
                  min(donor_stats$n_cells), max(donor_stats$n_cells),
                  median(donor_stats$n_cells)))

  # --------------------------------------------------------------------------
  # STEP B: Chunk Processing (Reading/saving ALL donors of target cell type)
  # --------------------------------------------------------------------------

  n_cells  <- as.numeric(adata$n_obs)
  n_chunks <- ceiling(n_cells / CHUNK_SIZE)
  gene_names <- rownames(adata$var)

  completed_chunks <- get_progress() %>%
    dplyr::filter(file == f_name) %>%
    dplyr::pull(chunk)

  for (i in seq_len(n_chunks)) {
    if (i %in% completed_chunks) {
      message(sprintf("      [SKIP] Chunk %d/%d (Already processed)", i, n_chunks))
      next
    }

    # Calculate Python indices (0-based)
    py_idx_start <- (i - 1) * CHUNK_SIZE
    py_idx_end   <- min(i * CHUNK_SIZE, n_cells)

    # Calculate R indices (1-based) for metadata subsetting
    r_idx_start <- py_idx_start + 1
    r_idx_end   <- r_idx_start + (py_idx_end - py_idx_start) - 1

    # Get metadata for this chunk
    chunk_meta_all <- full_meta_df[r_idx_start:r_idx_end, , drop = FALSE]

    # ----------------------------------------------------------------------
    # Filter Indices for this Chunk
    # Keep row IF: (Is Target Cell Type) — that's the ONLY filter
    # ----------------------------------------------------------------------
    keep_indices <- which(
      chunk_meta_all[[CELL_TYPE_COL]] == TARGET_CELL_TYPE &
        chunk_meta_all[[DONOR_ID_COL]] %in% valid_donor_ids
    )

    if (length(keep_indices) == 0) {
      if (i %% 10 == 0) message(sprintf("      [INFO] Chunk %d/%d: No valid cells.", i, n_chunks))
      update_progress(f_name, i)
      next
    }

    message(sprintf("      [PROC] Chunk %d/%d: Processing %d cells...",
                    i, n_chunks, length(keep_indices)))

    # Pull data from Python (Only load the chunk)
    chunk_slice <- adata[py_idx_start:(py_idx_end - 1)]

    # Transpose to Sparse Matrix (Genes x Cells) and subset
    counts_sparse <- t(chunk_slice$X)

    # Dimension safety check
    chunk_n_cells <- ncol(counts_sparse)
    meta_n_cells  <- nrow(chunk_meta_all)
    if (chunk_n_cells != meta_n_cells) {
      message(sprintf("      [WARN] Chunk %d: matrix cols (%d) != metadata rows (%d). Adjusting keep_indices.",
                      i, chunk_n_cells, meta_n_cells))
      keep_indices <- keep_indices[keep_indices <= chunk_n_cells]
      if (length(keep_indices) == 0) {
        message(sprintf("      [SKIP] Chunk %d: No valid indices after adjustment.", i))
        update_progress(f_name, i)
        rm(counts_sparse, chunk_slice); gc()
        next
      }
    }

    counts_sparse <- counts_sparse[, keep_indices, drop = FALSE]

    # Subset metadata
    chunk_meta_filtered <- chunk_meta_all[keep_indices, , drop = FALSE]

    rownames(counts_sparse) <- gene_names
    colnames(counts_sparse) <- rownames(chunk_meta_filtered)

    # Create temp SCE
    sce_small <- SingleCellExperiment(
      assays  = list(counts = counts_sparse),
      colData = chunk_meta_filtered
    )

    # ----------------------------------------------------------------------
    # Save/Append to Donor Shards
    # ----------------------------------------------------------------------
    donor_ids <- as.character(colData(sce_small)[[DONOR_ID_COL]])
    unique_chunk_donors <- unique(donor_ids)

    for (id in unique_chunk_donors) {
      safe_id <- gsub("[^[:alnum:]]", "_", id)
      donor_file <- file.path(SHARD_DIR, paste0("donor_", safe_id, ".qs"))

      donor_subset <- sce_small[, donor_ids == id]

      if (file.exists(donor_file)) {
        # Append to existing
        old_shard <- qs::qread(donor_file)
        common_genes <- intersect(rownames(old_shard), rownames(donor_subset))
        new_shard <- cbind(old_shard[common_genes, ], donor_subset[common_genes, ])
        qs::qsave(new_shard, donor_file)
      } else {
        # Create new
        qs::qsave(donor_subset, donor_file)
      }
    }

    # Cleanup & Progress
    update_progress(f_name, i)
    rm(sce_small, counts_sparse, chunk_slice, chunk_meta_filtered); gc()
  }

  adata$file$close()
  rm(full_meta_df); gc()
}

# ==============================================================================
# 4. FINAL MANIFEST GENERATION
# ==============================================================================
message("\n===> Building Master Manifest for Roussos_Full...")

all_shards <- list.files(SHARD_DIR, pattern = "^donor_.*\\.qs$", full.names = TRUE)

if (length(all_shards) == 0) {
  warning("No cells were saved. Check that the H5AD files exist and contain the target cell type.")
} else {

  manifest_df <- do.call(rbind, lapply(all_shards, function(f) {
    tmp <- qs::qread(f)
    data.frame(
      donor_id   = as.character(unique(colData(tmp)[[DONOR_ID_COL]])),
      cell_count = ncol(tmp),
      file_path  = f,
      stringsAsFactors = FALSE
    )
  }))

  sharding_final <- list(
    manifest    = manifest_df,
    config      = list(
      donor_col   = DONOR_ID_COL,
      n_donors    = nrow(manifest_df),
      total_cells = sum(manifest_df$cell_count),
      cell_filter = TARGET_CELL_TYPE,
      cell_type_short = CELL_TYPE_SHORT,
      filtered    = FALSE,
      description = "Unfiltered Roussos shards — no age or cell-count filters applied"
    )
  )

  qs::qsave(sharding_final, file.path(SHARD_DIR, "00_MASTER_MANIFEST.qs"))

  message("\n================================================================")
  message(" UNFILTERED SHARDING COMPLETE ")
  message(sprintf(" Cell Type:           %s (%s)", CELL_TYPE_SHORT, TARGET_CELL_TYPE))
  message(sprintf(" Filters Applied:     NONE"))
  message(sprintf(" Total Donors Saved:  %d", sharding_final$config$n_donors))
  message(sprintf(" Total Cells Saved:   %d", sharding_final$config$total_cells))
  message(sprintf(" Output:              %s", SHARD_DIR))
  message("================================================================")
}
