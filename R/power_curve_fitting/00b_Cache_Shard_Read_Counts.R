# ==============================================================================
# SCRIPT: Cache Shard Read Counts (Step 00b)
# PURPOSE: Pre-compute per-cell read counts (colSums(counts)) for all donor
#          shards and save as lightweight cache files. This avoids loading
#          full SCE objects (~50-200MB each) during metadata parsing in Step 01.
#
# v4.1.0: Loops over DATASET_CONFIGS to cache shards for all datasets
#          (e.g., Roussos + Roussos_DepthEq). Each dataset's cache is stored
#          under the model name (e.g., Shard_Cache/Astro_Shallow/ for DepthEq).
#
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Datalink/Power_Curve_Fitting/Shard_Cache/{ModelName}/donor_{ID}_reads.qs
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Matrix)
  library(SummarizedExperiment)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
    "config.R",
    "../config.R"
  )
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 00b: Cache Shard Read Counts (v4.1.0 — multi-dataset)")

# ==============================================================================
# 2. PROCESS EACH DATASET x CELL TYPE
# ==============================================================================

total_cached <- 0
total_skipped <- 0
total_failed <- 0

for (ds_config in DATASET_CONFIGS) {

  log_section(sprintf("Dataset: %s", ds_config$dataset))

  for (cell_type in ds_config$cell_types) {

    model_name <- get_dataset_model_name(cell_type, ds_config)
    log_section(sprintf("Cell Type: %s (model: %s, dataset: %s)",
                        cell_type, model_name, ds_config$dataset))

    paths <- get_paths_for_dataset(cell_type, ds_config)
    cache_dir <- paths$shard_cache_dir
    shard_dir <- paths$shard_dir

    safe_mkdir(cache_dir)

    # List all donor shard files
    shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
    # Exclude manifest
    shard_files <- shard_files[!grepl("MANIFEST", shard_files)]

    log_info("Found %d donor shards in %s", length(shard_files), shard_dir)

    if (length(shard_files) == 0) {
      log_warning("No shard files found for %s/%s", ds_config$dataset, cell_type)
      next
    }

    ct_cached <- 0
    ct_skipped <- 0
    ct_failed <- 0

    progress <- create_progress(length(shard_files),
                                label = sprintf("%s shard cache", model_name))

    for (shard_path in shard_files) {

      # Extract donor ID from filename: donor_{ID}.qs -> {ID}
      shard_basename <- basename(shard_path)
      donor_id <- sub("^donor_", "", sub("\\.qs$", "", shard_basename))

      cache_path <- file.path(cache_dir, paste0("donor_", donor_id, "_reads.qs"))

      # Skip if already cached (unless FORCE_RECALC)
      if (file.exists(cache_path) && !FORCE_RECALC) {
        ct_skipped <- ct_skipped + 1
        progress()
        next
      }

      # Load shard and compute per-cell read counts
      result <- tryCatch({
        sce <- qs::qread(shard_path)
        reads_per_cell <- Matrix::colSums(SummarizedExperiment::assay(sce, "counts"))
        qs::qsave(reads_per_cell, cache_path)
        rm(sce)
        "SUCCESS"
      }, error = function(e) {
        log_error("Failed to cache %s: %s", shard_basename, e$message)
        "FAILED"
      })

      if (result == "SUCCESS") {
        ct_cached <- ct_cached + 1
      } else {
        ct_failed <- ct_failed + 1
      }

      progress()

      # Periodic GC to manage memory
      if ((ct_cached + ct_skipped + ct_failed) %% 100 == 0) {
        gc(verbose = FALSE)
      }
    }

    log_info("%s: cached=%d, skipped=%d (already cached), failed=%d",
             model_name, ct_cached, ct_skipped, ct_failed)

    total_cached <- total_cached + ct_cached
    total_skipped <- total_skipped + ct_skipped
    total_failed <- total_failed + ct_failed
  }
}

# ==============================================================================
# 3. SUMMARY
# ==============================================================================

log_section("STEP 00b COMPLETE")
log_info("Total: cached=%d, skipped=%d, failed=%d",
         total_cached, total_skipped, total_failed)

if (total_failed > 0) {
  log_warning("%d shards failed to cache — check errors above", total_failed)
}

# Verify cache completeness across all datasets
for (ds_config in DATASET_CONFIGS) {
  for (cell_type in ds_config$cell_types) {
    paths <- get_paths_for_dataset(cell_type, ds_config)
    model_name <- get_dataset_model_name(cell_type, ds_config)
    n_shards <- length(list.files(paths$shard_dir, pattern = "^donor_.*\\.qs$")) - 1
    n_cached <- length(list.files(paths$shard_cache_dir, pattern = "^donor_.*_reads\\.qs$"))
    log_info("  %s (%s): %d/%d shards cached (%.1f%%)",
             model_name, ds_config$dataset,
             n_cached, n_shards, 100 * n_cached / max(n_shards, 1))
  }
}
