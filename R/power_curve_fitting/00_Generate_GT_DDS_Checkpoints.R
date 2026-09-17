# ==============================================================================
# SCRIPT: Generate Ground Truth DDS Checkpoints
# PURPOSE: Re-run DESeq2 on the GT half of each cartesian downsampling run
#          to extract per-gene dispersion estimates (not saved by original pipeline)
#
# v4.1.0: Loops over DATASET_CONFIGS to generate GT DDS for all datasets
#          (e.g., Roussos + Roussos_DepthEq). Each dataset's checkpoints are
#          stored under the model name (e.g., GT_DDS_Checkpoints/Astro_Shallow/).
#
# INPUT: GT sampling metadata + donor shards from cartesian downsampling
# OUTPUT: GT_DDS_Run_{N}.qs per model name (DESeqDataSet objects with dispersions)
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(Matrix)
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

log_section("STEP 0: Generate GT DDS Checkpoints (v4.1.0 — multi-dataset)")

# ==============================================================================
# 2. PROCESS EACH DATASET x CELL TYPE
# ==============================================================================

for (ds_config in DATASET_CONFIGS) {

  log_section(sprintf("Dataset: %s", ds_config$dataset))

  for (cell_type in ds_config$cell_types) {

    model_name <- get_dataset_model_name(cell_type, ds_config)
    log_section(sprintf("Cell Type: %s (model: %s, dataset: %s)",
                        cell_type, model_name, ds_config$dataset))

    # Get paths (mode doesn't matter for GT DDS — GT splits are identical across modes)
    paths <- get_paths_for_dataset(cell_type, ds_config, mode = "cells")

    # Create output directory
    safe_mkdir(paths$gt_dds_dir)

    for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {

      log_info("--- Run %d of %d ---", run_id, N_ROBUSTNESS_RUNS)

      # Output path
      dds_path <- file.path(paths$gt_dds_dir, sprintf("GT_DDS_Run_%d.qs", run_id))

      # Skip if already exists
      if (file.exists(dds_path) && !FORCE_RECALC) {
        log_info("GT DDS checkpoint already exists, skipping: %s", basename(dds_path))
        next
      }

      # Load GT sampling metadata to identify which donors and cells were used
      run_dir <- file.path(paths$cartesian_dir, paste0("Run_", run_id))
      gt_meta_path <- file.path(run_dir, "Ground_Truth_sampling_metadata.qs")

      if (!file.exists(gt_meta_path)) {
        log_warning("GT sampling metadata not found for run %d: %s", run_id, gt_meta_path)
        next
      }

      # The sampling metadata is a named list: donor_id -> cell indices
      gt_sampling <- qs::qread(gt_meta_path)
      gt_donor_ids <- names(gt_sampling)

      log_info("GT donors: %d | Shard dir: %s", length(gt_donor_ids), paths$shard_dir)

      # Aggregate pseudobulk from shards using the exact same cells as the original GT
      log_info("Aggregating pseudobulk counts (using original GT cell indices)...")
      agg_result <- aggregate_from_shards(
        donor_ids = gt_donor_ids,
        shard_dir = paths$shard_dir,
        cell_indices_map = gt_sampling,
        seed = SEED_BASE + run_id,
        verbose = FALSE
      )

      if (is.null(agg_result)) {
        log_error("Failed to aggregate GT donors for %s run %d", model_name, run_id)
        next
      }

      log_info("Pseudobulk: %d genes x %d donors",
               nrow(SummarizedExperiment::assay(agg_result$sce, "counts")),
               ncol(agg_result$sce))

      # Run DESeq2 and save DDS
      log_info("Running DESeq2 (saving DDS for dispersion extraction)...")
      status <- tryCatch({
        run_deseq_for_dds(
          agg_result = agg_result,
          save_path_dds = dds_path,
          save_path_results = NULL,
          min_gene_count = MIN_GENE_COUNT
        )
      }, error = function(e) {
        log_error("DESeq2 failed for %s run %d: %s", model_name, run_id, e$message)
        return(e$message)
      })

      log_info("Status: %s", status)

      # Clean up
      rm(agg_result)
      gc()
    }

    log_info("Completed GT DDS generation for %s (%s)", model_name, ds_config$dataset)
  }
}

log_section("STEP 0 COMPLETE")
log_info("GT DDS checkpoints generated for all datasets and cell types")
