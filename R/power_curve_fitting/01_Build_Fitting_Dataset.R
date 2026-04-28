# ==============================================================================
# SCRIPT: Build Fitting Dataset (v4.1.0)
# PURPOSE: Parse sampling metadata to extract actual n_cells and reads_per_cell
#          per iteration, merge cells + reads modes, and build unified fitting
#          datasets for logistic regression.
#
# v4.1.0: Loops over DATASET_CONFIGS to build fitting data for all datasets.
#          Each dataset's fitting data is stored under the model name
#          (e.g., Fitting_Data/Astro_Shallow/ for DepthEq).
#          Shared gene identification uses the primary dataset (first in
#          DATASET_CONFIGS with fit_combined=TRUE). Additional datasets
#          reuse this shared gene list for filtering.
#
# Produces per model name:
#   - {ModelName}/unified_fitting_data_per_iteration.qs  (full genes)
#   - {ModelName}/unified_fitting_data_per_iteration_shared.qs  (shared genes)
#
# Plus shared gene list (from primary dataset only):
#   - Fitting_Data/shared_truth_genes.qs
#
# INPUT: GT DESeqResults + GT DDS checkpoints + sampling metadata + shard cache
# OUTPUT: Unified fitting datasets with actual n_cells and reads_per_cell
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
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

log_section("STEP 1: Build Unified Fitting Dataset (v4.1.0 — multi-dataset)")

# ==============================================================================
# 2. PROCESS EACH DATASET
# ==============================================================================

# Store fitting data from the primary dataset for shared gene identification
primary_fitting_data <- list()

# Track which dataset is the primary (first with fit_combined=TRUE)
primary_ds_idx <- which(sapply(DATASET_CONFIGS, function(x) x$fit_combined))[1]
if (is.na(primary_ds_idx)) primary_ds_idx <- 1

for (ds_idx in seq_along(DATASET_CONFIGS)) {

  ds_config <- DATASET_CONFIGS[[ds_idx]]
  is_primary <- (ds_idx == primary_ds_idx)

  log_section(sprintf("Dataset: %s %s", ds_config$dataset,
                      if (is_primary) "(PRIMARY)" else ""))

  # ------------------------------------------------------------------
  # 2a. Load shard caches for this dataset
  # ------------------------------------------------------------------
  # Fast-skip: if every cell type in this dataset already has its per-iteration
  # fitting file on disk, we don't need to load any shard caches at all
  # (~1,300 qs::qread calls per cell type). Just load the per-iter files for
  # the primary dataset so shared-gene identification can proceed.

  all_per_iter_done <- all(vapply(ds_config$cell_types, function(ct) {
    p <- get_paths_for_dataset(ct, ds_config)
    file.exists(file.path(p$fitting_dir, "unified_fitting_data_per_iteration.qs"))
  }, logical(1)))

  if (all_per_iter_done && !FORCE_RECALC) {
    log_info("All per-iteration fitting files exist for dataset %s -- skipping shard cache load.",
             ds_config$dataset)
    if (is_primary) {
      for (cell_type in ds_config$cell_types) {
        paths <- get_paths_for_dataset(cell_type, ds_config)
        per_iter_path <- file.path(paths$fitting_dir,
                                   "unified_fitting_data_per_iteration.qs")
        primary_fitting_data[[cell_type]] <- qs::qread(per_iter_path)
        log_info("  Loaded per-iter fitting data for %s (primary)",
                 get_dataset_model_name(cell_type, ds_config))
      }
    }
    next
  }

  log_section("Loading shard read count caches")

  shard_caches <- list()

  for (cell_type in ds_config$cell_types) {
    model_name <- get_dataset_model_name(cell_type, ds_config)
    paths <- get_paths_for_dataset(cell_type, ds_config)
    cache_dir <- paths$shard_cache_dir

    if (!dir.exists(cache_dir)) {
      log_error("Shard cache directory not found for %s: %s", model_name, cache_dir)
      log_error("Run Step 00b first to generate shard caches")
      stop(sprintf("Missing shard cache for %s. Run Step 00b first.", model_name))
    }

    cache_files <- list.files(cache_dir, pattern = "^donor_.*_reads\\.qs$", full.names = TRUE)
    log_info("Loading %d shard caches for %s...", length(cache_files), model_name)

    ct_cache <- list()
    for (cf in cache_files) {
      donor_id <- sub("^donor_", "", sub("_reads\\.qs$", "", basename(cf)))
      ct_cache[[donor_id]] <- qs::qread(cf)
    }

    # Key by cell_type (not model_name) since sampling metadata uses cell_type donor IDs
    shard_caches[[cell_type]] <- ct_cache
    log_info("  Loaded %d donor caches for %s", length(ct_cache), model_name)
  }

  # ------------------------------------------------------------------
  # 2b. Process each cell type in this dataset
  # ------------------------------------------------------------------

  for (cell_type in ds_config$cell_types) {

    model_name <- get_dataset_model_name(cell_type, ds_config)
    log_section(sprintf("Cell Type: %s (model: %s, dataset: %s)",
                        cell_type, model_name, ds_config$dataset))

    paths <- get_paths_for_dataset(cell_type, ds_config)
    safe_mkdir(paths$fitting_dir)

    output_path_per_iter <- file.path(paths$fitting_dir,
                                       "unified_fitting_data_per_iteration.qs")

    if (file.exists(output_path_per_iter) && !FORCE_RECALC) {
      log_info("Fitting data already exists, loading: %s", basename(output_path_per_iter))
      per_iter_data <- qs::qread(output_path_per_iter)
      if (is_primary) {
        primary_fitting_data[[cell_type]] <- per_iter_data
      }
      next
    }

    shard_cache <- shard_caches[[cell_type]]

    # ------------------------------------------------------------------
    # Collect per-run GT gene characteristics + detection outcomes
    # across BOTH modes (cells + reads)
    # ------------------------------------------------------------------

    all_run_data <- list()

    # Check if this is GT-only mode
    is_gt_only <- !is.null(ds_config$special_mode) && ds_config$special_mode == "gt_only"
    
    if (is_gt_only) {
      log_info("GT-ONLY MODE: Loading GT from %s, iterations from %s",
               ds_config$cartesian_target, ds_config$cartesian_source)
    }

    for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {

      log_info("--- Processing Run %d ---", run_id)

      # Determine GT path based on mode
      if (is_gt_only) {
        # GT-only mode: Load GT from the target dataset (e.g., Roussos_ShallowFewer)
        gt_path <- file.path(
          CARTESIAN_BASE, cell_type, ds_config$cartesian_target, "cells",
          paste0("Run_", run_id), "Ground_Truth.qs"
        )
      } else {
        # Normal mode: Load GT from current dataset
        gt_path <- file.path(
          get_paths_for_dataset(cell_type, ds_config, "cells")$cartesian_dir,
          paste0("Run_", run_id), "Ground_Truth.qs"
        )
      }
      
      if (!file.exists(gt_path)) {
        log_warning("GT results not found for run %d, skipping", run_id)
        next
      }
      gt_res <- as.data.frame(qs::qread(gt_path))
      gt_res$gene_id <- rownames(gt_res)

      # Load GT DDS checkpoint (for dispersion)
      dds_path <- file.path(paths$gt_dds_dir, sprintf("GT_DDS_Run_%d.qs", run_id))
      if (!file.exists(dds_path)) {
        log_warning("GT DDS not found for %s run %d, skipping", model_name, run_id)
        next
      }
      gt_dds <- qs::qread(dds_path)
      gt_dispersions <- DESeq2::dispersions(gt_dds)
      names(gt_dispersions) <- rownames(gt_dds)

      # Define truth set: GT padj < threshold AND |GT LFC| > threshold
      truth_genes <- gt_res %>%
        filter(!is.na(padj),
               padj < GT_PADJ_THRESHOLD,
               abs(log2FoldChange) > GT_LFC_THRESHOLD) %>%
        mutate(
          gt_lfc = log2FoldChange,
          abs_gt_lfc = abs(log2FoldChange),
          gt_basemean = baseMean,
          gt_dispersion = gt_dispersions[gene_id]
        ) %>%
        filter(!is.na(gt_dispersion),
               gt_dispersion > 0,
               gt_basemean > 0) %>%
        select(gene_id, gt_lfc, abs_gt_lfc, gt_basemean, gt_dispersion)

      n_truth <- nrow(truth_genes)
      log_info("Truth set: %d genes (padj < %.2f, |LFC| > %.1f)",
               n_truth, GT_PADJ_THRESHOLD, GT_LFC_THRESHOLD)

      if (n_truth == 0) {
        log_warning("No truth genes for run %d, skipping", run_id)
        next
      }

      # Process BOTH modes for this run
      for (mode in MODES) {

        log_info("  Mode: %s", mode)

        # Determine iteration directory based on mode
        if (is_gt_only) {
          # GT-only mode: Load iterations from source dataset (e.g., Roussos_DepthEq)
          run_dir <- file.path(
            CARTESIAN_BASE, cell_type, ds_config$cartesian_source, mode,
            paste0("Run_", run_id)
          )
        } else {
          # Normal mode: Use current dataset
          run_dir <- file.path(
            get_paths_for_dataset(cell_type, ds_config, mode)$cartesian_dir,
            paste0("Run_", run_id)
          )
        }

        if (!dir.exists(run_dir)) {
          log_warning("Run directory not found: %s", run_dir)
          next
        }

        # List iteration result files (DE results)
        iter_files <- list_iteration_files(run_dir, mode)
        if (is.null(iter_files) || nrow(iter_files) == 0) {
          log_warning("No iteration files found for %s run %d", mode, run_id)
          next
        }

        # List sampling metadata files
        meta_files <- list_metadata_files(run_dir, mode)
        if (is.null(meta_files) || nrow(meta_files) == 0) {
          log_warning("No metadata files found for %s run %d", mode, run_id)
          next
        }

        log_info("  Found %d iteration files, %d metadata files",
                 nrow(iter_files), nrow(meta_files))

        # Get unique grid points (excluding NA y_value)
        grid_points <- iter_files %>%
          filter(!is.na(y_value)) %>%
          distinct(n_donors, y_value)

        log_info("  Grid points: %d", nrow(grid_points))

        # For each grid point, compute per-gene detection + actual depth
        for (gp_idx in seq_len(nrow(grid_points))) {
          gp <- grid_points[gp_idx, ]

          # Get iteration files for this grid point
          gp_iter_files <- iter_files %>%
            filter(n_donors == gp$n_donors, y_value == gp$y_value)

          # Get metadata files for this grid point
          gp_meta_files <- meta_files %>%
            filter(n_donors == gp$n_donors, y_value == gp$y_value)

          n_iters <- nrow(gp_iter_files)
          if (n_iters == 0) next

          # Process each iteration
          for (iter_idx in seq_len(n_iters)) {
            iter_num <- gp_iter_files$iteration[iter_idx]

            # Load DE results for this iteration
            iter_res <- tryCatch({
              as.data.frame(qs::qread(gp_iter_files$file_path[iter_idx]))
            }, error = function(e) {
              return(NULL)
            })
            if (is.null(iter_res)) next

            iter_res$gene_id <- rownames(iter_res)

            # Detection: padj < threshold
            detected_genes <- iter_res$gene_id[!is.na(iter_res$padj) &
                                                 iter_res$padj < DETECT_PADJ_THRESHOLD]

            # Per-gene detection for truth genes
            detected <- as.integer(truth_genes$gene_id %in% detected_genes)

            # Compute actual depth from metadata
            meta_match <- gp_meta_files %>%
              filter(iteration == iter_num)

            n_cells_mean <- NA_real_
            rpc_mean <- NA_real_

            if (nrow(meta_match) == 1) {
              depth_info <- compute_iteration_depth(meta_match$file_path[1], shard_cache)
              if (!is.null(depth_info)) {
                n_cells_mean <- depth_info$n_cells_mean
                rpc_mean <- depth_info$reads_per_cell_mean
              }
            }

            # Skip if depth computation failed
            if (is.na(n_cells_mean) || is.na(rpc_mean) ||
                n_cells_mean <= 0 || rpc_mean <= 0) {
              next
            }

            # Build per-iteration row for each truth gene
            iter_data <- data.frame(
              gene_id = truth_genes$gene_id,
              detected = detected,
              n_donors = gp$n_donors,
              n_cells = n_cells_mean,
              reads_per_cell = rpc_mean,
              mode = mode,
              run_id = run_id,
              iteration = iter_num,
              stringsAsFactors = FALSE
            )

            key <- sprintf("%s_%s_%s_R%d_N%d_Y%s_I%d",
                           model_name, mode, ds_config$dataset, run_id,
                           gp$n_donors, gp$y_value, iter_num)
            all_run_data[[key]] <- iter_data
          }
        }

        # Clean up
        rm(iter_files, meta_files)
      }

      # Merge truth gene characteristics (once per run, shared across modes)
      truth_key <- sprintf("truth_%s_%s_R%d", model_name, ds_config$dataset, run_id)
      all_run_data[[truth_key]] <- truth_genes %>%
        mutate(run_id = run_id)

      rm(gt_res, gt_dds, gt_dispersions, truth_genes)
      gc()
    }

    # ------------------------------------------------------------------
    # Combine all iteration data and merge with gene characteristics
    # ------------------------------------------------------------------

    log_info("Combining iteration data...")

    # Separate truth gene data from iteration data
    truth_keys <- grep("^truth_", names(all_run_data), value = TRUE)
    iter_keys <- setdiff(names(all_run_data), truth_keys)

    if (length(iter_keys) == 0) {
      log_warning("No iteration data collected for %s", model_name)
      next
    }

    # Combine truth genes across runs (take mean of gene characteristics)
    truth_combined <- bind_rows(all_run_data[truth_keys]) %>%
      group_by(gene_id) %>%
      summarise(
        gt_lfc = mean(gt_lfc, na.rm = TRUE),
        abs_gt_lfc = mean(abs_gt_lfc, na.rm = TRUE),
        gt_basemean = mean(gt_basemean, na.rm = TRUE),
        gt_dispersion = mean(gt_dispersion, na.rm = TRUE),
        .groups = "drop"
      )

    # Combine all iteration data
    iter_combined <- bind_rows(all_run_data[iter_keys])
    log_info("Combined iteration data: %d rows", nrow(iter_combined))

    # Merge with gene characteristics
    iter_combined <- iter_combined %>%
      left_join(truth_combined, by = "gene_id")

    # ------------------------------------------------------------------
    # Build per-iteration fitting dataset
    # ------------------------------------------------------------------

    per_iter_data <- iter_combined %>%
      mutate(
        cell_type = model_name,
        log_donors = log(n_donors),
        log_n_cells = log(n_cells),
        log_reads_per_cell = log(reads_per_cell),
        log_abs_gt_lfc = log(abs_gt_lfc),
        log_gt_basemean = log(gt_basemean),
        log_gt_dispersion = log(gt_dispersion),
        mean_detection = as.numeric(detected)
      )

    # Save per-iteration dataset
    qs::qsave(per_iter_data, output_path_per_iter)
    log_info("Saved per-iteration fitting data: %d rows, %d genes",
             nrow(per_iter_data), n_distinct(per_iter_data$gene_id))

    # Store for shared gene identification (primary dataset only)
    if (is_primary) {
      primary_fitting_data[[cell_type]] <- per_iter_data
    }

    # Report depth statistics
    log_info("Depth statistics (per-iteration):")
    depth_stats <- per_iter_data %>%
      group_by(mode) %>%
      summarise(
        n_cells_min = min(n_cells), n_cells_max = max(n_cells),
        n_cells_median = median(n_cells),
        rpc_min = min(reads_per_cell), rpc_max = max(reads_per_cell),
        rpc_median = median(reads_per_cell),
        rpc_cv = sd(reads_per_cell) / mean(reads_per_cell),
        .groups = "drop"
      )
    for (i in seq_len(nrow(depth_stats))) {
      ds <- depth_stats[i, ]
      log_info("  %s: n_cells=[%.0f, %.0f] (med=%.0f), rpc=[%.0f, %.0f] (med=%.0f, CV=%.3f)",
               ds$mode, ds$n_cells_min, ds$n_cells_max, ds$n_cells_median,
               ds$rpc_min, ds$rpc_max, ds$rpc_median, ds$rpc_cv)
    }

    rm(iter_combined, all_run_data, truth_combined)
    gc()
  }

  # Clean up shard caches for this dataset
  rm(shard_caches)
  gc()
}

# ==============================================================================
# 3. IDENTIFY SHARED GENE SET (from primary dataset only)
# ==============================================================================

log_section("Identifying shared truth gene set")

shared_genes_path <- file.path(OUTPUT_BASE, "Fitting_Data", "shared_truth_genes.qs")

if (length(primary_fitting_data) < length(DATASET_CONFIGS[[primary_ds_idx]]$cell_types)) {
  log_warning("Not all primary cell types have fitting data, skipping shared gene identification")

  # Try to load existing shared genes
  if (file.exists(shared_genes_path)) {
    shared_genes <- qs::qread(shared_genes_path)
    log_info("Loaded existing shared gene list: %d genes", length(shared_genes))
  } else {
    shared_genes <- character(0)
    log_warning("No existing shared gene list found")
  }

} else {

  # Identify shared genes across all cell types in the primary dataset
  shared_genes <- identify_shared_genes(primary_fitting_data)

  if (length(shared_genes) == 0) {
    log_warning("No shared genes found")
  } else {
    # Save shared gene list
    safe_mkdir(dirname(shared_genes_path))
    qs::qsave(shared_genes, shared_genes_path)
    log_info("Saved shared gene list (%d genes): %s",
             length(shared_genes), shared_genes_path)
  }
}

# ==============================================================================
# 4. SAVE SHARED FITTING DATA FOR ALL DATASETS
# ==============================================================================

log_section("Saving shared fitting data for all datasets")

if (length(shared_genes) > 0) {

  for (ds_config in DATASET_CONFIGS) {
    for (cell_type in ds_config$cell_types) {

      model_name <- get_dataset_model_name(cell_type, ds_config)
      paths <- get_paths_for_dataset(cell_type, ds_config)

      # Per-iteration shared
      per_iter_path <- file.path(paths$fitting_dir,
                                  "unified_fitting_data_per_iteration.qs")
      shared_per_iter_path <- file.path(paths$fitting_dir,
                                         "unified_fitting_data_per_iteration_shared.qs")

      if (file.exists(per_iter_path)) {
        # Skip if shared data already exists (unless FORCE_RECALC)
        if (file.exists(shared_per_iter_path) && !FORCE_RECALC) {
          log_info("Shared data already exists for %s, skipping", model_name)
          next
        }

        per_iter_df <- qs::qread(per_iter_path)
        shared_per_iter <- per_iter_df %>% filter(gene_id %in% shared_genes)

        if (nrow(shared_per_iter) > 0) {
          qs::qsave(shared_per_iter, shared_per_iter_path)
          log_info("Saved shared per-iteration data for %s: %d rows (%d genes)",
                   model_name, nrow(shared_per_iter),
                   n_distinct(shared_per_iter$gene_id))
        } else {
          log_warning("No shared genes found in %s fitting data", model_name)
        }
      } else {
        log_warning("Full fitting data not found for %s: %s",
                    model_name, per_iter_path)
      }
    }
  }
} else {
  log_warning("No shared genes available — skipping shared fitting data")
}

log_section("STEP 1 COMPLETE")
log_info("Unified fitting datasets built for all datasets and cell types")
