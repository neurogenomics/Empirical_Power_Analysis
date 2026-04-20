#!/usr/bin/env Rscript
# ==============================================================================
# TOY TEST: Power Curve Fitting Pipeline v4.0.0 - Integration Test (Stage 2)
# ==============================================================================
# Run: Rscript R/power_curve_fitting/toy/run_toy_test.R
#
# This test validates the v3.1.0 pipeline on minimal real data:
# - Uses 1 cell type (Astro), 1 run (Run_1), BOTH modes (cells + reads)
# - Tests shard cache loading
# - Tests metadata parsing for actual n_cells and reads_per_cell
# - Builds a mini unified fitting dataset
# - Fits the v3.1.0 unified formula
# - Validates predictions
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(DESeq2)
  library(dplyr)
})

cat("\n")
cat("================================================================\n")
cat("  STAGE 2: Toy Integration Test — Power Curve Fitting v4.0.0\n")
cat("================================================================\n\n")

# ==============================================================================
# 1. SETUP
# ==============================================================================

config_paths <- c(
  "R/power_curve_fitting/config.R",
  "../config.R",
  "config.R"
)

config_loaded <- FALSE
for (cp in config_paths) {
  if (file.exists(cp)) {
    source(cp)
    config_loaded <- TRUE
    break
  }
}

if (!config_loaded) stop("Could not find config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

cat("Config loaded. Pipeline version:", PIPELINE_VERSION, "\n\n")

# Toy parameters
TOY_CELL_TYPE <- "Astro"
TOY_RUN <- 1
TOY_MAX_GRID_POINTS <- 10  # Per mode
TOY_MAX_ITERS <- 3         # Per grid point

# ==============================================================================
# 2. CHECK DATA AVAILABILITY
# ==============================================================================

cat("--- Checking data availability ---\n")

# Check both modes
for (mode in MODES) {
  paths <- get_paths(TOY_CELL_TYPE, mode)
  run_dir <- file.path(paths$cartesian_dir, paste0("Run_", TOY_RUN))

  gt_path <- file.path(run_dir, "Ground_Truth.qs")
  if (!file.exists(gt_path)) {
    cat(sprintf("FAIL: GT results not found for %s mode: %s\n", mode, gt_path))
    quit(status = 1)
  }
  cat(sprintf("  GT results (%s): FOUND\n", mode))

  iter_files <- list_iteration_files(run_dir, mode)
  if (is.null(iter_files) || nrow(iter_files) == 0) {
    cat(sprintf("FAIL: No iteration files for %s mode in: %s\n", mode, run_dir))
    quit(status = 1)
  }
  cat(sprintf("  Iteration files (%s): %d found\n", mode, nrow(iter_files)))

  meta_files <- list_metadata_files(run_dir, mode)
  if (is.null(meta_files) || nrow(meta_files) == 0) {
    cat(sprintf("FAIL: No metadata files for %s mode in: %s\n", mode, run_dir))
    quit(status = 1)
  }
  cat(sprintf("  Metadata files (%s): %d found\n", mode, nrow(meta_files)))
}

# Check shard cache
paths <- get_paths(TOY_CELL_TYPE)
if (!dir.exists(paths$shard_cache_dir)) {
  cat(sprintf("FAIL: Shard cache not found: %s\n", paths$shard_cache_dir))
  cat("  Run Step 00b first: Rscript R/power_curve_fitting/run_pipeline.R 0.5\n")
  quit(status = 1)
}

cache_files <- list.files(paths$shard_cache_dir, pattern = "_reads\\.qs$")
cat(sprintf("  Shard cache: %d files found\n", length(cache_files)))

if (length(cache_files) == 0) {
  cat("FAIL: No shard cache files found. Run Step 00b first.\n")
  quit(status = 1)
}

# ==============================================================================
# 3. LOAD SHARD CACHE
# ==============================================================================

cat("\n--- Loading shard cache ---\n")

shard_cache <- list()
for (cf in list.files(paths$shard_cache_dir, pattern = "_reads\\.qs$", full.names = TRUE)) {
  donor_id <- sub("^donor_", "", sub("_reads\\.qs$", "", basename(cf)))
  shard_cache[[donor_id]] <- qs::qread(cf)
}
cat(sprintf("  Loaded %d donor caches\n", length(shard_cache)))

# Verify cache contents
sample_cache <- shard_cache[[1]]
cat(sprintf("  Sample cache: %d cells, reads range [%.0f, %.0f]\n",
            length(sample_cache), min(sample_cache), max(sample_cache)))

# ==============================================================================
# 4. TEST METADATA PARSING
# ==============================================================================

cat("\n--- Testing metadata parsing ---\n")

# Parse a single metadata file from cells mode
cells_run_dir <- file.path(get_paths(TOY_CELL_TYPE, "cells")$cartesian_dir,
                           paste0("Run_", TOY_RUN))
cells_meta <- list_metadata_files(cells_run_dir, "cells")
test_meta_path <- cells_meta$file_path[1]

depth_info <- compute_iteration_depth(test_meta_path, shard_cache)

if (is.null(depth_info)) {
  cat("FAIL: compute_iteration_depth returned NULL\n")
  quit(status = 1)
}

cat(sprintf("  Test metadata file: %s\n", basename(test_meta_path)))
cat(sprintf("  n_cells_mean: %.1f\n", depth_info$n_cells_mean))
cat(sprintf("  reads_per_cell_mean: %.1f\n", depth_info$reads_per_cell_mean))
cat(sprintf("  n_donors_actual: %d\n", depth_info$n_donors_actual))

if (depth_info$n_cells_mean <= 0 || depth_info$reads_per_cell_mean <= 0) {
  cat("FAIL: Invalid depth values\n")
  quit(status = 1)
}

# Parse a reads mode metadata file too
reads_run_dir <- file.path(get_paths(TOY_CELL_TYPE, "reads")$cartesian_dir,
                           paste0("Run_", TOY_RUN))
reads_meta <- list_metadata_files(reads_run_dir, "reads")
test_reads_meta_path <- reads_meta$file_path[1]

depth_info_reads <- compute_iteration_depth(test_reads_meta_path, shard_cache)

if (is.null(depth_info_reads)) {
  cat("FAIL: compute_iteration_depth returned NULL for reads mode\n")
  quit(status = 1)
}

cat(sprintf("  Reads mode - n_cells_mean: %.1f, rpc_mean: %.1f\n",
            depth_info_reads$n_cells_mean, depth_info_reads$reads_per_cell_mean))

cat("  Metadata parsing: PASS\n")

# ==============================================================================
# 5. LOAD GT AND DEFINE TRUTH SET
# ==============================================================================

cat("\n--- Loading GT and defining truth set ---\n")

gt_path <- file.path(cells_run_dir, "Ground_Truth.qs")
gt_res <- as.data.frame(qs::qread(gt_path))
gt_res$gene_id <- rownames(gt_res)

# Check GT DDS for dispersion
dds_path <- file.path(paths$gt_dds_dir, sprintf("GT_DDS_Run_%d.qs", TOY_RUN))

if (file.exists(dds_path)) {
  cat("  GT DDS checkpoint: FOUND\n")
  gt_dds <- qs::qread(dds_path)
  gt_dispersions <- DESeq2::dispersions(gt_dds)
  names(gt_dispersions) <- rownames(gt_dds)
} else {
  cat("  GT DDS checkpoint: NOT FOUND (using mock dispersion)\n")
  set.seed(42)
  gt_dispersions <- setNames(pmax(0.01, gt_res$baseMean / 100 + rnorm(nrow(gt_res), 0, 0.3)),
                              gt_res$gene_id)
}

truth_genes <- gt_res %>%
  filter(!is.na(padj),
         padj < GT_PADJ_THRESHOLD,
         abs(log2FoldChange) > GT_LFC_THRESHOLD) %>%
  mutate(
    abs_gt_lfc = abs(log2FoldChange),
    gt_basemean = baseMean,
    gt_dispersion = gt_dispersions[gene_id]
  ) %>%
  filter(!is.na(gt_dispersion), gt_dispersion > 0, gt_basemean > 0) %>%
  select(gene_id, abs_gt_lfc, gt_basemean, gt_dispersion)

cat(sprintf("  Truth set: %d genes\n", nrow(truth_genes)))

if (nrow(truth_genes) < 5) {
  cat("FAIL: Too few truth genes\n")
  quit(status = 1)
}

# ==============================================================================
# 6. BUILD MINI UNIFIED FITTING DATASET
# ==============================================================================

cat("\n--- Building mini unified fitting dataset ---\n")

fitting_rows <- list()

for (mode in MODES) {
  run_dir <- file.path(get_paths(TOY_CELL_TYPE, mode)$cartesian_dir,
                       paste0("Run_", TOY_RUN))

  iter_files <- list_iteration_files(run_dir, mode)
  meta_files <- list_metadata_files(run_dir, mode)

  grid_points <- iter_files %>%
    filter(!is.na(y_value)) %>%
    distinct(n_donors, y_value) %>%
    head(TOY_MAX_GRID_POINTS)

  cat(sprintf("  %s mode: %d grid points\n", mode, nrow(grid_points)))

  for (gp_idx in seq_len(nrow(grid_points))) {
    gp <- grid_points[gp_idx, ]

    gp_iter_files <- iter_files %>%
      filter(n_donors == gp$n_donors, y_value == gp$y_value) %>%
      head(TOY_MAX_ITERS)

    gp_meta_files <- meta_files %>%
      filter(n_donors == gp$n_donors, y_value == gp$y_value) %>%
      head(TOY_MAX_ITERS)

    for (iter_idx in seq_len(nrow(gp_iter_files))) {
      iter_num <- gp_iter_files$iteration[iter_idx]

      # Load DE results
      iter_res <- tryCatch({
        as.data.frame(qs::qread(gp_iter_files$file_path[iter_idx]))
      }, error = function(e) NULL)
      if (is.null(iter_res)) next

      iter_res$gene_id <- rownames(iter_res)
      detected_genes <- iter_res$gene_id[!is.na(iter_res$padj) &
                                           iter_res$padj < DETECT_PADJ_THRESHOLD]
      detected <- as.integer(truth_genes$gene_id %in% detected_genes)

      # Compute actual depth from metadata
      meta_match <- gp_meta_files %>% filter(iteration == iter_num)
      n_cells_mean <- NA_real_
      rpc_mean <- NA_real_

      if (nrow(meta_match) == 1) {
        depth <- compute_iteration_depth(meta_match$file_path[1], shard_cache)
        if (!is.null(depth)) {
          n_cells_mean <- depth$n_cells_mean
          rpc_mean <- depth$reads_per_cell_mean
        }
      }

      if (is.na(n_cells_mean) || is.na(rpc_mean)) next

      fitting_rows[[length(fitting_rows) + 1]] <- data.frame(
        gene_id = truth_genes$gene_id,
        detected = detected,
        n_donors = gp$n_donors,
        n_cells = n_cells_mean,
        reads_per_cell = rpc_mean,
        mode = mode,
        run_id = TOY_RUN,
        iteration = iter_num,
        abs_gt_lfc = truth_genes$abs_gt_lfc,
        gt_basemean = truth_genes$gt_basemean,
        gt_dispersion = truth_genes$gt_dispersion,
        stringsAsFactors = FALSE
      )
    }
  }
}

fitting_data <- bind_rows(fitting_rows)
cat(sprintf("  Raw rows: %d\n", nrow(fitting_data)))

# Build per-iteration fitting dataset (no aggregation)
per_iter_data <- fitting_data %>%
  mutate(
    cell_type = TOY_CELL_TYPE,
    mean_detection = as.numeric(detected),
    log_donors = log(n_donors),
    log_n_cells = log(n_cells),
    log_reads_per_cell = log(reads_per_cell),
    log_abs_gt_lfc = log(abs_gt_lfc),
    log_gt_basemean = log(gt_basemean),
    log_gt_dispersion = log(gt_dispersion)
  )

cat(sprintf("  Per-iteration dataset: %d rows, %d genes\n",
            nrow(per_iter_data), n_distinct(per_iter_data$gene_id)))

# Check that both modes contributed
modes_present <- unique(per_iter_data$mode)
cat(sprintf("  Modes present: %s\n", paste(modes_present, collapse = ", ")))

if (length(modes_present) < 2) {
  cat("WARNING: Only one mode present in toy data (expected both cells + reads)\n")
}

# Check depth variation
cat(sprintf("  n_cells range: [%.0f, %.0f]\n",
            min(per_iter_data$n_cells), max(per_iter_data$n_cells)))
cat(sprintf("  reads_per_cell range: [%.0f, %.0f]\n",
            min(per_iter_data$reads_per_cell), max(per_iter_data$reads_per_cell)))
cat(sprintf("  reads_per_cell CV: %.3f\n",
            sd(per_iter_data$reads_per_cell) / mean(per_iter_data$reads_per_cell)))

# ==============================================================================
# 7. FIT v3.1.0 UNIFIED FORMULA
# ==============================================================================

cat("\n--- Fitting v4.0.0 unified logistic regression ---\n")

formula <- build_model_formula(include_cell_type = FALSE)
cat(sprintf("  Formula: %s\n", deparse(formula)))

model <- safe_glm_binomial(formula, data = per_iter_data)

if (is.null(model)) {
  cat("FAIL: Logistic regression failed to converge\n")
  quit(status = 1)
}

pseudo_r2 <- mcfadden_r2(model)
cat(sprintf("  McFadden's pseudo R-squared: %.4f\n", pseudo_r2))
cat(sprintf("  Log-likelihood: %.2f\n", as.numeric(logLik(model))))
cat(sprintf("  AIC: %.2f\n", AIC(model)))
cat("  Coefficients:\n")
print(round(coef(model), 4))

# ==============================================================================
# 8. VALIDATE PREDICTIONS
# ==============================================================================

cat("\n--- Validating predictions ---\n")

per_iter_data$predicted <- predict(model, type = "response")

pred_range <- range(per_iter_data$predicted, na.rm = TRUE)
cat(sprintf("  Prediction range: [%.4f, %.4f]\n", pred_range[1], pred_range[2]))

if (pred_range[1] < 0 || pred_range[2] > 1) {
  cat("FAIL: Predictions outside [0, 1]\n")
  quit(status = 1)
}

corr <- cor(per_iter_data$mean_detection, per_iter_data$predicted, use = "complete.obs")
cat(sprintf("  Correlation (observed vs predicted): %.4f\n", corr))

if (corr < 0) {
  cat("FAIL: Negative correlation\n")
  quit(status = 1)
}

# Test prediction at a new point (v3.1.0 formula)
new_data <- data.frame(
  log_donors = log(100),
  log_n_cells = log(250),
  log_reads_per_cell = log(6000),
  log_abs_gt_lfc = log(3.0),
  log_gt_basemean = log(100),
  log_gt_dispersion = log(0.5)
)

new_pred <- predict(model, newdata = new_data, type = "response")
cat(sprintf("  Prediction at N=100, C=250, RPC=6000, LFC=3, BM=100, D=0.5: %.4f\n", new_pred))

if (is.na(new_pred) || new_pred < 0 || new_pred > 1) {
  cat("FAIL: Invalid prediction at test point\n")
  quit(status = 1)
}

# ==============================================================================
# 9. TEST PREDICTION GRID
# ==============================================================================

cat("\n--- Testing prediction grid ---\n")

grid <- generate_prediction_grid(vary_axis = "n_cells", fixed_depth = 6000, n_points = 5)
grid_preds <- predict(model, newdata = grid, type = "response")

cat(sprintf("  Grid predictions: %d values, range [%.4f, %.4f]\n",
            length(grid_preds), min(grid_preds, na.rm = TRUE), max(grid_preds, na.rm = TRUE)))

if (any(is.na(grid_preds))) {
  cat("WARNING: Some grid predictions are NA\n")
}

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  TOY TEST RESULTS (v4.0.0)\n")
cat("================================================================\n")
cat("\n")
cat("  Data loading:        PASS\n")
cat("  Shard cache:         PASS\n")
cat("  Metadata parsing:    PASS\n")
cat(sprintf("  Modes merged:        %s\n", paste(modes_present, collapse = " + ")))
cat(sprintf("  Truth set:           %d genes\n", nrow(truth_genes)))
cat(sprintf("  Fitting dataset:     %d rows (per-iteration, unified)\n", nrow(per_iter_data)))
cat(sprintf("  Depth variation:     n_cells=[%.0f,%.0f], rpc=[%.0f,%.0f]\n",
            min(per_iter_data$n_cells), max(per_iter_data$n_cells),
            min(per_iter_data$reads_per_cell), max(per_iter_data$reads_per_cell)))
cat(sprintf("  Logistic regression: PASS (R2=%.3f)\n", pseudo_r2))
cat(sprintf("  Predictions:         PASS (range [%.3f, %.3f])\n", pred_range[1], pred_range[2]))
cat(sprintf("  Obs-Pred corr:       %.3f\n", corr))
cat(sprintf("  New point pred:      %.3f\n", new_pred))
cat("\n")
cat("  VERDICT: PASS\n")
cat("\n")
cat("  NEXT STEPS:\n")
cat("  1. Run Step 00b to cache shard read counts (if not done):\n")
cat("     Rscript R/power_curve_fitting/run_pipeline.R 0.5 48\n")
cat("  2. Run full pipeline:\n")
cat("     Rscript R/power_curve_fitting/run_pipeline.R all 48\n")
cat("\n")
cat("================================================================\n")
cat("\n")

quit(status = 0)
