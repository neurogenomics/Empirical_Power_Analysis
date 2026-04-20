# ==============================================================================
# SCRIPT: Cross-Validate Modes (v4.0.0)
# PURPOSE: Validate that the unified model (trained on merged cells+reads data)
#          predicts well for each mode subset independently. Also compare
#          per-CT model predictions against the combined (no cell_type) model.
# INPUT: Fitted models from Step 2, fitting data from Step 1
# OUTPUT: Cross-validation plots and summary statistics
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(scales)
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

# Source shared prettify functions
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 4: Cross-Validate Modes (v4.0.0)")

cross_val_dir <- file.path(OUTPUT_BASE, "Plots", "Cross_Validation")
safe_mkdir(cross_val_dir)

# ==============================================================================
# 2. LOAD MODELS AND DATA (baseline: full_per_iteration)
# ==============================================================================

baseline_variant <- "full_per_iteration"
model_dir <- get_model_dir("full", "per_iteration")
checkpoint_path <- file.path(model_dir, "all_models_checkpoint.qs")

if (!file.exists(checkpoint_path)) {
  log_error("Model checkpoint not found: %s", checkpoint_path)
  stop("Missing model checkpoint for cross-validation")
}

checkpoint <- qs::qread(checkpoint_path)
per_ct_models <- checkpoint$per_cell_type
combined_model <- checkpoint$combined

# Load fitting data
all_data <- list()
for (ct in CELL_TYPES) {
  paths <- get_paths(ct)
  data_path <- file.path(paths$fitting_dir, "unified_fitting_data_per_iteration.qs")
  if (file.exists(data_path)) {
    all_data[[ct]] <- qs::qread(data_path)
  }
}
df_all <- bind_rows(all_data)

log_info("Loaded %d rows of fitting data across %d cell types",
         nrow(df_all), length(all_data))

# ==============================================================================
# 3. MODE-SPECIFIC VALIDATION
# ==============================================================================

# For each cell type, compare model predictions on cells-mode vs reads-mode subsets

cross_val_results <- list()

for (ct in CELL_TYPES) {

  log_info("=== Cross-validating: %s ===", ct)

  model <- per_ct_models[[ct]]
  if (is.null(model)) {
    log_warning("Missing model for %s, skipping", ct)
    next
  }

  df_ct <- df_all %>% filter(cell_type == ct)

  for (mode in MODES) {
    df_mode <- df_ct %>% filter(mode == !!mode)

    if (nrow(df_mode) == 0) {
      log_warning("No %s-mode data for %s", mode, ct)
      next
    }

    # Filter to valid rows
    df_mode <- df_mode %>%
      filter(!is.na(mean_detection))

    df_mode$predicted <- predict(model, newdata = df_mode, type = "response")

    corr <- cor(df_mode$mean_detection, df_mode$predicted, use = "complete.obs")
    rmse <- sqrt(mean((df_mode$mean_detection - df_mode$predicted)^2, na.rm = TRUE))
    mae <- mean(abs(df_mode$mean_detection - df_mode$predicted), na.rm = TRUE)

    log_info("  %s/%s: n=%d, corr=%.4f, RMSE=%.4f, MAE=%.4f",
             ct, mode, nrow(df_mode), corr, rmse, mae)

    cross_val_results[[paste(ct, mode, sep = "_")]] <- data.frame(
      cell_type = ct,
      mode = mode,
      n_obs = nrow(df_mode),
      correlation = corr,
      rmse = rmse,
      mae = mae,
      stringsAsFactors = FALSE
    )
  }
}

# ==============================================================================
# 4. COMBINED MODEL VALIDATION
# ==============================================================================

if (!is.null(combined_model)) {
  log_info("=== Validating combined (cell-type-agnostic) model ===")

  df_valid <- df_all %>%
    filter(!is.na(mean_detection))

  df_valid$pred_combined <- predict(combined_model, newdata = df_valid, type = "response")

  # Per-CT predictions from per-CT models
  df_valid$pred_per_ct <- NA_real_
  for (ct in CELL_TYPES) {
    model <- per_ct_models[[ct]]
    if (is.null(model)) next
    idx <- df_valid$cell_type == ct
    df_valid$pred_per_ct[idx] <- predict(model, newdata = df_valid[idx, ], type = "response")
  }

  # Compare combined vs per-CT predictions
  valid_both <- df_valid %>% filter(!is.na(pred_combined), !is.na(pred_per_ct))

  if (nrow(valid_both) > 0) {
    corr_models <- cor(valid_both$pred_combined, valid_both$pred_per_ct, use = "complete.obs")
    rmse_models <- sqrt(mean((valid_both$pred_combined - valid_both$pred_per_ct)^2, na.rm = TRUE))
    mae_models <- mean(abs(valid_both$pred_combined - valid_both$pred_per_ct), na.rm = TRUE)

    log_info("Combined vs per-CT: corr=%.4f, RMSE=%.4f, MAE=%.4f",
             corr_models, rmse_models, mae_models)

    # Scatter: combined vs per-CT predictions
    p_scatter <- ggplot(valid_both, aes(x = pred_per_ct, y = pred_combined, colour = cell_type)) +
      geom_point(alpha = 0.2, size = 0.8) +
      geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 1) +
      scale_colour_manual(values = CELL_TYPE_COLOURS, labels = prettify_labeller()) +
      scale_x_continuous(limits = c(0, 1), labels = percent) +
      scale_y_continuous(limits = c(0, 1), labels = percent) +
      labs(
        title = "Combined vs Per-CT Model Predictions",
        subtitle = sprintf("r=%.4f, RMSE=%.4f — Combined model has no cell_type term",
                           corr_models, rmse_models),
        x = "Per-CT Model Prediction",
        y = "Combined Model Prediction",
        colour = "Cell Type"
      ) +
      theme_publication(base_size = 16) +
      coord_equal()

    ggsave(file.path(cross_val_dir, "Combined_vs_PerCT.png"),
           p_scatter, width = PLOT_HEIGHT, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    # [REMOVED] Residual-by-mode plot removed in v4.0.0 cleanup.
    # Binary response makes raw residual scatter meaningless for logistic
    # regression (two deterministic diagonal lines). Replaced by summary
    # table in cross_validation_summary.csv.

    # Per-CT combined-vs-per-CT agreement metrics
    for (ct in CELL_TYPES) {
      df_ct <- valid_both %>% filter(cell_type == ct)
      if (nrow(df_ct) == 0) next

      ct_corr <- cor(df_ct$pred_combined, df_ct$pred_per_ct, use = "complete.obs")
      ct_rmse <- sqrt(mean((df_ct$pred_combined - df_ct$pred_per_ct)^2, na.rm = TRUE))
      ct_mae <- mean(abs(df_ct$pred_combined - df_ct$pred_per_ct), na.rm = TRUE)

      log_info("  Combined vs per-CT (%s): n=%d, r=%.4f, RMSE=%.4f, MAE=%.4f",
               ct, nrow(df_ct), ct_corr, ct_rmse, ct_mae)

      cross_val_results[[paste0("combined_vs_perCT_", ct)]] <- data.frame(
        cell_type = ct,
        mode = "combined_vs_perCT",
        n_obs = nrow(df_ct),
        correlation = ct_corr,
        rmse = ct_rmse,
        mae = ct_mae,
        stringsAsFactors = FALSE
      )
    }
  }
}

# ==============================================================================
# 5. SAVE SUMMARY
# ==============================================================================

if (length(cross_val_results) > 0) {
  summary_df <- bind_rows(cross_val_results)

  # Round to 3 decimal places for publication
  summary_df$correlation <- round(summary_df$correlation, 3)
  summary_df$rmse <- round(summary_df$rmse, 3)
  summary_df$mae <- round(summary_df$mae, 3)

  write.csv(summary_df, file.path(cross_val_dir, "cross_validation_summary.csv"),
            row.names = FALSE)
  log_info("Cross-validation summary saved: %d rows", nrow(summary_df))

  # Print summary table
  log_info("=== CROSS-VALIDATION SUMMARY ===")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    log_info("  %s/%s: n=%d, r=%.3f, RMSE=%.3f, MAE=%.3f",
             row$cell_type, row$mode, row$n_obs, row$correlation, row$rmse, row$mae)
  }
}

log_section("STEP 4 COMPLETE")
log_info("Cross-validation of unified model complete")
