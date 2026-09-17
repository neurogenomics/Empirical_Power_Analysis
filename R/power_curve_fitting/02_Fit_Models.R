# ==============================================================================
# SCRIPT: Fit Logistic Regression Models (v4.1.0)
# PURPOSE: Fit unified power equations using the v4.0.0 formula:
#
#   detection ~ log_donors + log_n_cells + log_reads_per_cell +
#               log_abs_gt_lfc + log_gt_basemean + log_gt_dispersion
#
# v4.1.0: Loops over DATASET_CONFIGS to fit models for all datasets.
#          Each dataset's per-CT models are saved under the model name
#          (e.g., Astro_Shallow_model.qs for DepthEq).
#          Combined/additive/interaction models are only fitted for datasets
#          with fit_combined=TRUE (requires 2+ cell types).
#          All coefficients are collected into a single all_coefficients.csv.
#
# INPUT: Unified fitting datasets from Step 1
# OUTPUT: Model objects, coefficient tables, comparison table
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(lmtest)
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

log_section("STEP 2: Fit Logistic Regression Models (v4.1.0 — multi-dataset)")

# ==============================================================================
# 2. MAIN LOOP: gene_set x agg_level
# ==============================================================================

# Accumulate comparison rows across all datasets
comparison_rows <- list()

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {

    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    log_section(sprintf("Variant: %s", variant_label))

    model_dir <- get_model_dir(gene_set, agg_level)
    safe_mkdir(model_dir)

    # ==================================================================
    # Collect ALL per-CT models and data across ALL datasets
    # ==================================================================

    # These accumulate across datasets for the final coefficient CSV
    all_per_ct_models <- list()
    all_per_ct_coefs <- list()
    all_per_ct_confints <- list()

    # For combined models: only use data from datasets with fit_combined=TRUE
    combined_data_list <- list()

    for (ds_config in DATASET_CONFIGS) {

      log_section(sprintf("Dataset: %s (variant: %s)", ds_config$dataset, variant_label))

      # ------------------------------------------------------------------
      # 2a. Load fitting data for this dataset's cell types
      # ------------------------------------------------------------------

      ds_data <- list()

      for (cell_type in ds_config$cell_types) {
        model_name <- get_dataset_model_name(cell_type, ds_config)
        paths <- get_paths_for_dataset(cell_type, ds_config)

        # Select correct dataset file
        shared_suffix <- if (gene_set == "shared") "_shared" else ""
        data_path <- file.path(paths$fitting_dir,
                                sprintf("unified_fitting_data_%s%s.qs",
                                        agg_level, shared_suffix))

        if (!file.exists(data_path)) {
          log_warning("Fitting data not found for %s/%s: %s",
                      model_name, variant_label, basename(data_path))
          next
        }

        df <- qs::qread(data_path)
        log_info("Loaded %s: %d rows, %d genes",
                 model_name, nrow(df), n_distinct(df$gene_id))

        # Ensure cell_type column uses model_name
        df$cell_type <- model_name

        ds_data[[model_name]] <- df
      }

      if (length(ds_data) == 0) {
        log_warning("No fitting data available for %s/%s", ds_config$dataset, variant_label)
        next
      }

      # ------------------------------------------------------------------
      # 2b. Fit per-cell-type models for this dataset
      # ------------------------------------------------------------------

      log_info("Fitting per-cell-type models for %s...", ds_config$dataset)

      formula_per_ct <- build_model_formula(include_cell_type = FALSE)

      for (model_name in names(ds_data)) {
        df_ct <- ds_data[[model_name]]

        if (nrow(df_ct) < 50) {
          log_warning("Too few observations for %s (%d rows), skipping",
                      model_name, nrow(df_ct))
          next
        }

        # Verify required columns exist
        required_cols <- c("mean_detection", "log_donors", "log_n_cells",
                           "log_reads_per_cell", "log_abs_gt_lfc", "log_gt_basemean",
                           "log_gt_dispersion")
        missing_cols <- setdiff(required_cols, colnames(df_ct))
        if (length(missing_cols) > 0) {
          log_error("Missing columns in %s: %s", model_name,
                    paste(missing_cols, collapse = ", "))
          next
        }

        # Remove rows with NA/Inf in predictors
        n_before <- nrow(df_ct)
        df_ct <- df_ct %>%
          filter(
            !is.na(mean_detection),
            is.finite(log_donors),
            is.finite(log_n_cells),
            is.finite(log_reads_per_cell),
            is.finite(log_abs_gt_lfc),
            is.finite(log_gt_basemean),
            is.finite(log_gt_dispersion)
          )
        n_after <- nrow(df_ct)
        if (n_before != n_after) {
          log_info("Removed %d rows with NA/Inf values (%d -> %d) for %s",
                   n_before - n_after, n_before, n_after, model_name)
        }

        log_info("Fitting %s (n=%d)...", model_name, nrow(df_ct))

        model <- safe_glm_binomial(formula_per_ct, data = df_ct)

        if (!is.null(model)) {
          all_per_ct_models[[model_name]] <- model
          all_per_ct_coefs[[model_name]] <- format_coefficients(model, model_name)

          # Extract confidence intervals
          ci <- tryCatch({
            confint.default(model, level = CI_LEVEL)
          }, error = function(e) {
            log_warning("CI extraction failed for %s: %s", model_name, e$message)
            NULL
          })
          if (!is.null(ci)) {
            all_per_ct_confints[[model_name]] <- ci
          }

          # Print summary
          pseudo_r2 <- mcfadden_r2(model)
          log_info("  McFadden's pseudo R-squared: %.4f", pseudo_r2)
          log_info("  Log-likelihood: %.2f", as.numeric(logLik(model)))
          log_info("  AIC: %.2f", AIC(model))

          # Save individual model
          qs::qsave(model, file.path(model_dir, paste0(model_name, "_model.qs")))
        }
      }

      # ------------------------------------------------------------------
      # 2c. Collect data for combined models (only from fit_combined datasets)
      # ------------------------------------------------------------------

      if (ds_config$fit_combined) {
        for (model_name in names(ds_data)) {
          combined_data_list[[model_name]] <- ds_data[[model_name]]
        }
      }

    }  # end ds_config loop

    # ==================================================================
    # 3. Fit combined models (from fit_combined datasets only)
    # ==================================================================

    # Identify which model names came from fit_combined datasets
    combined_model_names <- names(combined_data_list)
    model_combined <- NULL
    model_additive <- NULL
    model_interaction <- NULL

    if (length(combined_data_list) >= 2) {

      df_combined <- bind_rows(combined_data_list)

      # For combined models, use the original cell type names (Astro/Micro/Oligo)
      # which are already set in the cell_type column from Step 1
      df_combined$cell_type <- factor(df_combined$cell_type, levels = combined_model_names)

      log_info("Combined dataset: %d rows across %d cell types (%s), %d unique genes",
               nrow(df_combined), n_distinct(df_combined$cell_type),
               paste(combined_model_names, collapse = ", "),
               n_distinct(df_combined$gene_id))

      # Remove rows with NA/Inf in predictors
      n_before <- nrow(df_combined)
      df_combined <- df_combined %>%
        filter(
          !is.na(mean_detection),
          is.finite(log_donors),
          is.finite(log_n_cells),
          is.finite(log_reads_per_cell),
          is.finite(log_abs_gt_lfc),
          is.finite(log_gt_basemean),
          is.finite(log_gt_dispersion)
        )
      n_after <- nrow(df_combined)
      if (n_before != n_after) {
        log_info("Removed %d rows with NA/Inf values (%d -> %d)",
                 n_before - n_after, n_before, n_after)
      }

      # ------------------------------------------------------------------
      # 3a. Combined model (NO cell_type term)
      # ------------------------------------------------------------------

      log_info("Fitting combined model (cell-type-agnostic)...")
      formula_combined <- build_model_formula(include_cell_type = FALSE)
      log_info("Combined formula: %s", deparse(formula_combined))

      model_combined <- safe_glm_binomial(formula_combined, data = df_combined)

      if (!is.null(model_combined)) {
        pseudo_r2 <- mcfadden_r2(model_combined)
        log_info("Combined model: R²=%.4f, loglik=%.2f, AIC=%.2f",
                 pseudo_r2, as.numeric(logLik(model_combined)), AIC(model_combined))
        qs::qsave(model_combined, file.path(model_dir, "combined_model.qs"))
      }

      # ------------------------------------------------------------------
      # 3b. Additive: cell_type + predictors
      # ------------------------------------------------------------------

      log_info("Fitting combined models with cell_type for LRT...")

      formula_additive <- build_model_formula(include_cell_type = TRUE, interactions = FALSE)
      log_info("Additive formula: %s", deparse(formula_additive))
      model_additive <- safe_glm_binomial(formula_additive, data = df_combined)

      if (!is.null(model_additive)) {
        pseudo_r2 <- mcfadden_r2(model_additive)
        log_info("Additive model: R²=%.4f, loglik=%.2f, AIC=%.2f",
                 pseudo_r2, as.numeric(logLik(model_additive)), AIC(model_additive))
        qs::qsave(model_additive, file.path(model_dir, "combined_additive_model.qs"))
      }

      # ------------------------------------------------------------------
      # 3c. Full interaction: cell_type * predictors
      # ------------------------------------------------------------------

      formula_interaction <- build_model_formula(include_cell_type = TRUE, interactions = TRUE)
      log_info("Interaction formula: %s", deparse(formula_interaction))
      model_interaction <- safe_glm_binomial(formula_interaction, data = df_combined)

      if (!is.null(model_interaction)) {
        pseudo_r2 <- mcfadden_r2(model_interaction)
        log_info("Interaction model: R²=%.4f, loglik=%.2f, AIC=%.2f",
                 pseudo_r2, as.numeric(logLik(model_interaction)), AIC(model_interaction))
        qs::qsave(model_interaction, file.path(model_dir, "combined_interaction_model.qs"))
      }

      # ------------------------------------------------------------------
      # 3d. Likelihood ratio tests
      # ------------------------------------------------------------------

      lrt1_pval <- NA_real_; lrt1_chisq <- NA_real_; lrt1_df <- NA_integer_
      lrt2_pval <- NA_real_; lrt2_chisq <- NA_real_; lrt2_df <- NA_integer_

      if (!is.null(model_combined) && !is.null(model_additive)) {
        lrt1 <- tryCatch({
          lmtest::lrtest(model_combined, model_additive)
        }, error = function(e) {
          log_error("LRT1 (combined vs additive) failed: %s", e$message)
          NULL
        })
        if (!is.null(lrt1)) {
          lrt1_pval <- lrt1$`Pr(>Chisq)`[2]
          lrt1_chisq <- lrt1$Chisq[2]
          lrt1_df <- lrt1$Df[2]
          log_info("LRT1 (no CT vs +CT): Chi-sq=%.2f, df=%d, p=%.2e",
                   lrt1_chisq, lrt1_df, lrt1_pval)
        }
      }

      if (!is.null(model_additive) && !is.null(model_interaction)) {
        lrt2 <- tryCatch({
          lmtest::lrtest(model_additive, model_interaction)
        }, error = function(e) {
          log_error("LRT2 (additive vs interaction) failed: %s", e$message)
          NULL
        })
        if (!is.null(lrt2)) {
          lrt2_pval <- lrt2$`Pr(>Chisq)`[2]
          lrt2_chisq <- lrt2$Chisq[2]
          lrt2_df <- lrt2$Df[2]
          log_info("LRT2 (+CT vs CT*pred): Chi-sq=%.2f, df=%d, p=%.2e",
                   lrt2_chisq, lrt2_df, lrt2_pval)
        }
      }

    } else {
      log_info("Fewer than 2 cell types with fit_combined=TRUE — skipping combined models")
      lrt1_pval <- NA_real_; lrt1_chisq <- NA_real_; lrt1_df <- NA_integer_
      lrt2_pval <- NA_real_; lrt2_chisq <- NA_real_; lrt2_df <- NA_integer_
    }

    # ==================================================================
    # 4. Assess coefficient convergence across per-CT models
    # ==================================================================

    ci_overlap_all <- NA
    max_coef_diff <- NA_real_

    if (length(all_per_ct_confints) >= 2 && length(all_per_ct_coefs) >= 2) {
      common_terms <- Reduce(intersect, lapply(all_per_ct_coefs, function(x) x$term))
      compare_terms <- setdiff(common_terms, c("(Intercept)", "(phi)"))

      all_overlap <- TRUE
      max_diff <- 0

      for (term in compare_terms) {
        ct_estimates <- sapply(names(all_per_ct_coefs), function(ct) {
          all_per_ct_coefs[[ct]]$estimate[all_per_ct_coefs[[ct]]$term == term]
        })

        if (length(ct_estimates) >= 2) {
          pairwise_diffs <- abs(outer(ct_estimates, ct_estimates, "-"))
          max_diff <- max(max_diff, max(pairwise_diffs))
        }

        for (ct in names(all_per_ct_confints)) {
          ci <- all_per_ct_confints[[ct]]
          if (term %in% rownames(ci)) {
            ci_low <- ci[term, 1]
            ci_high <- ci[term, 2]

            for (ct2 in names(all_per_ct_confints)) {
              if (ct2 <= ct) next
              ci2 <- all_per_ct_confints[[ct2]]
              if (term %in% rownames(ci2)) {
                ci2_low <- ci2[term, 1]
                ci2_high <- ci2[term, 2]
                if (max(ci_low, ci2_low) > min(ci_high, ci2_high)) {
                  all_overlap <- FALSE
                }
              }
            }
          }
        }
      }

      ci_overlap_all <- all_overlap
      max_coef_diff <- max_diff
    }

    # ==================================================================
    # 5. Build comparison row
    # ==================================================================

    # Per-CT R² for all models (including Astro_Shallow etc.)
    all_model_names <- names(all_per_ct_models)
    per_ct_r2 <- sapply(all_model_names, function(mn) {
      mcfadden_r2(all_per_ct_models[[mn]])
    })

    combined_r2 <- mcfadden_r2(model_combined)
    additive_r2 <- mcfadden_r2(model_additive)
    interaction_r2 <- mcfadden_r2(model_interaction)

    # Build a comparison row with dynamic per-CT R² columns
    comparison_row <- data.frame(
      gene_set = gene_set,
      agg_level = agg_level,
      variant = variant_label,
      n_models = length(all_per_ct_models),
      r2_combined = combined_r2,
      r2_additive = additive_r2,
      r2_interaction = interaction_r2,
      lrt1_chisq = lrt1_chisq,
      lrt1_df = lrt1_df,
      lrt1_pval = lrt1_pval,
      lrt1_significant = !is.na(lrt1_pval) && lrt1_pval < LRT_ALPHA,
      lrt2_chisq = lrt2_chisq,
      lrt2_df = lrt2_df,
      lrt2_pval = lrt2_pval,
      lrt2_significant = !is.na(lrt2_pval) && lrt2_pval < LRT_ALPHA,
      ci_overlap_all = ci_overlap_all,
      max_coef_diff = max_coef_diff,
      stringsAsFactors = FALSE
    )

    # Add per-CT R² columns dynamically
    for (mn in all_model_names) {
      comparison_row[[paste0("r2_", mn)]] <- per_ct_r2[mn]
    }

    comparison_rows[[variant_label]] <- comparison_row

    # ==================================================================
    # 6. Save all results for this variant
    # ==================================================================

    # Save coefficient tables (all per-CT + combined models)
    all_coefs <- bind_rows(all_per_ct_coefs)
    if (!is.null(model_combined)) {
      all_coefs <- bind_rows(all_coefs, format_coefficients(model_combined, "combined"))
    }
    if (!is.null(model_additive)) {
      all_coefs <- bind_rows(all_coefs, format_coefficients(model_additive, "combined_additive"))
    }
    if (!is.null(model_interaction)) {
      all_coefs <- bind_rows(all_coefs, format_coefficients(model_interaction, "combined_interaction"))
    }
    write.csv(all_coefs, file.path(model_dir, "all_coefficients.csv"), row.names = FALSE)

    # Save model comparison table
    model_comparison <- data.frame(
      model = character(), n_obs = integer(), n_params = integer(),
      loglik = numeric(), aic = numeric(), bic = numeric(),
      pseudo_r2 = numeric(), stringsAsFactors = FALSE
    )

    for (mn in names(all_per_ct_models)) {
      m <- all_per_ct_models[[mn]]
      model_comparison <- rbind(model_comparison, data.frame(
        model = mn, n_obs = nobs(m), n_params = length(coef(m)),
        loglik = as.numeric(logLik(m)), aic = AIC(m), bic = BIC(m),
        pseudo_r2 = mcfadden_r2(m), stringsAsFactors = FALSE
      ))
    }

    for (name in c("combined", "combined_additive", "combined_interaction")) {
      m <- switch(name,
                  "combined" = model_combined,
                  "combined_additive" = model_additive,
                  "combined_interaction" = model_interaction)
      if (!is.null(m)) {
        model_comparison <- rbind(model_comparison, data.frame(
          model = name, n_obs = nobs(m), n_params = length(coef(m)),
          loglik = as.numeric(logLik(m)), aic = AIC(m), bic = BIC(m),
          pseudo_r2 = mcfadden_r2(m), stringsAsFactors = FALSE
        ))
      }
    }

    write.csv(model_comparison, file.path(model_dir, "model_comparison.csv"), row.names = FALSE)

    # Save all models as a single checkpoint
    qs::qsave(list(
      per_cell_type = all_per_ct_models,
      combined = model_combined,
      combined_additive = model_additive,
      combined_interaction = model_interaction,
      model_comparison = model_comparison,
      coefficients = all_coefs,
      confints = all_per_ct_confints,
      gene_set = gene_set,
      agg_level = agg_level,
      formula_per_ct = build_model_formula(include_cell_type = FALSE),
      formula_combined = build_model_formula(include_cell_type = FALSE),
      formula_additive = build_model_formula(include_cell_type = TRUE, interactions = FALSE),
      formula_interaction = build_model_formula(include_cell_type = TRUE, interactions = TRUE)
    ), file.path(model_dir, "all_models_checkpoint.qs"))

    log_info("Saved model checkpoint for variant: %s", variant_label)
  }
}

# ==============================================================================
# 7. SAVE MODEL COMPARISON TABLE
# ==============================================================================

log_section("Saving model comparison summary")

if (length(comparison_rows) > 0) {
  comparison_df <- bind_rows(comparison_rows)

  comparison_path <- file.path(OUTPUT_BASE, "Models", "model_comparison_v3.csv")
  safe_mkdir(dirname(comparison_path))
  write.csv(comparison_df, comparison_path, row.names = FALSE)
  log_info("Saved model comparison: %d variants", nrow(comparison_df))

  # Print summary
  log_info("=== MODEL COMPARISON SUMMARY ===")
  for (i in seq_len(nrow(comparison_df))) {
    row <- comparison_df[i, ]
    # Print per-CT R² values dynamically
    r2_cols <- grep("^r2_", names(row), value = TRUE)
    r2_cols <- setdiff(r2_cols, c("r2_combined", "r2_additive", "r2_interaction"))
    r2_strs <- sapply(r2_cols, function(col) {
      sprintf("%s=%.4f", sub("^r2_", "", col), row[[col]])
    })
    log_info("  %s: R²(combined)=%.4f, %s",
             row$variant, row$r2_combined, paste(r2_strs, collapse = ", "))
    if (!is.na(row$lrt1_pval)) {
      log_info("    LRT1 (no CT vs +CT): p=%.2e %s",
               row$lrt1_pval,
               ifelse(row$lrt1_significant, "[SIGNIFICANT]", "[not significant]"))
    }
    if (!is.na(row$lrt2_pval)) {
      log_info("    LRT2 (+CT vs CT*pred): p=%.2e %s",
               row$lrt2_pval,
               ifelse(row$lrt2_significant, "[SIGNIFICANT]", "[not significant]"))
    }
  }
}

log_section("STEP 2 COMPLETE")
log_info("Logistic regression models fitted for all datasets and variants")
