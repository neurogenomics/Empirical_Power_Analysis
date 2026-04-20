# ==============================================================================
# SCRIPT: Systematic Ordering Permutation Test (v1.0.0)
# PURPOSE: Test all 5^5 = 3,125 combinations of covariate quantiles to determine
#          how often the cell type power hierarchy (Micro>Astro>Oligo) holds vs
#          when it changes. Evaluates at 100 donors (representative mid-range).
# INPUT: Fitted models from power curve fitting
# OUTPUT: Ordering frequency table and edge case analysis in Summary directory
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
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
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("Systematic Ordering Permutation Test (v1.0.0)")

OUTPUT_DIR <- file.path(OUTPUT_BASE, "Summary")
safe_mkdir(OUTPUT_DIR)

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

ALL_CT_MODELS <- c("Astro", "Micro", "Oligo", "Astro_Shallow")
CT_COLOURS_4 <- c(Astro = "#E41A1C", Micro = "#377EB8", Oligo = "#4DAF4A", Astro_Shallow = "#E69F00")

# Evaluation point: 100 donors (representative mid-range)
EVAL_DONORS <- 100

# Quantile levels to test
QUANTILE_LEVELS <- c(0.05, 0.25, 0.50, 0.75, 0.95)
QUANTILE_NAMES <- c("Q5", "Q25", "Q50", "Q75", "Q95")

# Covariates to vary
COVARIATES <- c("n_cells", "reads_per_cell", "abs_gt_lfc", "gt_basemean", "gt_dispersion")

log_info("Testing %d^%d = %d combinations at %d donors", 
         length(QUANTILE_LEVELS), length(COVARIATES), 
         length(QUANTILE_LEVELS)^length(COVARIATES), EVAL_DONORS)

# ==============================================================================
# 3. LOAD MODELS AND COMPUTE QUANTILES
# ==============================================================================

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {
    
    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    log_section(sprintf("Testing variant: %s", variant_label))
    
    model_dir <- get_model_dir(gene_set, agg_level)
    
    # Load models
    checkpoint_path <- file.path(model_dir, "all_models_checkpoint.qs")
    if (!file.exists(checkpoint_path)) {
      log_warning("Model checkpoint not found for variant %s, skipping", variant_label)
      next
    }
    
    checkpoint <- qs::qread(checkpoint_path)
    per_ct_models <- checkpoint$per_cell_type
    
    # Check all models available
    missing_models <- setdiff(ALL_CT_MODELS, names(per_ct_models))
    if (length(missing_models) > 0) {
      log_warning("Missing models for variant %s: %s", variant_label, paste(missing_models, collapse = ", "))
      next
    }
    
    # Load fitting data to compute quantiles
    all_data_4 <- list()
    for (ct_model in ALL_CT_MODELS) {
      paths <- get_paths_for_dataset("Astro", list(dataset = "Roussos", model_names = NULL, fit_combined = TRUE))
      if (ct_model == "Astro_Shallow") {
        paths <- get_paths_for_dataset("Astro", list(dataset = "Roussos_DepthEq", model_names = c(Astro = "Astro_Shallow"), fit_combined = FALSE))
      } else if (ct_model != "Astro") {
        paths <- get_paths_for_dataset(ct_model, list(dataset = "Roussos", model_names = NULL, fit_combined = TRUE))
      }
      
      shared_suffix <- if (gene_set == "shared") "_shared" else ""
      data_path <- file.path(paths$fitting_dir,
                            sprintf("unified_fitting_data_%s%s.qs", agg_level, shared_suffix))
      if (file.exists(data_path)) {
        all_data_4[[ct_model]] <- qs::qread(data_path)
      }
    }
    df_4_combined <- bind_rows(all_data_4)
    
    if (nrow(df_4_combined) == 0) {
      log_warning("No fitting data for variant %s", variant_label)
      next
    }
    
    # Compute quantiles
    quantile_values <- list(
      n_cells = quantile(df_4_combined$n_cells, QUANTILE_LEVELS, na.rm = TRUE),
      reads_per_cell = quantile(df_4_combined$reads_per_cell, QUANTILE_LEVELS, na.rm = TRUE),
      abs_gt_lfc = quantile(df_4_combined$abs_gt_lfc, QUANTILE_LEVELS, na.rm = TRUE),
      gt_basemean = quantile(df_4_combined$gt_basemean, QUANTILE_LEVELS, na.rm = TRUE),
      gt_dispersion = quantile(df_4_combined$gt_dispersion, QUANTILE_LEVELS, na.rm = TRUE)
    )
    
    log_info("Quantile ranges:")
    for (cov in COVARIATES) {
      log_info("  %s: %.1f to %.1f", cov, min(quantile_values[[cov]]), max(quantile_values[[cov]]))
    }
    
    # ==============================================================================
    # 4. GENERATE ALL PERMUTATIONS AND TEST ORDERINGS
    # ==============================================================================
    
    log_info("Generating all %d permutations...", length(QUANTILE_LEVELS)^length(COVARIATES))
    
    # Create all combinations using expand.grid
    quantile_indices <- expand.grid(
      n_cells_idx = 1:5,
      reads_per_cell_idx = 1:5,
      abs_gt_lfc_idx = 1:5,
      gt_basemean_idx = 1:5,
      gt_dispersion_idx = 1:5
    )
    
    n_combinations <- nrow(quantile_indices)
    log_info("Testing %d combinations...", n_combinations)
    
    # Storage for results
    results <- list()
    
    # Progress tracking
    progress_interval <- max(1, floor(n_combinations / 20))
    
    for (i in 1:n_combinations) {
      if (i %% progress_interval == 0) {
        log_info("Progress: %d/%d (%.1f%%)", i, n_combinations, 100 * i / n_combinations)
      }
      
      idx <- quantile_indices[i, ]
      
      # Extract indices (idx is a 1-row data frame)
      n_cells_idx <- as.integer(idx[["n_cells_idx"]])
      rpc_idx <- as.integer(idx[["reads_per_cell_idx"]])
      lfc_idx <- as.integer(idx[["abs_gt_lfc_idx"]])
      basemean_idx <- as.integer(idx[["gt_basemean_idx"]])
      disp_idx <- as.integer(idx[["gt_dispersion_idx"]])
      
      # Build prediction data frame for this combination
      pred_df <- data.frame(
        n_donors = EVAL_DONORS,
        log_donors = log(EVAL_DONORS),
        log_n_cells = log(quantile_values$n_cells[n_cells_idx]),
        log_reads_per_cell = log(quantile_values$reads_per_cell[rpc_idx]),
        log_abs_gt_lfc = log(quantile_values$abs_gt_lfc[lfc_idx]),
        log_gt_basemean = log(quantile_values$gt_basemean[basemean_idx]),
        log_gt_dispersion = log(quantile_values$gt_dispersion[disp_idx])
      )
      
      # Predict for each model
      powers <- sapply(ALL_CT_MODELS, function(m) {
        pred_val <- tryCatch({
          p <- predict(per_ct_models[[m]], newdata = pred_df, type = "response")
          if (length(p) > 0) as.numeric(p[1]) else NA_real_
        }, error = function(e) NA_real_)
        return(pred_val)
      })
      
      # Skip if any NA predictions
      if (any(is.na(powers))) {
        if (i %% 100 == 0) {
          log_warning("Combination %d produced NA predictions, skipping", i)
        }
        next
      }
      
      # Determine ranking (highest power = rank 1)
      ranking <- rank(-powers, ties.method = "min")
      
      # Create ordering string (e.g., "Micro>Astro>Astro_Shallow>Oligo")
      model_order <- ALL_CT_MODELS[order(ranking)]
      ordering_str <- paste(model_order, collapse = " > ")
      
      # Store result
      results[[i]] <- data.frame(
        combination_id = i,
        n_cells_level = QUANTILE_NAMES[n_cells_idx],
        rpc_level = QUANTILE_NAMES[rpc_idx],
        lfc_level = QUANTILE_NAMES[lfc_idx],
        basemean_level = QUANTILE_NAMES[basemean_idx],
        dispersion_level = QUANTILE_NAMES[disp_idx],
        Astro_power = powers["Astro"],
        Micro_power = powers["Micro"],
        Oligo_power = powers["Oligo"],
        Astro_Shallow_power = powers["Astro_Shallow"],
        ordering = ordering_str,
        stringsAsFactors = FALSE
      )
    }
    
    results_df <- bind_rows(results)
    
    # ==============================================================================
    # 5. ANALYZE ORDERING FREQUENCIES
    # ==============================================================================
    
    log_section("Ordering Frequency Analysis")
    
    # Count frequencies
    ordering_counts <- results_df %>%
      group_by(ordering) %>%
      summarise(
        count = n(),
        proportion = n() / nrow(results_df),
        .groups = "drop"
      ) %>%
      arrange(desc(count))
    
    log_info("Unique orderings found: %d", nrow(ordering_counts))
    log_info("\nTop 10 most common orderings:")
    print(head(ordering_counts, 10))
    
    # Check if canonical ordering dominates
    canonical_ordering <- "Micro > Astro > Astro_Shallow > Oligo"
    canonical_count <- ordering_counts %>%
      filter(ordering == canonical_ordering) %>%
      pull(count)
    
    if (length(canonical_count) == 0) canonical_count <- 0
    
    log_info("\nCanonical ordering '%s': %d/%d (%.2f%%)", 
             canonical_ordering, canonical_count, n_combinations,
             100 * canonical_count / n_combinations)
    
    # Identify rare orderings (potential edge cases)
    rare_threshold <- 0.01  # Less than 1% of combinations
    rare_orderings <- ordering_counts %>%
      filter(proportion < rare_threshold) %>%
      arrange(proportion)
    
    if (nrow(rare_orderings) > 0) {
      log_info("\nRare orderings (<1%% of combinations) - potential edge cases:")
      print(rare_orderings)
    }
    
    # ==============================================================================
    # 6. EDGE CASE ANALYSIS: When does ordering change?
    # ==============================================================================
    
    log_section("Edge Case Analysis")
    
    # Find combinations where Oligo is not last (unusual)
    unusual_orderings <- results_df %>%
      filter(!grepl("Oligo$", ordering)) %>%
      select(combination_id, ordering, n_cells_level, rpc_level, lfc_level, 
             basemean_level, dispersion_level, Astro_power, Micro_power, Oligo_power)
    
    if (nrow(unusual_orderings) > 0) {
      log_info("Found %d combinations where Oligo is NOT last (unusual)", nrow(unusual_orderings))
      log_info("Sample of unusual cases:")
      print(head(unusual_orderings, 10))
    } else {
      log_info("Oligo is always last in all %d combinations - hierarchy is extremely stable", n_combinations)
    }
    
    # Find combinations where Astro > Micro (even more unusual)
    astro_beats_micro <- results_df %>%
      filter(Astro_power > Micro_power) %>%
      select(combination_id, ordering, n_cells_level, rpc_level, lfc_level,
             basemean_level, dispersion_level, Astro_power, Micro_power)
    
    if (nrow(astro_beats_micro) > 0) {
      log_info("Found %d combinations where Astro > Micro (very unusual)", nrow(astro_beats_micro))
      log_info("Sample:")
      print(head(astro_beats_micro, 5))
    } else {
      log_info("Micro always >= Astro in all %d combinations", n_combinations)
    }
    
    # ==============================================================================
    # 7. SAVE RESULTS
    # ==============================================================================
    
    # Save full results
    results_path <- file.path(OUTPUT_DIR, sprintf("ordering_permutation_results_%s.csv", variant_label))
    write.csv(results_df, results_path, row.names = FALSE)
    log_info("Saved full results: %s", results_path)
    
    # Save ordering frequency summary
    freq_path <- file.path(OUTPUT_DIR, sprintf("ordering_frequencies_%s.csv", variant_label))
    write.csv(ordering_counts, freq_path, row.names = FALSE)
    log_info("Saved frequency summary: %s", freq_path)
    
    # Save edge cases
    if (nrow(unusual_orderings) > 0) {
      edge_path <- file.path(OUTPUT_DIR, sprintf("edge_cases_unusual_%s.csv", variant_label))
      write.csv(unusual_orderings, edge_path, row.names = FALSE)
      log_info("Saved edge cases: %s", edge_path)
    }
    
    # ==============================================================================
    # 8. VISUALIZATION: Ordering Frequency Bar Chart
    # ==============================================================================
    
    log_info("Generating ordering frequency plot...")
    
    # Plot top orderings
    top_n <- min(15, nrow(ordering_counts))
    plot_data <- ordering_counts %>%
      head(top_n) %>%
      mutate(
        ordering_label = gsub(" > ", ">", ordering),
        ordering_label = factor(ordering_label, levels = rev(ordering_label))
      )
    
    p_freq <- ggplot(plot_data, aes(x = proportion, y = ordering_label, fill = proportion)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.1f%% (n=%d)", 100 * proportion, count)), 
                hjust = -0.1, size = 3.5) +
      scale_fill_gradient(low = "#E69F00", high = "#0072B2") +
      scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.3))) +
      labs(
        title = sprintf("Power Hierarchy Orderings (%s)", variant_label),
        subtitle = sprintf("All %d combinations of 5 covariates × 5 quantiles at %d donors", 
                           n_combinations, EVAL_DONORS),
        x = "Proportion of Combinations",
        y = "Ordering (Highest Power → Lowest)"
      ) +
      theme_publication(base_size = 14) +
      theme(legend.position = "none")
    
    freq_plot_path <- file.path(OUTPUT_DIR, sprintf("ordering_frequencies_%s.png", variant_label))
    ggsave(freq_plot_path, p_freq, width = 12, height = max(6, top_n * 0.4), dpi = 300)
    log_info("Saved frequency plot: %s", freq_plot_path)
    
    # ==============================================================================
    # 9. SUMMARY STATISTICS
    # ==============================================================================
    
    log_section("Summary Statistics")
    
    # Compute power ranges for each cell type
    power_summary <- results_df %>%
      summarise(
        Astro_min = min(Astro_power),
        Astro_max = max(Astro_power),
        Astro_median = median(Astro_power),
        Micro_min = min(Micro_power),
        Micro_max = max(Micro_power),
        Micro_median = median(Micro_power),
        Oligo_min = min(Oligo_power),
        Oligo_max = max(Oligo_power),
        Oligo_median = median(Oligo_power),
        Astro_Shallow_min = min(Astro_Shallow_power),
        Astro_Shallow_max = max(Astro_Shallow_power),
        Astro_Shallow_median = median(Astro_Shallow_power)
      )
    
    log_info("Power ranges at %d donors across all %d combinations:", EVAL_DONORS, n_combinations)
    print(power_summary)
    
    # Save summary
    summary_path <- file.path(OUTPUT_DIR, sprintf("ordering_test_summary_%s.txt", variant_label))
    sink(summary_path)
    cat(sprintf("Systematic Ordering Permutation Test Results\n"))
    cat(sprintf("Variant: %s\n", variant_label))
    cat(sprintf("Evaluation point: %d donors\n", EVAL_DONORS))
    cat(sprintf("Total combinations tested: %d (5^5)\n\n", n_combinations))
    cat(sprintf("Unique orderings found: %d\n\n", nrow(ordering_counts)))
    cat(sprintf("Top 5 orderings:\n"))
    print(head(ordering_counts, 5))
    cat(sprintf("\nCanonical ordering '%s': %.2f%% of combinations\n", 
                canonical_ordering, 100 * canonical_count / n_combinations))
    cat(sprintf("\nPower ranges:\n"))
    print(power_summary)
    if (nrow(unusual_orderings) > 0) {
      cat(sprintf("\nEdge cases (Oligo not last): %d combinations (%.2f%%)\n",
                  nrow(unusual_orderings), 100 * nrow(unusual_orderings) / n_combinations))
    }
    if (nrow(astro_beats_micro) > 0) {
      cat(sprintf("\nExtreme edge cases (Astro > Micro): %d combinations (%.2f%%)\n",
                  nrow(astro_beats_micro), 100 * nrow(astro_beats_micro) / n_combinations))
    }
    sink()
    
    log_info("Saved summary: %s", summary_path)
  }
}

log_section("Permutation Test Complete")
log_info("All results saved to: %s", OUTPUT_DIR)
