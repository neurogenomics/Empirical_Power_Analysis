# ==============================================================================
# SCRIPT: Visualise Power Curves and Surfaces (v4.0.0)
# PURPOSE: Generate power curve plots, cell type overlays, residual diagnostics,
#          and model comparison plots for the unified v3.1.0 models.
# INPUT: Fitted models from Step 2, fitting data from Step 1
# OUTPUT: PNG plots in Datalink/Power_Curve_Fitting/Plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
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

log_section("STEP 3: Visualise Power Curves (v4.0.0)")

# ==============================================================================
# 2. HELPER: Compute representative gene covariate values from data
# ==============================================================================

compute_representative_values <- function(df) {
  list(
    median_lfc = median(df$abs_gt_lfc, na.rm = TRUE),
    median_basemean = median(df$gt_basemean, na.rm = TRUE),
    median_dispersion = median(df$gt_dispersion, na.rm = TRUE),
    median_n_cells = median(df$n_cells, na.rm = TRUE),
    median_rpc = median(df$reads_per_cell, na.rm = TRUE),
    q25_lfc = quantile(df$abs_gt_lfc, 0.25, na.rm = TRUE),
    q75_lfc = quantile(df$abs_gt_lfc, 0.75, na.rm = TRUE),
    q25_basemean = quantile(df$gt_basemean, 0.25, na.rm = TRUE),
    q75_basemean = quantile(df$gt_basemean, 0.75, na.rm = TRUE)
  )
}

# ==============================================================================
# 3. VISUALISE EACH VARIANT (gene_set x agg_level)
# ==============================================================================

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {

    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    log_section(sprintf("Visualising variant: %s", variant_label))

    model_dir <- get_model_dir(gene_set, agg_level)
    plots_dir <- file.path(OUTPUT_BASE, "Plots", variant_label)
    safe_mkdir(plots_dir)

    # Load individual per-cell-type models, stripping training data to save memory
    # GLM objects store the full training frame (~700MB each); we only need coefficients for predict()
    per_ct_models <- list()
    ALL_CT_MODELS <- c("Astro", "Micro", "Oligo", "Astro_Shallow", "Astro_Shallow_Fewer")
    for (ct_model in ALL_CT_MODELS) {
      model_path <- file.path(model_dir, sprintf("%s_model.qs", ct_model))
      if (file.exists(model_path)) {
        full_model <- qs::qread(model_path)
        # Strip the large data components but keep what predict() needs
        full_model$data <- NULL
        full_model$model <- NULL
        full_model$residuals <- NULL
        full_model$fitted.values <- NULL
        full_model$effects <- NULL
        full_model$weights <- NULL
        full_model$prior.weights <- NULL
        full_model$linear.predictors <- NULL
        full_model$y <- NULL
        per_ct_models[[ct_model]] <- full_model
        rm(full_model); gc()
        log_info("Loaded model (slim): %s", ct_model)
      }
    }
    if (length(per_ct_models) == 0) {
      log_warning("No models found for variant %s, skipping", variant_label)
      next
    }

    # Load fitting data for all cell types
    all_data <- list()
    for (ct in CELL_TYPES) {
      paths <- get_paths(ct)
      shared_suffix <- if (gene_set == "shared") "_shared" else ""
      data_path <- file.path(paths$fitting_dir,
                              sprintf("unified_fitting_data_%s%s.qs",
                                      agg_level, shared_suffix))
      if (file.exists(data_path)) {
        all_data[[ct]] <- qs::qread(data_path)
      }
    }
    df_combined <- bind_rows(all_data)

    if (nrow(df_combined) == 0) {
      log_warning("No fitting data for variant %s", variant_label)
      next
    }

    # Compute representative values from data
    rep_vals <- compute_representative_values(df_combined)
    log_info("Representative values: LFC=%.2f, baseMean=%.1f, disp=%.3f, n_cells=%.0f, rpc=%.0f",
             rep_vals$median_lfc, rep_vals$median_basemean, rep_vals$median_dispersion,
             rep_vals$median_n_cells, rep_vals$median_rpc)

    # ------------------------------------------------------------------
    # 3a. Power vs Donors (per cell type, varying n_cells)
    # ------------------------------------------------------------------

    log_info("Generating power vs donors plots...")

    for (ct in names(per_ct_models)) {
      model <- per_ct_models[[ct]]
      if (is.null(model)) next

      # Select a few n_cells values to show as curves
      cell_steps_clean <- CELL_STEPS[!is.na(CELL_STEPS)]
      n_cells_display <- cell_steps_clean[seq(1, length(cell_steps_clean),
                                               length.out = min(5, length(cell_steps_clean)))]

      donor_seq <- exp(seq(log(5), log(600), length.out = 100))

      pred_data <- expand.grid(n_donors = donor_seq, n_cells_val = n_cells_display)
      pred_data$log_donors <- log(pred_data$n_donors)
      pred_data$log_n_cells <- log(pred_data$n_cells_val)
      pred_data$log_reads_per_cell <- log(rep_vals$median_rpc)
      pred_data$log_abs_gt_lfc <- log(rep_vals$median_lfc)
      pred_data$log_gt_basemean <- log(rep_vals$median_basemean)
      pred_data$log_gt_dispersion <- log(rep_vals$median_dispersion)

      pred_data$predicted <- predict(model, newdata = pred_data, type = "response")
      pred_data$n_cells_label <- factor(
        format(pred_data$n_cells_val, big.mark = ",", scientific = FALSE))

      p <- ggplot(pred_data, aes(x = n_donors, y = predicted, colour = n_cells_label)) +
        geom_line(linewidth = 1) +
        geom_hline(yintercept = POWER_TARGETS, linetype = "dashed", alpha = 0.3) +
        scale_x_log10(labels = comma) +
        scale_y_continuous(limits = c(0, 1), labels = percent) +
        labs(
          title = sprintf("%s: Power vs Donors", prettify_label(ct)),
          subtitle = sprintf("|LFC| = %.2f, base mean = %.0f, dispersion = %.2f, %s reads/cell",
                             rep_vals$median_lfc, rep_vals$median_basemean,
                             rep_vals$median_dispersion,
                             format(round(rep_vals$median_rpc), big.mark = ",")),
          x = "Number of Donors",
          y = "Predicted Detection Probability",
          colour = "Cells per Donor"
        ) +
        theme_publication(base_size = 16) +
        theme(legend.position = "right")

      ggsave(file.path(plots_dir, sprintf("Power_vs_Donors_%s.png", ct)),
             p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    }

    # ------------------------------------------------------------------
    # 3b. Cell Type Comparison Overlay (at matched covariates)
    # ------------------------------------------------------------------

    log_info("Generating cell type comparison plots...")

    if (length(per_ct_models) >= 2) {

      # Use median n_cells and rpc from the combined data
      n_cells_mid <- rep_vals$median_n_cells
      rpc_mid <- rep_vals$median_rpc

      donor_seq <- exp(seq(log(5), log(600), length.out = 100))

      overlay_data <- lapply(names(per_ct_models), function(ct) {
        model <- per_ct_models[[ct]]
        if (is.null(model)) return(NULL)

        pred_df <- data.frame(
          n_donors = donor_seq,
          log_donors = log(donor_seq),
          log_n_cells = log(n_cells_mid),
          log_reads_per_cell = log(rpc_mid),
          log_abs_gt_lfc = log(rep_vals$median_lfc),
          log_gt_basemean = log(rep_vals$median_basemean),
          log_gt_dispersion = log(rep_vals$median_dispersion)
        )

        pred_df$predicted <- predict(model, newdata = pred_df, type = "response")
        pred_df$cell_type <- ct
        pred_df
      })

      overlay_df <- bind_rows(overlay_data[!sapply(overlay_data, is.null)])

      p <- ggplot(overlay_df, aes(x = n_donors, y = predicted, colour = cell_type)) +
        geom_line(linewidth = 1.2) +
        geom_hline(yintercept = POWER_TARGETS, linetype = "dashed", alpha = 0.3) +
        scale_x_log10(labels = comma) +
        scale_y_continuous(limits = c(0, 1), labels = percent) +
        scale_colour_manual(values = CELL_TYPE_COLOURS, labels = prettify_labeller()) +
        labs(
          title = "Cell Type Comparison",
          subtitle = sprintf("%s cells/donor, %s reads/cell, |LFC| = %.2f, base mean = %.0f, dispersion = %.2f",
                             format(round(n_cells_mid), big.mark = ","),
                             format(round(rpc_mid), big.mark = ","),
                             rep_vals$median_lfc, rep_vals$median_basemean,
                             rep_vals$median_dispersion),
          x = "Number of Donors",
          y = "Predicted Detection Probability",
          colour = "Cell Type"
        ) +
        theme_publication(base_size = 16) +
        theme(legend.position = "right")

      ggsave(file.path(plots_dir, "CellType_Comparison_Donors.png"),
             p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

      # Also compare the combined model vs per-CT models
      if (!is.null(checkpoint$combined)) {
        combined_pred <- data.frame(
          n_donors = donor_seq,
          log_donors = log(donor_seq),
          log_n_cells = log(n_cells_mid),
          log_reads_per_cell = log(rpc_mid),
          log_abs_gt_lfc = log(rep_vals$median_lfc),
          log_gt_basemean = log(rep_vals$median_basemean),
          log_gt_dispersion = log(rep_vals$median_dispersion)
        )
        combined_pred$predicted <- predict(checkpoint$combined, newdata = combined_pred,
                                           type = "response")
        combined_pred$cell_type <- "Combined"

        overlay_with_combined <- bind_rows(overlay_df, combined_pred)
        overlay_with_combined$cell_type <- factor(overlay_with_combined$cell_type,
                                                   levels = c(CELL_TYPES, "Combined"))

        p2 <- ggplot(overlay_with_combined,
                     aes(x = n_donors, y = predicted, colour = cell_type,
                         linetype = cell_type)) +
          geom_line(linewidth = 1.2) +
          geom_hline(yintercept = POWER_TARGETS, linetype = "dashed", alpha = 0.3) +
          scale_x_log10(labels = comma) +
          scale_y_continuous(limits = c(0, 1), labels = percent) +
          scale_colour_manual(values = c(CELL_TYPE_COLOURS, Combined = "black"),
                              labels = prettify_labeller()) +
          scale_linetype_manual(values = c(Astro = "solid", Micro = "solid",
                                            Oligo = "solid", Combined = "dashed"),
                                labels = prettify_labeller()) +
          labs(
            title = "Per-Cell-Type vs Combined Model",
            subtitle = NULL,
            x = "Number of Donors",
            y = "Predicted Detection Probability",
            colour = "Cell Type",
            linetype = "Cell Type"
          ) +
          theme_publication(base_size = 16) +
          theme(legend.position = "right")

        ggsave(file.path(plots_dir, "PerCT_vs_Combined.png"),
               p2, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      }
    }

    # ------------------------------------------------------------------
    # 3c. [REMOVED] Residual diagnostics removed in v4.0.0 cleanup.
    #     The binary response (0/1) makes predicted-vs-observed scatter
    #     and normal QQ plots meaningless for logistic regression.
    #     Proper logistic diagnostics (calibration, binned residuals)
    #     may be added in a future version.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # 3d. Depth decomposition: n_cells vs reads_per_cell partial effects
    # ------------------------------------------------------------------

    log_info("Generating depth decomposition plots...")

    # Use the combined model (or first per-CT model) to show partial effects
    demo_model <- if (!is.null(checkpoint$combined)) {
      checkpoint$combined
    } else {
      per_ct_models[[1]]
    }

    if (!is.null(demo_model)) {
      # Partial effect of n_cells (holding rpc fixed)
      n_cells_seq <- exp(seq(log(5), log(3000), length.out = 100))
      partial_cells <- data.frame(
        n_cells = n_cells_seq,
        log_donors = log(100),
        log_n_cells = log(n_cells_seq),
        log_reads_per_cell = log(rep_vals$median_rpc),
        log_abs_gt_lfc = log(rep_vals$median_lfc),
        log_gt_basemean = log(rep_vals$median_basemean),
        log_gt_dispersion = log(rep_vals$median_dispersion)
      )
      partial_cells$predicted <- predict(demo_model, newdata = partial_cells, type = "response")

      # Partial effect of reads_per_cell (holding n_cells fixed)
      rpc_seq <- exp(seq(log(2000), log(15000), length.out = 100))
      partial_rpc <- data.frame(
        reads_per_cell = rpc_seq,
        log_donors = log(100),
        log_n_cells = log(rep_vals$median_n_cells),
        log_reads_per_cell = log(rpc_seq),
        log_abs_gt_lfc = log(rep_vals$median_lfc),
        log_gt_basemean = log(rep_vals$median_basemean),
        log_gt_dispersion = log(rep_vals$median_dispersion)
      )
      partial_rpc$predicted <- predict(demo_model, newdata = partial_rpc, type = "response")

      p_cells <- ggplot(partial_cells, aes(x = n_cells, y = predicted)) +
        geom_line(linewidth = 1.2, colour = "#377EB8") +
        scale_x_log10(labels = comma) +
        scale_y_continuous(limits = c(0, 1), labels = percent) +
        labs(
          title = "Partial Effect: Cells per Donor",
          subtitle = sprintf("100 donors, %s reads/cell, median gene covariates",
                             format(round(rep_vals$median_rpc), big.mark = ",")),
          x = "Cells per Donor",
          y = "Predicted Power"
        ) +
        theme_publication(base_size = 16)

      p_rpc <- ggplot(partial_rpc, aes(x = reads_per_cell, y = predicted)) +
        geom_line(linewidth = 1.2, colour = "#E41A1C") +
        scale_x_log10(labels = comma) +
        scale_y_continuous(limits = c(0, 1), labels = percent) +
        labs(
          title = "Partial Effect: Reads per Cell",
          subtitle = sprintf("100 donors, %s cells/donor, median gene covariates",
                             format(round(rep_vals$median_n_cells), big.mark = ",")),
          x = "Reads per Cell",
          y = "Predicted Power"
        ) +
        theme_publication(base_size = 16)

      p_depth <- cowplot::plot_grid(p_cells, p_rpc, ncol = 2, labels = c("A", "B"))

      ggsave(file.path(plots_dir, "Depth_Decomposition.png"),
             p_depth, width = PLOT_WIDTH * 1.5, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    }

    # ------------------------------------------------------------------
    # 3e. One-at-a-Time Sensitivity Grid (v4.2.0)
    #     Shows how the cell type hierarchy changes across covariate space.
    #     11 panels: 1 centre (all median) + 5 covariates × 2 levels (Q25/Q75)
    # ------------------------------------------------------------------

    log_info("Generating one-at-a-time sensitivity grid...")

    # Define the 5 models to display (including Astro_Shallow_Fewer)
    ALL_CT_MODELS <- c("Astro", "Micro", "Oligo", "Astro_Shallow", "Astro_Shallow_Fewer")
    CT_COLOURS_5 <- c(
      Astro = "#E41A1C",              # Red
      Micro = "#377EB8",              # Blue
      Oligo = "#4DAF4A",              # Green
      Astro_Shallow = "#E69F00",      # Orange (distinct from Astro)
      Astro_Shallow_Fewer = "#CC79A7" # Pink (distinct from all others)
    )
    CT_SHAPES_5 <- c(
      Astro = 16,              # Circle
      Micro = 17,              # Triangle
      Oligo = 15,              # Square
      Astro_Shallow = 1,       # Open circle
      Astro_Shallow_Fewer = 2  # Open triangle
    )

    # Compute quantiles from all 4 models' fitting data combined
    all_data_4 <- list()
    for (ct_model in ALL_CT_MODELS) {
      paths <- get_paths_for_dataset("Astro", list(dataset = "Roussos", model_names = NULL, fit_combined = TRUE))
      # Override model_name for non-Astro
      if (ct_model == "Astro_Shallow") {
        paths <- get_paths_for_dataset("Astro", list(dataset = "Roussos_DepthEq", model_names = c(Astro = "Astro_Shallow"), fit_combined = FALSE))
      } else if (ct_model == "Astro_Shallow_Fewer") {
        paths <- get_paths_for_dataset("Astro", list(dataset = "Roussos_ShallowFewer", model_names = c(Astro = "Astro_Shallow_Fewer"), fit_combined = FALSE))
      } else if (ct_model != "Astro") {
        paths <- get_paths_for_dataset(ct_model, list(dataset = "Roussos", model_names = NULL, fit_combined = TRUE))
      }

      shared_suffix <- if (gene_set == "shared") "_shared" else ""
      data_path <- file.path(paths$fitting_dir,
                              sprintf("unified_fitting_data_%s%s.qs",
                                      agg_level, shared_suffix))
      if (file.exists(data_path)) {
        all_data_4[[ct_model]] <- qs::qread(data_path)
      }
    }
    df_4_combined <- bind_rows(all_data_4)

    if (nrow(df_4_combined) > 0) {
      # Compute uniform intervals within the observed range for each covariate
      # 5 evenly spaced points from min to max (indices: 1=min, 2=25%, 3=mid, 4=75%, 5=max)
      uniform_seq <- function(x) {
        rng <- range(x, na.rm = TRUE)
        seq(rng[1], rng[2], length.out = 5)
      }
      quantiles <- list(
        n_cells = uniform_seq(df_4_combined$n_cells),
        reads_per_cell = uniform_seq(df_4_combined$reads_per_cell),
        abs_gt_lfc = uniform_seq(df_4_combined$abs_gt_lfc),
        gt_basemean = uniform_seq(df_4_combined$gt_basemean),
        gt_dispersion = uniform_seq(df_4_combined$gt_dispersion)
      )

      log_info("Uniform intervals (min/25%%/mid/75%%/max): n_cells=%.0f/%.0f/%.0f/%.0f/%.0f, rpc=%.0f/%.0f/%.0f/%.0f/%.0f",
               quantiles$n_cells[1], quantiles$n_cells[2], quantiles$n_cells[3], quantiles$n_cells[4], quantiles$n_cells[5],
               quantiles$reads_per_cell[1], quantiles$reads_per_cell[2], quantiles$reads_per_cell[3], quantiles$reads_per_cell[4], quantiles$reads_per_cell[5])

      # Build 21 prediction scenarios (1 centre + 5 covariates × 4 levels)
      donor_seq <- exp(seq(log(5), log(600), length.out = 100))

      # Centre panel: all at median (Q50 = index 3)
      centre_df <- expand.grid(n_donors = donor_seq, model = ALL_CT_MODELS)
      centre_df$log_donors <- log(centre_df$n_donors)
      centre_df$log_n_cells <- log(quantiles$n_cells[3])
      centre_df$log_reads_per_cell <- log(quantiles$reads_per_cell[3])
      centre_df$log_abs_gt_lfc <- log(quantiles$abs_gt_lfc[3])
      centre_df$log_gt_basemean <- log(quantiles$gt_basemean[3])
      centre_df$log_gt_dispersion <- log(quantiles$gt_dispersion[3])
      centre_df$panel_label <- "Centre (All Midpoint)"
      centre_df$row_label <- "Centre"
      centre_df$col_label <- ""

      # Function to create satellite panels
      make_satellite_df <- function(cov_name, level_idx, level_name) {
        df <- expand.grid(n_donors = donor_seq, model = ALL_CT_MODELS)
        df$log_donors <- log(df$n_donors)

        # All at median (index 3) except the varied covariate
        df$log_n_cells <- log(ifelse(cov_name == "n_cells", quantiles$n_cells[level_idx], quantiles$n_cells[3]))
        df$log_reads_per_cell <- log(ifelse(cov_name == "reads_per_cell", quantiles$reads_per_cell[level_idx], quantiles$reads_per_cell[3]))
        df$log_abs_gt_lfc <- log(ifelse(cov_name == "abs_gt_lfc", quantiles$abs_gt_lfc[level_idx], quantiles$abs_gt_lfc[3]))
        df$log_gt_basemean <- log(ifelse(cov_name == "gt_basemean", quantiles$gt_basemean[level_idx], quantiles$gt_basemean[3]))
        df$log_gt_dispersion <- log(ifelse(cov_name == "gt_dispersion", quantiles$gt_dispersion[level_idx], quantiles$gt_dispersion[3]))

        df$panel_label <- sprintf("%s %s", cov_name, level_name)
        df$row_label <- cov_name
        df$col_label <- level_name
        return(df)
      }

      # Build satellite panels (Min, 25%, 75%, Max for each covariate — uniform intervals)
      covariates <- c("n_cells", "reads_per_cell", "abs_gt_lfc", "gt_basemean", "gt_dispersion")
      satellite_dfs <- list()
      for (cov in covariates) {
        satellite_dfs[[paste0(cov, "_Q5")]] <- make_satellite_df(cov, 1, "Min")
        satellite_dfs[[paste0(cov, "_Q25")]] <- make_satellite_df(cov, 2, "25%")
        satellite_dfs[[paste0(cov, "_Q75")]] <- make_satellite_df(cov, 4, "75%")
        satellite_dfs[[paste0(cov, "_Q95")]] <- make_satellite_df(cov, 5, "Max")
      }

      # Predict for each model
      predict_for_model <- function(df_subset) {
        results <- list()
        for (m in ALL_CT_MODELS) {
          model <- per_ct_models[[m]]
          if (is.null(model)) next

          df_m <- df_subset %>% filter(model == m)
          if (nrow(df_m) == 0) next

          df_m$predicted <- predict(model, newdata = df_m, type = "response")
          results[[m]] <- df_m
        }
        bind_rows(results)
      }

      # Centre panel for separate save
      pred_centre <- predict_for_model(centre_df)

      # Satellite panels
      pred_satellites <- lapply(satellite_dfs, predict_for_model)

      # Create individual ggplot objects
      create_panel <- function(pred_df, title, show_legend = FALSE) {
        p <- ggplot(pred_df, aes(x = n_donors, y = predicted,
                                  colour = model, shape = model)) +
          geom_line(linewidth = 1, aes(linetype = model)) +
          geom_hline(yintercept = POWER_TARGETS, linetype = "dashed", alpha = 0.3) +
          scale_x_log10(labels = comma, limits = c(5, 600)) +
          scale_y_continuous(limits = c(0, 1), labels = percent) +
          scale_colour_manual(values = CT_COLOURS_5, labels = prettify_label) +
          scale_shape_manual(values = CT_SHAPES_5, labels = prettify_label) +
          scale_linetype_manual(values = c(Astro = "solid", Micro = "solid", Oligo = "solid", 
                                           Astro_Shallow = "dashed", Astro_Shallow_Fewer = "dotted"),
                               labels = prettify_label) +
          labs(
            title = title,
            subtitle = "All other covariates at median",
            x = "Number of Donors",
            y = "Predicted Power"
          ) +
          theme_publication(base_size = 14) +
          theme(
            legend.position = if (show_legend) "right" else "none",
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 11, colour = "grey40")
          )
        return(p)
      }

      # Build title strings with actual values (index 3 = Q50/median)
      centre_title <- sprintf("All Covariates at Midpoint\n%s cells/donor, %s reads/cell, |LFC| = %.2f, base mean = %.0f, disp. = %.2f",
                               format(round(quantiles$n_cells[3]), big.mark = ","),
                               format(round(quantiles$reads_per_cell[3]), big.mark = ","),
                               quantiles$abs_gt_lfc[3], quantiles$gt_basemean[3], quantiles$gt_dispersion[3])

      # Save centre panel separately with legend
      p_centre <- create_panel(pred_centre, centre_title, show_legend = TRUE)
      ggsave(file.path(plots_dir, "CellType_Comparison_Median.png"),
             p_centre, width = PLOT_WIDTH * 1.3, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      log_info("Saved median comparison: CellType_Comparison_Median.png")

      # Satellite panel titles with values
      pretty_cov_names <- c(
        n_cells = "Cells per Donor",
        reads_per_cell = "Reads per Cell",
        abs_gt_lfc = "|LFC|",
        gt_basemean = "Base Mean Expression",
        gt_dispersion = "Dispersion"
      )

      p_satellites <- list()
      for (cov in covariates) {
        q5_val <- quantiles[[cov]][1]
        q25_val <- quantiles[[cov]][2]
        q75_val <- quantiles[[cov]][4]
        q95_val <- quantiles[[cov]][5]
        if (cov %in% c("n_cells", "reads_per_cell", "gt_basemean")) {
          q5_str <- format(round(q5_val), big.mark = ",")
          q25_str <- format(round(q25_val), big.mark = ",")
          q75_str <- format(round(q75_val), big.mark = ",")
          q95_str <- format(round(q95_val), big.mark = ",")
        } else {
          q5_str <- sprintf("%.2f", q5_val)
          q25_str <- sprintf("%.2f", q25_val)
          q75_str <- sprintf("%.2f", q75_val)
          q95_str <- sprintf("%.2f", q95_val)
        }

        title_q5 <- sprintf("%s = %s\n(Min)", pretty_cov_names[cov], q5_str)
        title_q25 <- sprintf("%s = %s\n(25%%)", pretty_cov_names[cov], q25_str)
        title_q75 <- sprintf("%s = %s\n(75%%)", pretty_cov_names[cov], q75_str)
        title_q95 <- sprintf("%s = %s\n(Max)", pretty_cov_names[cov], q95_str)

        p_satellites[[paste0(cov, "_Q5")]] <- create_panel(pred_satellites[[paste0(cov, "_Q5")]], title_q5)
        p_satellites[[paste0(cov, "_Q25")]] <- create_panel(pred_satellites[[paste0(cov, "_Q25")]], title_q25)
        p_satellites[[paste0(cov, "_Q75")]] <- create_panel(pred_satellites[[paste0(cov, "_Q75")]], title_q75)
        p_satellites[[paste0(cov, "_Q95")]] <- create_panel(pred_satellites[[paste0(cov, "_Q95")]], title_q95)
      }

      # Assemble grid: just the satellite rows with 4 columns each
      # Create a dummy plot with legend to extract for the grid
      p_dummy <- create_panel(pred_satellites[[1]], "", show_legend = TRUE)
      get_legend <- function(p) {
        g <- ggplotGrob(p)
        legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
        return(legend)
      }
      legend_grob <- get_legend(p_dummy)

      # Build rows: 4 columns per covariate (Q5 | Q25 | Q75 | Q95)
      row_list <- list()
      for (cov in covariates) {
        p_q5 <- p_satellites[[paste0(cov, "_Q5")]]
        p_q25 <- p_satellites[[paste0(cov, "_Q25")]]
        p_q75 <- p_satellites[[paste0(cov, "_Q75")]]
        p_q95 <- p_satellites[[paste0(cov, "_Q95")]]
        # 4-column row for this covariate
        cov_row <- cowplot::plot_grid(p_q5, p_q25, p_q75, p_q95, ncol = 4)
        row_list[[length(row_list) + 1]] <- cov_row
      }

      # Combine all covariate rows
      p_grid <- cowplot::plot_grid(plotlist = row_list, ncol = 1, rel_heights = rep(1, length(covariates)))

      # Add legend on the right side
      p_grid_with_legend <- cowplot::plot_grid(p_grid, legend_grob, ncol = 2, rel_widths = c(4, 0.8))

      # Add main title
      p_grid_title <- cowplot::ggdraw() +
        cowplot::draw_label("Sensitivity Analysis Across Covariate Range",
                            x = 0.5, y = 0.95, size = 15, fontface = "bold") +
        cowplot::draw_label("Each row varies one covariate at uniform intervals (Min, 25%, 75%, Max); all others held at midpoint",
                            x = 0.5, y = 0.88, size = 10, colour = "grey40")

      p_final <- cowplot::plot_grid(p_grid_title, p_grid_with_legend, ncol = 1, rel_heights = c(0.10, 1))

      # 4-column layout: wider figure to accommodate
      ggsave(file.path(plots_dir, "CellType_Sensitivity_Grid.png"),
             p_final, width = PLOT_WIDTH * 2.4, height = PLOT_HEIGHT * 3.5, dpi = PLOT_DPI, limitsize = FALSE)

      log_info("Saved sensitivity grid: CellType_Sensitivity_Grid.png")
    } else {
      log_warning("No fitting data for sensitivity grid, skipping")
    }

    log_info("Visualisation complete for variant: %s", variant_label)
  }
}

# ==============================================================================
# 4. MODEL COMPARISON PLOTS (across variants)
# ==============================================================================

log_section("Generating model comparison plots")

comparison_plots_dir <- file.path(OUTPUT_BASE, "Plots", "Comparison")
safe_mkdir(comparison_plots_dir)

# Load model comparison table
comparison_path <- file.path(OUTPUT_BASE, "Models", "model_comparison_v3.csv")
if (!file.exists(comparison_path)) {
  log_warning("model_comparison_v3.csv not found, skipping comparison plots")
} else {

  comp_df <- read.csv(comparison_path, stringsAsFactors = FALSE)

  # ------------------------------------------------------------------
  # 4a. R² comparison: per-CT vs combined
  # ------------------------------------------------------------------

  log_info("Generating R² comparison plot...")

  r2_long <- comp_df %>%
    select(variant, r2_Astro, r2_Micro, r2_Oligo, r2_combined,
           r2_additive, r2_interaction) %>%
    pivot_longer(cols = starts_with("r2_"),
                 names_to = "model", values_to = "r2") %>%
    mutate(model = sub("r2_", "", model),
           model = factor(model, levels = c("Astro", "Micro", "Oligo",
                                            "combined", "additive", "interaction")))

  p_r2 <- ggplot(r2_long, aes(x = variant, y = r2, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(CELL_TYPE_COLOURS,
                                 combined = "grey40",
                                 additive = "grey60",
                                 interaction = "grey80"),
                      labels = prettify_labeller()) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      title = "Pseudo R² Comparison Across Model Variants",
      subtitle = "Combined = no cell type term; Additive = +cell type; Interaction = cell type × predictors",
      x = "Model Variant (gene_set_aggregation)",
      y = "Pseudo R²",
      fill = "Model"
    ) +
    theme_publication(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(comparison_plots_dir, "R2_Comparison.png"),
         p_r2, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

  # ------------------------------------------------------------------
  # 4b. LRT Summary
  # ------------------------------------------------------------------

  log_info("Generating LRT summary plot...")

  # Build long-form data with chi-squared values, df, and p-values
  lrt_df <- comp_df %>%
    select(variant, lrt1_chisq, lrt1_df, lrt1_pval, lrt1_significant,
           lrt2_chisq, lrt2_df, lrt2_pval, lrt2_significant) %>%
    pivot_longer(
      cols = c(lrt1_chisq, lrt2_chisq),
      names_to = "test", values_to = "chisq"
    ) %>%
    mutate(
      test_label = ifelse(grepl("lrt1", test),
                          "LRT1: No CT vs +CT (df=2)",
                          "LRT2: +CT vs CT*pred (df=12)"),
      df = ifelse(grepl("lrt1", test), lrt1_df, lrt2_df),
      pval = ifelse(grepl("lrt1", test), lrt1_pval, lrt2_pval),
      significant = ifelse(grepl("lrt1", test), lrt1_significant, lrt2_significant),
      # Format p-value label: show "p < 10^-300" for p = 0
      p_label = ifelse(pval == 0, "p < 10^-300",
                       sprintf("p = %.1e", pval)),
      chisq_label = format(round(chisq), big.mark = ",")
    )

  p_lrt <- ggplot(lrt_df, aes(x = variant, y = chisq, fill = test_label)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = chisq_label),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 4.5) +
    geom_text(aes(label = p_label, y = chisq * 0.5),
              position = position_dodge(width = 0.7),
              vjust = 0.5, size = 4, colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("LRT1: No CT vs +CT (df=2)" = "#E41A1C",
                                 "LRT2: +CT vs CT*pred (df=12)" = "#377EB8")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "LRT Results: Does Cell Type Matter?",
      subtitle = "All tests overwhelmingly significant (p < 10^-300)",
      x = "Model Variant",
      y = expression(chi^2 ~ "statistic"),
      fill = "LRT Test"
    ) +
    theme_publication(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(comparison_plots_dir, "LRT_Summary.png"),
         p_lrt, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

  # ------------------------------------------------------------------
  # 4c. Coefficient Forest Plots
  #     Shows per-CT models (coloured circles) + combined and additive
  #     models (black diamonds: filled = combined, open = additive).
  #     Interaction model excluded (base terms = Astro reference level).
  #     Interaction terms (cell_type:predictor) excluded.
  # ------------------------------------------------------------------

  log_info("Generating coefficient forest plots...")

  all_variant_coefs <- list()

  # Models to include: all 5 cell types (including sensitivity variants) + combined + additive
  # Exclude combined_interaction (base terms duplicate Astro reference level)
  ALL_CELL_TYPE_MODELS <- c(CELL_TYPES, "Astro_Shallow", "Astro_Shallow_Fewer")
  FOREST_MODELS <- c(ALL_CELL_TYPE_MODELS, "combined", "combined_additive")

  for (gene_set in GENE_SETS) {
    for (agg_level in AGGREGATION_LEVELS) {
      model_dir <- get_model_dir(gene_set, agg_level)
      coef_path <- file.path(model_dir, "all_coefficients.csv")
      if (!file.exists(coef_path)) next

      coefs <- read.csv(coef_path, stringsAsFactors = FALSE)
      coefs <- coefs %>%
        # Keep only the 5 target models
        filter(model %in% FOREST_MODELS) %>%
        # Exclude intercepts, phi, cell_type intercepts, and interaction terms
        filter(!grepl("^\\(|^cell_type|:", term)) %>%
        mutate(
          gene_set = gene_set,
          agg_level = agg_level,
          variant = paste0(gene_set, "_", agg_level),
          ci_low = estimate - qnorm(1 - (1 - CI_LEVEL) / 2) * std_error,
          ci_high = estimate + qnorm(1 - (1 - CI_LEVEL) / 2) * std_error,
          # Display labels for legend
          model_label = case_when(
            model %in% ALL_CELL_TYPE_MODELS ~ model,
            model == "combined" ~ "Combined",
            model == "combined_additive" ~ "Additive",
            TRUE ~ model
          )
        )
      all_variant_coefs[[paste(gene_set, agg_level, sep = "_")]] <- coefs
    }
  }

  if (length(all_variant_coefs) > 0) {
    coef_df <- bind_rows(all_variant_coefs)

    # Set factor order for legend: all cell types first, then pooled models
    model_order <- c(ALL_CELL_TYPE_MODELS, "Combined", "Additive")
    coef_df$model_label <- factor(coef_df$model_label, levels = model_order)

    # Colour scale: all cell types use CELL_TYPE_COLOURS, pooled models use black
    forest_colours <- c(CELL_TYPE_COLOURS, Combined = "black", Additive = "black")

    # Shape scale: Astro = circle (16), Micro = triangle (17), Oligo = square (15),
    # Astro_Shallow = open circle (1), Astro_Shallow_Fewer = open triangle (2),
    # combined = filled diamond (18), additive = open diamond (5)
    forest_shapes <- c(
      Astro = 16, Micro = 17, Oligo = 15,
      Astro_Shallow = 1, Astro_Shallow_Fewer = 2,
      Combined = 18, Additive = 5
    )

    plot_terms <- unique(coef_df$term)
    plot_terms <- plot_terms[!plot_terms %in% c("(Intercept)", "(phi)")]

    for (term_name in plot_terms) {
      term_data <- coef_df %>% filter(term == term_name)
      if (nrow(term_data) == 0) next

      p_forest <- ggplot(term_data, aes(x = estimate, y = variant,
                                          colour = model_label,
                                          shape = model_label)) +
        geom_point(size = 3, position = position_dodge(width = 0.6)) +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                        height = 0.2, position = position_dodge(width = 0.6)) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_colour_manual(values = forest_colours, labels = prettify_labeller()) +
        scale_shape_manual(values = forest_shapes, labels = prettify_labeller()) +
        labs(
          title = sprintf("Coefficient: %s", term_name),
          subtitle = "Per-Cell-Type models (circles) vs pooled models (diamonds)",
          x = "Estimate (95% CI)",
          y = "Model Variant",
          colour = "Model",
          shape = "Model"
        ) +
        theme_publication(base_size = 16)

      safe_term <- gsub("[^[:alnum:]]", "_", term_name)
      ggsave(file.path(comparison_plots_dir, sprintf("Forest_%s.png", safe_term)),
             p_forest, width = PLOT_WIDTH, height = PLOT_HEIGHT * 0.7, dpi = PLOT_DPI)
    }
  }
}

log_section("STEP 3 COMPLETE")
log_info("All visualisations generated")
