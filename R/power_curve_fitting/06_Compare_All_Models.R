# ==============================================================================
# SCRIPT: Compare All Models Including Astro Shallow Variants (v4.2.0)
# PURPOSE: Generate comparison visualisations across all fitted models
#          (Astro, Micro, Oligo, Astro_Shallow, Astro_Shallow_Fewer) to assess 
#          depth equalization and cell count matching effects.
# INPUT: Fitted models from Step 2 (both gene sets)
# OUTPUT: Comparison plots in Datalink/Power_Curve_Fitting/Plots/with_shallow/
# 
# v4.2.0: Added Astro_Shallow_Fewer (depth + cell count matched to Microglia)
#         alongside Astro_Shallow (depth-only matched) for comprehensive
#         sensitivity analysis of depth vs cell count contributions.
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(cowplot)
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

log_section("STEP 6: Compare All Models (v4.1.0)")

# ==============================================================================
# 2. OUTPUT DIRECTORY SETUP (Parallel Structure)
# ==============================================================================

# Create parallel output structure that preserves original outputs
OUTPUT_BASE_SHALLOW <- project_path("Datalink", "Power_Curve_Fitting")
PLOTS_DIR_SHALLOW <- file.path(OUTPUT_BASE_SHALLOW, "Plots", "with_shallow")
SUMMARY_DIR_SHALLOW <- file.path(OUTPUT_BASE_SHALLOW, "Summary", "with_shallow")

safe_mkdir(OUTPUT_BASE_SHALLOW)
safe_mkdir(PLOTS_DIR_SHALLOW)
safe_mkdir(SUMMARY_DIR_SHALLOW)

log_info("Output base: %s", OUTPUT_BASE_SHALLOW)

# ==============================================================================
# 3. CONFIGURATION
# ==============================================================================

# Models to compare (all 5 per-CT models including both Astro shallow variants)
ALL_MODELS <- c("Astro", "Micro", "Oligo", "Astro_Shallow", "Astro_Shallow_Fewer")

# Colour palette: Distinct colors for all models
# Astro variants use distinct red-family colors for clear bar chart distinction
MODEL_COLOURS <- c(
  Astro = "#E41A1C",                # Red (native astrocytes)
  Micro = "#377EB8",                # Blue
  Oligo = "#4DAF4A",                # Green
  Astro_Shallow = "#FF7F0E",        # Orange (depth-equalized astrocytes)
  Astro_Shallow_Fewer = "#9467BD"  # Purple (depth + cell count matched astrocytes)
)

MODEL_LINETYPES <- c(
  Astro = "solid",
  Micro = "solid",
  Oligo = "solid",
  Astro_Shallow = "dashed",
  Astro_Shallow_Fewer = "dotted"   # Distinct linetype for fewer cells variant
)

MODEL_SHAPES <- c(
  Astro = 16,               # Filled circle
  Micro = 17,               # Filled triangle
  Oligo = 15,               # Filled square
  Astro_Shallow = 1,        # Open circle
  Astro_Shallow_Fewer = 2   # Open triangle (distinct from Astro_Shallow)
)

# Power targets for required sample size comparison
POWER_TARGETS_COMPARE <- c(0.50, 0.80, 0.90)

# ==============================================================================
# 4. HELPER: Load All Models for a Gene Set
# ==============================================================================

load_all_models <- function(gene_set, agg_level = "per_iteration") {
  variant_label <- sprintf("%s_%s", gene_set, agg_level)
  model_dir <- get_model_dir(gene_set, agg_level)
  
  models <- list()
  for (model_name in ALL_MODELS) {
    model_path <- file.path(model_dir, sprintf("%s_model.qs", model_name))
    if (file.exists(model_path)) {
      models[[model_name]] <- qs::qread(model_path)
      log_info("Loaded model: %s (%s)", model_name, variant_label)
    } else {
      log_warning("Model not found: %s (%s)", model_name, variant_label)
    }
  }
  return(models)
}

# ==============================================================================
# 5. HELPER: Load Fitting Data for Representative Values
# ==============================================================================

load_all_fitting_data <- function(gene_set, agg_level = "per_iteration") {
  shared_suffix <- if (gene_set == "shared") "_shared" else ""
  
  all_data <- list()
  
  # Load from all model output directories
  for (model_name in ALL_MODELS) {
    fitting_dir <- file.path(OUTPUT_BASE, "Fitting_Data", model_name)
    data_path <- file.path(fitting_dir,
                           sprintf("unified_fitting_data_%s%s.qs",
                                   agg_level, shared_suffix))
    if (file.exists(data_path)) {
      all_data[[model_name]] <- qs::qread(data_path)
    }
  }
  
  if (length(all_data) == 0) return(NULL)
  return(bind_rows(all_data))
}

compute_representative_values <- function(df) {
  list(
    median_lfc = median(df$abs_gt_lfc, na.rm = TRUE),
    median_basemean = median(df$gt_basemean, na.rm = TRUE),
    median_dispersion = median(df$gt_dispersion, na.rm = TRUE),
    median_n_cells = median(df$n_cells, na.rm = TRUE),
    median_rpc = median(df$reads_per_cell, na.rm = TRUE)
  )
}

# ==============================================================================
# 6. HELPER: Find Required Donors for Target Power
# ==============================================================================

find_required_donors <- function(model, target_power, rep_vals,
                                 max_donors = 1000, tolerance = 0.01) {
  if (is.null(model)) return(NA_real_)
  
  # Binary search for required donor count
  low <- 5
  high <- max_donors
  
  # Check if target is achievable
  pred_high <- predict(model, newdata = data.frame(
    log_donors = log(high),
    log_n_cells = log(rep_vals$median_n_cells),
    log_reads_per_cell = log(rep_vals$median_rpc),
    log_abs_gt_lfc = log(rep_vals$median_lfc),
    log_gt_basemean = log(rep_vals$median_basemean),
    log_gt_dispersion = log(rep_vals$median_dispersion)
  ), type = "response")
  
  if (pred_high < target_power) return(NA_real_)  # Cannot achieve target
  
  # Binary search
  while (high - low > 1) {
    mid <- (low + high) / 2
    pred_mid <- predict(model, newdata = data.frame(
      log_donors = log(mid),
      log_n_cells = log(rep_vals$median_n_cells),
      log_reads_per_cell = log(rep_vals$median_rpc),
      log_abs_gt_lfc = log(rep_vals$median_lfc),
      log_gt_basemean = log(rep_vals$median_basemean),
      log_gt_dispersion = log(rep_vals$median_dispersion)
    ), type = "response")
    
    if (pred_mid < target_power) {
      low <- mid
    } else {
      high <- mid
    }
  }
  
  return(ceiling(high))
}

# ==============================================================================
# 7. GENERATE PLOTS FOR EACH GENE SET
# ==============================================================================

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {
    
    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    log_section(sprintf("Comparing models for variant: %s", variant_label))
    
    # Create output directory for this variant
    variant_plots_dir <- file.path(PLOTS_DIR_SHALLOW, variant_label)
    safe_mkdir(variant_plots_dir)
    
    # Load models
    models <- load_all_models(gene_set, agg_level)
    if (length(models) == 0) {
      log_warning("No models found for variant %s, skipping", variant_label)
      next
    }
    
    # Load fitting data for representative values
    fitting_data <- load_all_fitting_data(gene_set, agg_level)
    if (is.null(fitting_data) || nrow(fitting_data) == 0) {
      log_warning("No fitting data for variant %s, using defaults", variant_label)
      rep_vals <- list(
        median_lfc = DEFAULT_REPRESENTATIVE_LFC,
        median_basemean = DEFAULT_REPRESENTATIVE_BASEMEAN,
        median_dispersion = DEFAULT_REPRESENTATIVE_DISPERSION,
        median_n_cells = 500,
        median_rpc = 6000
      )
    } else {
      rep_vals <- compute_representative_values(fitting_data)
    }
    
    log_info("Representative values: LFC=%.2f, baseMean=%.1f, disp=%.3f, n_cells=%.0f, rpc=%.0f",
             rep_vals$median_lfc, rep_vals$median_basemean, rep_vals$median_dispersion,
             rep_vals$median_n_cells, rep_vals$median_rpc)
    
    available_models <- names(models)
    log_info("Available models: %s", paste(available_models, collapse = ", "))
    
    # ------------------------------------------------------------------
    # 7a. PLOT 1: Power Curve Overlays (All 4 Models)
    # ------------------------------------------------------------------
    
    log_info("Generating power curve overlay plot...")
    
    donor_seq <- exp(seq(log(5), log(600), length.out = 100))
    
    overlay_data <- lapply(available_models, function(model_name) {
      model <- models[[model_name]]
      if (is.null(model)) return(NULL)
      
      pred_df <- data.frame(
        n_donors = donor_seq,
        log_donors = log(donor_seq),
        log_n_cells = log(rep_vals$median_n_cells),
        log_reads_per_cell = log(rep_vals$median_rpc),
        log_abs_gt_lfc = log(rep_vals$median_lfc),
        log_gt_basemean = log(rep_vals$median_basemean),
        log_gt_dispersion = log(rep_vals$median_dispersion)
      )
      
      pred_df$predicted <- predict(model, newdata = pred_df, type = "response")
      pred_df$model <- model_name
      pred_df
    })
    
    overlay_df <- bind_rows(overlay_data[!sapply(overlay_data, is.null)])
    overlay_df$model <- factor(overlay_df$model, levels = ALL_MODELS)
    
    p_overlay <- ggplot(overlay_df, aes(x = n_donors, y = predicted,
                                         colour = model, linetype = model)) +
      geom_line(linewidth = 1.2) +
      # Only the 80% reference line — the default 25% interval gridlines
      # supplied by theme_publication carry the rest.
      geom_hline(yintercept = 0.8, linetype = "dotted", alpha = 0.6) +
      scale_x_log10(labels = comma) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent) +
      scale_colour_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
      scale_linetype_manual(values = MODEL_LINETYPES, labels = prettify_labeller()) +
      labs(
        title = "Power Curve Comparison",
        subtitle = sprintf("|LFC| = %.2f, base mean = %.0f, dispersion = %.2f, %s cells/donor, %s reads/cell",
                           rep_vals$median_lfc, rep_vals$median_basemean,
                           rep_vals$median_dispersion,
                           format(round(rep_vals$median_n_cells), big.mark = ","),
                           format(round(rep_vals$median_rpc), big.mark = ",")),
        x = "Number of Donors",
        y = "Predicted Detection Probability",
        colour = "Cell Type",
        linetype = "Cell Type"
      ) +
      theme_publication(base_size = 16) +
      theme(legend.position = "right")
    
    ggsave(file.path(variant_plots_dir, "Power_Curve_Overlay.png"),
           p_overlay, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    # No-legend variant for the FigPower A/B composite — the bottom panel
    # (CellType_Sensitivity_Simple) carries the shared legend, so this top
    # panel must not duplicate it.
    p_overlay_nolegend <- p_overlay + theme(legend.position = "none")
    ggsave(file.path(variant_plots_dir, "Power_Curve_Overlay_NoLegend.png"),
           p_overlay_nolegend, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    # ------------------------------------------------------------------
    # 7a-LFC. PLOT 1 variants: Power Curve Overlays at fixed LFC values
    # ------------------------------------------------------------------

    for (fixed_lfc in c(1, 2)) {
      log_info("Generating power curve overlay plot at LFC=%d...", fixed_lfc)

      overlay_lfc_data <- lapply(available_models, function(model_name) {
        model <- models[[model_name]]
        if (is.null(model)) return(NULL)

        pred_df <- data.frame(
          n_donors = donor_seq,
          log_donors = log(donor_seq),
          log_n_cells = log(rep_vals$median_n_cells),
          log_reads_per_cell = log(rep_vals$median_rpc),
          log_abs_gt_lfc = log(fixed_lfc),
          log_gt_basemean = log(rep_vals$median_basemean),
          log_gt_dispersion = log(rep_vals$median_dispersion)
        )

        pred_df$predicted <- predict(model, newdata = pred_df, type = "response")
        pred_df$model <- model_name
        pred_df
      })

      overlay_lfc_df <- bind_rows(overlay_lfc_data[!sapply(overlay_lfc_data, is.null)])
      overlay_lfc_df$model <- factor(overlay_lfc_df$model, levels = ALL_MODELS)

      p_overlay_lfc <- ggplot(overlay_lfc_df, aes(x = n_donors, y = predicted,
                                                     colour = model, linetype = model)) +
        geom_line(linewidth = 1.2) +
        geom_hline(yintercept = POWER_TARGETS, linetype = "dotted", alpha = 0.4) +
        scale_x_log10(labels = comma) +
        scale_y_continuous(limits = c(0, 1), labels = percent) +
        scale_colour_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
        scale_linetype_manual(values = MODEL_LINETYPES, labels = prettify_labeller()) +
        labs(
          title = sprintf("Power Curve Comparison (|LFC| = %d)", fixed_lfc),
          subtitle = sprintf("base mean = %.0f, dispersion = %.2f, %s cells/donor, %s reads/cell",
                             rep_vals$median_basemean,
                             rep_vals$median_dispersion,
                             format(round(rep_vals$median_n_cells), big.mark = ","),
                             format(round(rep_vals$median_rpc), big.mark = ",")),
          x = "Number of Donors",
          y = "Predicted Detection Probability",
          colour = "Cell Type",
          linetype = "Cell Type"
        ) +
        theme_publication(base_size = 16) +
        theme(legend.position = "right")

      ggsave(file.path(variant_plots_dir, sprintf("Power_Curve_Overlay_LFC%d.png", fixed_lfc)),
             p_overlay_lfc, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    }

    # ------------------------------------------------------------------
    # 7a2. PLOT 1b: Simplified Power Curve (Practical Reference Point)
    #      For main Figure 3B — clean, no subtitle, publication-ready
    # ------------------------------------------------------------------

    log_info("Generating simplified sensitivity plot...")

    # Use practical reference point values
    PRACTICAL_LFC <- 1.5
    PRACTICAL_BASEMEAN <- 50
    PRACTICAL_DISPERSION <- 0.5
    PRACTICAL_N_CELLS <- 250
    PRACTICAL_RPC <- rep_vals$median_rpc

    simple_data <- lapply(available_models, function(model_name) {
      model <- models[[model_name]]
      if (is.null(model)) return(NULL)

      pred_df <- data.frame(
        n_donors = donor_seq,
        log_donors = log(donor_seq),
        log_n_cells = log(PRACTICAL_N_CELLS),
        log_reads_per_cell = log(PRACTICAL_RPC),
        log_abs_gt_lfc = log(PRACTICAL_LFC),
        log_gt_basemean = log(PRACTICAL_BASEMEAN),
        log_gt_dispersion = log(PRACTICAL_DISPERSION)
      )

      pred_df$predicted <- predict(model, newdata = pred_df, type = "response")
      pred_df$model <- model_name
      pred_df
    })

    simple_df <- bind_rows(simple_data[!sapply(simple_data, is.null)])
    simple_df$model <- factor(simple_df$model, levels = ALL_MODELS)

    p_simple <- ggplot(simple_df, aes(x = n_donors, y = predicted,
                                       colour = model, linetype = model)) +
      geom_line(linewidth = 1.2) +
      # Only the 80% reference line — match panel A.
      geom_hline(yintercept = 0.8, linetype = "dotted", alpha = 0.6) +
      scale_x_log10(labels = comma) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent) +
      scale_colour_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
      scale_linetype_manual(values = MODEL_LINETYPES, labels = prettify_labeller()) +
      labs(
        x = "Number of Donors",
        y = "Predicted Detection Probability",
        colour = "Cell Type",
        linetype = "Cell Type"
      ) +
      theme_publication(base_size = 16) +
      theme(legend.position = "bottom",
            legend.box = "horizontal") +
      guides(colour = guide_legend(nrow = 2),
             linetype = guide_legend(nrow = 2))

    ggsave(file.path(variant_plots_dir, "CellType_Sensitivity_Simple.png"),
           p_simple, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    # ------------------------------------------------------------------
    # 7b. PLOT 2: Coefficient Forest Plots
    # ------------------------------------------------------------------
    
    log_info("Generating coefficient forest plots...")
    
    # Load coefficient data
    model_dir <- get_model_dir(gene_set, agg_level)
    coef_path <- file.path(model_dir, "all_coefficients.csv")
    
    if (file.exists(coef_path)) {
      coef_df <- read.csv(coef_path, stringsAsFactors = FALSE)
      
      # Filter to available models only, exclude intercepts and interaction terms
      coef_df <- coef_df %>%
        filter(model %in% available_models) %>%
        filter(!grepl("^\\(|^cell_type|:", term)) %>%
        mutate(
          model = factor(model, levels = ALL_MODELS),
          ci_low = estimate - qnorm(1 - (1 - CI_LEVEL) / 2) * std_error,
          ci_high = estimate + qnorm(1 - (1 - CI_LEVEL) / 2) * std_error
        )
      
      predictors <- unique(coef_df$term)
      
      for (term_name in predictors) {
        term_data <- coef_df %>% filter(term == term_name)
        if (nrow(term_data) == 0) next
        
        p_forest <- ggplot(term_data, aes(x = estimate, y = model,
                                           colour = model, shape = model)) +
          geom_point(size = 4) +
          geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          scale_colour_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
          scale_shape_manual(values = MODEL_SHAPES, labels = prettify_labeller()) +
          labs(
            title = sprintf("Coefficient: %s", prettify_label(term_name)),
            subtitle = "Error bars show 95% CI",
            x = "Estimate",
            y = "Model"
          ) +
          theme_publication(base_size = 16) +
          theme(legend.position = "none")
        
        safe_term <- gsub("[^[:alnum:]]", "_", term_name)
        ggsave(file.path(variant_plots_dir, sprintf("Forest_%s.png", safe_term)),
               p_forest, width = PLOT_WIDTH * 0.8, height = PLOT_HEIGHT * 0.6, dpi = PLOT_DPI)
      }
      
      # Combined forest plot (all coefficients)
      coef_df$term_pretty <- prettify_label(coef_df$term)
      p_forest_all <- ggplot(coef_df, aes(x = estimate, y = term_pretty,
                                            colour = model, shape = model)) +
        geom_point(size = 3, position = position_dodge(width = 0.6)) +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2,
                       position = position_dodge(width = 0.6)) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_colour_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
        scale_shape_manual(values = MODEL_SHAPES, labels = prettify_labeller()) +
        labs(
          title = "Logistic Regression Coefficients",
          subtitle = NULL,
          x = "Coefficient Estimate (95% CI)",
          y = "Predictor",
          colour = "Cell Type",
          shape = "Cell Type"
        ) +
        theme_publication(base_size = 16) +
        theme(legend.position = "right")
      
      ggsave(file.path(variant_plots_dir, "Forest_All_Coefficients.png"),
             p_forest_all, width = PLOT_WIDTH * 1.2, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    }
    
    # ------------------------------------------------------------------
    # 7c. PLOT 3: R² Comparison Bar Chart
    # ------------------------------------------------------------------
    
    log_info("Generating R² comparison plot...")
    
    r2_values <- sapply(available_models, function(model_name) {
      model <- models[[model_name]]
      if (is.null(model)) return(NA_real_)
      mcfadden_r2(model)
    })
    
    r2_df <- data.frame(
      model = factor(names(r2_values), levels = ALL_MODELS),
      r2 = as.numeric(r2_values)
    )
    
    p_r2 <- ggplot(r2_df, aes(x = model, y = r2, fill = model)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.3f", r2)), vjust = -0.5, size = 5) +
      scale_fill_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
      scale_x_discrete(labels = prettify_labeller()) +
      scale_y_continuous(limits = c(0, max(r2_df$r2, na.rm = TRUE) * 1.15),
                         labels = scales::percent) +
      labs(
        title = "Pseudo-R² Comparison",
        subtitle = "McFadden's R² — higher is better",
        x = "Model",
        y = "Pseudo-R²"
      ) +
      theme_publication(base_size = 16) +
      theme(legend.position = "none")
    
    ggsave(file.path(variant_plots_dir, "R2_Comparison.png"),
           p_r2, width = PLOT_WIDTH * 0.8, height = PLOT_HEIGHT * 0.8, dpi = PLOT_DPI)
    
    # ------------------------------------------------------------------
    # 7d. PLOT 4: Required Sample Size Comparison
    # ------------------------------------------------------------------
    
    log_info("Generating required sample size comparison...")
    
    required_n <- expand.grid(
      model = available_models,
      target_power = POWER_TARGETS_COMPARE,
      stringsAsFactors = FALSE
    )
    
    required_n$required_donors <- mapply(function(m, tp) {
      find_required_donors(models[[m]], tp, rep_vals)
    }, required_n$model, required_n$target_power)
    
    required_n$model <- factor(required_n$model, levels = ALL_MODELS)
    required_n$target_label <- sprintf("%.0f%% Power", required_n$target_power * 100)
    
    p_required <- ggplot(required_n, aes(x = target_label, y = required_donors, fill = model)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = ifelse(is.na(required_donors), ">1000", as.character(required_donors))),
                position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
      scale_fill_manual(values = MODEL_COLOURS, labels = prettify_labeller()) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Required Donors to Achieve Target Power",
        subtitle = sprintf("|LFC| = %.1f, %s cells/donor, %s reads/cell",
                           rep_vals$median_lfc,
                           format(round(rep_vals$median_n_cells), big.mark = ","),
                           format(round(rep_vals$median_rpc), big.mark = ",")),
        x = "Target Power",
        y = "Required Number of Donors",
        fill = "Model"
      ) +
      theme_publication(base_size = 16) +
      theme(legend.position = "right")
    
    ggsave(file.path(variant_plots_dir, "Required_Sample_Size.png"),
           p_required, width = PLOT_WIDTH, height = PLOT_HEIGHT * 0.8, dpi = PLOT_DPI)
    
    # ------------------------------------------------------------------
    # 7e. PLOT 5: Power Surface Difference Heatmap (Astro - Astro_Shallow)
    # ------------------------------------------------------------------
    
    if ("Astro" %in% available_models && "Astro_Shallow" %in% available_models) {
      log_info("Generating power difference heatmap...")
      
      # Create grid over donors × cells
      donor_grid <- DONOR_STEPS
      cell_grid <- CELL_STEPS[!is.na(CELL_STEPS)]
      
      grid_df <- expand.grid(
        n_donors = donor_grid,
        n_cells = cell_grid
      )
      
      grid_df$log_donors <- log(grid_df$n_donors)
      grid_df$log_n_cells <- log(grid_df$n_cells)
      grid_df$log_reads_per_cell <- log(rep_vals$median_rpc)
      grid_df$log_abs_gt_lfc <- log(rep_vals$median_lfc)
      grid_df$log_gt_basemean <- log(rep_vals$median_basemean)
      grid_df$log_gt_dispersion <- log(rep_vals$median_dispersion)
      
      # Predict with both models
      grid_df$power_astro <- predict(models[["Astro"]], newdata = grid_df, type = "response")
      grid_df$power_shallow <- predict(models[["Astro_Shallow"]], newdata = grid_df, type = "response")
      grid_df$power_diff <- grid_df$power_astro - grid_df$power_shallow
      
      # Format labels
      grid_df$donors_label <- factor(grid_df$n_donors)
      grid_df$cells_label <- factor(grid_df$n_cells)
      
      p_heatmap <- ggplot(grid_df, aes(x = donors_label, y = cells_label, fill = power_diff)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = sprintf("%.2f", power_diff)), size = 3.5, colour = "black") +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                              midpoint = 0, limits = c(-0.3, 0.3),
                              oob = scales::squish,
                              labels = scales::percent) +
        labs(
          title = "Power Difference: Astrocytes - Astro (Depth-Eq.)",
          subtitle = sprintf("Positive = Astrocytes higher power; |LFC| = %.1f, %s reads/cell",
                             rep_vals$median_lfc,
                             format(round(rep_vals$median_rpc), big.mark = ",")),
          x = "Number of Donors",
          y = "Cells per Donor",
          fill = "Power\nDifference"
        ) +
        theme_publication(base_size = 16) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        )
      
      ggsave(file.path(variant_plots_dir, "Power_Difference_Heatmap.png"),
             p_heatmap, width = PLOT_WIDTH * 1.2, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    }
    
    log_info("Visualisation complete for variant: %s", variant_label)
  }
}

# ==============================================================================
# 8. GENERATE SUMMARY STATISTICS TABLE
# ==============================================================================

log_section("Generating summary statistics")

summary_list <- list()

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {
    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    
    models <- load_all_models(gene_set, agg_level)
    fitting_data <- load_all_fitting_data(gene_set, agg_level)
    
    if (length(models) == 0) next
    
    rep_vals <- if (!is.null(fitting_data) && nrow(fitting_data) > 0) {
      compute_representative_values(fitting_data)
    } else {
      list(
        median_lfc = DEFAULT_REPRESENTATIVE_LFC,
        median_basemean = DEFAULT_REPRESENTATIVE_BASEMEAN,
        median_dispersion = DEFAULT_REPRESENTATIVE_DISPERSION,
        median_n_cells = 500,
        median_rpc = 6000
      )
    }
    
    for (model_name in names(models)) {
      model <- models[[model_name]]
      if (is.null(model)) next
      
      r2 <- mcfadden_r2(model)
      
      # Required donors for each power target
      req_50 <- find_required_donors(model, 0.50, rep_vals)
      req_80 <- find_required_donors(model, 0.80, rep_vals)
      req_90 <- find_required_donors(model, 0.90, rep_vals)
      
      summary_list[[length(summary_list) + 1]] <- data.frame(
        variant = variant_label,
        gene_set = gene_set,
        model = model_name,
        r2 = r2,
        req_donors_50pct = req_50,
        req_donors_80pct = req_80,
        req_donors_90pct = req_90,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  
  # Add Astro vs Astro_Shallow and Astro_Shallow_Fewer comparisons
  summary_wide <- summary_df %>%
    select(variant, model, r2, req_donors_80pct) %>%
    pivot_wider(names_from = model, values_from = c(r2, req_donors_80pct))
  
  # Astro vs Astro_Shallow (depth-equalized)
  if ("r2_Astro" %in% names(summary_wide) && "r2_Astro_Shallow" %in% names(summary_wide)) {
    summary_wide <- summary_wide %>%
      mutate(
        r2_diff_shallow = r2_Astro - r2_Astro_Shallow,
        req_donors_diff_shallow = req_donors_80pct_Astro - req_donors_80pct_Astro_Shallow
      )
  }
  
  # Astro vs Astro_Shallow_Fewer (depth + cell count matched)
  if ("r2_Astro" %in% names(summary_wide) && "r2_Astro_Shallow_Fewer" %in% names(summary_wide)) {
    summary_wide <- summary_wide %>%
      mutate(
        r2_diff_shallow_fewer = r2_Astro - r2_Astro_Shallow_Fewer,
        req_donors_diff_shallow_fewer = req_donors_80pct_Astro - req_donors_80pct_Astro_Shallow_Fewer
      )
  }
  
  # Astro_Shallow vs Astro_Shallow_Fewer (depth-only vs depth + count matched)
  if ("r2_Astro_Shallow" %in% names(summary_wide) && "r2_Astro_Shallow_Fewer" %in% names(summary_wide)) {
    summary_wide <- summary_wide %>%
      mutate(
        r2_diff_shallow_variants = r2_Astro_Shallow - r2_Astro_Shallow_Fewer,
        req_donors_diff_shallow_variants = req_donors_80pct_Astro_Shallow - req_donors_80pct_Astro_Shallow_Fewer
      )
  }
  
  # Save summaries
  write.csv(summary_df, file.path(SUMMARY_DIR_SHALLOW, "model_summary.csv"), row.names = FALSE)
  write.csv(summary_wide, file.path(SUMMARY_DIR_SHALLOW, "model_comparison_wide.csv"), row.names = FALSE)
  
  log_info("Saved summary tables to %s", SUMMARY_DIR_SHALLOW)
  
  # Print summary
  log_info("\n--- Model Summary ---")
  print(summary_df)
}

# ==============================================================================
# 9. COEFFICIENT DIFFERENCES TABLE (Astro vs shallow variants)
# ==============================================================================

log_section("Computing coefficient differences")

# Helper function to compute coefficient differences between two models
compute_coef_diff <- function(coef_df, model1, model2, variant_label) {
  coefs1 <- coef_df %>% filter(model == model1) %>% filter(!grepl("^\\(", term))
  coefs2 <- coef_df %>% filter(model == model2) %>% filter(!grepl("^\\(", term))
  
  if (nrow(coefs1) == 0 || nrow(coefs2) == 0) return(NULL)
  
  inner_join(
    coefs1 %>% select(term, estimate_1 = estimate, se_1 = std_error),
    coefs2 %>% select(term, estimate_2 = estimate, se_2 = std_error),
    by = "term"
  ) %>%
    mutate(
      variant = variant_label,
      comparison = sprintf("%s vs %s", model2, model1),
      diff = estimate_2 - estimate_1,
      pooled_se = sqrt(se_1^2 + se_2^2),
      z_diff = diff / pooled_se,
      p_diff = 2 * pnorm(-abs(z_diff)),
      significant = p_diff < 0.05
    )
}

coef_diff_list <- list()

for (gene_set in GENE_SETS) {
  for (agg_level in AGGREGATION_LEVELS) {
    variant_label <- sprintf("%s_%s", gene_set, agg_level)
    model_dir <- get_model_dir(gene_set, agg_level)
    coef_path <- file.path(model_dir, "all_coefficients.csv")
    
    if (!file.exists(coef_path)) next
    
    coef_df <- read.csv(coef_path, stringsAsFactors = FALSE)
    
    # Compare Astro vs Astro_Shallow
    diff_shallow <- compute_coef_diff(coef_df, "Astro", "Astro_Shallow", variant_label)
    if (!is.null(diff_shallow)) {
      coef_diff_list[[paste(variant_label, "shallow", sep = "_")]] <- diff_shallow
    }
    
    # Compare Astro vs Astro_Shallow_Fewer
    diff_fewer <- compute_coef_diff(coef_df, "Astro", "Astro_Shallow_Fewer", variant_label)
    if (!is.null(diff_fewer)) {
      coef_diff_list[[paste(variant_label, "fewer", sep = "_")]] <- diff_fewer
    }
    
    # Compare Astro_Shallow vs Astro_Shallow_Fewer
    diff_variants <- compute_coef_diff(coef_df, "Astro_Shallow", "Astro_Shallow_Fewer", variant_label)
    if (!is.null(diff_variants)) {
      coef_diff_list[[paste(variant_label, "variants", sep = "_")]] <- diff_variants
    }
  }
}

if (length(coef_diff_list) > 0) {
  coef_diff_df <- bind_rows(coef_diff_list)
  write.csv(coef_diff_df, file.path(SUMMARY_DIR_SHALLOW, "coefficient_differences.csv"), row.names = FALSE)
  
  log_info("Saved coefficient differences to %s/coefficient_differences.csv", SUMMARY_DIR_SHALLOW)
  
  # Print key differences by comparison type
  for (comp in unique(coef_diff_df$comparison)) {
    log_info("\n--- Coefficient Differences (%s) ---", comp)
    comp_data <- coef_diff_df %>% 
      filter(comparison == comp) %>%
      select(variant, term, diff, z_diff, significant)
    print(comp_data)
  }
}

# ==============================================================================
# 10. COMPLETION
# ==============================================================================

log_section("STEP 6 COMPLETE")
log_info("All comparison plots generated in: %s", PLOTS_DIR_SHALLOW)
log_info("Summary tables generated in: %s", SUMMARY_DIR_SHALLOW)
