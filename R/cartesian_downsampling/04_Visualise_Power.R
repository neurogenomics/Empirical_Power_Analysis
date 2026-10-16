# ==============================================================================
# SCRIPT: Visualize Power Analysis Results
# PURPOSE: Generate heatmaps for precision, recall, and other metrics
# INPUT: Aggregated summary data from Step 3
# OUTPUT: PNG heatmaps in plots directory
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required. Install: install.packages('ggnewscale')")
  }
  library(ggnewscale)
  
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",  # From project root
    "config.R",                            # From script directory
    "../config.R"                          # From tests/ or toy/
  )
  
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  
  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

# Source shared prettify functions
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 4: Visualize Power")

# Mode-aware configuration
y_col <- get_y_col_name()
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Counts per Donor" else "Cells per Donor"
log_info("Mode: %s | Y-column: %s", DOWNSAMPLE_MODE, y_col)

# Mode-aware factor function
clean_y_factor <- if (DOWNSAMPLE_MODE == "reads") clean_reads_factor else clean_cells_factor

# ==============================================================================
# 2. PLOTTING FUNCTION
# ==============================================================================

#' Create metric heatmap with optional feasibility overlay
#' 
#' @param data Data frame with Donors, Y-axis column, and metric columns
#' @param val_col Column name for the metric to plot
#' @param title_str Plot title
#' @param fix_scale Whether to fix scale to 0-1
#' @param feas_data Optional feasibility data for overlay (cells mode only)
#' @param y_col_name Name of the Y-axis column ("Cells" or "Target_Reads")
#' @param y_axis_label Label for Y-axis
#' @return ggplot object
plot_metric_heatmap <- function(data, val_col, title_str,
                                 fix_scale = TRUE, feas_data = NULL,
                                 y_col_name = y_col, y_axis_label = y_label,
                                 facet_col = "LFC_Bin",
                                 facet_intervals = LFC_INTERVALS,
                                 is_proportion = TRUE,
                                 label_digits = 0) {

  # Facet labels are the bare bin names. Gene counts are intentionally omitted:
  # paper convention is NO "(n = x genes)" annotations on these grid heatmaps.
  data <- data %>%
    mutate(Facet_Display_Label = as.character(!!sym(facet_col)))

  # Order facets by bin definition order
  unique_bins <- rev(names(facet_intervals))
  ordered_labels <- unique(data$Facet_Display_Label[order(match(data[[facet_col]], unique_bins))])
  data$Facet_Display_Label <- factor(data$Facet_Display_Label, levels = ordered_labels)

  # Grey out entire facets that carry essentially no ground-truth genes (mean
  # gene count below threshold): metrics for a near-empty bin (e.g. High
  # expression within a small-LFC bin) are misleading, so blank the whole facet
  # to a uniform grey rather than showing spurious values.
  if ("Ref_GT_Count" %in% names(data)) {
    facet_gene_means <- tapply(data$Ref_GT_Count, data$Facet_Display_Label,
                               function(v) mean(v, na.rm = TRUE))
    empty_facets <- names(facet_gene_means)[is.nan(facet_gene_means) |
                                              facet_gene_means < EMPTY_FACET_GENE_THRESHOLD]
    if (length(empty_facets) > 0) {
      data[[val_col]][as.character(data$Facet_Display_Label) %in% empty_facets] <- NA
    }
  }

  # Join feasibility data if provided (cells mode only — feas_data has "Cells" column)
  if (!is.null(feas_data) && y_col_name == "Cells") {
    data <- data %>% left_join(feas_data, by = c("Donors", "Cells"))
  }

  # Base plot (dynamic Y-axis column)
  p <- ggplot(data, aes(x = Donors, y = !!sym(y_col_name))) +
    geom_tile(aes(fill = !!sym(val_col)), color = "white") +
    facet_wrap(~Facet_Display_Label, ncol = 2) +
    labs(title = title_str, x = "Donors", y = y_axis_label) +
    theme_publication(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 13),
      panel.spacing = unit(0.1, "lines"),
      # Harmonised with 07_Visualise_Continuous_Grid.R so the grid + surface
      # composites in R/figure_composition have matching legend widths.
      legend.position = "right",
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(1.3, "cm"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  # Derive a clean legend title from the metric column name
  legend_title <- gsub("^Mean_", "Mean ", val_col)
  legend_title <- gsub("_Harm$", " (Harmonised)", legend_title)
  legend_title <- gsub("_", " ", legend_title)

  # Color scale (order = 1 ensures gradient legend appears above threshold legend)
  if (fix_scale) {
    p <- p + scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                                   na.value = "grey50", name = legend_title,
                                   guide = guide_colorbar(order = 1))
  } else {
    p <- p + scale_fill_viridis_c(option = "magma", na.value = "grey50", name = legend_title,
                                   guide = guide_colorbar(order = 1))
  }
  
  # Threshold and feasibility overlays only apply to proportion metrics
  if (is_proportion) {
    # Determine metric label for threshold legend (Recall or Precision)
    threshold_metric_label <- if (grepl("Recall", val_col)) "Recall" else "Precision"

    # Add threshold factor column for legend integration
    data <- data %>%
      mutate(
        Threshold_Met = factor(
          ifelse(!is.na(!!sym(val_col)) & !!sym(val_col) >= 0.8,
                 paste0("\u2265 80% ", threshold_metric_label),
                 paste0("< 80% ", threshold_metric_label)),
          levels = c(paste0("\u2265 80% ", threshold_metric_label),
                     paste0("< 80% ", threshold_metric_label))
        )
      )

    # If no cells cross the 80% threshold, skip the threshold legend entirely
    any_threshold_met <- any(!is.na(data[[val_col]]) & data[[val_col]] >= 0.8)

    # Overlay tiles with colour mapped to threshold status for legend
    p <- p +
      new_scale_color() +
      geom_tile(data = data, aes(color = Threshold_Met), fill = NA, linewidth = 0.6,
                width = 0.85, height = 0.85) +
      scale_color_manual(
        name = "Threshold",
        values = setNames(
          c("#00FFFF", "grey70"),
          c(paste0("\u2265 80% ", threshold_metric_label),
            paste0("< 80% ", threshold_metric_label))
        ),
        limits = c(paste0("\u2265 80% ", threshold_metric_label),
                   paste0("< 80% ", threshold_metric_label)),
        drop = FALSE,
        guide = if (any_threshold_met) {
          guide_legend(order = 2, override.aes = list(fill = c(NA, "grey85"), color = c("#00FFFF", "grey70"), linewidth = 1.2))
        } else {
          "none"
        }
      )

    # Add feasibility overlay for recall plots (cells mode only)
    if (grepl("Mean_Recall", val_col) && !is.null(feas_data) && y_col_name == "Cells") {
      highlight_subset <- data %>%
        filter(!is.na(Feasibility_Count) & Feasibility_Count > 0,
               !is.na(!!sym(val_col)) & !!sym(val_col) >= 0.8)

      if (nrow(highlight_subset) > 0) {
        highlight_subset$Feasibility_Count <- factor(highlight_subset$Feasibility_Count)
        discrete_blues <- c("1" = "#00FFFF", "2" = "#33CCFF", "3" = "#3399FF",
                            "4" = "#0066FF", "5" = "#0000FF")

        p <- p +
          new_scale_color() +
          geom_tile(data = highlight_subset,
                    aes(x = Donors, y = !!sym(y_col_name), color = Feasibility_Count),
                    fill = NA, linewidth = 1.2) +
          scale_color_manual(name = "Datasets\nAvailable", values = discrete_blues)
      }
    }
  }
  
  # Add value labels with adaptive text colour (dark on light tiles, white on dark)
  # Threshold ~0.45 on viridis/magma: values above are dark enough for white text
  if (is_proportion) {
    data$.text_col <- ifelse(is.na(data[[val_col]]) | data[[val_col]] > 0.45, "black", "white")
  } else {
    # For count data, normalise to 0-1 range for colour threshold
    val_range <- range(data[[val_col]], na.rm = TRUE)
    val_norm <- if (diff(val_range) > 0) {
      (data[[val_col]] - val_range[1]) / diff(val_range)
    } else {
      rep(0.5, nrow(data))
    }
    data$.text_col <- ifelse(is.na(data[[val_col]]) | val_norm > 0.45, "black", "white")
  }
  for (tc in unique(data$.text_col)) {
    sub <- data[data$.text_col == tc, ]
    if (nrow(sub) == 0) next
    if (is_proportion) {
      p <- p + geom_text(data = sub,
                         aes(label = ifelse(is.na(!!sym(val_col)), "",
                                            round(!!sym(val_col) * 100))),
                         color = tc, size = 3.0, show.legend = FALSE)
    } else {
      p <- p + geom_text(data = sub,
                         aes(label = ifelse(is.na(!!sym(val_col)), "",
                                            round(!!sym(val_col), !!label_digits))),
                         color = tc, size = 3.0, show.legend = FALSE)
    }
  }
  
  return(p)
}

# ==============================================================================
# 3. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

# Store read stats for global normalization
global_read_stats <- list()

for (cell_type in cell_types) {

  # Collect per-method pval_q summaries for the aggregated across-method plot
  pval_q_method_parts <- list()

  for (method in get_viz_methods()) {

  log_section(sprintf("Cell Type: %s | Method: %s", cell_type, method))

  # Get paths
  paths <- get_paths(cell_type, method = method)

  # Skip if visualisation already complete (no-overwrite mode)
  viz_sentinel <- file.path(paths$plots_dir, ".step4_visualisation_complete")
  if (!FORCE_RECALC && file.exists(viz_sentinel)) {
    log_info("Skipping %s — visualisation already complete (sentinel exists)", cell_type)
    next
  }

  # Load checkpoint
  checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
  
  if (!file.exists(checkpoint_file)) {
    log_warning("Checkpoint not found for %s, run Step 3 first", cell_type)
    next
  }
  
  checkpoint_data <- qs::qread(checkpoint_file)
  df_summary <- checkpoint_data$df_summary
  df_xist_summary <- checkpoint_data$df_xist_summary
  df_reads_summary <- checkpoint_data$df_reads_summary
  df_bm_summary <- checkpoint_data$df_bm_summary
  df_cart_summary <- checkpoint_data$df_cart_summary
  df_exprfilt_summary <- checkpoint_data$df_exprfilt_summary
  df_pval_q_summary <- checkpoint_data$df_pval_q_summary

  log_info("Loaded summary data: %d rows", nrow(df_summary))
  if (!is.null(df_bm_summary)) log_info("BaseMean summary: %d rows", nrow(df_bm_summary))
  if (!is.null(df_cart_summary)) log_info("Grid summary: %d rows", nrow(df_cart_summary))
  if (!is.null(df_exprfilt_summary)) log_info("Expr-filtered LFC summary: %d rows", nrow(df_exprfilt_summary))
  if (!is.null(df_pval_q_summary)) log_info("P-value quartile summary: %d rows", nrow(df_pval_q_summary))

  # Check if Sub_Count data is available (requires step 3 re-run after v2.3 update)
  has_sub_count <- "Mean_Sub_Count" %in% names(df_summary)
  if (!has_sub_count) {
    log_warning("Mean_Sub_Count not found in checkpoint — re-run Step 3 to generate DEG count heatmaps")
  }
  has_sub_count_lowexpr <- "Mean_Sub_Count_LowExpr" %in% names(df_summary)
  if (!has_sub_count_lowexpr) {
    log_warning("Mean_Sub_Count_LowExpr not found in checkpoint — re-run Step 3 to generate low-expression DEG count heatmaps")
  }

  # Derive F1 score and log DEG count columns on all summary dataframes
  add_derived_metrics <- function(df, has_sub) {
    if (is.null(df)) return(df)
    # F1 = harmonic mean of precision and recall
    df$Mean_F1 <- ifelse(
      !is.na(df$Mean_Precision) & !is.na(df$Mean_Recall) &
        (df$Mean_Precision + df$Mean_Recall) > 0,
      2 * df$Mean_Precision * df$Mean_Recall / (df$Mean_Precision + df$Mean_Recall),
      NA_real_
    )
    if ("Mean_Precision_Harm" %in% names(df) && "Mean_Recall_Harm" %in% names(df)) {
      df$Mean_F1_Harm <- ifelse(
        !is.na(df$Mean_Precision_Harm) & !is.na(df$Mean_Recall_Harm) &
          (df$Mean_Precision_Harm + df$Mean_Recall_Harm) > 0,
        2 * df$Mean_Precision_Harm * df$Mean_Recall_Harm / (df$Mean_Precision_Harm + df$Mean_Recall_Harm),
        NA_real_
      )
    }
    if (has_sub && "Mean_Sub_Count" %in% names(df)) {
      df$Log_Sub_Count <- ifelse(
        !is.na(df$Mean_Sub_Count) & df$Mean_Sub_Count > 0,
        log10(df$Mean_Sub_Count),
        NA_real_
      )
    }
    if ("Mean_Sub_Count_LowExpr" %in% names(df)) {
      df$Log_Sub_Count_LowExpr <- ifelse(
        !is.na(df$Mean_Sub_Count_LowExpr) & df$Mean_Sub_Count_LowExpr > 0,
        log10(df$Mean_Sub_Count_LowExpr),
        NA_real_
      )
    }
    df
  }

  df_summary <- add_derived_metrics(df_summary, has_sub_count)
  df_bm_summary <- add_derived_metrics(df_bm_summary, has_sub_count)
  df_cart_summary <- add_derived_metrics(df_cart_summary, has_sub_count)
  df_exprfilt_summary <- add_derived_metrics(df_exprfilt_summary, has_sub_count)
  df_pval_q_summary <- add_derived_metrics(df_pval_q_summary, has_sub_count)
  
  # Load feasibility data (optional, cells mode only — feasibility data uses cell counts)
  feas_map <- if (DOWNSAMPLE_MODE == "cells") load_feasibility_map(cell_type) else NULL
  
  # ==============================================================================
  # 4. PREPARE DATA FOR PLOTTING
  # ==============================================================================
  
  # Convert to factors (mode-aware)
  df_summary$Donors <- clean_donors_factor(df_summary$Donors)
  df_summary[[y_col]] <- clean_y_factor(df_summary[[y_col]])
  
  # Complete the grid (fill missing combinations with NA)
  df_summary <- df_summary %>%
    complete(Donors, !!sym(y_col), LFC_Bin, Gene_Category)
  
  df_xist_summary$Donors <- clean_donors_factor(df_xist_summary$Donors)
  df_xist_summary[[y_col]] <- clean_y_factor(df_xist_summary[[y_col]])
  df_xist_summary <- df_xist_summary %>%
    complete(Donors, !!sym(y_col), LFC_Bin)
  
  # Store read stats for later
  global_read_stats[[paste(cell_type, method, sep = "_")]] <- list(
    data = df_reads_summary,
    plot_dir = paths$plots_dir,
    cell_type_name = sprintf("%s (%s)", cell_type, method)
  )
  
  # ==============================================================================
  # 5. GENERATE HEATMAPS
  # ==============================================================================
  
  metrics_to_plot <- list(
    list(col = "Mean_Precision", desc = "Mean Precision", fix = TRUE, is_prop = TRUE, digits = 0),
    list(col = "Mean_Recall", desc = "Mean Recall", fix = TRUE, is_prop = TRUE, digits = 0),
    list(col = "Mean_F1", desc = "Mean F1 Score", fix = TRUE, is_prop = TRUE, digits = 0)
  )
  if (has_sub_count) {
    metrics_to_plot <- c(metrics_to_plot, list(
      list(col = "Mean_Sub_Count", desc = "Mean DEGs Found", fix = FALSE, is_prop = FALSE, digits = 0),
      list(col = "Log_Sub_Count", desc = "Mean DEGs Found (log10)", fix = FALSE, is_prop = FALSE, digits = 1)
    ))
  }
  if (has_sub_count_lowexpr) {
    metrics_to_plot <- c(metrics_to_plot, list(
      list(col = "Mean_Sub_Count_LowExpr", desc = "Mean DEGs Found (baseMean<50)", fix = FALSE, is_prop = FALSE, digits = 0),
      list(col = "Log_Sub_Count_LowExpr", desc = "Mean DEGs Found (baseMean<50, log10)", fix = FALSE, is_prop = FALSE, digits = 1)
    ))
  }
  
  for (cat in CHR_CATEGORIES) {
    
    df_cat <- df_summary %>% filter(Gene_Category == cat)
    
    if (nrow(df_cat) == 0) {
      log_warning("No data for category %s", cat)
      next
    }
    
    for (m in metrics_to_plot) {
      
      fname <- paste0(m$col, "_", cat, ".png")
      fpath <- file.path(paths$plots_dir, fname)
      
      # Build clean title: omit "Total" category, show others
      cat_suffix <- if (cat == "Total") "" else paste0(" — ", cat)
      plot_title <- paste0(m$desc, ": ", prettify_label(cell_type), cat_suffix)

      p <- plot_metric_heatmap(
        data = df_cat,
        val_col = m$col,
        title_str = plot_title,
        fix_scale = m$fix,
        feas_data = feas_map,
        is_proportion = if (!is.null(m$is_prop)) m$is_prop else TRUE,
        label_digits = if (!is.null(m$digits)) m$digits else 0
      )
      
      # Dimensions harmonised with 07 semi-log output (14x10) so composite
      # panels line up in R/figure_composition.
      ggsave(fpath, p, width = 14, height = 10, dpi = 150)
      log_info("Saved: %s", fname)
    }
  }

  # ==============================================================================
  # 5b. GENERATE INDIVIDUAL LFC BIN PLOTS (Total only, with 0.8 threshold line)
  # ==============================================================================
  
  log_info("Generating individual LFC bin plots for Total category...")
  
  # Create subdirectory for individual bin plots
  individual_bins_dir <- file.path(paths$plots_dir, "Individual_Bins")
  if (!dir.exists(individual_bins_dir)) {
    dir.create(individual_bins_dir, recursive = TRUE)
    log_info("Created individual bins directory: %s", individual_bins_dir)
  }
  
  # Only Total category, only Recall and Precision (no harmonized, no autosomal/sexchr)
  simple_metrics <- list(
    list(col = "Mean_Recall", desc = "Mean Recall"),
    list(col = "Mean_Precision", desc = "Mean Precision")
  )
  
  df_total <- df_summary %>% filter(Gene_Category == "Total")
  
  if (nrow(df_total) > 0) {
    
    for (m in simple_metrics) {
      
      val_col <- m$col
      
      # Generate separate plot for each LFC bin
      for (bin_name in unique(df_total$LFC_Bin)) {
        
        df_bin <- df_total %>% filter(LFC_Bin == bin_name)
        
        if (nrow(df_bin) == 0) next
        
        # Prepare data for UNIFORM DISCRETE heatmap (same as faceted plots)
        # Keep as factors for uniform grid spacing, not numeric for axis-matched spacing
        df_plot <- df_bin %>%
          filter(!is.na(Donors), 
                 !is.na(!!sym(y_col)),
                 as.character(!!sym(y_col)) != "Max",
                 !is.na(!!sym(val_col)))
        
        # Get GT count info
        bin_info <- df_bin %>%
          summarise(
            Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1),
            Str_GT = dplyr::first(GT_Count_String),
            .groups = "drop"
          )
        
        # Build individual bin plot with uniform discrete grid (like main plots)
        # Using geom_tile with factors creates uniform spacing
        p_bin <- ggplot(df_plot, aes(x = Donors, y = !!sym(y_col), fill = !!sym(val_col))) +
          geom_tile(color = "white") +
          scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                               na.value = "grey50", name = m$desc) +
          coord_cartesian(clip = "off") +  # Prevent edge clipping
          labs(
            title = paste0(m$desc, ": ", prettify_label(cell_type)),
            x = "Number of Donors",
            y = y_label
          ) +
          theme_publication(base_size = 20) +  # Larger text for readability
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            panel.grid = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            plot.title = element_text(face = "bold", size = 22)
          )
        
        # Adaptive text colour for individual bin plots
        df_plot$.text_col <- ifelse(is.na(df_plot[[val_col]]) | df_plot[[val_col]] > 0.45, "black", "white")
        for (.tc in unique(df_plot$.text_col)) {
          .sub <- df_plot[df_plot$.text_col == .tc, ]
          if (nrow(.sub) > 0)
            p_bin <- p_bin + geom_text(data = .sub,
              aes(label = round(!!sym(val_col) * 100)),
              color = .tc, size = 4.5, show.legend = FALSE)
        }

        # Save - clean filename: remove +, spaces, parens, dots, and other special chars
        bin_suffix <- gsub("[^a-zA-Z0-9_-]", "", gsub(" ", "_", bin_name))
        fname <- paste0(val_col, "_Total_", bin_suffix, "_Individual_Discrete.png")
        fpath <- file.path(individual_bins_dir, fname)

        ggsave(fpath, p_bin, width = 14, height = 12, dpi = 200)  # Wider for readability
        log_info("Saved individual discrete plot: %s", fname)
      }
    }
  }
  
  log_info("Completed individual discrete plots for %s", cell_type)

  # ==============================================================================
  # 5c. BASEMEAN-BINNED HEATMAPS
  # ==============================================================================

  if (!is.null(df_bm_summary)) {
    log_info("Generating BaseMean-binned heatmaps...")

    df_bm_summary$Donors <- clean_donors_factor(df_bm_summary$Donors)
    df_bm_summary[[y_col]] <- clean_y_factor(df_bm_summary[[y_col]])
    df_bm_summary <- df_bm_summary %>%
      complete(Donors, !!sym(y_col), BaseMean_Bin, Gene_Category)

    for (cat in CHR_CATEGORIES) {
      df_bm_cat <- df_bm_summary %>% filter(Gene_Category == cat)
      if (nrow(df_bm_cat) == 0) next

      for (m in metrics_to_plot) {
        cat_suffix <- if (cat == "Total") "" else paste0(" - ", cat)
        plot_title <- paste0(m$desc, " by Gene Counts: ", prettify_label(cell_type), cat_suffix)

        p <- plot_metric_heatmap(
          data = df_bm_cat,
          val_col = m$col,
          title_str = plot_title,
            fix_scale = m$fix,
          feas_data = feas_map,
          facet_col = "BaseMean_Bin",
          facet_intervals = BASEMEAN_INTERVALS,
          is_proportion = if (!is.null(m$is_prop)) m$is_prop else TRUE,
          label_digits = if (!is.null(m$digits)) m$digits else 0
        )

        fname <- paste0(m$col, "_BaseMean_", cat, ".png")
        ggsave(file.path(paths$plots_dir, fname), p, width = 14, height = 10, dpi = 150)
        log_info("Saved: %s", fname)
      }
    }
  }

  # ==============================================================================
  # 5d. GRID (LFC x BASEMEAN) HEATMAPS
  # ==============================================================================

  if (!is.null(df_cart_summary)) {
    log_info("Generating grid (LFC x BaseMean) heatmaps...")

    # Build ordered interval list for facet ordering
    cart_intervals <- list()
    for (lfc_n in names(LFC_COARSE)) {
      for (bm_n in names(BASEMEAN_COARSE)) {
        label <- paste0(lfc_n, " / ", bm_n)
        cart_intervals[[label]] <- c(0, Inf)
      }
    }

    df_cart_summary$Donors <- clean_donors_factor(df_cart_summary$Donors)
    df_cart_summary[[y_col]] <- clean_y_factor(df_cart_summary[[y_col]])
    df_cart_summary <- df_cart_summary %>%
      complete(Donors, !!sym(y_col), Cartesian_Bin, Gene_Category)

    for (cat in CHR_CATEGORIES) {
      df_cart_cat <- df_cart_summary %>% filter(Gene_Category == cat)
      if (nrow(df_cart_cat) == 0) next

      for (m in metrics_to_plot) {
        cat_suffix <- if (cat == "Total") "" else paste0(" - ", cat)
        plot_title <- paste0(m$desc, " (LFC x Gene Counts): ", prettify_label(cell_type), cat_suffix)

        p <- plot_metric_heatmap(
          data = df_cart_cat,
          val_col = m$col,
          title_str = plot_title,
            fix_scale = m$fix,
          feas_data = feas_map,
          facet_col = "Cartesian_Bin",
          facet_intervals = cart_intervals,
          is_proportion = if (!is.null(m$is_prop)) m$is_prop else TRUE,
          label_digits = if (!is.null(m$digits)) m$digits else 0
        )

        fname <- paste0(m$col, "_Cartesian_", cat, ".png")
        ggsave(file.path(paths$plots_dir, fname), p, width = 14, height = 10, dpi = 150)
        log_info("Saved: %s", fname)
      }
    }
  }

  # ==============================================================================
  # 5e. EXPRESSION-FILTERED LFC HEATMAPS (baseMean > EXPR_FILTER_THRESHOLD)
  # ==============================================================================

  if (!is.null(df_exprfilt_summary)) {
    log_info("Generating expression-filtered LFC heatmaps (baseMean > %d)...", EXPR_FILTER_THRESHOLD)

    df_exprfilt_summary$Donors <- clean_donors_factor(df_exprfilt_summary$Donors)
    df_exprfilt_summary[[y_col]] <- clean_y_factor(df_exprfilt_summary[[y_col]])
    df_exprfilt_summary <- df_exprfilt_summary %>%
      complete(Donors, !!sym(y_col), LFC_Bin, Gene_Category)

    for (cat in CHR_CATEGORIES) {
      df_ef_cat <- df_exprfilt_summary %>% filter(Gene_Category == cat)
      if (nrow(df_ef_cat) == 0) next

      for (m in metrics_to_plot) {
        cat_suffix <- if (cat == "Total") "" else paste0(" - ", cat)
        plot_title <- paste0(m$desc, " (Expr > ", EXPR_FILTER_THRESHOLD, "): ",
                             prettify_label(cell_type), cat_suffix)

        p <- plot_metric_heatmap(
          data = df_ef_cat,
          val_col = m$col,
          title_str = plot_title,
            fix_scale = m$fix,
          feas_data = feas_map,
          facet_col = "LFC_Bin",
          facet_intervals = LFC_INTERVALS,
          is_proportion = if (!is.null(m$is_prop)) m$is_prop else TRUE,
          label_digits = if (!is.null(m$digits)) m$digits else 0
        )

        fname <- paste0(m$col, "_ExprFilt_", cat, ".png")
        ggsave(file.path(paths$plots_dir, fname), p, width = 14, height = 10, dpi = 150)
        log_info("Saved: %s", fname)
      }
    }
  }

  # ==============================================================================
  # 5f. P-VALUE QUARTILE HEATMAPS (sensitivity check vs LFC binning)
  # Bins Q1-Q4 are defined by GT padj quartiles; Q1 = most significant 25%.
  # Default: Recall + Precision. Set PVAL_Q_FULL_PLOTS=true for F1 too.
  # ==============================================================================
  if (!is.null(df_pval_q_summary)) {
    log_info("Generating p-value quartile heatmaps...")

    df_pval_q_summary$Donors <- clean_donors_factor(df_pval_q_summary$Donors)
    df_pval_q_summary[[y_col]] <- clean_y_factor(df_pval_q_summary[[y_col]])
    # Complete the grid so incomplete regions render as grey (NA fill) tiles,
    # matching every other grid in the pipeline (Chris, 2026-08-17).
    df_pval_q_summary <- df_pval_q_summary %>%
      complete(Donors, !!sym(y_col), PVal_Quartile)

    # Preserve base quartile names before relabelling (needed for cross-method aggregation).
    df_pval_q_summary$PVal_Quartile_Base <- df_pval_q_summary$PVal_Quartile

    # Build padj range labels using mean bounds stored per-quartile in the summary.
    # Replaces bare "Q1 (most sig)" etc. with e.g. "Q1 (most sig) (padj <= 1.2e-10)".
    if (all(c("PQ_Lo", "PQ_Hi") %in% colnames(df_pval_q_summary))) {
      pq_bounds_tbl <- df_pval_q_summary %>%
        group_by(PVal_Quartile) %>%
        summarise(PQ_Lo = mean(PQ_Lo, na.rm = TRUE),
                  PQ_Hi = mean(PQ_Hi, na.rm = TRUE), .groups = "drop")
      # Label using only the inner quartile boundaries (Q1's upper = q25, Q4's
      # lower = q75). Q1 uses "padj <= hi", Q4 uses "padj >= lo", Q2/Q3 the range.
      pq_label_vec <- setNames(
        mapply(function(nm, lo, hi) {
          if (grepl("^Q1", nm)) {
            sprintf("%s (padj \u2264 %.1e)", nm, hi)
          } else if (grepl("^Q4", nm)) {
            sprintf("%s (padj \u2265 %.1e)", nm, lo)
          } else {
            sprintf("%s (%.1e\u2013%.1e)", nm, lo, hi)
          }
        }, pq_bounds_tbl$PVal_Quartile, pq_bounds_tbl$PQ_Lo, pq_bounds_tbl$PQ_Hi),
        pq_bounds_tbl$PVal_Quartile
      )
      df_pval_q_summary$PVal_Quartile <- pq_label_vec[df_pval_q_summary$PVal_Quartile]
    }

    # Build display order from new labels — Q4 first (top-left), Q1 last (bottom-right)
    pval_q_new_labels <- unique(df_pval_q_summary$PVal_Quartile)
    pval_q_display_order <- setNames(
      replicate(length(pval_q_new_labels), c(0, 1), simplify = FALSE),
      sort(pval_q_new_labels, decreasing = TRUE)
    )

    pq_metrics <- list(
      list(col = "Mean_Recall",    desc = "Mean Recall"),
      list(col = "Mean_Precision", desc = "Mean Precision")
    )
    if (PVAL_Q_FULL_PLOTS) {
      pq_metrics <- c(pq_metrics, list(list(col = "Mean_F1", desc = "Mean F1")))
    }

    for (m in pq_metrics) {
      p <- plot_metric_heatmap(
        data = df_pval_q_summary,
        val_col = m$col,
        title_str = paste0(m$desc, " by padj Quartile: ", prettify_label(cell_type), " (", method, ")"),
        fix_scale = TRUE,
        feas_data = feas_map,
        facet_col = "PVal_Quartile",
        facet_intervals = pval_q_display_order,
        is_proportion = TRUE,
        label_digits = 0
      )
      fname <- paste0(m$col, "_PValQuartile_", method, ".png")
      ggsave(file.path(paths$plots_dir, fname), p, width = 14, height = 8, dpi = 150)
      log_info("Saved: %s", fname)
    }

    # Store for cross-method aggregated plot (generated after method loop).
    # Use PVal_Quartile_Base (Q1/Q2/Q3/Q4) not the range-labelled string, because
    # bounds differ slightly between methods and would create 8-12 spurious facet groups.
    pval_q_method_parts[[method]] <- df_pval_q_summary %>%
      mutate(Method = method,
             PVal_Quartile = PVal_Quartile_Base)
  }

  # ==============================================================================
  # 6. XIST SANITY CHECK PLOT
  # ==============================================================================
  
  p_xist <- ggplot(df_xist_summary %>% filter(!is.na(Donors), !is.na(!!sym(y_col))),
                   aes(x = Donors, y = !!sym(y_col), fill = Mean_Detection)) +
    geom_tile(color = "white") +
    facet_wrap(~LFC_Bin, ncol = 2) +
    labs(title = paste0(prettify_label(cell_type), ": XIST Detection"),
         x = "Donors", y = y_label, fill = "Detection\nRate") +
    theme_publication(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "left",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50") +
    geom_text(aes(label = ifelse(is.na(Mean_Detection), "", round(Mean_Detection * 100)),
                  color = ifelse(is.na(Mean_Detection) | Mean_Detection > 0.45, "black", "white")),
              size = 3.0, show.legend = FALSE) +
    scale_color_identity()
  
  ggsave(file.path(paths$plots_dir, "XIST_Sanity_Check.png"), p_xist, width = 12, height = 8)
  log_info("Saved: XIST_Sanity_Check.png")
  
  # Write sentinel so subsequent runs skip this cell type
  writeLines(as.character(Sys.time()), viz_sentinel)
  log_info("Completed visualization for %s/%s (sentinel written)", cell_type, method)
  }

  # ---------------------------------------------------------------------------
  # Aggregated p-value quartile plot across all methods (saved in first method's dir)
  # ---------------------------------------------------------------------------
  if (length(pval_q_method_parts) >= 1) {
    df_pval_q_all <- bind_rows(pval_q_method_parts)
    df_pval_q_agg <- df_pval_q_all %>%
      group_by(Donors, !!sym(y_col), PVal_Quartile) %>%
      summarise(
        Mean_Recall    = mean(Mean_Recall,    na.rm = TRUE),
        Mean_Precision = mean(Mean_Precision, na.rm = TRUE),
        Ref_GT_Count   = mean(Ref_GT_Count,   na.rm = TRUE),
        PQ_Lo          = mean(PQ_Lo,          na.rm = TRUE),
        PQ_Hi          = mean(PQ_Hi,          na.rm = TRUE),
        .groups = "drop"
      )

    # Relabel base quartiles (Q1-Q4) with the mean padj range across methods,
    # using only the inner quartile boundaries (Q1 "padj <= hi", Q4 "padj >= lo",
    # Q2/Q3 the range) to match the per-method plots.
    if (all(c("PQ_Lo", "PQ_Hi") %in% colnames(df_pval_q_agg))) {
      agg_bounds_tbl <- df_pval_q_agg %>%
        group_by(PVal_Quartile) %>%
        summarise(PQ_Lo = mean(PQ_Lo, na.rm = TRUE),
                  PQ_Hi = mean(PQ_Hi, na.rm = TRUE), .groups = "drop")
      agg_label_vec <- setNames(
        mapply(function(nm, lo, hi) {
          if (grepl("^Q1", nm)) {
            sprintf("%s (padj \u2264 %.1e)", nm, hi)
          } else if (grepl("^Q4", nm)) {
            sprintf("%s (padj \u2265 %.1e)", nm, lo)
          } else {
            sprintf("%s (%.1e\u2013%.1e)", nm, lo, hi)
          }
        }, agg_bounds_tbl$PVal_Quartile, agg_bounds_tbl$PQ_Lo, agg_bounds_tbl$PQ_Hi),
        agg_bounds_tbl$PVal_Quartile
      )
      df_pval_q_agg$PVal_Quartile <- agg_label_vec[df_pval_q_agg$PVal_Quartile]
    }

    agg_plots_dir <- get_paths(cell_type, method = CROSS_METHOD_LABEL)$plots_dir
    dir.create(agg_plots_dir, recursive = TRUE, showWarnings = FALSE)
    pval_q_agg_labels <- unique(df_pval_q_agg$PVal_Quartile)
    pval_q_display_order <- setNames(
      replicate(length(pval_q_agg_labels), c(0, 1), simplify = FALSE),
      sort(pval_q_agg_labels, decreasing = TRUE)
    )

    for (agg_m in list(
      list(col = "Mean_Recall",    desc = "Mean Recall",    fname = "Mean_Recall_PValQuartile_AllMethods.png"),
      list(col = "Mean_Precision", desc = "Mean Precision", fname = "Mean_Precision_PValQuartile_AllMethods.png")
    )) {
      p_agg <- plot_metric_heatmap(
        data = df_pval_q_agg,
        val_col = agg_m$col,
        title_str = paste0(agg_m$desc, " by padj Quartile: ", prettify_label(cell_type), " (All Methods)"),
        fix_scale = TRUE,
        feas_data = NULL,
        facet_col = "PVal_Quartile",
        facet_intervals = pval_q_display_order,
        is_proportion = TRUE,
        label_digits = 0
      )
      ggsave(file.path(agg_plots_dir, agg_m$fname), p_agg, width = 14, height = 8, dpi = 150)
      log_info("Saved aggregated pval_q %s plot for %s", agg_m$desc, cell_type)
    }
  }

}

# ==============================================================================
# 7. STANDARDIZED READ DEPTH PLOTS
# ==============================================================================

if (DOWNSAMPLE_MODE == "reads") {
  log_section("Skipping Read Depth Diagnostic Plots (reads mode — read depth is the independent variable)")
} else {
  log_section("Generating Standardized Read Depth Plots")
  
  if (length(global_read_stats) > 0) {
    
    # Calculate global min/max for consistent scale
    all_reads <- unlist(lapply(global_read_stats, function(x) x$data$Mean_Total_Reads))
    global_max_read <- max(all_reads, na.rm = TRUE)
    global_min_read <- min(all_reads[all_reads > 0], na.rm = TRUE)
    if (is.infinite(global_min_read)) global_min_read <- 1
    
    log_info("Global read range: %.0f - %.0f", global_min_read, global_max_read)
    
    for (ct in names(global_read_stats)) {
      
      item <- global_read_stats[[ct]]
      df_reads <- item$data
      plot_dir <- item$plot_dir
      cell_name <- item$cell_type_name
      
      df_reads$Donors <- clean_donors_factor(df_reads$Donors)
      df_reads[[y_col]] <- clean_y_factor(df_reads[[y_col]])
      df_reads <- df_reads %>% complete(Donors, !!sym(y_col))
      
      p_reads <- ggplot(df_reads, aes(x = Donors, y = !!sym(y_col), fill = Mean_Total_Reads)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c(option = "plasma", trans = "log10",
                             limits = c(global_min_read, global_max_read),
                             na.value = "grey90", name = "Total Counts") +
        geom_text(aes(label = ifelse(is.na(Mean_Total_Reads), "",
                                     scales::scientific(Mean_Total_Reads, digits = 2)),
                      color = ifelse(is.na(Mean_Total_Reads) |
                        (Mean_Total_Reads - global_min_read) / (global_max_read - global_min_read) > 0.45,
                        "black", "white")),
                  size = 3.0, show.legend = FALSE) +
        scale_color_identity() +
        labs(title = paste0(prettify_label(cell_name), ": Total Counts per Matrix"),
             x = "Number of Donors", y = y_label) +
        theme_publication(base_size = 16) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "left",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0)
        )
      
      ggsave(file.path(plot_dir, "Total_Reads_Per_Matrix.png"), p_reads, width = 8, height = 6)
      log_info("Saved: Total_Reads_Per_Matrix.png for %s", ct)
    }
  }
}

log_section("STEP 4 COMPLETE")
log_info("All visualizations generated")
