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
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Reads per Donor" else "Cells per Donor"
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
#' @param sub_str Plot subtitle
#' @param fix_scale Whether to fix scale to 0-1
#' @param feas_data Optional feasibility data for overlay (cells mode only)
#' @param y_col_name Name of the Y-axis column ("Cells" or "Target_Reads")
#' @param y_axis_label Label for Y-axis
#' @return ggplot object
plot_metric_heatmap <- function(data, val_col, title_str, sub_str, 
                                 fix_scale = TRUE, feas_data = NULL,
                                 y_col_name = y_col, y_axis_label = y_label) {
  
  # Add header labels with CV information
  bin_labels_df <- data %>%
    group_by(LFC_Bin) %>%
    summarise(
      Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 2),
      Mean_Harm = round(mean(Ref_Harm_Count, na.rm = TRUE), 2),
      CV_GT_Val = mean(CV_GT, na.rm = TRUE),
      CV_Harm_Val = mean(CV_Harm, na.rm = TRUE),
      Str_GT = dplyr::first(GT_Count_String),
      Str_Harm = dplyr::first(Harm_Count_String),
      .groups = "drop"
    )
  
  data <- data %>%
    dplyr::select(-any_of(c("Mean_GT", "Mean_Harm", "CV_GT_Val", "CV_Harm_Val", "Str_GT", "Str_Harm"))) %>%
    left_join(bin_labels_df, by = "LFC_Bin")
  
  is_harm_plot <- grepl("_Harm", val_col)
  
  data <- data %>%
    mutate(LFC_Display_Label = if (is_harm_plot) {
      sprintf("%s (n = %.0f genes)", LFC_Bin, Mean_Harm)
    } else {
      sprintf("%s (n = %.0f genes)", LFC_Bin, Mean_GT)
    })
  
  # Order facets by LFC bin
  unique_bins <- rev(names(LFC_INTERVALS))
  ordered_labels <- unique(data$LFC_Display_Label[order(match(data$LFC_Bin, unique_bins))])
  data$LFC_Display_Label <- factor(data$LFC_Display_Label, levels = ordered_labels)
  
  # Join feasibility data if provided (cells mode only — feas_data has "Cells" column)
  if (!is.null(feas_data) && y_col_name == "Cells") {
    data <- data %>% left_join(feas_data, by = c("Donors", "Cells"))
  }
  
  # Base plot (dynamic Y-axis column)
  p <- ggplot(data, aes(x = Donors, y = !!sym(y_col_name))) +
    geom_tile(aes(fill = !!sym(val_col)), color = "white") +
    facet_wrap(~LFC_Display_Label, ncol = 2) +
    labs(title = title_str, subtitle = sub_str, x = "Donors", y = y_axis_label) +
    theme_publication(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "left",
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

  # Overlay tiles with colour mapped to threshold status for legend
  p <- p +
    new_scale_color() +
    geom_tile(data = data, aes(color = Threshold_Met), fill = NA, linewidth = 0.6,
              width = 0.85, height = 0.85) +
    scale_color_manual(
      name = "Threshold",
      values = setNames(
        c("red", "grey70"),
        c(paste0("\u2265 80% ", threshold_metric_label),
          paste0("< 80% ", threshold_metric_label))
      ),
      drop = FALSE,
      guide = guide_legend(order = 2, override.aes = list(fill = c(NA, "grey85"), linewidth = 1.2))
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
  
  # Add value labels with adaptive text colour (dark on light tiles, white on dark)
  # Threshold ~0.45 on viridis/magma: values above are dark enough for white text
  data$.text_col <- ifelse(is.na(data[[val_col]]) | data[[val_col]] > 0.45, "black", "white")
  for (tc in unique(data$.text_col)) {
    sub <- data[data$.text_col == tc, ]
    if (nrow(sub) == 0) next
    p <- p + geom_text(data = sub,
                       aes(label = ifelse(is.na(!!sym(val_col)), "",
                                          round(!!sym(val_col) * 100))),
                       color = tc, size = 3.0, show.legend = FALSE)
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
  
  log_section(sprintf("Cell Type: %s", cell_type))
  
  # Get paths
  paths <- get_paths(cell_type)
  
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
  
  log_info("Loaded summary data: %d rows", nrow(df_summary))
  
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
  global_read_stats[[cell_type]] <- list(
    data = df_reads_summary,
    plot_dir = paths$plots_dir,
    cell_type_name = cell_type
  )
  
  # ==============================================================================
  # 5. GENERATE HEATMAPS
  # ==============================================================================
  
  metrics_to_plot <- list(
    list(col = "Mean_Precision", desc = "Mean Precision", fix = TRUE),
    list(col = "Mean_Recall", desc = "Mean Recall", fix = TRUE),
    list(col = "Mean_Recall_Harm", desc = "Mean Recall (Harmonised)", fix = TRUE)
  )
  
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
        sub_str = NULL,
        fix_scale = m$fix,
        feas_data = feas_map
      )
      
      ggsave(fpath, p, width = 12, height = 8)
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
            subtitle = paste0(bin_name, " (n = ", round(bin_info$Mean_GT), " genes)"),
            x = "Number of Donors",
            y = y_label
          ) +
          theme_publication(base_size = 20) +  # Larger text for readability
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            panel.grid = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 22),
            plot.subtitle = element_text(size = 16)
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
  
  log_info("Completed visualization for %s", cell_type)
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
                             na.value = "grey90", name = "Total Reads") +
        geom_text(aes(label = ifelse(is.na(Mean_Total_Reads), "",
                                     scales::scientific(Mean_Total_Reads, digits = 2)),
                      color = ifelse(is.na(Mean_Total_Reads) |
                        (Mean_Total_Reads - global_min_read) / (global_max_read - global_min_read) > 0.45,
                        "black", "white")),
                  size = 3.0, show.legend = FALSE) +
        scale_color_identity() +
        labs(title = paste0(prettify_label(cell_name), ": Total Reads per Matrix"),
             subtitle = NULL,
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
