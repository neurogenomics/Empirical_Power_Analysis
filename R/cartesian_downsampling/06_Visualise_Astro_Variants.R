# ==============================================================================
# SCRIPT: Visualize Astro Variant Datasets (Depth-Equalized and Fewer Cells)
# PURPOSE: Generate heatmaps for Astrocyte variant datasets and create
#          side-by-side comparison plots across standard, depth-equalized,
#          and fewer-cells configurations.
# INPUT: Aggregated summary data from cartesian downsampling for Astro variants
# OUTPUT: Individual heatmaps + comparison plots in Astro/Plots/Variants/
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required. Install: install.packages('ggnewscale')")
  }
  library(ggnewscale)
  
  if (!requireNamespace("shadowtext", quietly = TRUE)) {
    stop("Package 'shadowtext' is required. Install: install.packages('shadowtext')")
  }
  library(shadowtext)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",
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
  
  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 6: Visualize Astro Variant Datasets")

# Mode-aware configuration
y_col <- get_y_col_name()
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Reads per Donor" else "Cells per Donor"
clean_y_factor <- if (DOWNSAMPLE_MODE == "reads") clean_reads_factor else clean_cells_factor
log_info("Mode: %s | Y-column: %s", DOWNSAMPLE_MODE, y_col)

# ==============================================================================
# 2. ASTRO VARIANT CONFIGURATION
# ==============================================================================

# Define Astro variant datasets with display labels
ASTRO_VARIANTS <- list(
  list(
    dataset = "Roussos",
    label = "Standard",
    color = "#2E86AB",
    description = "Standard (Deep Sequencing)"
  ),
  list(
    dataset = "Roussos_DepthEq",
    label = "Depth-Equalized",
    color = "#A23B72",
    description = "Depth-Equalized to Micro (~4,426 reads/cell)"
  ),
  list(
    dataset = "Roussos_ShallowFewer",
    label = "Depth-Equalized + Fewer Cells",
    color = "#F18F01",
    description = "Depth-Equalized + Cell Count Matched to Micro"
  )
)

# Extract dataset names for processing
variant_datasets <- sapply(ASTRO_VARIANTS, function(x) x$dataset)
log_info("Processing %d Astro variants: %s", length(variant_datasets), 
         paste(variant_datasets, collapse = ", "))

# ==============================================================================
# 3. PLOTTING FUNCTIONS (adapted from 04_Visualise_Power.R)
# ==============================================================================

#' Create metric heatmap for variant comparison
#' 
#' @param data Data frame with Donors, Y-axis column, and metric columns
#' @param val_col Column name for the metric to plot
#' @param title_str Plot title
#' @param sub_str Plot subtitle
#' @param fix_scale Whether to fix scale to 0-1
#' @param y_col_name Name of the Y-axis column
#' @param y_axis_label Label for Y-axis
#' @return ggplot object
plot_variant_heatmap <- function(data, val_col, title_str, sub_str, 
                                  fix_scale = TRUE,
                                  y_col_name = y_col, 
                                  y_axis_label = y_label) {
  
  # Build facet labels with GT count info
  bin_labels_df <- data %>%
    group_by(LFC_Bin) %>%
    summarise(
      Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1),
      Mean_Harm = round(mean(Ref_Harm_Count, na.rm = TRUE), 1),
      .groups = "drop"
    )
  
  data <- data %>%
    dplyr::select(-any_of(c("Mean_GT", "Mean_Harm"))) %>%
    left_join(bin_labels_df, by = "LFC_Bin")
  
  is_harm_plot <- grepl("_Harm", val_col)
  
  data <- data %>%
    mutate(LFC_Display_Label = if (is_harm_plot) {
      sprintf("%s\nMean GT: %.1f", LFC_Bin, Mean_Harm)
    } else {
      sprintf("%s\nMean GT: %.1f", LFC_Bin, Mean_GT)
    })
  
  # Order facets by LFC bin
  unique_bins <- rev(names(LFC_INTERVALS))
  ordered_labels <- unique(data$LFC_Display_Label[order(match(data$LFC_Bin, unique_bins))])
  data$LFC_Display_Label <- factor(data$LFC_Display_Label, levels = ordered_labels)
  
  # Base plot
  p <- ggplot(data, aes(x = Donors, y = !!sym(y_col_name))) +
    geom_tile(aes(fill = !!sym(val_col)), color = "white") +
    facet_wrap(~LFC_Display_Label) +
    labs(title = title_str, subtitle = sub_str, x = "Donors", y = y_axis_label) +
    theme_publication(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
  # Color scale
  if (fix_scale) {
    p <- p + scale_fill_viridis_c(option = "magma", limits = c(0, 1), 
                                   na.value = "grey50", name = "Value")
  } else {
    p <- p + scale_fill_viridis_c(option = "magma", na.value = "grey50", name = "Value")
  }
  
  # Red border for high values (>= 0.8)
  high_val_subset <- data %>% filter(!is.na(!!sym(val_col)) & !!sym(val_col) >= 0.8)
  
  if (nrow(high_val_subset) > 0) {
    p <- p + geom_tile(data = high_val_subset, aes(x = Donors, y = !!sym(y_col_name)),
                       fill = NA, color = "red", linewidth = 0.6, 
                       width = 0.85, height = 0.85)
  }
  
  # Value labels
  p <- p + geom_shadowtext(aes(label = ifelse(is.na(!!sym(val_col)), "", 
                                                sprintf("%.2f", !!sym(val_col)))),
                            color = "white", bg.colour = "black", bg.r = 0.08, size = 2.5)
  
  return(p)
}

# ==============================================================================
# 4. LOAD DATA FOR ALL VARIANTS
# ==============================================================================

log_section("Loading Data for All Astro Variants")

variant_data <- list()

for (variant in ASTRO_VARIANTS) {
  dataset_name <- variant$dataset
  log_info("Loading data for %s (%s)...", dataset_name, variant$label)
  
  # Get paths using existing utility but with custom dataset
  paths <- get_paths("Astro", dataset = dataset_name)
  checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
  
  if (!file.exists(checkpoint_file)) {
    log_warning("Checkpoint not found for %s: %s", dataset_name, checkpoint_file)
    log_warning("Run Steps 1-3 for %s first", dataset_name)
    next
  }
  
  checkpoint_data <- qs::qread(checkpoint_file)
  
  variant_data[[dataset_name]] <- list(
    label = variant$label,
    color = variant$color,
    description = variant$description,
    summary = checkpoint_data$df_summary,
    xist = checkpoint_data$df_xist_summary,
    reads = checkpoint_data$df_reads_summary,
    paths = paths
  )
  
  log_info("Loaded %s: %d summary rows", dataset_name, 
           nrow(checkpoint_data$df_summary))
}

if (length(variant_data) == 0) {
  log_error("No variant data found. Ensure Steps 1-3 have been run for Astro variants.")
  stop("No data to visualize")
}

log_info("Successfully loaded %d/%d variants", length(variant_data), length(ASTRO_VARIANTS))

# ==============================================================================
# 5. CREATE OUTPUT DIRECTORY
# ==============================================================================

# Output directory for variant comparison plots
variant_plots_dir <- file.path(OUTPUT_BASE, "Astro", "Variants", DOWNSAMPLE_MODE, "Plots")
if (!dir.exists(variant_plots_dir)) {
  dir.create(variant_plots_dir, recursive = TRUE)
  log_info("Created variant plots directory: %s", variant_plots_dir)
}

# ==============================================================================
# 6. GENERATE INDIVIDUAL HEATMAPS FOR EACH VARIANT
# ==============================================================================

log_section("Generating Individual Variant Heatmaps")

for (dataset_name in names(variant_data)) {
  vdata <- variant_data[[dataset_name]]
  log_info("Processing heatmaps for %s...", dataset_name)
  
  df_summary <- vdata$summary
  df_xist <- vdata$xist
  
  # Convert to factors
  df_summary$Donors <- clean_donors_factor(df_summary$Donors)
  df_summary[[y_col]] <- clean_y_factor(df_summary[[y_col]])
  df_summary <- df_summary %>%
    complete(Donors, !!sym(y_col), LFC_Bin, Gene_Category)
  
  df_xist$Donors <- clean_donors_factor(df_xist$Donors)
  df_xist[[y_col]] <- clean_y_factor(df_xist[[y_col]])
  df_xist <- df_xist %>%
    complete(Donors, !!sym(y_col), LFC_Bin)
  
  # Generate heatmaps for each chromosome category
  for (cat in CHR_CATEGORIES) {
    df_cat <- df_summary %>% filter(Gene_Category == cat)
    if (nrow(df_cat) == 0) next
    
    # Mean Recall
    p_recall <- plot_variant_heatmap(
      data = df_cat,
      val_col = "Mean_Recall",
      title_str = sprintf("%s - %s", vdata$label, cat),
      sub_str = sprintf("Mean Recall (%s)", vdata$description),
      fix_scale = TRUE
    )
    
    fname <- sprintf("%s_Recall_%s.png", dataset_name, cat)
    ggsave(file.path(variant_plots_dir, fname), p_recall, width = 12, height = 8)
    log_info("  Saved: %s", fname)
    
    # Mean Precision
    p_precision <- plot_variant_heatmap(
      data = df_cat,
      val_col = "Mean_Precision",
      title_str = sprintf("%s - %s", vdata$label, cat),
      sub_str = sprintf("Mean Precision (%s)", vdata$description),
      fix_scale = TRUE
    )
    
    fname <- sprintf("%s_Precision_%s.png", dataset_name, cat)
    ggsave(file.path(variant_plots_dir, fname), p_precision, width = 12, height = 8)
    log_info("  Saved: %s", fname)
  }
  
  # XIST sanity check
  p_xist <- ggplot(df_xist %>% filter(!is.na(Donors), !is.na(!!sym(y_col))),
                   aes(x = Donors, y = !!sym(y_col), fill = Mean_Detection)) +
    geom_tile(color = "white") +
    facet_wrap(~LFC_Bin) +
    labs(title = sprintf("XIST Detection: %s", vdata$label),
         subtitle = vdata$description,
         x = "Donors", y = y_label, fill = "Detection\nRate") +
    theme_publication(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50") +
    geom_shadowtext(aes(label = ifelse(is.na(Mean_Detection), "", sprintf("%.2f", Mean_Detection))),
                    color = "white", bg.colour = "black", bg.r = 0.08, size = 2.5)
  
  ggsave(file.path(variant_plots_dir, sprintf("%s_XIST.png", dataset_name)), 
         p_xist, width = 12, height = 8)
  log_info("  Saved: %s_XIST.png", dataset_name)
}

# ==============================================================================
# 7. GENERATE SIDE-BY-SIDE COMPARISON PLOTS
# ==============================================================================

log_section("Generating Side-by-Side Comparison Plots")

# Combine data from all variants for comparison
combined_data <- list()

for (dataset_name in names(variant_data)) {
  vdata <- variant_data[[dataset_name]]
  df <- vdata$summary
  df$Dataset <- dataset_name
  df$Dataset_Label <- vdata$label
  df$Dataset_Color <- vdata$color
  combined_data[[dataset_name]] <- df
}

df_combined <- do.call(rbind, combined_data)
rownames(df_combined) <- NULL

log_info("Combined data: %d rows across %d variants", nrow(df_combined), length(combined_data))

# Filter to Total category for main comparison
df_total <- df_combined %>% filter(Gene_Category == "Total")

# Convert to factors for plotting
df_total$Donors <- clean_donors_factor(df_total$Donors)
df_total[[y_col]] <- clean_y_factor(df_total[[y_col]])
df_total$Dataset_Label <- factor(df_total$Dataset_Label, 
                                  levels = sapply(ASTRO_VARIANTS, function(x) x$label))

# Create comparison plots for each metric
metrics_to_compare <- list(
  list(col = "Mean_Recall", name = "Recall", desc = "Mean Recall"),
  list(col = "Mean_Precision", name = "Precision", desc = "Mean Precision"),
  list(col = "Mean_Recall_Harm", name = "Recall_Harm", desc = "Mean Recall (Harmonized)")
)

for (metric in metrics_to_compare) {
  if (!(metric$col %in% colnames(df_total))) {
    log_warning("Metric %s not found, skipping comparison plot", metric$col)
    next
  }
  
  log_info("Creating %s comparison plot...", metric$name)
  
  # Create faceted comparison
  p_compare <- ggplot(df_total, aes(x = Donors, y = !!sym(y_col), fill = !!sym(metric$col))) +
    geom_tile(color = "white") +
    facet_grid(LFC_Bin ~ Dataset_Label, switch = "y") +
    labs(title = sprintf("Astro Variants Comparison — %s", metric$desc),
         subtitle = sprintf("Mode: %s | Gene Category: Total", DOWNSAMPLE_MODE),
         x = "Number of Donors", y = y_label, fill = metric$desc) +
    theme_publication(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 11, color = "white"),
      strip.text.y = element_text(face = "bold", size = 10),
      strip.background.x = element_rect(fill = "#333333", color = NA),
      strip.background.y = element_rect(fill = "#666666", color = NA),
      panel.spacing = unit(0.5, "lines")
    ) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50")
  
  # Add red borders for high values
  high_vals <- df_total %>% filter(!is.na(!!sym(metric$col)) & !!sym(metric$col) >= 0.8)
  if (nrow(high_vals) > 0) {
    p_compare <- p_compare + 
      geom_tile(data = high_vals, fill = NA, color = "red", linewidth = 0.5)
  }
  
  fname <- sprintf("Comparison_%s_Total.png", metric$name)
  ggsave(file.path(variant_plots_dir, fname), p_compare, width = 16, height = 10)
  log_info("  Saved: %s", fname)
}

# ==============================================================================
# 8. GENERATE DIFFERENCE HEATMAPS (Variant vs Standard)
# ==============================================================================

log_section("Generating Difference Heatmaps (Variant vs Standard)")

# Only compute differences if we have the standard dataset
if ("Roussos" %in% names(variant_data) && length(variant_data) > 1) {
  standard_data <- variant_data[["Roussos"]]$summary %>%
    filter(Gene_Category == "Total")
  
  for (dataset_name in names(variant_data)) {
    if (dataset_name == "Roussos") next  # Skip standard
    
    vdata <- variant_data[[dataset_name]]
    variant_summary <- vdata$summary %>% filter(Gene_Category == "Total")
    
    # Join with standard data
    df_diff <- variant_summary %>%
      inner_join(standard_data, 
                 by = c("Donors", y_col, "LFC_Bin", "Gene_Category"),
                 suffix = c("_var", "_std")) %>%
      mutate(
        Recall_Diff = Mean_Recall_var - Mean_Recall_std,
        Precision_Diff = Mean_Precision_var - Mean_Precision_std
      )
    
    # Plot recall difference
    df_diff$Donors <- clean_donors_factor(df_diff$Donors)
    df_diff[[y_col]] <- clean_y_factor(df_diff[[y_col]])
    
    p_diff <- ggplot(df_diff, aes(x = Donors, y = !!sym(y_col), fill = Recall_Diff)) +
      geom_tile(color = "white") +
      facet_wrap(~LFC_Bin) +
      labs(title = sprintf("%s vs Standard — Recall Difference", vdata$label),
           subtitle = sprintf("Positive = Higher recall in %s", vdata$label),
           x = "Donors", y = y_label, fill = "Recall\nDifference") +
      theme_publication(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           limits = c(-1, 1), na.value = "grey50")
    
    fname <- sprintf("Difference_%s_vs_Standard_Recall.png", dataset_name)
    ggsave(file.path(variant_plots_dir, fname), p_diff, width = 12, height = 8)
    log_info("  Saved: %s", fname)
  }
} else {
  log_warning("Standard dataset (Roussos) not available, skipping difference heatmaps")
}

# ==============================================================================
# 9. GENERATE SUMMARY TABLE
# ==============================================================================

log_section("Generating Summary Table")

summary_rows <- list()

for (dataset_name in names(variant_data)) {
  vdata <- variant_data[[dataset_name]]
  df <- vdata$summary %>% filter(Gene_Category == "Total")
  
  # Calculate summary statistics
  for (lfc_bin in unique(df$LFC_Bin)) {
    df_bin <- df %>% filter(LFC_Bin == lfc_bin)
    
    # Find grid point with highest recall
    best_idx <- which.max(df_bin$Mean_Recall)
    best_row <- df_bin[best_idx, ]
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Dataset = dataset_name,
      Dataset_Label = vdata$label,
      LFC_Bin = lfc_bin,
      Max_Recall = max(df_bin$Mean_Recall, na.rm = TRUE),
      Max_Recall_Donors = best_row$Donors,
      Max_Recall_Y = best_row[[y_col]],
      Mean_Recall_All = mean(df_bin$Mean_Recall, na.rm = TRUE),
      Median_Recall_All = median(df_bin$Mean_Recall, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

df_summary_table <- do.call(rbind, summary_rows)
rownames(df_summary_table) <- NULL

# Save summary table
summary_csv <- file.path(variant_plots_dir, "Variant_Summary_Table.csv")
write.csv(df_summary_table, summary_csv, row.names = FALSE)
log_info("Saved summary table: %s", summary_csv)

# Print summary
log_info("\n=== Astro Variants Summary ===")
for (dataset_name in names(variant_data)) {
  vdata <- variant_data[[dataset_name]]
  ds_summary <- df_summary_table %>% filter(Dataset == dataset_name)
  log_info("%s (%s):", dataset_name, vdata$label)
  for (i in 1:nrow(ds_summary)) {
    log_info("  %s: Max Recall = %.3f at (Donors=%s, %s=%s)",
             ds_summary$LFC_Bin[i],
             ds_summary$Max_Recall[i],
             ds_summary$Max_Recall_Donors[i],
             y_col,
             ds_summary$Max_Recall_Y[i])
  }
}

# ==============================================================================
# 10. SAVE COMBINED CHECKPOINT
# ==============================================================================

checkpoint_path <- file.path(variant_plots_dir, "Astro_Variants_Combined_Data.qs")
qs::qsave(list(
  variant_data = variant_data,
  combined_summary = df_combined,
  summary_table = df_summary_table,
  astro_variants = ASTRO_VARIANTS,
  mode = DOWNSAMPLE_MODE,
  timestamp = Sys.time()
), checkpoint_path)
log_info("Saved combined checkpoint: %s", checkpoint_path)

log_section("STEP 6 COMPLETE")
log_info("All Astro variant visualizations generated in: %s", variant_plots_dir)
log_info("Total plots generated: %d individual + %d comparison + %d difference",
         length(variant_data) * (length(CHR_CATEGORIES) * 2 + 1),  # Individual
         length(metrics_to_compare),  # Comparison
         max(0, length(variant_data) - 1))  # Difference
