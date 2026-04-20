# ==============================================================================
# SCRIPT: Cross-Cell-Type Aggregated Power Heatmaps
# PURPOSE: Pool raw results across all cell types and runs to produce
#          aggregated recall/precision heatmaps. Each (gene, cell_type,
#          run, iteration) is an independent observation — the same gene
#          appearing as a DEG in 3 cell types contributes 3 observations.
# INPUT: Aggregated_Summary_Data.qs checkpoints from Step 3 (per cell type)
# OUTPUT: Aggregated heatmaps in Global/{Dataset}/Plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  
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

log_section("STEP 5: Cross-Cell-Type Aggregated Heatmaps")

# Mode-aware configuration
y_col <- get_y_col_name()
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Reads per Donor" else "Cells per Donor"
clean_y_factor <- if (DOWNSAMPLE_MODE == "reads") clean_reads_factor else clean_cells_factor
log_info("Mode: %s | Y-column: %s", DOWNSAMPLE_MODE, y_col)

# ==============================================================================
# 2. AGGREGATION FUNCTIONS
# ==============================================================================

#' Load raw results from all cell type checkpoints
#'
#' @param cell_types Character vector of cell type names
#' @param output_base Base output directory
#' @param dataset Dataset name (e.g., "Roussos")
#' @return Data frame with all raw results tagged by cell_type, or NULL if none found
load_all_raw_results <- function(cell_types, output_base = OUTPUT_BASE, 
                                   dataset = SHARD_DATASET,
                                   mode = DOWNSAMPLE_MODE) {
  
  all_raw <- list()
  
  for (ct in cell_types) {
    checkpoint_file <- file.path(output_base, ct, dataset, mode, "Plots", 
                                  "Aggregated_Summary_Data.qs")
    
    if (!file.exists(checkpoint_file)) {
      log_warning("Checkpoint not found for %s: %s", ct, checkpoint_file)
      next
    }
    
    checkpoint_data <- qs::qread(checkpoint_file)
    
    if (is.null(checkpoint_data$df_raw) || nrow(checkpoint_data$df_raw) == 0) {
      log_warning("No raw data in checkpoint for %s", ct)
      next
    }
    
    df <- checkpoint_data$df_raw
    df$Cell_Type <- ct
    all_raw[[ct]] <- df
    
    log_info("Loaded %d raw rows from %s", nrow(df), ct)
  }
  
  if (length(all_raw) == 0) {
    log_error("No raw data found for any cell type")
    return(NULL)
  }
  
  combined <- do.call(rbind, all_raw)
  rownames(combined) <- NULL
  
  log_info("Combined %d total rows across %d cell types", 
           nrow(combined), length(all_raw))
  
  return(combined)
}

#' Compute pooled aggregation across cell types and runs
#'
#' For each (Donors, Cells, LFC_Bin, Gene_Category), pools all observations
#' across cell types, runs, and iterations. Each row in df_raw represents one
#' observation with its own Recall, Precision, and GT_Count.
#'
#' Pooled Recall = sum(Recall_i * GT_Count_i) / sum(GT_Count_i)
#'   This equals total_recovered / total_GT_observations across all cell types.
#'
#' Pooled Precision = sum(Precision_i * Called_i) / sum(Called_i)
#'   where Called_i is inferred from the relationship:
#'   Precision = TP / (TP + FP), and TP = Recall * GT_Count
#'   So Called = TP / Precision (when Precision > 0)
#'
#' @param df_combined Combined raw data frame from load_all_raw_results()
#' @return Data frame with pooled metrics per grid point
compute_pooled_aggregation <- function(df_combined, y_col_name = y_col) {
  
  pooled <- df_combined %>%
    mutate(
      # Recover TP from recall: TP = Recall * GT_Count
      TP_Recall = Recall * GT_Count,
      # Recover TP from precision: TP_prec = Precision * (TP_prec + FP)
      # We know TP_prec / Precision = TP_prec + FP (when Precision > 0)
      # But we need the denominator (total called significant in bin)
      # Since Precision = TP / Called, and TP = Recall * GT_Count (approximately),
      # Called = TP / Precision (when Precision > 0)
      Called_Sig = ifelse(!is.na(Precision) & Precision > 0, 
                          (Recall * GT_Count) / Precision, 
                          NA_real_),
      TP_Precision = ifelse(!is.na(Precision), Recall * GT_Count, NA_real_),
      # Harmonized versions
      TP_Recall_Harm = if ("Recall_Harm" %in% colnames(df_combined) && 
                           "GT_Harm_Count" %in% colnames(df_combined)) {
        Recall_Harm * GT_Harm_Count
      } else {
        rep(0, n())
      },
      GT_Harm_Count_Safe = ifelse(is.na(GT_Harm_Count), 0, GT_Harm_Count)
    ) %>%
    group_by(Donors, !!sym(y_col_name), LFC_Bin, Gene_Category) %>%
    summarise(
      # Pooled recall: total recovered / total GT observations
      Pooled_Recall = ifelse(sum(GT_Count, na.rm = TRUE) > 0,
                              sum(TP_Recall, na.rm = TRUE) / sum(GT_Count, na.rm = TRUE),
                              0),
      # Pooled precision: total TP / total called
      Pooled_Precision = ifelse(sum(Called_Sig, na.rm = TRUE) > 0,
                                 sum(TP_Precision, na.rm = TRUE) / sum(Called_Sig, na.rm = TRUE),
                                 0),
      # Pooled harmonized recall
      Pooled_Recall_Harm = ifelse(sum(GT_Harm_Count_Safe, na.rm = TRUE) > 0,
                                   sum(TP_Recall_Harm, na.rm = TRUE) / sum(GT_Harm_Count_Safe, na.rm = TRUE),
                                   0),
      # Metadata for headers
      Total_GT_Observations = sum(GT_Count, na.rm = TRUE),
      Mean_GT_Per_Obs = mean(GT_Count, na.rm = TRUE),
      N_Cell_Types = n_distinct(Cell_Type),
      N_Observations = n(),
      Cell_Types_Present = paste(sort(unique(Cell_Type)), collapse = "+"),
      .groups = "drop"
    )
  
  return(pooled)
}

# ==============================================================================
# 3. PLOTTING FUNCTION
# ==============================================================================

#' Create aggregated metric heatmap
#'
#' @param data Data frame from compute_pooled_aggregation()
#' @param val_col Column name for the metric to plot
#' @param title_str Plot title
#' @param sub_str Plot subtitle
#' @return ggplot object
plot_aggregated_heatmap <- function(data, val_col, title_str, sub_str,
                                     y_col_name = y_col, y_axis_label = y_label) {
  
  # Build facet labels with pooled GT count info
  bin_labels_df <- data %>%
    group_by(LFC_Bin) %>%
    summarise(
      Mean_Pooled_GT = round(mean(Total_GT_Observations, na.rm = TRUE), 0),
      CTs = dplyr::first(Cell_Types_Present),
      N_Obs = round(mean(N_Observations, na.rm = TRUE), 0),
      .groups = "drop"
    )
  
  data <- data %>%
    dplyr::select(-any_of(c("Mean_Pooled_GT", "CTs", "N_Obs"))) %>%
    left_join(bin_labels_df, by = "LFC_Bin")
  
  data <- data %>%
    mutate(LFC_Display_Label = sprintf(
      "%s\nPooled GT obs: %s\n[%s]",
      LFC_Bin, 
      scales::comma(Mean_Pooled_GT),
      CTs
    ))
  
  # Order facets by LFC bin
  unique_bins <- rev(names(LFC_INTERVALS))
  ordered_labels <- unique(data$LFC_Display_Label[order(match(data$LFC_Bin, unique_bins))])
  data$LFC_Display_Label <- factor(data$LFC_Display_Label, levels = ordered_labels)
  
  # Base plot (dynamic Y-axis column)
  p <- ggplot(data, aes(x = Donors, y = !!sym(y_col_name))) +
    geom_tile(aes(fill = !!sym(val_col)), color = "white") +
    facet_wrap(~LFC_Display_Label) +
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
  
  # Color scale (0-1, magma)
  p <- p + scale_fill_viridis_c(option = "magma", limits = c(0, 1), 
                                  na.value = "grey50", name = "Value")
  
  # Red border for high values (>= 0.8)
  high_val_subset <- data %>% filter(!is.na(!!sym(val_col)) & !!sym(val_col) >= 0.8)
  
  if (nrow(high_val_subset) > 0) {
    p <- p + geom_tile(data = high_val_subset, aes(x = Donors, y = !!sym(y_col_name)),
                       fill = NA, color = "red", linewidth = 0.6,
                       width = 0.85, height = 0.85)
  }
  
  # Value labels (white text with black outline for readability on all tiles)
  p <- p + geom_shadowtext(aes(label = ifelse(is.na(!!sym(val_col)), "",
                                               sprintf("%.2f", !!sym(val_col)))),
                           color = "white", bg.colour = "black", bg.r = 0.08, size = 2.5)
  
  return(p)
}

# ==============================================================================
# 4. LOAD AND AGGREGATE DATA
# ==============================================================================

cell_types <- get_cell_types()
log_info("Detected %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

if (length(cell_types) < 2) {
  log_warning("Only %d cell type found. Cross-cell-type aggregation requires >= 2. Skipping.", 
              length(cell_types))
} else {
  
  # Load raw data from all cell types
  df_combined <- load_all_raw_results(cell_types)
  
  if (is.null(df_combined)) {
    log_error("Failed to load raw data. Ensure Step 3 has been run for all cell types.")
  } else {
    
    # Filter to Total chromosome category only
    df_total <- df_combined %>% filter(Gene_Category == "Total")
    
    log_info("Filtered to Total category: %d rows", nrow(df_total))
    
    # Compute pooled aggregation
    pooled <- compute_pooled_aggregation(df_total)
    
    log_info("Pooled aggregation: %d grid points", nrow(pooled))
    
    # ==============================================================================
    # 5. CREATE OUTPUT DIRECTORY AND SAVE CHECKPOINT
    # ==============================================================================
    
    global_plots_dir <- file.path(OUTPUT_BASE, "Global", SHARD_DATASET, DOWNSAMPLE_MODE, "Plots")
    if (!dir.exists(global_plots_dir)) {
      dir.create(global_plots_dir, recursive = TRUE)
      log_info("Created output directory: %s", global_plots_dir)
    }
    
    # Save pooled checkpoint
    checkpoint_path <- file.path(global_plots_dir, "Pooled_Cross_CellType_Data.qs")
    qs::qsave(list(
      pooled = pooled,
      df_combined_total = df_total,
      cell_types = cell_types,
      timestamp = Sys.time()
    ), checkpoint_path)
    log_info("Saved pooled checkpoint: %s", checkpoint_path)
    
    # ==============================================================================
    # 6. GENERATE HEATMAPS
    # ==============================================================================
    
    log_section("Generating Aggregated Heatmaps")
    
    # Convert to factors for plotting (mode-aware)
    pooled$Donors <- clean_donors_factor(pooled$Donors)
    pooled[[y_col]] <- clean_y_factor(pooled[[y_col]])
    
    # Complete the grid (fill missing combinations with NA)
    pooled <- pooled %>%
      complete(Donors, !!sym(y_col), LFC_Bin, Gene_Category)
    
    ct_label <- paste(sapply(cell_types, prettify_label), collapse = " + ")
    
    metrics_to_plot <- list(
      list(
        col = "Pooled_Recall",
        desc = "Pooled Recall",
        sub = sprintf("Aggregated across %s — all runs and iterations", ct_label),
        fname = "Pooled_Recall_Total.png"
      ),
      list(
        col = "Pooled_Precision",
        desc = "Pooled Precision",
        sub = sprintf("Aggregated across %s — all runs and iterations", ct_label),
        fname = "Pooled_Precision_Total.png"
      ),
      list(
        col = "Pooled_Recall_Harm",
        desc = "Pooled Recall (Harmonised)",
        sub = sprintf("Aggregated across %s — all runs and iterations", ct_label),
        fname = "Pooled_Recall_Harm_Total.png"
      )
    )
    
    for (m in metrics_to_plot) {
      
      # Check if the column exists and has non-zero values
      if (!(m$col %in% colnames(pooled))) {
        log_warning("Column %s not found in pooled data, skipping", m$col)
        next
      }
      
      fpath <- file.path(global_plots_dir, m$fname)
      
      p <- plot_aggregated_heatmap(
        data = pooled,
        val_col = m$col,
        title_str = m$desc,
        sub_str = m$sub
      )
      
      ggsave(fpath, p, width = 18, height = 12)
      log_info("Saved: %s", m$fname)
    }
    
    log_info("All aggregated heatmaps generated in %s", global_plots_dir)
  }
}

log_section("STEP 5 COMPLETE")
log_info("Cross-cell-type aggregation finished")
