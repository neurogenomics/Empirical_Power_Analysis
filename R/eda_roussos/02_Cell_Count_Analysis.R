# ==============================================================================
# SCRIPT: Cell Count Analysis - Cells Per Donor Distribution
# PURPOSE: Generate histograms showing distribution of cells per donor
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Histograms in Datalink/Exploratory_Analysis/{CellType}/{Dataset}/Histograms/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Load shared utilities for prettified labels
# Use hardcoded fallback path for Rscript compatibility
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Load config - skip if SCRIPT_DIR already set (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({ dirname(sys.frame(1)$ofile) }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "eda_roussos")
  })
  source(file.path(script_path, "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### STEP 02: CELL COUNT ANALYSIS ###\n")

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

process_cell_type_counts <- function(cell_type, current_dataset, current_label) {
  
  message(sprintf("\n=== Processing: %s [%s] ===", cell_type, current_dataset))
  
  # Setup paths
  shard_dir <- file.path(SHARD_BASE, cell_type, current_dataset)
  output_paths <- get_output_paths(cell_type, current_dataset)
  create_output_dirs(output_paths)
  
  # Check if shard directory exists
  if (!dir.exists(shard_dir)) {
    warning(sprintf("Shard directory not found: %s", shard_dir))
    return(NULL)
  }
  
  # Load manifest
  manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
  if (!file.exists(manifest_path)) {
    warning(sprintf("Manifest not found: %s", manifest_path))
    return(NULL)
  }
  
  manifest_df <- qs::qread(manifest_path)$manifest
  message(sprintf("  Found %d donor shards", nrow(manifest_df)))
  
  # Extract cell counts per donor
  message("  Counting cells per donor...")
  
  # Use cell_count column from manifest if available
  if ("cell_count" %in% colnames(manifest_df)) {
    cell_counts <- manifest_df$cell_count
  } else {
    # Fallback: load from actual shard files
    cell_counts <- sapply(manifest_df$donor_id, function(donor_id) {
      safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
      shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
      
      if (!file.exists(shard_path)) return(NA_integer_)
      
      sce <- qs::qread(shard_path)
      return(ncol(sce))
    })
    cell_counts <- cell_counts[!is.na(cell_counts)]
  }
  
  if (length(cell_counts) == 0) {
    warning(sprintf("No valid cell counts for %s", cell_type))
    return(NULL)
  }
  
  message(sprintf("  Cell counts: %d donors", length(cell_counts)))
  
  # Calculate summary statistics
  stats <- data.frame(
    cell_type = cell_type,
    n_donors = length(cell_counts),
    mean_cells = mean(cell_counts),
    median_cells = median(cell_counts),
    min_cells = min(cell_counts),
    max_cells = max(cell_counts),
    q25_cells = quantile(cell_counts, 0.25),
    q75_cells = quantile(cell_counts, 0.75),
    sd_cells = sd(cell_counts)
  )
  
  message(sprintf("  Mean: %.1f, Median: %.1f, Range: %d-%d", 
                  stats$mean_cells, stats$median_cells, 
                  stats$min_cells, stats$max_cells))
  
  # Create histogram
  message("  Generating histogram...")
  
  plot_df <- data.frame(cells = cell_counts)
  pretty_ct <- prettify_label(cell_type)
  
  p <- ggplot(plot_df, aes(x = cells)) +
    geom_histogram(bins = HISTOGRAM_BINS, fill = "#4E79A7", 
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = stats$median_cells, linetype = "dashed", 
               color = "red", linewidth = 1) +
    theme_publication(base_size = 16) +
    labs(
      title = sprintf("%s — Cells Per Donor", pretty_ct),
      subtitle = sprintf("Median: %s | Mean: %s | n = %d",
                         comma(round(stats$median_cells)),
                         comma(round(stats$mean_cells)),
                         stats$n_donors),
      x = "Cells per Donor",
      y = "Number of Donors"
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  # Save plot
  output_file <- file.path(output_paths$histograms, 
                           "Histogram_Cells_Per_Donor.png")
  ggsave(output_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  message(sprintf("  Saved: %s", basename(output_file)))
  
  # Save summary stats
  stats_file <- file.path(output_paths$tables, 
                          "Cell_Count_Summary.csv")
  write.csv(stats, stats_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(stats_file)))
  
  message(sprintf("  ✓ Completed cell count analysis for %s", cell_type))
  return(stats)
}

# ==============================================================================
# 3. RUN FOR ALL DATASETS × CELL TYPES
# ==============================================================================

for (current_dataset in SHARD_DATASETS) {
  
  current_label <- COHORT_LABELS[current_dataset]
  message(sprintf("\n--- Dataset: %s (%s) ---", current_dataset, current_label))
  
  # Skip-if-exists: check whether histogram PNGs already exist for ALL cell types
  all_exist <- all(sapply(CELL_TYPES, function(ct) {
    hist_dir <- file.path(OUTPUT_BASE, ct, current_dataset, "Histograms")
    file.exists(file.path(hist_dir, "Histogram_Cells_Per_Donor.png"))
  }))
  
  if (all_exist && !FORCE_RECALC) {
    message(sprintf("[SKIP] Cell count outputs already exist for %s (set FORCE_RECALC=TRUE to regenerate)",
                    current_dataset))
    next
  }
  
  message(sprintf("Processing %d cell types: %s", 
                  length(CELL_TYPES), paste(CELL_TYPES, collapse = ", ")))
  
  all_stats <- lapply(CELL_TYPES, function(ct) {
    tryCatch({
      process_cell_type_counts(ct, current_dataset, current_label)
    }, error = function(e) {
      warning(sprintf("Critical failure for %s: %s", ct, e$message))
      return(NULL)
    })
  })
  
  # Combine stats
  all_stats <- do.call(rbind, all_stats[!sapply(all_stats, is.null)])
  
  if (nrow(all_stats) > 0) {
    combined_file <- file.path(OUTPUT_BASE, 
                               paste0("Cell_Count_Summary_All_CellTypes_", current_dataset, ".csv"))
    write.csv(all_stats, combined_file, row.names = FALSE)
    message(sprintf("\nSaved combined summary: %s", basename(combined_file)))
    
    # ========================================================================
    # COMBINED OVERLAY HISTOGRAM (all cell types on one plot)
    # ========================================================================
    message("\n  Generating combined overlay histogram...")
    
    # Collect cell counts from all cell types
    all_cell_counts_list <- list()
    for (ct in CELL_TYPES) {
      shard_dir <- file.path(SHARD_BASE, ct, current_dataset)
      manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
      if (!file.exists(manifest_path)) next
      
      manifest_df <- qs::qread(manifest_path)$manifest
      if ("cell_count" %in% colnames(manifest_df)) {
        counts <- manifest_df$cell_count
      } else {
        counts <- sapply(manifest_df$donor_id, function(donor_id) {
          safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
          shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
          if (!file.exists(shard_path)) return(NA_integer_)
          sce <- qs::qread(shard_path)
          return(ncol(sce))
        })
        counts <- counts[!is.na(counts)]
      }
      if (length(counts) > 0) {
        all_cell_counts_list[[ct]] <- data.frame(
          cell_type = ct,
          cells = counts,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(all_cell_counts_list) > 0) {
      combined_df <- bind_rows(all_cell_counts_list)
      combined_df$cell_type_pretty <- prettify_label(combined_df$cell_type)
      combined_df$cell_type_pretty <- factor(combined_df$cell_type_pretty,
        levels = prettify_label(CELL_TYPES))
      
      # Compute global x-axis max for shared scaling
      global_max <- max(combined_df$cells, na.rm = TRUE) * 1.05
      
      # Combined overlay histogram (linear scale)
      p_combined <- ggplot(combined_df, aes(x = cells, fill = cell_type_pretty)) +
        geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.45, 
                       position = "identity", color = NA) +
        theme_publication(base_size = 16) +
        scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
        coord_cartesian(xlim = c(0, global_max)) +
        labs(
          title = "Cells Per Donor — All Cell Types",
          x = "Cells per Donor",
          y = "Number of Donors",
          fill = "Cell Type"
        ) +
        scale_x_continuous(labels = comma) +
        scale_y_continuous(labels = comma)
      
      combined_plot_path <- file.path(OUTPUT_BASE, 
        paste0("Histogram_Cells_Per_Donor_Combined_", current_dataset, ".png"))
      ggsave(combined_plot_path, p_combined, width = 10, height = 6)
      message(sprintf("  Saved combined histogram: %s", basename(combined_plot_path)))
      
      # Combined overlay histogram (log10 x-axis)
      p_combined_log <- ggplot(combined_df, aes(x = cells, fill = cell_type_pretty)) +
        geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.45,
                       position = "identity", color = NA) +
        theme_publication(base_size = 16) +
        scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
        scale_x_log10(labels = comma) +
        labs(
          title = "Cells Per Donor — All Cell Types (Log Scale)",
          x = "Cells per Donor (log10)",
          y = "Number of Donors",
          fill = "Cell Type"
        ) +
        scale_y_continuous(labels = comma)
      
      combined_log_path <- file.path(OUTPUT_BASE,
        paste0("Histogram_Cells_Per_Donor_Combined_Log10_", current_dataset, ".png"))
      ggsave(combined_log_path, p_combined_log, width = 10, height = 6)
      message(sprintf("  Saved combined log10 histogram: %s", basename(combined_log_path)))
    }
  }
  
  message(sprintf("\n  %s: %d/%d cell types processed successfully", 
                  current_dataset, nrow(all_stats), length(CELL_TYPES)))
}

message("\n### STEP 02 COMPLETE ###")
