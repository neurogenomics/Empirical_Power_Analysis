# ==============================================================================
# SCRIPT: Read Depth Analysis - Counts Per Cell Distribution
# PURPOSE: Generate histograms showing distribution of library size (UMIs) per cell
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Histograms in Datalink/Exploratory_Analysis/{CellType}/{Dataset}/Histograms/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(ggplot2)
  library(Matrix)
  library(scales)
  library(dplyr)
  library(BiocParallel)
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

message("\n### STEP 03: READ DEPTH ANALYSIS ###\n")

# Setup parallel backend
BPPARAM <- MulticoreParam(workers = N_CORES, progressbar = FALSE)

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

process_cell_type_reads <- function(cell_type, current_dataset, current_label) {
  
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
  
  # Extract library sizes in parallel
  message("  Calculating library sizes per cell...")
  reads_list <- bplapply(manifest_df$donor_id, function(donor_id) {
    safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
    
    if (!file.exists(shard_path)) return(NULL)
    
    sce <- qs::qread(shard_path)
    lib_sizes <- Matrix::colSums(assay(sce, "counts"))
    return(as.numeric(lib_sizes))
  }, BPPARAM = BPPARAM)
  
  # Combine all reads
  all_reads <- unlist(reads_list, use.names = FALSE)
  
  if (length(all_reads) == 0) {
    warning(sprintf("No valid read counts for %s", cell_type))
    return(NULL)
  }
  
  message(sprintf("  Library sizes: %s cells", comma(length(all_reads))))
  
  # Calculate summary statistics
  stats <- data.frame(
    cell_type = cell_type,
    n_cells = length(all_reads),
    mean_reads = mean(all_reads),
    median_reads = median(all_reads),
    min_reads = min(all_reads),
    max_reads = max(all_reads),
    q25_reads = quantile(all_reads, 0.25),
    q75_reads = quantile(all_reads, 0.75),
    sd_reads = sd(all_reads)
  )
  
  message(sprintf("  Mean: %s, Median: %s, Range: %s-%s", 
                  comma(stats$mean_reads), comma(stats$median_reads),
                  comma(stats$min_reads), comma(stats$max_reads)))
  
  # Create histogram
  message("  Generating histogram...")
  
  plot_df <- data.frame(reads = all_reads)
  pretty_ct <- prettify_label(cell_type)
  
  p <- ggplot(plot_df, aes(x = reads)) +
    geom_histogram(bins = HISTOGRAM_BINS, fill = "#E15759", 
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = stats$median_reads, linetype = "dashed",
               color = "blue", linewidth = 1) +
    theme_publication(base_size = 16) +
    labs(
      title = sprintf("%s — Counts Per Cell", pretty_ct),
      subtitle = sprintf("Median: %s | Mean: %s | n = %s cells",
                         comma(round(stats$median_reads)),
                         comma(round(stats$mean_reads)),
                         comma(stats$n_cells)),
      x = "Counts (UMIs) per Cell",
      y = "Number of Cells"
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  # Save plot
  output_file <- file.path(output_paths$histograms,
                           "Histogram_Reads_Per_Cell.png")
  ggsave(output_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  message(sprintf("  Saved: %s", basename(output_file)))
  
  # Save summary stats
  stats_file <- file.path(output_paths$tables,
                          "Read_Depth_Summary.csv")
  write.csv(stats, stats_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(stats_file)))
  
  message(sprintf("  ✓ Completed read depth analysis for %s", cell_type))
  return(stats)
}

# ==============================================================================
# 3. RUN FOR ALL DATASETS × CELL TYPES
# ==============================================================================

for (current_dataset in SHARD_DATASETS) {
  
  current_label <- COHORT_LABELS[current_dataset]
  message(sprintf("\n--- Dataset: %s (%s) ---", current_dataset, current_label))
  
  # Skip-if-exists: check whether read depth histogram PNGs already exist for ALL cell types
  all_exist <- all(sapply(CELL_TYPES, function(ct) {
    hist_dir <- file.path(OUTPUT_BASE, ct, current_dataset, "Histograms")
    file.exists(file.path(hist_dir, "Histogram_Reads_Per_Cell.png"))
  }))
  
  if (all_exist && !FORCE_RECALC) {
    message(sprintf("[SKIP] Read depth outputs already exist for %s (set FORCE_RECALC=TRUE to regenerate)",
                    current_dataset))
    next
  }
  
  message(sprintf("Processing %d cell types: %s",
                  length(CELL_TYPES), paste(CELL_TYPES, collapse = ", ")))
  
  all_stats <- lapply(CELL_TYPES, function(ct) {
    tryCatch({
      process_cell_type_reads(ct, current_dataset, current_label)
    }, error = function(e) {
      warning(sprintf("Critical failure for %s: %s", ct, e$message))
      return(NULL)
    })
  })
  
  # Combine stats
  all_stats <- do.call(rbind, all_stats[!sapply(all_stats, is.null)])
  
  if (nrow(all_stats) > 0) {
    combined_file <- file.path(OUTPUT_BASE, 
                               paste0("Read_Depth_Summary_All_CellTypes_", current_dataset, ".csv"))
    write.csv(all_stats, combined_file, row.names = FALSE)
    message(sprintf("\nSaved combined summary: %s", basename(combined_file)))
    
    # ========================================================================
    # COMBINED OVERLAY HISTOGRAM (all cell types on one plot)
    # ========================================================================
    message("\n  Generating combined overlay histogram...")
    
    # Collect library sizes from all cell types (sample for performance)
    MAX_CELLS_PER_CT <- 100000  # Sample to keep plot manageable
    all_reads_list <- list()
    
    for (ct in CELL_TYPES) {
      shard_dir <- file.path(SHARD_BASE, ct, current_dataset)
      manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
      if (!file.exists(manifest_path)) next
      
      manifest_df <- qs::qread(manifest_path)$manifest
      
      reads_list <- bplapply(manifest_df$donor_id, function(donor_id) {
        safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
        shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
        if (!file.exists(shard_path)) return(NULL)
        sce <- qs::qread(shard_path)
        lib_sizes <- Matrix::colSums(assay(sce, "counts"))
        return(as.numeric(lib_sizes))
      }, BPPARAM = BPPARAM)
      
      ct_reads <- unlist(reads_list, use.names = FALSE)
      
      # Subsample if too large
      if (length(ct_reads) > MAX_CELLS_PER_CT) {
        set.seed(42)
        ct_reads <- sample(ct_reads, MAX_CELLS_PER_CT)
      }
      
      if (length(ct_reads) > 0) {
        all_reads_list[[ct]] <- data.frame(
          cell_type = ct,
          reads = ct_reads,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(all_reads_list) > 0) {
      combined_df <- bind_rows(all_reads_list)
      combined_df$cell_type_pretty <- prettify_label(combined_df$cell_type)
      combined_df$cell_type_pretty <- factor(combined_df$cell_type_pretty,
        levels = prettify_label(CELL_TYPES))
      
      # Compute global x-axis max for shared scaling
      global_max <- max(combined_df$reads, na.rm = TRUE) * 1.05
      
      # Combined overlay histogram (linear scale)
      p_combined <- ggplot(combined_df, aes(x = reads, fill = cell_type_pretty)) +
        geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.45,
                       position = "identity", color = NA) +
        theme_publication(base_size = 16) +
        scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
        coord_cartesian(xlim = c(0, global_max)) +
        labs(
          title = "Counts Per Cell — All Cell Types",
          x = "Counts (UMIs) per Cell",
          y = "Number of Cells",
          fill = "Cell Type"
        ) +
        scale_x_continuous(labels = comma) +
        scale_y_continuous(labels = comma)
      
      combined_plot_path <- file.path(OUTPUT_BASE,
        paste0("Histogram_Reads_Per_Cell_Combined_", current_dataset, ".png"))
      ggsave(combined_plot_path, p_combined, width = 10, height = 6)
      message(sprintf("  Saved combined histogram: %s", basename(combined_plot_path)))
      
      # Combined overlay histogram (log10 x-axis)
      p_combined_log <- ggplot(combined_df, aes(x = reads, fill = cell_type_pretty)) +
        geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.45,
                       position = "identity", color = NA) +
        theme_publication(base_size = 16) +
        scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
        scale_x_log10(labels = comma) +
        labs(
          title = "Counts Per Cell — All Cell Types (Log Scale)",
          x = "Counts (UMIs) per Cell (log10)",
          y = "Number of Cells",
          fill = "Cell Type"
        ) +
        scale_y_continuous(labels = comma)
      
      combined_log_path <- file.path(OUTPUT_BASE,
        paste0("Histogram_Reads_Per_Cell_Combined_Log10_", current_dataset, ".png"))
      ggsave(combined_log_path, p_combined_log, width = 10, height = 6)
      message(sprintf("  Saved combined log10 histogram: %s", basename(combined_log_path)))
    }
  }
  
  message(sprintf("\n  %s: %d/%d cell types processed successfully", 
                  current_dataset, nrow(all_stats), length(CELL_TYPES)))
}

message("\n### STEP 03 COMPLETE ###")
