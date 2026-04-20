# ==============================================================================
# SCRIPT: QC Metrics - Cell-Level Quality Control
# PURPOSE: Calculate and visualize QC metrics (MT%, library size, genes detected)
# INPUT: Donor shards from Datalink/Shards/{CellType}/Roussos_Full/
# OUTPUT: QC plots and tables in Datalink/Exploratory_Analysis/{CellType}/Roussos_Full/QC/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(ggplot2)
  library(Matrix)
  library(scales)
  library(BiocParallel)
})

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

message("\n### STEP 05: QC METRICS ###\n")

# Setup parallel backend
BPPARAM <- MulticoreParam(workers = N_CORES, progressbar = FALSE)

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

process_cell_type_qc <- function(cell_type) {
  
  message(sprintf("\n=== Processing: %s ===", cell_type))
  
  # Setup paths
  shard_dir <- file.path(SHARD_BASE, cell_type, SHARD_DATASET)
  output_paths <- get_output_paths(cell_type)
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
  
  # Calculate QC metrics per shard in parallel
  message("  Calculating QC metrics...")
  qc_list <- bplapply(manifest_df$donor_id, function(donor_id) {
    safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
    
    if (!file.exists(shard_path)) return(NULL)
    
    sce <- qs::qread(shard_path)
    
    list(
      mt_percent = calculate_mt_percent(sce, MT_PATTERN),
      genes_detected = calculate_genes_detected(sce),
      library_size = calculate_library_size(sce)
    )
  }, BPPARAM = BPPARAM)
  
  qc_list <- qc_list[!sapply(qc_list, is.null)]
  
  if (length(qc_list) == 0) {
    warning(sprintf("No valid QC metrics for %s", cell_type))
    return(NULL)
  }
  
  # Combine metrics
  mt_percent <- unlist(lapply(qc_list, `[[`, "mt_percent"))
  genes_detected <- unlist(lapply(qc_list, `[[`, "genes_detected"))
  library_size <- unlist(lapply(qc_list, `[[`, "library_size"))
  
  message(sprintf("  QC metrics: %s cells", comma(length(mt_percent))))
  
  # Calculate summary statistics
  qc_stats <- data.frame(
    cell_type = cell_type,
    n_cells = length(mt_percent),
    mean_mt_percent = mean(mt_percent, na.rm = TRUE),
    median_mt_percent = median(mt_percent, na.rm = TRUE),
    mean_genes_detected = mean(genes_detected),
    median_genes_detected = median(genes_detected),
    mean_library_size = mean(library_size),
    median_library_size = median(library_size)
  )
  
  message(sprintf("  MT%%: %.2f%% | Genes: %.0f | Library: %s",
                  qc_stats$median_mt_percent,
                  qc_stats$median_genes_detected,
                  comma(qc_stats$median_library_size)))
  
  # Plot 1: Mitochondrial Percent
  message("  Generating QC plots...")
  
  p1 <- ggplot(data.frame(mt = mt_percent), aes(x = mt)) +
    geom_histogram(bins = HISTOGRAM_BINS, fill = "#F28E2B",
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = qc_stats$median_mt_percent, linetype = "dashed",
               color = "red", linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0(cell_type, " - Mitochondrial Percentage"),
      subtitle = sprintf("Median: %.2f%% | Mean: %.2f%%",
                         qc_stats$median_mt_percent,
                         qc_stats$mean_mt_percent),
      x = "Mitochondrial %",
      y = "Number of Cells"
    ) +
    scale_y_continuous(labels = comma)
  
  ggsave(file.path(output_paths$qc, "QC_Mitochondrial_Percent.png"),
         p1, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  
  # Plot 2: Genes Detected
  p2 <- ggplot(data.frame(genes = genes_detected), aes(x = genes)) +
    geom_histogram(bins = HISTOGRAM_BINS, fill = "#59A14F",
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = qc_stats$median_genes_detected, linetype = "dashed",
               color = "red", linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0(cell_type, " - Genes Detected Per Cell"),
      subtitle = sprintf("Median: %s | Mean: %s",
                         comma(qc_stats$median_genes_detected),
                         comma(qc_stats$mean_genes_detected)),
      x = "Genes Detected",
      y = "Number of Cells"
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  ggsave(file.path(output_paths$qc, "QC_Genes_Detected.png"),
         p2, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  
  # Plot 3: Library Size
  p3 <- ggplot(data.frame(lib = library_size), aes(x = lib)) +
    geom_histogram(bins = HISTOGRAM_BINS, fill = "#B07AA1",
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = qc_stats$median_library_size, linetype = "dashed",
               color = "red", linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0(cell_type, " - Library Size Per Cell"),
      subtitle = sprintf("Median: %s | Mean: %s",
                         comma(qc_stats$median_library_size),
                         comma(qc_stats$mean_library_size)),
      x = "Library Size (Total UMIs)",
      y = "Number of Cells"
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  ggsave(file.path(output_paths$qc, "QC_Library_Size.png"),
         p3, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  
  message("  Saved 3 QC plots")
  
  # Save summary stats
  stats_file <- file.path(output_paths$tables, "QC_Summary_Stats.csv")
  write.csv(qc_stats, stats_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(stats_file)))
  
  message(sprintf("  ✓ Completed QC metrics for %s", cell_type))
  return(qc_stats)
}

# ==============================================================================
# 3. RUN FOR ALL CELL TYPES
# ==============================================================================

message(sprintf("Processing %d cell types: %s",
                length(CELL_TYPES), paste(CELL_TYPES, collapse = ", ")))

all_stats <- lapply(CELL_TYPES, function(ct) {
  tryCatch({
    process_cell_type_qc(ct)
  }, error = function(e) {
    warning(sprintf("Critical failure for %s: %s", ct, e$message))
    return(NULL)
  })
})

# Combine stats
all_stats <- do.call(rbind, all_stats[!sapply(all_stats, is.null)])

if (nrow(all_stats) > 0) {
  combined_file <- file.path(OUTPUT_BASE, "QC_Summary_All_CellTypes.csv")
  write.csv(all_stats, combined_file, row.names = FALSE)
  message(sprintf("\nSaved combined summary: %s", combined_file))
}

message(sprintf("\n### STEP 05 COMPLETE: %d/%d cell types processed successfully ###",
                nrow(all_stats), length(CELL_TYPES)))
