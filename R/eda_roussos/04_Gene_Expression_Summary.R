# ==============================================================================
# SCRIPT: Gene Expression Summary - Gene-Level Statistics
# PURPOSE: Calculate mean expression, detection rate, and identify zero-inflated genes
# INPUT: Donor shards from Datalink/Shards/{CellType}/Roussos_Full/
# OUTPUT: Tables in Datalink/Exploratory_Analysis/{CellType}/Roussos_Full/Tables/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(Matrix)
  library(dplyr)
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

message("\n### STEP 04: GENE EXPRESSION SUMMARY ###\n")

# Setup parallel backend
BPPARAM <- MulticoreParam(workers = N_CORES, progressbar = FALSE)

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

process_cell_type_genes <- function(cell_type) {
  
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
  
  # Load all shards and combine counts
  message("  Loading all shards...")
  shard_list <- bplapply(manifest_df$donor_id, function(donor_id) {
    safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
    
    if (!file.exists(shard_path)) return(NULL)
    qs::qread(shard_path)
  }, BPPARAM = BPPARAM)
  
  shard_list <- shard_list[!sapply(shard_list, is.null)]
  
  if (length(shard_list) == 0) {
    warning(sprintf("No valid shards for %s", cell_type))
    return(NULL)
  }
  
  # Combine into single SCE
  message("  Combining shards...")
  sce_combined <- do.call(cbind, shard_list)
  message(sprintf("  Combined: %d genes × %s cells", 
                  nrow(sce_combined), format(ncol(sce_combined), big.mark = ",")))
  
  # Get counts matrix
  counts_mat <- assay(sce_combined, "counts")
  
  # Calculate statistics
  message("  Calculating gene statistics...")
  
  mean_expr <- calculate_mean_expression(counts_mat)
  detection_rate <- calculate_detection_rate(counts_mat)
  zero_inflated <- identify_zero_inflated(counts_mat, ZERO_INFLATED_THRESHOLD)
  
  # Create summary table
  gene_summary <- data.frame(
    gene = rownames(counts_mat),
    mean_expression = mean_expr,
    detection_rate = detection_rate,
    zero_inflated = zero_inflated,
    stringsAsFactors = FALSE
  )
  
  # Sort by mean expression
  gene_summary <- gene_summary[order(gene_summary$mean_expression, decreasing = TRUE), ]
  
  message(sprintf("  Zero-inflated genes: %d/%d (%.1f%%)",
                  sum(zero_inflated), length(zero_inflated),
                  100 * sum(zero_inflated) / length(zero_inflated)))
  
  # Save full summary
  summary_file <- file.path(output_paths$tables, "Gene_Expression_Summary.csv")
  write.csv(gene_summary, summary_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(summary_file)))
  
  # Get top expressed genes
  top_genes <- get_top_expressed_genes(counts_mat, TOP_GENES_N)
  
  # Save top genes
  top_file <- file.path(output_paths$tables, "Top_Expressed_Genes.csv")
  write.csv(top_genes, top_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(top_file)))
  
  message(sprintf("  ✓ Completed gene expression summary for %s", cell_type))
  
  return(list(
    n_genes = nrow(gene_summary),
    n_zero_inflated = sum(zero_inflated),
    top_gene = top_genes$gene[1],
    top_gene_expr = top_genes$mean_expression[1]
  ))
}

# ==============================================================================
# 3. RUN FOR ALL CELL TYPES
# ==============================================================================

message(sprintf("Processing %d cell types: %s",
                length(CELL_TYPES), paste(CELL_TYPES, collapse = ", ")))

results <- lapply(CELL_TYPES, function(ct) {
  tryCatch({
    process_cell_type_genes(ct)
  }, error = function(e) {
    warning(sprintf("Critical failure for %s: %s", ct, e$message))
    return(NULL)
  })
})

# Summary
n_success <- sum(!sapply(results, is.null))
message(sprintf("\n### STEP 04 COMPLETE: %d/%d cell types processed successfully ###",
                n_success, length(CELL_TYPES)))
