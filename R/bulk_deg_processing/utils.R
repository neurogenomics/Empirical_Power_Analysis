# ==============================================================================
# BULK DEG PROCESSING PIPELINE - UTILITY FUNCTIONS
# ==============================================================================
# Shared functions used across pipeline steps
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(biomaRt)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
})

# ==============================================================================
# GENE ID CLEANING
# ==============================================================================

#' Clean Ensembl Gene IDs
#' 
#' Removes version suffixes from Ensembl gene IDs
#' Example: ENSG00000123456.12 -> ENSG00000123456
#' 
#' @param gene_ids Character vector of Ensembl gene IDs
#' @return Character vector of cleaned gene IDs
clean_gene_ids <- function(gene_ids) {
  if (length(gene_ids) == 0) return(character(0))
  cleaned <- sub("\\..*", "", gene_ids)
  return(cleaned)
}

# ==============================================================================
# TISSUE COUNTING
# ==============================================================================

#' Count Tissues Per Gene
#' 
#' Calculate how many unique tissues each gene appears in
#' 
#' @param dt data.table with columns 'gene' and 'tissue'
#' @return data.table with columns 'gene' and 'tissue_count'
count_tissues_per_gene <- function(dt) {
  if (nrow(dt) == 0) {
    return(data.table(gene = character(0), tissue_count = integer(0)))
  }
  
  if (!all(c("gene", "tissue") %in% names(dt))) {
    stop("Input must have 'gene' and 'tissue' columns")
  }
  
  gene_counts <- dt[, .(tissue_count = uniqueN(tissue)), by = gene]
  return(gene_counts)
}

# ==============================================================================
# CONSISTENCY FILTERING
# ==============================================================================

#' Filter Genes by Tissue Consistency
#' 
#' Keep only genes that appear in at least a threshold proportion of tissues
#' 
#' @param dt data.table with columns 'gene', 'tissue', and significance columns
#' @param p_threshold P-value threshold (lfsr)
#' @param consistency_threshold Proportion of tissues required (0-1)
#' @return data.table with filtered genes and averaged statistics
filter_consistent_genes <- function(dt, p_threshold = 0.05, consistency_threshold = 0.90) {
  if (nrow(dt) == 0) {
    return(data.table(
      gene = character(0),
      log2FoldChange = numeric(0),
      lfcSE = numeric(0),
      padj = numeric(0)
    ))
  }
  
  # Validate required columns
  required_cols <- c("gene", "tissue", "effsize", "effsize_se", "lfsr")
  if (!all(required_cols %in% names(dt))) {
    stop(sprintf("Input must have columns: %s", paste(required_cols, collapse = ", ")))
  }
  
  # Clean gene IDs
  dt[, clean_gene := clean_gene_ids(gene)]
  
  # Calculate tissue threshold
  total_tissues <- uniqueN(dt$tissue)
  min_tissues <- ceiling(total_tissues * consistency_threshold)
  
  message(sprintf("Total tissues: %d", total_tissues))
  message(sprintf("Required tissues for %.0f%% consistency: %d", 
                  consistency_threshold * 100, min_tissues))
  
  # Filter for significance and calculate averages
  consistent_genes <- dt[lfsr < p_threshold, .(
    n_count = .N,
    log2FoldChange = mean(effsize, na.rm = TRUE),
    lfcSE = mean(effsize_se, na.rm = TRUE),
    padj = mean(lfsr, na.rm = TRUE)
  ), by = .(gene = clean_gene)]
  
  # Apply consistency threshold
  final_list <- consistent_genes[n_count >= min_tissues & abs(log2FoldChange) > 0]
  
  # Remove helper column
  final_list[, n_count := NULL]
  
  message(sprintf("Genes passing filter: %d", nrow(final_list)))
  
  return(final_list)
}

# ==============================================================================
# CHROMOSOME ANNOTATION
# ==============================================================================

#' Annotate Genes with Chromosome Location
#' 
#' Map Ensembl gene IDs to chromosome names using org.Hs.eg.db
#' 
#' @param gene_ids Character vector of Ensembl gene IDs (without version suffixes)
#' @return data.table with columns 'gene', 'chromosome', 'is_sex_chr'
annotate_chromosomes <- function(gene_ids) {
  if (length(gene_ids) == 0) {
    return(data.table(
      gene = character(0),
      chromosome = character(0),
      is_sex_chr = logical(0)
    ))
  }
  
  message(sprintf("Annotating %d genes with chromosome locations...", length(gene_ids)))
  
  # Map to chromosomes
  chr_map <- mapIds(org.Hs.eg.db,
                    keys = gene_ids,
                    column = "CHR",
                    keytype = "ENSEMBL",
                    multiVals = "first")
  
  # Create result table
  result <- data.table(
    gene = gene_ids,
    chromosome = as.character(chr_map[gene_ids]),
    is_sex_chr = as.character(chr_map[gene_ids]) %in% c("X", "Y")
  )
  
  # Count unmapped genes
  n_unmapped <- sum(is.na(result$chromosome))
  if (n_unmapped > 0) {
    message(sprintf("Warning: %d genes could not be mapped to chromosomes", n_unmapped))
  }
  
  return(result)
}

# ==============================================================================
# INTERVAL BINNING
# ==============================================================================

#' Assign Genes to LFC Intervals
#' 
#' Bin genes by absolute log fold-change into specified intervals
#' 
#' @param lfc_values Numeric vector of log fold-change values
#' @param intervals Named list of intervals (e.g., list("0.5-1.0" = c(0.5, 1.0)))
#' @return Character vector of interval labels (same length as lfc_values)
bin_by_lfc_interval <- function(lfc_values, intervals) {
  if (length(lfc_values) == 0) return(character(0))
  
  abs_lfc <- abs(lfc_values)
  interval_labels <- rep(NA_character_, length(abs_lfc))
  
  for (label in names(intervals)) {
    range_vals <- intervals[[label]]
    min_val <- range_vals[1]
    max_val <- range_vals[2]
    
    # Assign to interval: min <= LFC < max
    in_interval <- abs_lfc >= min_val & abs_lfc < max_val
    interval_labels[in_interval] <- label
  }
  
  return(interval_labels)
}

# ==============================================================================
# BIOMART QUERYING
# ==============================================================================

#' Query BioMart for Gene Annotations
#' 
#' Fetch gene symbols and descriptions from Ensembl BioMart
#' 
#' @param gene_ids Character vector of Ensembl gene IDs (without version suffixes)
#' @param dataset BioMart dataset name (default: "hsapiens_gene_ensembl")
#' @param attributes Character vector of attributes to fetch
#' @return data.frame with gene annotations
query_biomart_safe <- function(gene_ids, 
                                dataset = "hsapiens_gene_ensembl",
                                attributes = c("ensembl_gene_id", "external_gene_name", "description")) {
  if (length(gene_ids) == 0) {
    return(data.frame(
      ensembl_gene_id = character(0),
      external_gene_name = character(0),
      description = character(0)
    ))
  }
  
  message(sprintf("Querying BioMart for %d genes...", length(gene_ids)))
  
  # Connect to Ensembl
  ensembl <- tryCatch({
    useEnsembl(biomart = "genes", dataset = dataset)
  }, error = function(e) {
    stop(sprintf("Failed to connect to BioMart: %s", e$message))
  })
  
  # Fetch annotations
  gene_info <- tryCatch({
    getBM(attributes = attributes,
          filters = "ensembl_gene_id",
          values = gene_ids,
          mart = ensembl)
  }, error = function(e) {
    stop(sprintf("BioMart query failed: %s", e$message))
  })
  
  message(sprintf("Retrieved annotations for %d genes", nrow(gene_info)))
  
  return(gene_info)
}

# ==============================================================================
# CHECKPOINT VALIDATION
# ==============================================================================

#' Validate Checkpoint Data
#' 
#' Check if cached data meets current requirements
#' 
#' @param data Data object to validate
#' @param type Type of checkpoint ("bulk_all", "filtered", "annotated")
#' @return Logical indicating if data is valid
validate_checkpoint <- function(data, type = "bulk_all") {
  if (is.null(data)) {
    message("  ✗ Data is NULL")
    return(FALSE)
  }
  
  if (type == "bulk_all") {
    # Validate bulk_all checkpoint
    required_cols <- c("gene", "tissue", "effsize", "effsize_se", "lfsr")
    if (!all(required_cols %in% names(data))) {
      message(sprintf("  ✗ Missing required columns: %s", 
                      paste(setdiff(required_cols, names(data)), collapse = ", ")))
      return(FALSE)
    }
    
    if (nrow(data) == 0) {
      message("  ✗ Empty data")
      return(FALSE)
    }
    
  } else if (type == "filtered") {
    # Validate filtered checkpoint
    required_cols <- c("gene", "log2FoldChange", "lfcSE", "padj")
    if (!all(required_cols %in% names(data))) {
      message(sprintf("  ✗ Missing required columns: %s",
                      paste(setdiff(required_cols, names(data)), collapse = ", ")))
      return(FALSE)
    }
    
    if (nrow(data) == 0) {
      message("  ✗ Empty data")
      return(FALSE)
    }
    
    # Check for NA in critical columns
    if (any(is.na(data$gene))) {
      message("  ✗ NA values in gene column")
      return(FALSE)
    }
    
  } else if (type == "annotated") {
    # Validate annotated checkpoint
    required_cols <- c("gene", "external_gene_name", "description", "log2FoldChange")
    if (!all(required_cols %in% names(data))) {
      message(sprintf("  ✗ Missing required columns: %s",
                      paste(setdiff(required_cols, names(data)), collapse = ", ")))
      return(FALSE)
    }
    
    if (nrow(data) == 0) {
      message("  ✗ Empty data")
      return(FALSE)
    }
  }
  
  message("  ✓ Checkpoint validation passed")
  return(TRUE)
}

# ==============================================================================
# SANITY CHECKS
# ==============================================================================

#' Sanity Check Helper
#' 
#' Run a set of sanity checks on data and report results
#' 
#' @param data Data object to check
#' @param description Description of what's being checked
#' @param checks Named list of check functions
#' @return Logical indicating if all checks passed
sanity_check <- function(data, description, checks) {
  cat(sprintf("\n[SANITY CHECK] %s\n", description))
  issues <- c()
  
  for (check_name in names(checks)) {
    check_fn <- checks[[check_name]]
    result <- check_fn(data)
    if (!isTRUE(result)) {
      issues <- c(issues, sprintf("⚠️  %s: %s", check_name, result))
    }
  }
  
  if (length(issues) == 0) {
    cat("  ✅ PASS - All sanity checks passed\n")
    return(TRUE)
  } else {
    cat("  🔴 SANITY CHECK FAILURES:\n")
    for (issue in issues) cat(sprintf("     %s\n", issue))
    cat("\n  ⚠️  These indicate fundamental data issues.\n")
    cat("     Investigate before proceeding.\n")
    return(FALSE)
  }
}
