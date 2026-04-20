# ==============================================================================
# UTILS: EDA Roussos Pipeline - Shared Utility Functions
# ==============================================================================
# This file contains shared helper functions used across multiple pipeline steps.
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
})

# ==============================================================================
# 1. AGE PARSING FUNCTIONS
# ==============================================================================

#' Check if a string represents an unknown/missing value
#' 
#' @param x Character vector to check
#' @return Logical vector indicating unknown values
is_unknown_string <- function(x) {
  s <- tolower(trimws(as.character(x)))
  grepl("^(unknown|not available|na|n/a|nan|)$", s)
}

#' Check if a column looks like age stage format
#' 
#' @param x Character or factor vector
#' @return Logical - TRUE if majority of values match age stage pattern
looks_like_age_stage <- function(x) {
  if (!(is.character(x) || is.factor(x))) return(FALSE)
  s <- tolower(trimws(as.character(x)))
  matches <- grepl("(year|month|week|day)[ -]?old", s)
  return(sum(matches, na.rm = TRUE) > (0.5 * length(matches)))
}

#' Parse age stage strings to numeric years
#' 
#' Handles formats like "65 year-old stage", "6 month-old stage", etc.
#' Uses the same logic as the cartesian downsampling pipeline.
#' 
#' @param x Character vector of age strings
#' @return Numeric vector of ages in years
#' @examples
#' parse_age_stage_to_years(c("80 year-old stage", "6 month-old stage"))
parse_age_stage_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
  
  # Extract numeric value
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))
  
  # Detect time units
  is_year    <- grepl("year[ -]?old", s)
  is_month   <- grepl("month[ -]?old", s)
  is_week    <- grepl("week[ -]?old", s)
  is_day     <- grepl("day[ -]?old", s)
  is_unknown <- is_unknown_string(s)
  
  # Convert to years
  y <- rep(NA_real_, length(s))
  y[is_year]    <- n[is_year]
  y[is_month]   <- n[is_month] / 12
  y[is_week]    <- n[is_week] / 52.1429
  y[is_day]     <- n[is_day] / 365.25
  y[is_unknown] <- NA_real_
  
  # Fallback: if it's just a number with no unit, assume years
  just_number <- !is_year & !is_month & !is_week & !is_day & !is_unknown & !is.na(n)
  y[just_number] <- n[just_number]
  
  return(y)
}

# ==============================================================================
# 2. PSEUDOBULK AGGREGATION FUNCTIONS
# ==============================================================================

#' Aggregate pseudobulk from donor shards
#' 
#' @param donor_ids Character vector of donor IDs
#' @param label Character label for this cohort (e.g., "Full_Cohort")
#' @param shard_dir Path to donor shards directory
#' @param id_col Name of donor ID column in metadata
#' @return SingleCellExperiment with pseudobulk counts
get_pb_from_list <- function(donor_ids, label, shard_dir, id_col) {
  
  pb_list <- lapply(donor_ids, function(id) {
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
    
    if (!file.exists(shard_path)) return(NULL)
    
    sce_shard <- qs::qread(shard_path)
    counts_vec <- Matrix::rowSums(assay(sce_shard, "counts"))
    
    # Capture cell count
    n_cells_donor <- ncol(sce_shard)
    
    # Extract metadata row
    meta_row <- as.data.frame(colData(sce_shard)[1, , drop = FALSE])
    meta_row$n_cells <- n_cells_donor
    
    return(list(counts = counts_vec, meta = meta_row))
  })
  
  pb_list <- pb_list[!sapply(pb_list, is.null)]
  if (length(pb_list) == 0) return(NULL)
  
  # Align genes
  all_genes <- names(pb_list[[1]]$counts)
  combined_counts <- do.call(cbind, lapply(pb_list, function(x) {
    as.matrix(x$counts[match(all_genes, names(x$counts))])
  }))
  
  # Align metadata
  combined_meta <- do.call(rbind, lapply(pb_list, `[[`, "meta"))
  combined_meta$Cohort <- label
  
  rownames(combined_counts) <- all_genes
  colnames(combined_counts) <- as.character(combined_meta[[id_col]])
  rownames(combined_meta) <- colnames(combined_counts)
  
  return(SingleCellExperiment(
    assays = list(counts = combined_counts),
    colData = DataFrame(combined_meta)
  ))
}

# ==============================================================================
# 3. GENE STATISTICS FUNCTIONS
# ==============================================================================

#' Calculate mean expression per gene
#' 
#' @param counts_matrix Sparse or dense matrix (genes x cells)
#' @return Numeric vector of mean counts per gene
calculate_mean_expression <- function(counts_matrix) {
  Matrix::rowMeans(counts_matrix)
}

#' Calculate detection rate per gene
#' 
#' @param counts_matrix Sparse or dense matrix (genes x cells)
#' @return Numeric vector of detection rates (0-1)
calculate_detection_rate <- function(counts_matrix) {
  Matrix::rowMeans(counts_matrix > 0)
}

#' Identify zero-inflated genes
#' 
#' @param counts_matrix Sparse or dense matrix (genes x cells)
#' @param threshold Proportion of zeros to flag as zero-inflated
#' @return Logical vector - TRUE for zero-inflated genes
identify_zero_inflated <- function(counts_matrix, threshold = 0.90) {
  zero_prop <- Matrix::rowMeans(counts_matrix == 0)
  return(zero_prop > threshold)
}

#' Get top expressed genes
#' 
#' @param counts_matrix Sparse or dense matrix (genes x cells)
#' @param n Number of top genes to return
#' @return Data frame with gene names and mean expression
get_top_expressed_genes <- function(counts_matrix, n = 50) {
  mean_expr <- calculate_mean_expression(counts_matrix)
  top_idx <- order(mean_expr, decreasing = TRUE)[1:min(n, length(mean_expr))]
  
  data.frame(
    gene = rownames(counts_matrix)[top_idx],
    mean_expression = mean_expr[top_idx],
    rank = 1:length(top_idx),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 4. QC METRIC FUNCTIONS
# ==============================================================================

#' Calculate mitochondrial percentage per cell
#' 
#' @param sce SingleCellExperiment object
#' @param mt_pattern Regex pattern for mitochondrial genes
#' @return Numeric vector of MT% per cell
calculate_mt_percent <- function(sce, mt_pattern = "^MT-") {
  counts_mat <- assay(sce, "counts")
  mt_genes <- grepl(mt_pattern, rownames(sce), ignore.case = TRUE)
  
  if (sum(mt_genes) == 0) {
    warning("No mitochondrial genes found with pattern: ", mt_pattern)
    return(rep(NA_real_, ncol(sce)))
  }
  
  mt_counts <- Matrix::colSums(counts_mat[mt_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts_mat)
  
  return(100 * mt_counts / total_counts)
}

#' Calculate number of genes detected per cell
#' 
#' @param sce SingleCellExperiment object
#' @return Integer vector of genes detected per cell
calculate_genes_detected <- function(sce) {
  counts_mat <- assay(sce, "counts")
  Matrix::colSums(counts_mat > 0)
}

#' Calculate library size per cell
#' 
#' @param sce SingleCellExperiment object
#' @return Numeric vector of total counts per cell
calculate_library_size <- function(sce) {
  counts_mat <- assay(sce, "counts")
  Matrix::colSums(counts_mat)
}

# ==============================================================================
# 5. DIRECTORY MANAGEMENT
# ==============================================================================

#' Create output directories if they don't exist
#' 
#' @param paths List of paths from get_output_paths()
#' @return NULL (creates directories as side effect)
create_output_dirs <- function(paths) {
  for (path in paths) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      message(sprintf("[UTILS] Created directory: %s", path))
    }
  }
  invisible(NULL)
}

# ==============================================================================
# 6. MANIFEST LOADING
# ==============================================================================

#' Load donor IDs from master manifest
#' 
#' @param shard_dir Path to donor shards directory
#' @param id_col Name of donor ID column
#' @return Character vector of donor IDs
load_donor_ids_from_manifest <- function(shard_dir, id_col) {
  manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
  
  if (!file.exists(manifest_path)) {
    stop(sprintf("Manifest not found: %s", manifest_path))
  }
  
  manifest_obj <- qs::qread(manifest_path)
  manifest_df <- manifest_obj$manifest
  
  # Use donor_id column from manifest instead of loading each file
  if (id_col %in% colnames(manifest_df)) {
    all_donors <- unique(as.character(manifest_df[[id_col]]))
  } else {
    # Fallback: try to load from files in current shard_dir
    all_donors <- unique(sapply(manifest_df[[id_col]], function(donor_id) {
      safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
      shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
      
      if (file.exists(shard_path)) {
        return(donor_id)
      } else {
        return(NA_character_)
      }
    }))
    all_donors <- all_donors[!is.na(all_donors)]
  }
  
  return(all_donors)
}

# ==============================================================================
# VERSION
# ==============================================================================

UTILS_VERSION <- "1.0.0"
message(sprintf("[UTILS] EDA Roussos Utils v%s loaded", UTILS_VERSION))
