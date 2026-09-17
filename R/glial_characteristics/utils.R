# ==============================================================================
# UTILS: Glial Characteristics Pipeline - Shared Utility Functions
# ==============================================================================
# This file contains shared helper functions used across all pipeline steps.
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(qs)
  library(dplyr)
  library(BiocParallel)
})

# Load shared prettify_labels for consistent cell type formatting
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. DIRECTORY MANAGEMENT
# ==============================================================================

#' Create output directories if they don't exist
#' 
#' @param paths Character vector of paths to create
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
# 2. MANIFEST LOADING
# ==============================================================================

#' Load manifest for a dataset/cell type combination
#' 
#' @param cell_type Cell type (Astro, Micro, Oligo)
#' @param dataset Dataset name (Roussos, DeJager)
#' @param shard_base Base path to shards
#' @param datasets_config List of dataset configurations
#' @return List with manifest data.frame and metadata
load_manifest <- function(cell_type, dataset, shard_base, datasets_config) {
  
  ds_config <- datasets_config[[dataset]]
  shard_dir <- file.path(shard_base, cell_type, ds_config$shard_dir_suffix)
  manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
  
  if (!file.exists(manifest_path)) {
    stop(sprintf("Manifest not found: %s", manifest_path))
  }
  
  manifest_obj <- qs::qread(manifest_path)
  
  return(list(
    manifest = manifest_obj$manifest,
    config = manifest_obj$config,
    shard_dir = shard_dir,
    ds_config = ds_config
  ))
}

#' Get donor IDs from manifest
#' 
#' @param manifest_obj Object returned by load_manifest()
#' @return Character vector of donor IDs
get_donor_ids <- function(manifest_obj) {
  unique(as.character(manifest_obj$manifest$donor_id))
}

#' Get cell counts from manifest
#' 
#' @param manifest_obj Object returned by load_manifest()
#' @return Named numeric vector: cell count per donor
get_cell_counts_from_manifest <- function(manifest_obj) {
  df <- manifest_obj$manifest
  counts <- df$cell_count
  names(counts) <- df$donor_id
  return(counts)
}

# ==============================================================================
# 3. SHARD LOADING
# ==============================================================================

#' Load a single donor shard
#' 
#' @param donor_id Donor ID
#' @param shard_dir Path to shard directory
#' @return SingleCellExperiment or NULL if not found
load_shard <- function(donor_id, shard_dir) {
  safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
  shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
  
  if (!file.exists(shard_path)) {
    return(NULL)
  }
  
  qs::qread(shard_path)
}

#' Extract per-cell metrics from a shard
#' 
#' @param sce SingleCellExperiment object
#' @return Data frame with per-cell metrics
extract_cell_metrics <- function(sce) {
  counts_mat <- assay(sce, "counts")
  
  data.frame(
    umis = as.numeric(Matrix::colSums(counts_mat)),
    genes = as.integer(Matrix::colSums(counts_mat > 0)),
    stringsAsFactors = FALSE
  )
}

#' Extract per-cell metrics from all donors in parallel
#' 
#' @param donor_ids Character vector of donor IDs
#' @param shard_dir Path to shard directory
#' @param n_cores Number of cores for parallel processing
#' @return List with: cell_metrics (data.frame), donor_summaries (data.frame)
extract_all_cell_metrics <- function(donor_ids, shard_dir, n_cores = 12) {
  
  BPPARAM <- MulticoreParam(workers = n_cores, progressbar = FALSE)
  
  results <- bplapply(donor_ids, function(donor_id) {
    sce <- load_shard(donor_id, shard_dir)
    if (is.null(sce)) return(NULL)
    
    metrics <- extract_cell_metrics(sce)
    metrics$donor_id <- donor_id
    metrics$n_cells <- ncol(sce)
    
    return(metrics)
  }, BPPARAM = BPPARAM)
  
  # Filter out NULLs
  results <- results[!sapply(results, is.null)]
  
  if (length(results) == 0) {
    stop("No valid shards found")
  }
  
  # Combine all cell-level metrics
  cell_metrics <- bind_rows(results)
  
  # Compute per-donor summaries
  donor_summaries <- cell_metrics %>%
    group_by(donor_id) %>%
    summarise(
      n_cells = first(n_cells),
      median_umis = median(umis),
      mean_umis = mean(umis),
      median_genes = median(genes),
      mean_genes = mean(genes),
      .groups = "drop"
    )
  
  return(list(
    cell_metrics = cell_metrics,
    donor_summaries = donor_summaries
  ))
}

# ==============================================================================
# 4. SUMMARY STATISTICS
# ==============================================================================

#' Compute summary statistics for a numeric vector
#' 
#' @param x Numeric vector
#' @param name Name of the metric
#' @return Data frame with summary statistics
compute_summary_stats <- function(x, name = "value") {
  data.frame(
    metric = name,
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Compute inter-cell-type ratios
#' 
#' @param stats_df Data frame with summary stats (must have cell_type, median columns)
#' @param reference_ct Reference cell type for denominator
#' @param metric_name Name of the metric being compared
#' @return Data frame with ratios
compute_celltype_ratios <- function(stats_df, reference_ct = "Micro", metric_name = "value") {
  
  # Get reference median
  ref_median <- stats_df$median[stats_df$cell_type == reference_ct]
  
  if (length(ref_median) != 1) {
    stop(sprintf("Reference cell type '%s' not found or duplicated", reference_ct))
  }
  
  # Compute ratios for each cell type
  ratios <- stats_df %>%
    mutate(
      ratio_vs_ref = median / ref_median,
      reference_ct = reference_ct,
      metric = metric_name
    ) %>%
    select(cell_type, median, ratio_vs_ref, reference_ct, metric)
  
  return(ratios)
}

# ==============================================================================
# 5. BOOTSTRAP FUNCTIONS
# ==============================================================================

#' Bootstrap the median of a numeric vector
#' 
#' @param x Numeric vector
#' @param n_boot Number of bootstrap resamples
#' @param seed Random seed
#' @return Numeric vector of bootstrap medians
bootstrap_median <- function(x, n_boot = 10000, seed = 42) {
  set.seed(seed)
  n <- length(x)
  
  boot_medians <- replicate(n_boot, {
    boot_idx <- sample.int(n, n, replace = TRUE)
    median(x[boot_idx])
  })
  
  return(boot_medians)
}

#' Bootstrap the ratio of two medians (from paired data)
#' 
#' @param x1 Numeric vector for numerator
#' @param x2 Numeric vector for denominator (same length)
#' @param n_boot Number of bootstrap resamples
#' @param seed Random seed
#' @return Numeric vector of bootstrap ratios
bootstrap_median_ratio_paired <- function(x1, x2, n_boot = 10000, seed = 42) {
  
  if (length(x1) != length(x2)) {
    stop("x1 and x2 must have the same length for paired bootstrap")
  }
  
  set.seed(seed)
  n <- length(x1)
  
  boot_ratios <- replicate(n_boot, {
    boot_idx <- sample.int(n, n, replace = TRUE)
    median(x1[boot_idx]) / median(x2[boot_idx])
  })
  
  return(boot_ratios)
}

#' Bootstrap the ratio of two medians (from unpaired data)
#' 
#' @param x1 Numeric vector for numerator
#' @param x2 Numeric vector for denominator
#' @param n_boot Number of bootstrap resamples
#' @param seed Random seed
#' @return Numeric vector of bootstrap ratios
bootstrap_median_ratio <- function(x1, x2, n_boot = 10000, seed = 42) {
  
  set.seed(seed)
  n1 <- length(x1)
  n2 <- length(x2)
  
  boot_ratios <- replicate(n_boot, {
    boot_idx1 <- sample.int(n1, n1, replace = TRUE)
    boot_idx2 <- sample.int(n2, n2, replace = TRUE)
    median(x1[boot_idx1]) / median(x2[boot_idx2])
  })
  
  return(boot_ratios)
}

#' Compute bootstrap confidence interval
#' 
#' @param boot_values Numeric vector of bootstrap statistics
#' @param ci_level Confidence level (default 0.95)
#' @param method CI method: "percentile" or "bca" 
#' @return Named numeric vector: lower, median, upper
bootstrap_ci <- function(boot_values, ci_level = 0.95, method = "percentile") {
  
  alpha <- 1 - ci_level
  
  if (method == "percentile") {
    ci <- quantile(boot_values, probs = c(alpha/2, 0.5, 1 - alpha/2), na.rm = TRUE)
    names(ci) <- c("lower", "median", "upper")
  } else {
    stop("Only 'percentile' method currently implemented")
  }
  
  return(ci)
}

# ==============================================================================
# 6. SUBSAMPLE UTILITY
# ==============================================================================

#' Subsample a data frame to a maximum number of rows
#' 
#' @param df Data frame to subsample
#' @param max_rows Maximum number of rows
#' @param seed Random seed
#' @return Subsampled data frame
subsample_df <- function(df, max_rows, seed = 42) {
  if (nrow(df) <= max_rows) {
    return(df)
  }
  
  set.seed(seed)
  idx <- sample.int(nrow(df), max_rows)
  df[idx, , drop = FALSE]
}

# ==============================================================================
# 7. SEX HARMONIZATION
# ==============================================================================

#' Harmonise sex encoding to standard "male"/"female"
#' 
#' @param sex_values Character vector of sex values
#' @param encoding Named vector: c(male = "...", female = "...")
#' @return Character vector with standardised values
harmonise_sex <- function(sex_values, encoding) {
  sex_std <- rep(NA_character_, length(sex_values))
  sex_std[sex_values == encoding["male"]] <- "male"
  sex_std[sex_values == encoding["female"]] <- "female"
  return(sex_std)
}

# ==============================================================================
# VERSION
# ==============================================================================

UTILS_VERSION <- "1.0.0"
message(sprintf("[UTILS] Glial Characteristics Utils v%s loaded", UTILS_VERSION))
