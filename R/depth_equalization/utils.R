# ==============================================================================
# UTILS: Depth Equalization Pipeline
# ==============================================================================
# Utility functions for multinomial downsampling of SingleCellExperiment shards.
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

# ==============================================================================
# 1. CORE DOWNSAMPLING FUNCTION
# ==============================================================================

#' Multinomially downsample a count matrix by a fixed fraction
#'
#' For each cell (column), the total UMI count is reduced to
#' round(original_total * fraction) by drawing from a multinomial
#' distribution with probabilities proportional to the original gene counts.
#' This preserves the gene-level proportions within each cell while
#' uniformly reducing sequencing depth.
#'
#' @param counts_matrix A sparse or dense matrix (genes x cells) of integer counts
#' @param fraction Numeric in (0, 1). The fraction of reads to retain per cell.
#' @param seed Integer seed for reproducibility. Set before the entire matrix operation.
#' @return A dgCMatrix with downsampled counts (same dimensions as input)
#'
#' @details
#' - Cells with 0 total counts are left unchanged (all zeros).
#' - The multinomial draw ensures the downsampled cell sums to exactly
#'   round(original_total * fraction).
#' - Gene proportions are preserved in expectation; stochastic variation
#'   increases for cells with very low counts.
downsample_counts_matrix <- function(counts_matrix, fraction, seed) {
  
  if (fraction <= 0 || fraction >= 1) {
    stop(sprintf("fraction must be in (0, 1), got %.4f", fraction))
  }
  
  set.seed(seed)
  
  n_genes <- nrow(counts_matrix)
  n_cells <- ncol(counts_matrix)
  
  # Convert to regular matrix for column-wise operations if sparse
  # (rmultinom needs numeric vectors, not sparse columns)
  is_sparse <- is(counts_matrix, "sparseMatrix")
  
  # Pre-compute column sums
  col_totals <- colSums(counts_matrix)
  
  # Build new matrix column by column
  # Use a list of sparse vectors for memory efficiency
  new_cols <- vector("list", n_cells)
  
  for (j in seq_len(n_cells)) {
    original_total <- col_totals[j]
    
    if (original_total == 0) {
      # Empty cell — leave as zeros
      new_cols[[j]] <- rep(0L, n_genes)
      next
    }
    
    new_total <- round(original_total * fraction)
    
    if (new_total == 0) {
      # Fraction too small for this cell — result is all zeros
      new_cols[[j]] <- rep(0L, n_genes)
      next
    }
    
    # Extract column as dense vector
    if (is_sparse) {
      cell_counts <- as.numeric(counts_matrix[, j])
    } else {
      cell_counts <- counts_matrix[, j]
    }
    
    # Compute probabilities (proportional to counts)
    probs <- cell_counts / original_total
    
    # Multinomial draw
    new_cols[[j]] <- as.integer(rmultinom(1, size = new_total, prob = probs))
  }
  
  # Assemble into sparse matrix
  new_matrix <- do.call(cbind, lapply(new_cols, function(x) {
    Matrix::Matrix(x, ncol = 1, sparse = TRUE)
  }))
  
  # Preserve dimension names
  rownames(new_matrix) <- rownames(counts_matrix)
  colnames(new_matrix) <- colnames(counts_matrix)
  
  return(new_matrix)
}

# ==============================================================================
# 2. SHARD-LEVEL DOWNSAMPLING
# ==============================================================================

#' Downsample a single donor shard
#'
#' Loads a SingleCellExperiment shard, multinomially downsamples the counts
#' matrix, and returns the modified SCE with all metadata preserved.
#'
#' @param shard_path Path to the input .qs shard file
#' @param fraction Downsampling fraction (0, 1)
#' @param seed Seed for reproducibility
#' @return A SingleCellExperiment with downsampled counts
downsample_shard <- function(shard_path, fraction, seed) {
  
  sce <- qs::qread(shard_path)
  
  # Downsample the counts assay
  original_counts <- assay(sce, "counts")
  new_counts <- downsample_counts_matrix(original_counts, fraction, seed)
  
  # Replace counts in the SCE
  assay(sce, "counts") <- new_counts
  
  return(sce)
}

# ==============================================================================
# 3. MANIFEST GENERATION
# ==============================================================================

#' Generate a manifest for the depth-equalized shard directory
#'
#' Mirrors the format of the original sharding pipeline's 00_MASTER_MANIFEST.qs:
#' - $manifest: data.frame with donor_id, cell_count, file_path
#' - $feasibility: named list mapping cell-step thresholds to eligible donor IDs
#' - $config: list with metadata about the equalization
#'
#' @param shard_dir Path to the depth-equalized shard directory
#' @param source_config List with source pipeline metadata
#' @param cell_steps Vector of cell count thresholds for feasibility map
#' @return List with manifest, feasibility, and config components
generate_manifest <- function(shard_dir, source_config, 
                              cell_steps = c(5, 10, 20, 50, 75, 100, 250, 500, 
                                             750, 1000, 1500, 2000, 3000, NA)) {
  
  shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  
  if (length(shard_files) == 0) {
    stop("No shard files found in: ", shard_dir)
  }
  
  log_info("Building manifest from %d shard files...", length(shard_files))
  
  manifest_df <- do.call(rbind, lapply(shard_files, function(f) {
    tmp <- qs::qread(f)
    data.frame(
      donor_id   = as.character(unique(colData(tmp)[["donor_id"]])),
      cell_count = ncol(tmp),
      file_path  = f,
      stringsAsFactors = FALSE
    )
  }))
  
  # Feasibility map (same format as sharding pipeline)
  feasibility_map <- lapply(cell_steps, function(n) {
    if (is.na(n)) return(manifest_df$donor_id)
    manifest_df$donor_id[manifest_df$cell_count >= n]
  })
  names(feasibility_map) <- paste0("step_", ifelse(is.na(cell_steps), "Max", cell_steps))
  
  manifest_final <- list(
    manifest    = manifest_df,
    feasibility = feasibility_map,
    config      = list(
      donor_col           = "donor_id",
      n_donors            = nrow(manifest_df),
      total_cells         = sum(manifest_df$cell_count),
      cell_filter         = source_config$cell_type,
      min_age             = source_config$min_age,
      min_cells           = source_config$min_cells,
      depth_equalized     = TRUE,
      source_dataset      = source_config$source_dataset,
      target_mean_rpc     = source_config$target_mean_rpc,
      downsample_fraction = source_config$downsample_fraction,
      equalization_date   = Sys.time()
    )
  )
  
  return(manifest_final)
}

# ==============================================================================
# 4. VALIDATION HELPERS
# ==============================================================================

#' Compute summary statistics for a shard's read depths
#'
#' @param shard_path Path to a .qs shard file
#' @return Named list with n_cells, mean_rpc, median_rpc, min_rpc, max_rpc, total_reads
summarise_shard_depth <- function(shard_path) {
  sce <- qs::qread(shard_path)
  reads_per_cell <- Matrix::colSums(assay(sce, "counts"))
  
  list(
    n_cells     = length(reads_per_cell),
    mean_rpc    = mean(reads_per_cell),
    median_rpc  = median(reads_per_cell),
    min_rpc     = min(reads_per_cell),
    max_rpc     = max(reads_per_cell),
    total_reads = sum(reads_per_cell)
  )
}

# ==============================================================================
# 5. CELL COUNT MATCHING (Astro_Shallow_Fewer Extension)
# ==============================================================================

#' Compute cell count matching statistics between source and target
#'
#' Calculates the scaling factor needed to match the source cell count
#' distribution to the target using mean ratio scaling.
#'
#' @param source_summary Data frame with donor_id, n_cells for source
#' @param target_summary Data frame with donor_id, n_cells for target
#' @return List with scaling_factor, target_mean, source_mean, and donor counts
compute_scaling_factor <- function(source_summary, target_summary) {
  target_mean <- mean(target_summary$n_cells, na.rm = TRUE)
  source_mean <- mean(source_summary$n_cells, na.rm = TRUE)
  scaling_factor <- target_mean / source_mean
  
  list(
    scaling_factor = scaling_factor,
    target_mean = target_mean,
    source_mean = source_mean,
    target_donors = nrow(target_summary),
    source_donors = nrow(source_summary),
    target_median = median(target_summary$n_cells, na.rm = TRUE),
    source_median = median(source_summary$n_cells, na.rm = TRUE)
  )
}

#' Load manifest and extract cell count summary
#'
#' @param manifest_path Path to 00_MASTER_MANIFEST.qs file
#' @return Data frame with donor_id and n_cells
load_cell_count_summary <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(sprintf("Manifest not found: %s", manifest_path))
  }
  
  manifest <- qs::qread(manifest_path)
  
  data.frame(
    donor_id = manifest$manifest$donor_id,
    n_cells = manifest$manifest$cell_count,
    stringsAsFactors = FALSE
  )
}

#' Apply cell count matching to shards via mean scaling
#'
#' For each donor shard, randomly subsample cells to match the target
#' distribution via mean scaling. Donors with insufficient cells are dropped.
#'
#' @param source_shard_dir Path to source shards (e.g., Roussos_DepthEq)
#' @param target_shard_dir Path to output shards (e.g., Roussos_ShallowFewer)
#' @param matching_stats From compute_scaling_factor()
#' @param min_cells_threshold Minimum cells to retain donor
#' @param seed_base Base seed for reproducibility (each donor uses seed_base + index)
#' @return Data frame with per-donor results
match_cell_counts <- function(source_shard_dir, target_shard_dir,
                              matching_stats, min_cells_threshold = 10,
                              seed_base = 9998) {
  
  # Get source shard files
  source_files <- list.files(source_shard_dir, pattern = "^donor_.*\\.qs$", 
                              full.names = TRUE)
  source_files <- sort(source_files)  # Deterministic ordering
  
  n_total <- length(source_files)
  log_info("Processing %d shards for cell count matching...", n_total)
  log_info("Scaling factor: %.4f (target mean: %.1f, source mean: %.1f)",
           matching_stats$scaling_factor, matching_stats$target_mean, 
           matching_stats$source_mean)
  
  # Create output directory
  if (!dir.exists(target_shard_dir)) {
    dir.create(target_shard_dir, recursive = TRUE)
  }
  
  # Process each shard
  results <- vector("list", n_total)
  
  for (i in seq_len(n_total)) {
    source_path <- source_files[i]
    shard_basename <- basename(source_path)
    donor_id <- sub("^donor_", "", sub("\\.qs$", "", shard_basename))
    
    # Load source shard
    sce <- tryCatch({
      qs::qread(source_path)
    }, error = function(e) {
      log_error("Failed to load shard %s: %s", shard_basename, e$message)
      return(NULL)
    })
    
    if (is.null(sce)) {
      results[[i]] <- list(
        donor_id = donor_id,
        success = FALSE,
        error = "Failed to load shard",
        n_cells_before = NA,
        n_cells_after = NA,
        dropped = TRUE
      )
      next
    }
    
    n_cells_before <- ncol(sce)
    
    # Calculate target cell count
    n_cells_target <- round(n_cells_before * matching_stats$scaling_factor)
    
    # Check if donor should be dropped
    if (n_cells_target < min_cells_threshold) {
      results[[i]] <- list(
        donor_id = donor_id,
        success = FALSE,
        error = sprintf("Insufficient cells after scaling: %d < %d",
                        n_cells_target, min_cells_threshold),
        n_cells_before = n_cells_before,
        n_cells_after = n_cells_target,
        dropped = TRUE
      )
      log_info("  %s: DROPPED (%d cells -> %d cells < threshold %d)",
               donor_id, n_cells_before, n_cells_target, min_cells_threshold)
      next
    }
    
    # Randomly subsample cells
    set.seed(seed_base + i)
    chosen_idx <- sort(sample(seq_len(n_cells_before), n_cells_target))
    
    # Subset the SCE
    sce_matched <- sce[, chosen_idx]
    
    # Save to target directory
    target_path <- file.path(target_shard_dir, shard_basename)
    tryCatch({
      qs::qsave(sce_matched, target_path)
      
      results[[i]] <- list(
        donor_id = donor_id,
        success = TRUE,
        error = NA,
        n_cells_before = n_cells_before,
        n_cells_after = n_cells_target,
        dropped = FALSE
      )
      
      if (i %% 100 == 0 || i == n_total) {
        log_info("  Processed %d/%d shards...", i, n_total)
      }
      
    }, error = function(e) {
      results[[i]] <<- list(
        donor_id = donor_id,
        success = FALSE,
        error = sprintf("Failed to save: %s", e$message),
        n_cells_before = n_cells_before,
        n_cells_after = n_cells_target,
        dropped = FALSE
      )
    })
  }
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      donor_id = r$donor_id,
      success = r$success,
      error = ifelse(is.na(r$error), "", r$error),
      n_cells_before = r$n_cells_before,
      n_cells_after = r$n_cells_after,
      dropped = r$dropped,
      stringsAsFactors = FALSE
    )
  }))
  
  return(results_df)
}

#' Generate manifest for cell-count-matched shards
#'
#' Similar to generate_manifest() but for the matched dataset
#'
#' @param shard_dir Path to the matched shard directory
#' @param source_config List with source pipeline metadata
#' @param matching_stats From compute_scaling_factor()
#' @param match_results Data frame from match_cell_counts()
#' @param cell_steps Vector of cell count thresholds for feasibility map
#' @return List with manifest, feasibility, and config components
generate_matched_manifest <- function(shard_dir, source_config, matching_stats,
                                      match_results,
                                      cell_steps = c(5, 10, 20, 50, 75, 100, 250, 500,
                                                     750, 1000, 1500, 2000, 3000, NA)) {
  
  # Filter to successful donors only
  successful_results <- match_results[match_results$success, ]
  
  if (nrow(successful_results) == 0) {
    stop("No successful donors for manifest generation")
  }
  
  # Build manifest data frame
  manifest_df <- data.frame(
    donor_id = successful_results$donor_id,
    cell_count = successful_results$n_cells_after,
    file_path = file.path(shard_dir, paste0("donor_", successful_results$donor_id, ".qs")),
    stringsAsFactors = FALSE
  )
  
  # Feasibility map
  feasibility_map <- lapply(cell_steps, function(n) {
    if (is.na(n)) return(manifest_df$donor_id)
    manifest_df$donor_id[manifest_df$cell_count >= n]
  })
  names(feasibility_map) <- paste0("step_", ifelse(is.na(cell_steps), "Max", cell_steps))
  
  # Config with matching info
  manifest_final <- list(
    manifest = manifest_df,
    feasibility = feasibility_map,
    config = list(
      donor_col = "donor_id",
      n_donors = nrow(manifest_df),
      total_cells = sum(manifest_df$cell_count),
      cell_filter = source_config$cell_type,
      min_age = source_config$min_age,
      min_cells = source_config$min_cells,
      depth_equalized = TRUE,
      cell_count_matched = TRUE,
      source_dataset = source_config$source_dataset,
      target_dataset = FEWER_EXTENSION$target_dataset,
      target_mean_rpc = source_config$target_mean_rpc,
      downsample_fraction = source_config$downsample_fraction,
      cell_count_scaling_factor = matching_stats$scaling_factor,
      target_mean_cells = matching_stats$target_mean,
      source_mean_cells = matching_stats$source_mean,
      donors_retained = sum(match_results$success),
      donors_dropped = sum(!match_results$success),
      equalization_date = Sys.time()
    )
  )
  
  return(manifest_final)
}
