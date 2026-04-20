# ==============================================================================
# SCRIPT: QC Validation Pipeline - Utility Functions
# PURPOSE: Shared functions for computing and evaluating QC metrics
# INPUT: SingleCellExperiment shard files (.qs)
# OUTPUT: QC metric data frames
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(Matrix)
})

# ==============================================================================
# MITOCHONDRIAL GENE IDENTIFICATION (via Ensembl/biomaRt)
# ==============================================================================

#' Fetch mitochondrial gene identifiers from Ensembl via biomaRt
#'
#' Queries Ensembl for all genes on chromosome MT and returns both
#' Ensembl gene IDs and HGNC symbols. Results are cached to avoid
#' repeated queries.
#'
#' @param cache_file Path to cache file (.qs). If exists, loads from cache.
#' @param dataset Ensembl dataset (default: "hsapiens_gene_ensembl")
#' @return Named list with $ensembl_ids (character vector) and $symbols (character vector)
get_mito_genes <- function(cache_file = NULL, dataset = "hsapiens_gene_ensembl") {
  # Check cache first
  if (!is.null(cache_file) && file.exists(cache_file)) {
    message("[MITO] Loading cached mitochondrial gene list from: ", cache_file)
    cached <- qs::qread(cache_file)
    # Handle legacy cache format (character vector of symbols only)
    if (is.character(cached)) {
      message("[MITO] Legacy cache format (symbols only). Re-querying for Ensembl IDs...")
    } else if (is.list(cached) && all(c("ensembl_ids", "symbols") %in% names(cached))) {
      message(sprintf("[MITO] Loaded %d Ensembl IDs and %d symbols from cache",
        length(cached$ensembl_ids), length(cached$symbols)))
      return(cached)
    }
  }

  message("[MITO] Querying Ensembl for mitochondrial genes...")
  suppressPackageStartupMessages(library(biomaRt))

  ensembl <- tryCatch({
    biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  }, error = function(e) {
    message(sprintf("[ERROR] biomaRt connection failed: %s", e$message))
    message("[MITO] Falling back to ^MT- prefix matching")
    return(NULL)
  })

  if (is.null(ensembl)) {
    return(NULL)  # Caller should handle fallback
  }

  mito_df <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters = "chromosome_name",
    values = "MT",
    mart = ensembl
  )

  mito_result <- list(
    ensembl_ids = unique(mito_df$ensembl_gene_id[mito_df$ensembl_gene_id != ""]),
    symbols = unique(mito_df$hgnc_symbol[mito_df$hgnc_symbol != ""])
  )
  message(sprintf("[MITO] Found %d Ensembl IDs and %d symbols for MT chromosome",
    length(mito_result$ensembl_ids), length(mito_result$symbols)))

  # Cache results
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    qs::qsave(mito_result, cache_file)
    message("[MITO] Cached to: ", cache_file)
  }

  return(mito_result)
}

#' Identify mitochondrial genes in a gene set
#'
#' Tries matching against Ensembl IDs first, then HGNC symbols, then ^MT- prefix.
#' Handles both the new list format from get_mito_genes() and legacy character vectors.
#'
#' @param gene_names Character vector of gene names (rownames of SCE)
#' @param mito_genes Either a named list with $ensembl_ids and $symbols (from get_mito_genes),
#'   or a character vector of gene symbols (legacy format), or NULL for ^MT- fallback
#' @return Logical vector, TRUE for mitochondrial genes
identify_mito_genes <- function(gene_names, mito_genes = NULL) {
  if (is.null(mito_genes)) {
    message("[MITO] No Ensembl gene list provided. Using ^MT- prefix fallback.")
    return(grepl("^MT-", gene_names))
  }

  # Handle new list format (with both Ensembl IDs and symbols)
  if (is.list(mito_genes) && "ensembl_ids" %in% names(mito_genes)) {
    # Try Ensembl IDs first (for ENSG-style gene names)
    is_mito <- gene_names %in% mito_genes$ensembl_ids
    if (sum(is_mito) > 0) {
      return(is_mito)
    }
    # Try symbols (for MT-CO1-style gene names)
    is_mito <- gene_names %in% mito_genes$symbols
    if (sum(is_mito) > 0) {
      return(is_mito)
    }
    # Fall back to ^MT- prefix
    message("[WARN] No Ensembl IDs or symbols matched. Falling back to ^MT- prefix.")
    return(grepl("^MT-", gene_names))
  }

  # Handle legacy character vector format (symbols only)
  if (is.character(mito_genes) && length(mito_genes) > 0) {
    is_mito <- gene_names %in% mito_genes
    if (sum(is_mito) > 0) {
      return(is_mito)
    }
    message("[WARN] No mito gene symbols matched. Falling back to ^MT- prefix.")
    return(grepl("^MT-", gene_names))
  }

  # Fallback
  message("[MITO] Unrecognised mito_genes format. Using ^MT- prefix fallback.")
  return(grepl("^MT-", gene_names))
}

# ==============================================================================
# PER-CELL QC METRIC COMPUTATION
# ==============================================================================

#' Compute per-cell QC metrics from a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment object with a "counts" assay
#' @param mito_genes Character vector of mitochondrial gene symbols
#' @return data.frame with columns: n_umis, n_genes, pct_mito
compute_cell_qc_metrics <- function(sce, mito_genes = NULL) {
  counts_mat <- assay(sce, "counts")
  gene_names <- rownames(sce)

  # Total UMIs per cell
  n_umis <- Matrix::colSums(counts_mat)

  # Genes detected per cell (count > 0)
  n_genes <- Matrix::colSums(counts_mat > 0)

  # Mitochondrial percentage
  is_mito <- identify_mito_genes(gene_names, mito_genes)
  if (any(is_mito)) {
    mito_counts <- Matrix::colSums(counts_mat[is_mito, , drop = FALSE])
    pct_mito <- ifelse(n_umis > 0, (mito_counts / n_umis) * 100, 0)
  } else {
    message("[WARN] No mitochondrial genes found in this SCE. Setting pct_mito = NA.")
    pct_mito <- rep(NA_real_, ncol(sce))
  }

  data.frame(
    n_umis = as.numeric(n_umis),
    n_genes = as.numeric(n_genes),
    pct_mito = as.numeric(pct_mito),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# THRESHOLD CHECKING
# ==============================================================================

#' Check per-cell QC metrics against Fullard et al. 2025 thresholds
#'
#' @param qc_df data.frame from compute_cell_qc_metrics
#' @param min_umis Minimum UMI count (default: 1500)
#' @param max_umis Maximum UMI count (default: 110000)
#' @param min_genes Minimum gene count (default: 1100)
#' @param max_genes Maximum gene count (default: 12500)
#' @param max_mito_pct Maximum mitochondrial percentage (default: 5.0)
#' @return data.frame with original columns plus logical violation flags
check_thresholds <- function(qc_df,
                             min_umis = 1500, max_umis = 110000,
                             min_genes = 1100, max_genes = 12500,
                             max_mito_pct = 5.0) {
  qc_df$fail_low_umis  <- qc_df$n_umis < min_umis
  qc_df$fail_high_umis <- qc_df$n_umis > max_umis
  qc_df$fail_low_genes <- qc_df$n_genes < min_genes
  qc_df$fail_high_genes <- qc_df$n_genes > max_genes

  if (all(is.na(qc_df$pct_mito))) {
    qc_df$fail_high_mito <- NA
  } else {
    qc_df$fail_high_mito <- !is.na(qc_df$pct_mito) & qc_df$pct_mito >= max_mito_pct
  }

  # Any violation
  fail_cols <- c("fail_low_umis", "fail_high_umis", "fail_low_genes",
                 "fail_high_genes", "fail_high_mito")
  qc_df$fail_any <- rowSums(qc_df[, fail_cols], na.rm = TRUE) > 0

  return(qc_df)
}

#' Summarise threshold violations
#'
#' @param qc_df data.frame from check_thresholds (with fail_ columns)
#' @return Named list with counts and percentages for each violation type
summarise_violations <- function(qc_df) {
  n_total <- nrow(qc_df)

  violation_cols <- c("fail_low_umis", "fail_high_umis", "fail_low_genes",
                      "fail_high_genes", "fail_high_mito", "fail_any")

  counts <- sapply(violation_cols, function(col) {
    sum(qc_df[[col]], na.rm = TRUE)
  })

  pcts <- round(counts / n_total * 100, 4)

  list(
    n_total = n_total,
    counts = counts,
    percentages = pcts
  )
}

# ==============================================================================
# PER-DONOR SUMMARY
# ==============================================================================

#' Summarise QC metrics for a single donor
#'
#' @param qc_df data.frame from check_thresholds for one donor
#' @param donor_id Character string identifying the donor
#' @return Single-row data.frame with summary statistics
summarise_donor <- function(qc_df, donor_id) {
  violations <- summarise_violations(qc_df)

  data.frame(
    donor_id = donor_id,
    n_cells = violations$n_total,
    median_umis = median(qc_df$n_umis),
    median_genes = median(qc_df$n_genes),
    median_mito_pct = median(qc_df$pct_mito, na.rm = TRUE),
    mean_umis = mean(qc_df$n_umis),
    mean_genes = mean(qc_df$n_genes),
    mean_mito_pct = mean(qc_df$pct_mito, na.rm = TRUE),
    min_umis = min(qc_df$n_umis),
    max_umis = max(qc_df$n_umis),
    min_genes = min(qc_df$n_genes),
    max_genes = max(qc_df$n_genes),
    max_mito_pct = max(qc_df$pct_mito, na.rm = TRUE),
    n_fail_low_umis = violations$counts[["fail_low_umis"]],
    n_fail_high_umis = violations$counts[["fail_high_umis"]],
    n_fail_low_genes = violations$counts[["fail_low_genes"]],
    n_fail_high_genes = violations$counts[["fail_high_genes"]],
    n_fail_high_mito = violations$counts[["fail_high_mito"]],
    n_fail_any = violations$counts[["fail_any"]],
    pct_fail_any = violations$percentages[["fail_any"]],
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# FAIL-EARLY LOGIC
# ==============================================================================

#' Check whether fail-early conditions are met
#'
#' @param all_qc_df Combined data.frame of QC metrics processed so far
#' @param n_donors_processed Number of donors processed so far
#' @param threshold_10 Violation rate threshold after 10 donors (default: 0.50)
#' @param threshold_50 Violation rate threshold after 50 donors (default: 0.20)
#' @return List with $should_stop (logical) and $reason (character or NULL)
check_fail_early <- function(all_qc_df, n_donors_processed,
                             threshold_10 = 0.50, threshold_50 = 0.20) {
  if (nrow(all_qc_df) == 0) {
    return(list(should_stop = FALSE, reason = NULL))
  }

  overall_violation_rate <- mean(all_qc_df$fail_any, na.rm = TRUE)

  if (n_donors_processed >= 10 && n_donors_processed < 50) {
    if (overall_violation_rate > threshold_10) {
      reason <- sprintf(
        "FAIL-EARLY: After %d donors, %.1f%% of cells violate QC thresholds (limit: %.0f%%). DATA IS LIKELY RAW/PRE-QC.",
        n_donors_processed, overall_violation_rate * 100, threshold_10 * 100
      )
      return(list(should_stop = TRUE, reason = reason))
    }
  }

  if (n_donors_processed >= 50) {
    if (overall_violation_rate > threshold_50) {
      reason <- sprintf(
        "FAIL-EARLY: After %d donors, %.1f%% of cells violate QC thresholds (limit: %.0f%%). DATA IS LIKELY RAW/PRE-QC.",
        n_donors_processed, overall_violation_rate * 100, threshold_50 * 100
      )
      return(list(should_stop = TRUE, reason = reason))
    }
  }

  return(list(should_stop = FALSE, reason = NULL))
}

# ==============================================================================
# SHARD FILE UTILITIES
# ==============================================================================

#' List donor shard files in a directory
#'
#' @param shard_dir Path to shard directory
#' @return Character vector of shard file paths (excluding manifest/progress files)
list_donor_shards <- function(shard_dir) {
  all_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  # Exclude manifest and progress files
  all_files <- all_files[!grepl("MANIFEST|progress", basename(all_files))]
  return(sort(all_files))
}

#' Extract donor ID from shard filename
#'
#' @param shard_path Path to shard file
#' @return Character string donor ID
extract_donor_id <- function(shard_path) {
  fname <- basename(shard_path)
  sub("^donor_", "", sub("\\.qs$", "", fname))
}

# ==============================================================================
# SANITY CHECKS (per three-stage-pipeline skill)
# ==============================================================================

#' Run sanity checks on a count matrix
#'
#' @param counts_mat A sparse or dense count matrix
#' @param description Character description for logging
#' @return Logical TRUE if all checks pass, FALSE otherwise
sanity_check_counts <- function(counts_mat, description = "Count matrix") {
  cat(sprintf("\n[SANITY CHECK] %s\n", description))
  issues <- c()

  # Check 1: All integers
  sample_vals <- as.numeric(counts_mat[1:min(1000, nrow(counts_mat)),
                                        1:min(100, ncol(counts_mat))])
  sample_vals <- sample_vals[!is.na(sample_vals)]
  non_int_pct <- mean(sample_vals != floor(sample_vals)) * 100
  if (non_int_pct > 0) {
    issues <- c(issues, sprintf("Non-integer values: %.2f%%", non_int_pct))
  }

  # Check 2: No negatives
  if (any(sample_vals < 0)) {
    issues <- c(issues, "CRITICAL: Negative counts detected")
  }

  # Check 3: Max value reasonable
  max_val <- max(sample_vals)
  if (max_val <= 1) {
    issues <- c(issues, sprintf("Max value = %.2f (expect >1 for raw counts — may be normalized)", max_val))
  }

  # Check 4: Sparsity
  zero_pct <- mean(sample_vals == 0) * 100
  if (zero_pct < 50 || zero_pct > 99.9) {
    issues <- c(issues, sprintf("Unusual sparsity: %.1f%% zeros (expect 70-95%%)", zero_pct))
  }

  if (length(issues) == 0) {
    cat("  ✅ PASS - All sanity checks passed\n")
    return(TRUE)
  } else {
    cat("  🔴 SANITY CHECK FAILURES:\n")
    for (issue in issues) cat(sprintf("     ⚠️  %s\n", issue))
    return(FALSE)
  }
}
