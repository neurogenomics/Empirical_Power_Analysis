# ==============================================================================
# SCRIPT: Unit Tests for QC Validation Utilities
# PURPOSE: Stage 1 validation — test all QC metric functions with mock data
# RUN: Rscript R/qc_validation/tests/test_qc_utils.R
# ==============================================================================

library(testthat)

# Source config and utils (hardcoded path for reliability)
source(file.path(PROJECT_ROOT, "R", "qc_validation", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

cat("\n========================================\n")
cat("Stage 1: Unit Tests for QC Validation\n")
cat("========================================\n\n")

# ==============================================================================
# Helper: Create mock SingleCellExperiment
# ==============================================================================

create_mock_sce <- function(n_cells = 100, n_genes = 500,
                            mean_umis = 5000, include_mito = TRUE,
                            mito_pct_target = 2.0) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
  })

  # Create sparse count matrix
  gene_names <- paste0("Gene_", seq_len(n_genes))
  cell_names <- paste0("Cell_", seq_len(n_cells))

  # Generate counts with realistic sparsity (~85%)
  set.seed(42)
  n_nonzero <- round(n_genes * n_cells * 0.15)
  rows <- sample(n_genes, n_nonzero, replace = TRUE)
  cols <- sample(n_cells, n_nonzero, replace = TRUE)
  vals <- rpois(n_nonzero, lambda = 3) + 1L

  counts_mat <- sparseMatrix(
    i = rows, j = cols, x = vals,
    dims = c(n_genes, n_cells),
    dimnames = list(gene_names, cell_names)
  )

  # Add mitochondrial genes
  if (include_mito) {
    mito_names <- c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6")
    n_mito <- length(mito_names)
    gene_names[1:n_mito] <- mito_names
    rownames(counts_mat)[1:n_mito] <- mito_names

    # Adjust mito counts to approximate target percentage
    total_per_cell <- Matrix::colSums(counts_mat)
    target_mito_per_cell <- total_per_cell * (mito_pct_target / 100)
    for (m in 1:n_mito) {
      counts_mat[m, ] <- rpois(n_cells, lambda = target_mito_per_cell / n_mito)
    }
  }

  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat)
  )

  return(sce)
}

# ==============================================================================
# Test: identify_mito_genes
# ==============================================================================

test_that("identify_mito_genes works with symbol list (new format)", {
  cat("[TEST] identify_mito_genes with symbol list (new format)...\n")

  gene_names <- c("MT-CO1", "MT-ND1", "GAPDH", "ACTB", "MT-ATP6")
  mito_list <- list(
    ensembl_ids = c("ENSG00000198804", "ENSG00000198712"),
    symbols = c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6", "MT-CYB")
  )

  result <- identify_mito_genes(gene_names, mito_list)

  expect_equal(result, c(TRUE, TRUE, FALSE, FALSE, TRUE))
  cat("  ✓ Correctly identified 3/5 genes as mitochondrial via symbols\n")
})

test_that("identify_mito_genes works with Ensembl IDs (new format)", {
  cat("[TEST] identify_mito_genes with Ensembl IDs (new format)...\n")

  gene_names <- c("ENSG00000198804", "ENSG00000198712", "ENSG00000000003", "ENSG00000000005")
  mito_list <- list(
    ensembl_ids = c("ENSG00000198804", "ENSG00000198712", "ENSG00000198886"),
    symbols = c("MT-CO1", "MT-CO2", "MT-ND4")
  )

  result <- identify_mito_genes(gene_names, mito_list)

  expect_equal(result, c(TRUE, TRUE, FALSE, FALSE))
  cat("  ✓ Correctly identified 2/4 genes as mitochondrial via Ensembl IDs\n")
})

test_that("identify_mito_genes works with legacy character vector", {
  cat("[TEST] identify_mito_genes with legacy character vector...\n")

  gene_names <- c("MT-CO1", "MT-ND1", "GAPDH", "ACTB", "MT-ATP6")
  mito_list <- c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6", "MT-CYB")

  result <- identify_mito_genes(gene_names, mito_list)

  expect_equal(result, c(TRUE, TRUE, FALSE, FALSE, TRUE))
  cat("  ✓ Correctly identified 3/5 genes via legacy format\n")
})

test_that("identify_mito_genes falls back to ^MT- prefix", {
  cat("[TEST] identify_mito_genes with ^MT- fallback...\n")

  gene_names <- c("MT-CO1", "MT-ND1", "GAPDH", "ACTB", "MT-ATP6")

  result <- identify_mito_genes(gene_names, mito_genes = NULL)

  expect_equal(result, c(TRUE, TRUE, FALSE, FALSE, TRUE))
  cat("  ✓ Correctly fell back to ^MT- prefix matching\n")
})

test_that("identify_mito_genes handles no matches gracefully", {
  cat("[TEST] identify_mito_genes with no matches...\n")

  gene_names <- c("GAPDH", "ACTB", "TP53")
  mito_list <- list(
    ensembl_ids = c("ENSG00000198804"),
    symbols = c("MT-CO1", "MT-ND1")
  )

  # Neither Ensembl IDs nor symbols match → falls back to ^MT- prefix (also nothing)
  result <- identify_mito_genes(gene_names, mito_list)

  expect_equal(sum(result), 0)
  cat("  ✓ Correctly returned no matches\n")
})

# ==============================================================================
# Test: compute_cell_qc_metrics
# ==============================================================================

test_that("compute_cell_qc_metrics returns correct structure", {
  cat("[TEST] compute_cell_qc_metrics structure...\n")

  sce <- create_mock_sce(n_cells = 50, n_genes = 200)
  mito_list <- list(
    ensembl_ids = character(0),
    symbols = c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6")
  )

  result <- compute_cell_qc_metrics(sce, mito_list)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 50)
  expect_true(all(c("n_umis", "n_genes", "pct_mito") %in% names(result)))
  expect_true(all(result$n_umis >= 0))
  expect_true(all(result$n_genes >= 0))
  expect_true(all(result$pct_mito >= 0 & result$pct_mito <= 100, na.rm = TRUE))

  cat(sprintf("  ✓ Returned %d rows with correct columns\n", nrow(result)))
  cat(sprintf("  ✓ UMI range: [%d, %d]\n", min(result$n_umis), max(result$n_umis)))
  cat(sprintf("  ✓ Gene range: [%d, %d]\n", min(result$n_genes), max(result$n_genes)))
  cat(sprintf("  ✓ Mito%% range: [%.2f, %.2f]\n", min(result$pct_mito), max(result$pct_mito)))
})

test_that("compute_cell_qc_metrics handles no mito genes", {
  cat("[TEST] compute_cell_qc_metrics with no mito genes...\n")

  sce <- create_mock_sce(n_cells = 20, n_genes = 100, include_mito = FALSE)

  result <- compute_cell_qc_metrics(sce, mito_genes = list(
    ensembl_ids = character(0),
    symbols = c("MT-CO1", "MT-ND1")
  ))

  expect_true(is.data.frame(result))
  # Should fall back to ^MT- prefix, which also finds nothing → pct_mito = NA
  # (unless mock data happens to have MT- prefixed genes)
  cat(sprintf("  ✓ Handled gracefully, pct_mito NAs: %d/%d\n",
    sum(is.na(result$pct_mito)), nrow(result)))
})

# ==============================================================================
# Test: check_thresholds
# ==============================================================================

test_that("check_thresholds correctly flags violations", {
  cat("[TEST] check_thresholds...\n")

  qc_df <- data.frame(
    n_umis = c(500, 2000, 50000, 120000, 5000),
    n_genes = c(800, 2000, 5000, 8000, 15000),
    pct_mito = c(1.0, 3.0, 4.9, 6.0, 0.5)
  )

  result <- check_thresholds(qc_df)

  # Cell 1: low UMIs (500 < 1500), low genes (800 < 1100)
  expect_true(result$fail_low_umis[1])
  expect_true(result$fail_low_genes[1])
  expect_false(result$fail_high_mito[1])

  # Cell 2: all pass
  expect_false(result$fail_any[2])

  # Cell 3: all pass (4.9 < 5.0)
  expect_false(result$fail_any[3])

  # Cell 4: high UMIs (120000 > 110000), high mito (6.0 >= 5.0)
  expect_true(result$fail_high_umis[4])
  expect_true(result$fail_high_mito[4])

  # Cell 5: high genes (15000 > 12500)
  expect_true(result$fail_high_genes[5])

  cat("  ✓ All threshold checks correct\n")
})

# ==============================================================================
# Test: summarise_violations
# ==============================================================================

test_that("summarise_violations returns correct counts", {
  cat("[TEST] summarise_violations...\n")

  qc_df <- data.frame(
    n_umis = c(500, 2000, 5000),
    n_genes = c(800, 2000, 5000),
    pct_mito = c(1.0, 3.0, 0.5)
  )
  qc_df <- check_thresholds(qc_df)
  result <- summarise_violations(qc_df)

  expect_equal(result$n_total, 3)
  expect_equal(result$counts[["fail_low_umis"]], 1)
  expect_equal(result$counts[["fail_low_genes"]], 1)
  expect_equal(result$counts[["fail_any"]], 1)  # Only cell 1 fails

  cat(sprintf("  ✓ Correctly counted %d violations out of %d cells\n",
    result$counts[["fail_any"]], result$n_total))
})

# ==============================================================================
# Test: summarise_donor
# ==============================================================================

test_that("summarise_donor returns single-row summary", {
  cat("[TEST] summarise_donor...\n")

  qc_df <- data.frame(
    n_umis = c(2000, 3000, 4000),
    n_genes = c(1500, 2000, 2500),
    pct_mito = c(1.0, 2.0, 3.0)
  )
  qc_df <- check_thresholds(qc_df)
  result <- summarise_donor(qc_df, "Donor_Test")

  expect_equal(nrow(result), 1)
  expect_equal(result$donor_id, "Donor_Test")
  expect_equal(result$n_cells, 3)
  expect_equal(result$median_umis, 3000)
  expect_equal(result$median_genes, 2000)
  expect_equal(result$n_fail_any, 0)

  cat("  ✓ Single-row summary with correct values\n")
})

# ==============================================================================
# Test: check_fail_early
# ==============================================================================

test_that("check_fail_early triggers at 10 donors with >50% violations", {
  cat("[TEST] check_fail_early at 10 donors...\n")

  # Create data where 60% of cells fail
  qc_df <- data.frame(
    fail_any = c(rep(TRUE, 60), rep(FALSE, 40))
  )

  result <- check_fail_early(qc_df, n_donors_processed = 10)
  expect_true(result$should_stop)
  expect_true(grepl("FAIL-EARLY", result$reason))

  cat("  ✓ Correctly triggered fail-early at 10 donors\n")
})

test_that("check_fail_early does NOT trigger with low violations", {
  cat("[TEST] check_fail_early with low violations...\n")

  qc_df <- data.frame(
    fail_any = c(rep(TRUE, 5), rep(FALSE, 95))
  )

  result <- check_fail_early(qc_df, n_donors_processed = 10)
  expect_false(result$should_stop)

  result <- check_fail_early(qc_df, n_donors_processed = 50)
  expect_false(result$should_stop)

  cat("  ✓ Correctly did NOT trigger fail-early\n")
})

test_that("check_fail_early triggers at 50 donors with >20% violations", {
  cat("[TEST] check_fail_early at 50 donors...\n")

  qc_df <- data.frame(
    fail_any = c(rep(TRUE, 25), rep(FALSE, 75))
  )

  result <- check_fail_early(qc_df, n_donors_processed = 50)
  expect_true(result$should_stop)

  cat("  ✓ Correctly triggered fail-early at 50 donors\n")
})

# ==============================================================================
# Test: list_donor_shards and extract_donor_id
# ==============================================================================

test_that("extract_donor_id parses filenames correctly", {
  cat("[TEST] extract_donor_id...\n")

  expect_equal(extract_donor_id("/path/to/donor_Donor_123.qs"), "Donor_123")
  expect_equal(extract_donor_id("donor_Donor_1.qs"), "Donor_1")
  expect_equal(extract_donor_id("/a/b/donor_Donor_1000.qs"), "Donor_1000")

  cat("  ✓ All donor ID extractions correct\n")
})

# ==============================================================================
# Test: sanity_check_counts
# ==============================================================================

test_that("sanity_check_counts passes for valid count matrix", {
  cat("[TEST] sanity_check_counts with valid data...\n")

  set.seed(42)
  mat <- Matrix::sparseMatrix(
    i = sample(100, 500, replace = TRUE),
    j = sample(50, 500, replace = TRUE),
    x = rpois(500, 5) + 1L,
    dims = c(100, 50)
  )

  result <- sanity_check_counts(mat, "Mock valid counts")
  expect_true(result)
  cat("  ✓ Passed for valid integer count matrix\n")
})

test_that("sanity_check_counts fails for normalized data", {
  cat("[TEST] sanity_check_counts with normalized data...\n")

  set.seed(42)
  mat <- matrix(rnorm(500, mean = 0, sd = 1), nrow = 100, ncol = 50)

  result <- sanity_check_counts(mat, "Mock normalized data")
  expect_false(result)
  cat("  ✓ Correctly failed for non-integer normalized data\n")
})

# ==============================================================================
# Summary
# ==============================================================================

cat("\n========================================\n")
cat("Stage 1: All unit tests PASSED\n")
cat("========================================\n")
