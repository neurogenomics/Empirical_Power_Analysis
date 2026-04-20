#!/usr/bin/env Rscript
# ==============================================================================
# TESTS: Depth Equalization Pipeline - Unit Tests
# ==============================================================================
# Run with: Rscript R/depth_equalization/tests/test_depth_equalization.R
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(SingleCellExperiment)
  library(Matrix)
})

# Source utils
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "depth_equalization")
source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

cat("\n========================================\n")
cat("  Depth Equalization Unit Tests\n")
cat("========================================\n\n")

# ==============================================================================
# Helper: Create a toy SCE for testing
# ==============================================================================

make_toy_sce <- function(n_genes = 100, n_cells = 50, mean_reads = 10000, seed = 42) {
  set.seed(seed)
  
  # Generate realistic-ish count matrix
  # Each cell gets a random total, then counts are distributed across genes
  cell_totals <- rpois(n_cells, lambda = mean_reads)
  
  counts_list <- lapply(cell_totals, function(total) {
    # Random gene probabilities
    probs <- rgamma(n_genes, shape = 0.5, rate = 1)
    probs <- probs / sum(probs)
    rmultinom(1, size = total, prob = probs)[, 1]
  })
  
  counts_mat <- do.call(cbind, counts_list)
  counts_mat <- Matrix::Matrix(counts_mat, sparse = TRUE)
  
  rownames(counts_mat) <- paste0("ENSG", sprintf("%011d", seq_len(n_genes)))
  colnames(counts_mat) <- paste0("cell_", seq_len(n_cells))
  
  col_data <- DataFrame(
    donor_id = rep("Donor_TEST", n_cells),
    sex = rep("female", n_cells),
    development_stage = rep("65-year-old human stage", n_cells),
    Source = rep("HBCC", n_cells),
    self_reported_ethnicity = rep("European", n_cells),
    disease = rep("normal", n_cells),
    cell_type = rep("Astro", n_cells)
  )
  
  SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = col_data
  )
}

# ==============================================================================
# TEST 1: Downsampled totals match expected fraction
# ==============================================================================

test_that("Downsampled cell totals match expected fraction", {
  set.seed(123)
  n_genes <- 200
  n_cells <- 100
  mean_reads <- 10000
  fraction <- 0.4679
  
  sce <- make_toy_sce(n_genes = n_genes, n_cells = n_cells, mean_reads = mean_reads)
  original_counts <- assay(sce, "counts")
  original_totals <- colSums(original_counts)
  
  ds_counts <- downsample_counts_matrix(original_counts, fraction, seed = 42)
  ds_totals <- colSums(ds_counts)
  
  # Each cell's new total should be round(original * fraction)
  expected_totals <- round(original_totals * fraction)
  
  expect_equal(ds_totals, expected_totals,
               info = "Downsampled totals should exactly match round(original * fraction)")
})

# ==============================================================================
# TEST 2: Dimensions are preserved
# ==============================================================================

test_that("Matrix dimensions are preserved after downsampling", {
  sce <- make_toy_sce(n_genes = 150, n_cells = 30)
  original_counts <- assay(sce, "counts")
  
  ds_counts <- downsample_counts_matrix(original_counts, 0.5, seed = 99)
  
  expect_equal(nrow(ds_counts), nrow(original_counts),
               info = "Gene count should be preserved")
  expect_equal(ncol(ds_counts), ncol(original_counts),
               info = "Cell count should be preserved")
  expect_equal(rownames(ds_counts), rownames(original_counts),
               info = "Gene names should be preserved")
  expect_equal(colnames(ds_counts), colnames(original_counts),
               info = "Cell names should be preserved")
})

# ==============================================================================
# TEST 3: Gene proportions are approximately preserved
# ==============================================================================

test_that("Gene proportions are approximately preserved", {
  # Use a realistic gene expression distribution: few highly expressed genes,
  # many lowly expressed genes. High dynamic range ensures correlation is meaningful.
  set.seed(77)
  n_genes <- 100
  # Negative binomial-like: high variance, wide dynamic range
  gene_means <- c(rep(5000, 5), rep(1000, 15), rep(100, 30), rep(10, 50))
  counts_vec <- rpois(n_genes, lambda = gene_means)
  total_reads <- sum(counts_vec)  # ~50,000 total
  
  counts_mat <- Matrix::Matrix(matrix(counts_vec, ncol = 1), sparse = TRUE)
  rownames(counts_mat) <- paste0("gene_", seq_len(n_genes))
  colnames(counts_mat) <- "cell_1"
  
  original_props <- counts_vec / total_reads
  
  ds_mat <- downsample_counts_matrix(counts_mat, 0.4679, seed = 42)
  ds_vec <- as.numeric(ds_mat[, 1])
  ds_props <- ds_vec / sum(ds_vec)
  
  # Correlation between original and downsampled proportions should be very high
  # With high dynamic range and ~50K reads, multinomial preserves structure well
  prop_cor <- cor(original_props, ds_props)
  expect_gt(prop_cor, 0.99,
            label = sprintf("Gene proportion correlation (%.4f) should be > 0.99", prop_cor))
})

# ==============================================================================
# TEST 4: Reproducibility with same seed
# ==============================================================================

test_that("Same seed produces identical results", {
  sce <- make_toy_sce(n_genes = 100, n_cells = 20)
  counts <- assay(sce, "counts")
  
  ds1 <- downsample_counts_matrix(counts, 0.4679, seed = 42)
  ds2 <- downsample_counts_matrix(counts, 0.4679, seed = 42)
  
  expect_identical(as.matrix(ds1), as.matrix(ds2),
                   info = "Same seed should produce identical downsampled matrices")
})

# ==============================================================================
# TEST 5: Different seeds produce different results
# ==============================================================================

test_that("Different seeds produce different results", {
  sce <- make_toy_sce(n_genes = 100, n_cells = 20)
  counts <- assay(sce, "counts")
  
  ds1 <- downsample_counts_matrix(counts, 0.4679, seed = 42)
  ds2 <- downsample_counts_matrix(counts, 0.4679, seed = 43)
  
  # They should differ (extremely unlikely to be identical by chance)
  expect_false(identical(as.matrix(ds1), as.matrix(ds2)),
               info = "Different seeds should produce different results")
})

# ==============================================================================
# TEST 6: No negative counts
# ==============================================================================

test_that("No negative counts in downsampled matrix", {
  sce <- make_toy_sce(n_genes = 100, n_cells = 50)
  counts <- assay(sce, "counts")
  
  ds <- downsample_counts_matrix(counts, 0.4679, seed = 42)
  
  expect_true(all(ds >= 0),
              info = "All counts should be non-negative")
})

# ==============================================================================
# TEST 7: All counts are integers
# ==============================================================================

test_that("All downsampled counts are integers", {
  sce <- make_toy_sce(n_genes = 100, n_cells = 50)
  counts <- assay(sce, "counts")
  
  ds <- downsample_counts_matrix(counts, 0.4679, seed = 42)
  ds_dense <- as.matrix(ds)
  
  expect_true(all(ds_dense == floor(ds_dense)),
              info = "All counts should be integers")
})

# ==============================================================================
# TEST 8: Shard-level downsampling preserves metadata
# ==============================================================================

test_that("downsample_shard preserves all SCE metadata", {
  # Create and save a toy shard
  sce <- make_toy_sce(n_genes = 100, n_cells = 30)
  tmp_path <- tempfile(fileext = ".qs")
  qs::qsave(sce, tmp_path)
  
  ds_sce <- downsample_shard(tmp_path, 0.4679, seed = 42)
  
  # Check metadata preservation
  expect_equal(ncol(ds_sce), ncol(sce), info = "Cell count preserved")
  expect_equal(nrow(ds_sce), nrow(sce), info = "Gene count preserved")
  expect_equal(colData(ds_sce), colData(sce), info = "colData preserved")
  expect_equal(rownames(ds_sce), rownames(sce), info = "Gene names preserved")
  expect_equal(colnames(ds_sce), colnames(sce), info = "Cell names preserved")
  
  # Check counts are actually different (downsampled)
  expect_true(sum(assay(ds_sce, "counts")) < sum(assay(sce, "counts")),
              info = "Total reads should be reduced")
  
  # Cleanup
  unlink(tmp_path)
})

# ==============================================================================
# TEST 9: Fraction validation
# ==============================================================================

test_that("Invalid fractions are rejected", {
  sce <- make_toy_sce(n_genes = 10, n_cells = 5)
  counts <- assay(sce, "counts")
  
  expect_error(downsample_counts_matrix(counts, 0, seed = 42),
               "fraction must be in")
  expect_error(downsample_counts_matrix(counts, 1, seed = 42),
               "fraction must be in")
  expect_error(downsample_counts_matrix(counts, -0.5, seed = 42),
               "fraction must be in")
  expect_error(downsample_counts_matrix(counts, 1.5, seed = 42),
               "fraction must be in")
})

# ==============================================================================
# TEST 10: Empty cells handled correctly
# ==============================================================================

test_that("Cells with zero counts remain zero", {
  # Create matrix with one empty cell
  counts_mat <- Matrix::Matrix(
    cbind(c(10, 20, 30), c(0, 0, 0), c(5, 15, 25)),
    sparse = TRUE
  )
  rownames(counts_mat) <- paste0("gene_", 1:3)
  colnames(counts_mat) <- paste0("cell_", 1:3)
  
  ds <- downsample_counts_matrix(counts_mat, 0.5, seed = 42)
  
  # Empty cell should remain empty
  expect_equal(sum(ds[, 2]), 0,
               info = "Empty cell should remain all zeros after downsampling")
  
  # Non-empty cells should have reduced counts
  expect_lt(sum(ds[, 1]), sum(counts_mat[, 1]))
  expect_lt(sum(ds[, 3]), sum(counts_mat[, 3]))
})

# ==============================================================================
# TEST 11: Global mean RPC matches target (statistical test)
# ==============================================================================

test_that("Global mean RPC approximately matches target after downsampling", {
  # Simulate a population of cells with Astro-like depth distribution
  set.seed(42)
  n_cells <- 500
  n_genes <- 50
  
  # Astro-like: mean ~9456, right-skewed
  cell_totals <- rlnorm(n_cells, meanlog = log(8000), sdlog = 0.5)
  cell_totals <- round(cell_totals)
  
  counts_list <- lapply(cell_totals, function(total) {
    probs <- rgamma(n_genes, shape = 0.3, rate = 1)
    probs <- probs / sum(probs)
    rmultinom(1, size = max(1, total), prob = probs)[, 1]
  })
  
  counts_mat <- Matrix::Matrix(do.call(cbind, counts_list), sparse = TRUE)
  rownames(counts_mat) <- paste0("gene_", seq_len(n_genes))
  colnames(counts_mat) <- paste0("cell_", seq_len(n_cells))
  
  original_mean_rpc <- mean(colSums(counts_mat))
  fraction <- 4426 / original_mean_rpc
  
  ds <- downsample_counts_matrix(counts_mat, fraction, seed = 42)
  ds_mean_rpc <- mean(colSums(ds))
  
  # Should be within 2% of target (rounding introduces small deviations)
  relative_error <- abs(ds_mean_rpc - 4426) / 4426
  expect_lt(relative_error, 0.02,
            label = sprintf("Mean RPC %.1f should be within 2%% of 4426 (error: %.4f)", 
                            ds_mean_rpc, relative_error))
})

# ==============================================================================
# TEST 12: Manifest generation format
# ==============================================================================

test_that("Generated manifest has correct structure", {
  # Create a temp directory with a few toy shards
  tmp_dir <- tempdir()
  shard_dir <- file.path(tmp_dir, "test_shards")
  dir.create(shard_dir, showWarnings = FALSE)
  
  for (i in 1:3) {
    sce <- make_toy_sce(n_genes = 50, n_cells = 10 * i, seed = i)
    qs::qsave(sce, file.path(shard_dir, sprintf("donor_Donor_%d.qs", i)))
  }
  
  config <- list(
    cell_type = "Astro",
    min_age = 40,
    min_cells = 30,
    source_dataset = "Roussos",
    target_mean_rpc = 4426,
    downsample_fraction = 0.4679
  )
  
  manifest <- generate_manifest(shard_dir, config)
  
  # Check structure
  expect_true("manifest" %in% names(manifest))
  expect_true("feasibility" %in% names(manifest))
  expect_true("config" %in% names(manifest))
  
  # Check manifest data frame
  expect_equal(nrow(manifest$manifest), 3)
  expect_true(all(c("donor_id", "cell_count", "file_path") %in% names(manifest$manifest)))
  
  # Check config
  expect_true(manifest$config$depth_equalized)
  expect_equal(manifest$config$source_dataset, "Roussos")
  expect_equal(manifest$config$target_mean_rpc, 4426)
  expect_equal(manifest$config$n_donors, 3)
  
  # Cleanup
  unlink(shard_dir, recursive = TRUE)
})

# ==============================================================================
# RUN ALL TESTS
# ==============================================================================

cat("\n========================================\n")
cat("  All tests completed!\n")
cat("========================================\n")
