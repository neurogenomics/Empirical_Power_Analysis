# ==============================================================================
# UNIT TESTS: EDA Roussos Utils - All Functions
# PURPOSE: Test all utility functions from utils.R
# USAGE: Rscript R/eda_roussos/tests/test_all_utils.R
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(SingleCellExperiment)
  library(Matrix)
})

# Load utils - use commandArgs for Rscript compatibility
script_path <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # Fallback for Rscript: use current working directory
  file.path(getwd(), "R", "eda_roussos", "tests")
})

source(file.path(script_path, "..", "utils.R"))

message("\n### TESTING EDA ROUSSOS UTILS ###\n")

# ==============================================================================
# 1. AGE PARSING TESTS
# ==============================================================================

test_that("Age parsing handles year-old stage format", {
  input <- c("80 year-old stage", "65 year-old stage", "90 year-old stage")
  result <- parse_age_stage_to_years(input)
  expect_equal(result, c(80, 65, 90))
})

test_that("Age parsing handles month-old stage format", {
  input <- c("6 month-old stage", "12 month-old stage", "18 month-old stage")
  result <- parse_age_stage_to_years(input)
  expect_equal(result, c(0.5, 1, 1.5))
})

test_that("Age parsing handles week-old stage format", {
  input <- c("52 week-old stage", "104 week-old stage")
  result <- parse_age_stage_to_years(input)
  expect_equal(result, c(52 / 52.1429, 104 / 52.1429), tolerance = 0.01)
})

test_that("Age parsing handles day-old stage format", {
  input <- c("365 day-old stage", "730 day-old stage")
  result <- parse_age_stage_to_years(input)
  expect_equal(result, c(365 / 365.25, 730 / 365.25), tolerance = 0.01)
})

test_that("Age parsing handles unknown/NA values", {
  input <- c("unknown", "not available", "NA", "n/a")
  result <- parse_age_stage_to_years(input)
  expect_true(all(is.na(result)))
})

test_that("looks_like_age_stage detects age format", {
  valid <- c("80 year-old stage", "6 month-old stage")
  invalid <- c("80", "6 months", "adult")
  
  expect_true(looks_like_age_stage(valid))
  expect_false(looks_like_age_stage(invalid))
})

# ==============================================================================
# 2. GENE STATISTICS TESTS
# ==============================================================================

test_that("calculate_mean_expression works correctly", {
  mat <- Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, sparse = TRUE)
  rownames(mat) <- c("Gene1", "Gene2")
  
  result <- calculate_mean_expression(mat)
  expect_equal(as.numeric(result), c(3, 4))  # Gene1: (1+3+5)/3=3, Gene2: (2+4+6)/3=4
})

test_that("calculate_detection_rate works correctly", {
  mat <- Matrix(c(0, 1, 2, 0, 3, 0), nrow = 2, ncol = 3, sparse = TRUE)
  rownames(mat) <- c("Gene1", "Gene2")
  
  result <- calculate_detection_rate(mat)
  expect_equal(as.numeric(result), c(2/3, 1/3))  # Gene1: 2/3 cells, Gene2: 1/3 cells
})

test_that("identify_zero_inflated works correctly", {
  mat <- Matrix(c(0, 0, 0, 1, 0, 1), nrow = 2, ncol = 3, sparse = TRUE)
  rownames(mat) <- c("Gene1", "Gene2")
  
  result <- identify_zero_inflated(mat, threshold = 0.5)
  expect_equal(as.logical(result), c(TRUE, FALSE))  # Gene1: 100% zeros, Gene2: 33% zeros
})

test_that("get_top_expressed_genes returns correct format", {
  mat <- Matrix(c(1, 5, 2, 6, 3, 7), nrow = 3, ncol = 2, sparse = TRUE)
  rownames(mat) <- c("GeneA", "GeneB", "GeneC")
  
  result <- get_top_expressed_genes(mat, n = 2)
  
  expect_equal(nrow(result), 2)
  expect_equal(colnames(result), c("gene", "mean_expression", "rank"))
  expect_equal(result$gene[1], "GeneC")  # Highest mean: (3+7)/2 = 5
  expect_equal(result$rank, c(1, 2))
})

# ==============================================================================
# 3. QC METRIC TESTS
# ==============================================================================

test_that("calculate_mt_percent works correctly", {
  counts <- Matrix(c(10, 90, 20, 80), nrow = 2, ncol = 2, sparse = TRUE)
  rownames(counts) <- c("MT-CO1", "ACTB")
  
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  result <- calculate_mt_percent(sce, mt_pattern = "^MT-")
  expect_equal(result, c(10, 20))  # Cell1: 10%, Cell2: 20%
})

test_that("calculate_mt_percent handles no MT genes", {
  counts <- Matrix(c(10, 90, 20, 80), nrow = 2, ncol = 2, sparse = TRUE)
  rownames(counts) <- c("ACTB", "GAPDH")
  
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  expect_warning(result <- calculate_mt_percent(sce, mt_pattern = "^MT-"))
  expect_true(all(is.na(result)))
})

test_that("calculate_genes_detected works correctly", {
  counts <- Matrix(c(0, 1, 2, 0, 3, 4), nrow = 3, ncol = 2, sparse = TRUE)
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  result <- calculate_genes_detected(sce)
  expect_equal(result, c(2, 2))  # Both cells have 2 genes detected
})

test_that("calculate_library_size works correctly", {
  counts <- Matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, ncol = 2, sparse = TRUE)
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  result <- calculate_library_size(sce)
  expect_equal(result, c(60, 150))  # Cell1: 10+20+30=60, Cell2: 40+50+60=150
})

# ==============================================================================
# 4. SUMMARY
# ==============================================================================

message("\n### TEST SUMMARY ###")
message("✓ All tests passed!")
message("✓ Age parsing functions (6 tests)")
message("✓ Gene statistics functions (4 tests)")
message("✓ QC metric functions (4 tests)")
message("\n✓ EDA Roussos utils are working correctly!\n")
