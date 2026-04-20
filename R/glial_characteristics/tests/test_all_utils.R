# ==============================================================================
# TESTS: Glial Characteristics Pipeline - Unit Tests for Utility Functions
# ==============================================================================
# Run with: Rscript tests/test_all_utils.R
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
})

# Load utilities
source(file.path(dirname(dirname(sys.frame(1)$ofile %||% ".")), "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### RUNNING UTILS TESTS ###\n")

# ==============================================================================
# TEST: Summary Statistics
# ==============================================================================

test_that("compute_summary_stats returns correct structure", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  result <- compute_summary_stats(x, "test_metric")
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$n, 10)
  expect_equal(result$mean, 5.5)
  expect_equal(result$median, 5.5)
  expect_equal(result$min, 1)
  expect_equal(result$max, 10)
  expect_equal(result$metric, "test_metric")
})

# ==============================================================================
# TEST: Bootstrap Functions
# ==============================================================================

test_that("bootstrap_median returns correct number of samples", {
  x <- rnorm(100)
  boot_vals <- bootstrap_median(x, n_boot = 100, seed = 42)
  
  expect_length(boot_vals, 100)
  expect_true(all(is.numeric(boot_vals)))
})

test_that("bootstrap_median is reproducible with same seed", {
  x <- rnorm(100)
  boot1 <- bootstrap_median(x, n_boot = 100, seed = 42)
  boot2 <- bootstrap_median(x, n_boot = 100, seed = 42)
  
  expect_equal(boot1, boot2)
})

test_that("bootstrap_median_ratio returns correct number of samples", {
  x1 <- rnorm(100, mean = 10)
  x2 <- rnorm(100, mean = 5)
  boot_ratios <- bootstrap_median_ratio(x1, x2, n_boot = 100, seed = 42)
  
  expect_length(boot_ratios, 100)
  expect_true(all(is.numeric(boot_ratios)))
  # Expect ratios to be roughly around 2 (10/5)
  expect_true(mean(boot_ratios) > 1.5 && mean(boot_ratios) < 2.5)
})

test_that("bootstrap_ci returns correct structure", {
  boot_vals <- rnorm(1000, mean = 5, sd = 1)
  ci <- bootstrap_ci(boot_vals, ci_level = 0.95)
  
  expect_named(ci, c("lower", "median", "upper"))
  expect_true(ci["lower"] < ci["median"])
  expect_true(ci["median"] < ci["upper"])
  # 95% CI should span roughly 4 SDs for normal
  expect_true(ci["upper"] - ci["lower"] < 4)
})

# ==============================================================================
# TEST: Subsample Function
# ==============================================================================

test_that("subsample_df respects max_rows", {
  df <- data.frame(x = 1:1000, y = rnorm(1000))
  
  result <- subsample_df(df, max_rows = 100, seed = 42)
  expect_equal(nrow(result), 100)
  
  # If df smaller than max_rows, return unchanged
  small_df <- data.frame(x = 1:50, y = rnorm(50))
  result_small <- subsample_df(small_df, max_rows = 100, seed = 42)
  expect_equal(nrow(result_small), 50)
})

test_that("subsample_df is reproducible", {
  df <- data.frame(x = 1:1000, y = rnorm(1000))
  
  result1 <- subsample_df(df, max_rows = 100, seed = 42)
  result2 <- subsample_df(df, max_rows = 100, seed = 42)
  

  expect_equal(result1$x, result2$x)
})

# ==============================================================================
# TEST: Sex Harmonization
# ==============================================================================

test_that("harmonise_sex works for Roussos encoding", {
  sex_vals <- c("male", "female", "male", "female")
  encoding <- c(male = "male", female = "female")
  
  result <- harmonise_sex(sex_vals, encoding)
  expect_equal(result, sex_vals)
})

test_that("harmonise_sex works for DeJager encoding", {
  sex_vals <- c("1", "0", "1", "0")
  encoding <- c(male = "1", female = "0")
  
  result <- harmonise_sex(sex_vals, encoding)
  expect_equal(result, c("male", "female", "male", "female"))
})

test_that("harmonise_sex handles unknown values", {
  sex_vals <- c("1", "0", "unknown", "NA")
  encoding <- c(male = "1", female = "0")
  
  result <- harmonise_sex(sex_vals, encoding)
  expect_equal(result[1:2], c("male", "female"))
  expect_true(is.na(result[3]))
  expect_true(is.na(result[4]))
})

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n### ALL TESTS PASSED ###\n")
