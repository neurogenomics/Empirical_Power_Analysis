# ==============================================================================
# UNIT TEST: Interval Binning
# ==============================================================================
# Tests the bin_by_lfc_interval() function
# ==============================================================================

# Source config and utils
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/bulk_deg_processing/config.R",  # From project root
    "../config.R",                      # From tests/ directory
    "config.R"                          # From script directory
  )
  
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  
  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# TEST CASES
# ==============================================================================

test_interval_binning <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: Interval Binning\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Define test intervals
  intervals <- list(
    "1.0 to 2.0" = c(1.0, 2.0),
    "2.0 to 3.0" = c(2.0, 3.0),
    "> 3.0"      = c(3.0, Inf)
  )
  
  # Test 1: Values in first interval
  cat("[TEST 1] Values in first interval (1.0-2.0)\n")
  lfc_values <- c(1.0, 1.5, 1.9)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- rep("1.0 to 2.0", 3)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 2: Values in second interval
  cat("\n[TEST 2] Values in second interval (2.0-3.0)\n")
  lfc_values <- c(2.0, 2.5, 2.9)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- rep("2.0 to 3.0", 3)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 3: Values in third interval (>3.0)
  cat("\n[TEST 3] Values in third interval (>3.0)\n")
  lfc_values <- c(3.0, 4.5, 10.0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- rep("> 3.0", 3)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 4: Boundary values (exact boundaries)
  cat("\n[TEST 4] Boundary values\n")
  lfc_values <- c(1.0, 2.0, 3.0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- c("1.0 to 2.0", "2.0 to 3.0", "> 3.0")
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 5: Negative values (should use absolute value)
  cat("\n[TEST 5] Negative values (uses absolute value)\n")
  lfc_values <- c(-1.5, -2.5, -4.0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- c("1.0 to 2.0", "2.0 to 3.0", "> 3.0")
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 6: Mixed positive and negative
  cat("\n[TEST 6] Mixed positive and negative\n")
  lfc_values <- c(1.5, -2.5, 4.0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- c("1.0 to 2.0", "2.0 to 3.0", "> 3.0")
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 7: Values below minimum (should be NA)
  cat("\n[TEST 7] Values below minimum interval\n")
  lfc_values <- c(0.1, 0.5, 0.9)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- rep(NA_character_, 3)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected all NA, got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 8: Zero value (should be NA)
  cat("\n[TEST 8] Zero value\n")
  lfc_values <- c(0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- NA_character_
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected NA, got:", result, "\n")
    all_passed <- FALSE
  }
  
  # Test 9: Empty input
  cat("\n[TEST 9] Empty input\n")
  lfc_values <- numeric(0)
  result <- bin_by_lfc_interval(lfc_values, intervals)
  expected <- character(0)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Final summary
  cat("\n========================================\n")
  if (all_passed) {
    cat("✓ ALL TESTS PASSED\n")
  } else {
    cat("✗ SOME TESTS FAILED\n")
  }
  cat("========================================\n\n")
  
  return(all_passed)
}

# Run tests
test_result <- test_interval_binning()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
