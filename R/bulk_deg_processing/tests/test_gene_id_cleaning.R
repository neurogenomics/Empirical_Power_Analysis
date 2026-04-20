# ==============================================================================
# UNIT TEST: Gene ID Cleaning
# ==============================================================================
# Tests the clean_gene_ids() function
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

test_gene_id_cleaning <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: Gene ID Cleaning\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Test 1: Standard Ensembl IDs with versions
  cat("[TEST 1] Standard Ensembl IDs with versions\n")
  input <- c("ENSG00000123456.12", "ENSG00000234567.5", "ENSG00000345678.1")
  expected <- c("ENSG00000123456", "ENSG00000234567", "ENSG00000345678")
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 2: IDs without versions (should remain unchanged)
  cat("\n[TEST 2] IDs without versions\n")
  input <- c("ENSG00000123456", "ENSG00000234567")
  expected <- c("ENSG00000123456", "ENSG00000234567")
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 3: Mixed (some with versions, some without)
  cat("\n[TEST 3] Mixed IDs\n")
  input <- c("ENSG00000123456.12", "ENSG00000234567", "ENSG00000345678.3")
  expected <- c("ENSG00000123456", "ENSG00000234567", "ENSG00000345678")
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", paste(expected, collapse = ", "), "\n")
    cat("    Got:", paste(result, collapse = ", "), "\n")
    all_passed <- FALSE
  }
  
  # Test 4: Empty input
  cat("\n[TEST 4] Empty input\n")
  input <- character(0)
  expected <- character(0)
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Test 5: Single ID
  cat("\n[TEST 5] Single ID\n")
  input <- "ENSG00000123456.12"
  expected <- "ENSG00000123456"
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", expected, "\n")
    cat("    Got:", result, "\n")
    all_passed <- FALSE
  }
  
  # Test 6: Multiple dots (should only remove first suffix)
  cat("\n[TEST 6] Multiple dots in ID\n")
  input <- "ENSG00000123456.12.extra"
  expected <- "ENSG00000123456"  # Should remove everything after first dot
  result <- clean_gene_ids(input)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:", expected, "\n")
    cat("    Got:", result, "\n")
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
test_result <- test_gene_id_cleaning()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
