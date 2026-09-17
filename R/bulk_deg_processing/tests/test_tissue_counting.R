# ==============================================================================
# UNIT TEST: Tissue Counting
# ==============================================================================
# Tests the count_tissues_per_gene() function
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

test_tissue_counting <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: Tissue Counting\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Test 1: Gene in multiple tissues
  cat("[TEST 1] Gene in multiple tissues\n")
  dt <- data.table(
    gene = c("GENE1", "GENE1", "GENE1", "GENE2", "GENE2"),
    tissue = c("Brain", "Heart", "Liver", "Brain", "Heart")
  )
  result <- count_tissues_per_gene(dt)
  expected <- data.table(
    gene = c("GENE1", "GENE2"),
    tissue_count = c(3L, 2L)
  )
  setkey(result, gene)
  setkey(expected, gene)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:\n")
    print(expected)
    cat("    Got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 2: Gene in single tissue
  cat("\n[TEST 2] Gene in single tissue\n")
  dt <- data.table(
    gene = c("GENE1", "GENE2"),
    tissue = c("Brain", "Heart")
  )
  result <- count_tissues_per_gene(dt)
  expected <- data.table(
    gene = c("GENE1", "GENE2"),
    tissue_count = c(1L, 1L)
  )
  setkey(result, gene)
  setkey(expected, gene)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Test 3: Duplicate gene-tissue pairs (should count unique tissues only)
  cat("\n[TEST 3] Duplicate gene-tissue pairs\n")
  dt <- data.table(
    gene = c("GENE1", "GENE1", "GENE1", "GENE1"),
    tissue = c("Brain", "Brain", "Heart", "Heart")
  )
  result <- count_tissues_per_gene(dt)
  expected <- data.table(
    gene = "GENE1",
    tissue_count = 2L
  )
  setkey(result, gene)
  setkey(expected, gene)
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected:\n")
    print(expected)
    cat("    Got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 4: Empty input
  cat("\n[TEST 4] Empty input\n")
  dt <- data.table(gene = character(0), tissue = character(0))
  result <- count_tissues_per_gene(dt)
  expected <- data.table(gene = character(0), tissue_count = integer(0))
  
  if (identical(result, expected)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Test 5: Missing required columns (should error)
  cat("\n[TEST 5] Missing required columns (should error)\n")
  dt <- data.table(gene = c("GENE1", "GENE2"))
  error_caught <- FALSE
  tryCatch({
    result <- count_tissues_per_gene(dt)
  }, error = function(e) {
    error_caught <<- TRUE
  })
  
  if (error_caught) {
    cat("  ✓ PASS (error correctly thrown)\n")
  } else {
    cat("  ✗ FAIL (should have thrown error)\n")
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
test_result <- test_tissue_counting()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
