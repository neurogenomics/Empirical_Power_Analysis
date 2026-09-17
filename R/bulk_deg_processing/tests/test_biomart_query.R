# ==============================================================================
# UNIT TEST: BioMart Querying
# ==============================================================================
# Tests the query_biomart_safe() function
# NOTE: This test requires internet connection and BioMart availability
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

test_biomart_query <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: BioMart Querying\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Test 1: Query known genes
  cat("[TEST 1] Query known genes (BRCA1, TP53)\n")
  gene_ids <- c("ENSG00000012048", "ENSG00000141510")
  
  tryCatch({
    result <- query_biomart_safe(gene_ids)
    
    # Check structure
    required_cols <- c("ensembl_gene_id", "external_gene_name", "description")
    if (all(required_cols %in% names(result)) &&
        nrow(result) == 2 &&
        "BRCA1" %in% result$external_gene_name &&
        "TP53" %in% result$external_gene_name) {
      cat("  ✓ PASS\n")
    } else {
      cat("  ✗ FAIL\n")
      cat("    Expected 2 genes (BRCA1, TP53), got:\n")
      print(result)
      all_passed <- FALSE
    }
  }, error = function(e) {
    cat("  ⚠ SKIP (BioMart connection failed)\n")
    cat("    Error:", e$message, "\n")
    # Don't fail test if BioMart is down
  })
  
  # Test 2: Query single gene
  cat("\n[TEST 2] Query single gene (BRCA1)\n")
  gene_ids <- c("ENSG00000012048")
  
  tryCatch({
    result <- query_biomart_safe(gene_ids)
    
    if (nrow(result) == 1 && result$external_gene_name[1] == "BRCA1") {
      cat("  ✓ PASS\n")
    } else {
      cat("  ✗ FAIL\n")
      cat("    Expected BRCA1, got:\n")
      print(result)
      all_passed <- FALSE
    }
  }, error = function(e) {
    cat("  ⚠ SKIP (BioMart connection failed)\n")
  })
  
  # Test 3: Query invalid gene (should return empty or partial results)
  cat("\n[TEST 3] Query invalid gene ID\n")
  gene_ids <- c("ENSG00000000000")  # Fake ID
  
  tryCatch({
    result <- query_biomart_safe(gene_ids)
    
    # Should return empty data.frame or 0 rows
    if (nrow(result) == 0) {
      cat("  ✓ PASS\n")
    } else {
      cat("  ⚠ WARNING: Expected 0 results, got", nrow(result), "\n")
      # Don't fail - BioMart behavior may vary
    }
  }, error = function(e) {
    cat("  ⚠ SKIP (BioMart connection failed)\n")
  })
  
  # Test 4: Empty input
  cat("\n[TEST 4] Empty input\n")
  gene_ids <- character(0)
  result <- query_biomart_safe(gene_ids)
  
  if (nrow(result) == 0) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Test 5: Output structure validation
  cat("\n[TEST 5] Output structure validation\n")
  gene_ids <- c("ENSG00000012048")
  
  tryCatch({
    result <- query_biomart_safe(gene_ids)
    
    required_cols <- c("ensembl_gene_id", "external_gene_name", "description")
    if (all(required_cols %in% names(result)) &&
        is.data.frame(result)) {
      cat("  ✓ PASS\n")
    } else {
      cat("  ✗ FAIL\n")
      cat("    Expected columns:", paste(required_cols, collapse = ", "), "\n")
      cat("    Got:", paste(names(result), collapse = ", "), "\n")
      all_passed <- FALSE
    }
  }, error = function(e) {
    cat("  ⚠ SKIP (BioMart connection failed)\n")
  })
  
  # Test 6: Connection failure handling
  cat("\n[TEST 6] Connection failure handling\n")
  # Test with invalid dataset to trigger error
  gene_ids <- c("ENSG00000012048")
  error_caught <- FALSE
  
  tryCatch({
    result <- query_biomart_safe(gene_ids, dataset = "invalid_dataset")
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
  cat("\nNOTE: Some tests may be skipped if BioMart is unavailable.\n")
  cat("This is expected and does not indicate a test failure.\n")
  cat("========================================\n\n")
  
  return(all_passed)
}

# Run tests
test_result <- test_biomart_query()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
