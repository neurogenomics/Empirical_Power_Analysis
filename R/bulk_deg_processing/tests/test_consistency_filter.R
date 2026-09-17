# ==============================================================================
# UNIT TEST: Consistency Filtering
# ==============================================================================
# Tests the filter_consistent_genes() function
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

test_consistency_filter <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: Consistency Filtering\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Test 1: Gene meeting 90% threshold (3/3 tissues)
  cat("[TEST 1] Gene meeting 90% threshold\n")
  dt <- data.table(
    gene = rep("ENSG00000123456.1", 3),
    tissue = c("Brain", "Heart", "Liver"),
    effsize = c(1.5, 1.6, 1.4),
    effsize_se = c(0.1, 0.1, 0.1),
    lfsr = c(0.01, 0.02, 0.01)
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  if (nrow(result) == 1 && result$gene[1] == "ENSG00000123456") {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 1 gene (ENSG00000123456), got", nrow(result), "\n")
    all_passed <- FALSE
  }
  
  # Test 2: Gene NOT meeting threshold (1/3 tissues significant)
  cat("\n[TEST 2] Gene NOT meeting threshold\n")
  dt <- data.table(
    gene = rep("ENSG00000123456.1", 3),
    tissue = c("Brain", "Heart", "Liver"),
    effsize = c(1.5, 1.6, 1.4),
    effsize_se = c(0.1, 0.1, 0.1),
    lfsr = c(0.01, 0.10, 0.10)  # Only 1/3 significant
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  if (nrow(result) == 0) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 0 genes, got", nrow(result), "\n")
    all_passed <- FALSE
  }
  
  # Test 3: Exactly at threshold (90% = 9/10 tissues)
  cat("\n[TEST 3] Exactly at threshold (9/10 tissues)\n")
  dt <- data.table(
    gene = rep("ENSG00000123456.1", 10),
    tissue = paste0("Tissue", 1:10),
    effsize = rep(1.5, 10),
    effsize_se = rep(0.1, 10),
    lfsr = c(rep(0.01, 9), 0.10)  # 9/10 significant
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  if (nrow(result) == 1) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 1 gene, got", nrow(result), "\n")
    all_passed <- FALSE
  }
  
  # Test 4: Averaging across tissues
  cat("\n[TEST 4] Averaging across tissues\n")
  dt <- data.table(
    gene = rep("ENSG00000123456.1", 3),
    tissue = c("Brain", "Heart", "Liver"),
    effsize = c(1.0, 2.0, 3.0),
    effsize_se = c(0.1, 0.2, 0.3),
    lfsr = c(0.01, 0.02, 0.03)
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  expected_lfc <- mean(c(1.0, 2.0, 3.0))
  expected_se <- mean(c(0.1, 0.2, 0.3))
  expected_padj <- mean(c(0.01, 0.02, 0.03))
  
  if (nrow(result) == 1 &&
      abs(result$log2FoldChange[1] - expected_lfc) < 0.001 &&
      abs(result$lfcSE[1] - expected_se) < 0.001 &&
      abs(result$padj[1] - expected_padj) < 0.001) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected LFC:", expected_lfc, "Got:", result$log2FoldChange[1], "\n")
    cat("    Expected SE:", expected_se, "Got:", result$lfcSE[1], "\n")
    cat("    Expected padj:", expected_padj, "Got:", result$padj[1], "\n")
    all_passed <- FALSE
  }
  
  # Test 5: Gene ID cleaning (version removal)
  cat("\n[TEST 5] Gene ID cleaning\n")
  dt <- data.table(
    gene = rep("ENSG00000123456.12", 3),
    tissue = c("Brain", "Heart", "Liver"),
    effsize = c(1.5, 1.6, 1.4),
    effsize_se = c(0.1, 0.1, 0.1),
    lfsr = c(0.01, 0.02, 0.01)
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  if (nrow(result) == 1 && result$gene[1] == "ENSG00000123456") {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected gene ID: ENSG00000123456, got:", result$gene[1], "\n")
    all_passed <- FALSE
  }
  
  # Test 6: Empty input
  cat("\n[TEST 6] Empty input\n")
  dt <- data.table(
    gene = character(0),
    tissue = character(0),
    effsize = numeric(0),
    effsize_se = numeric(0),
    lfsr = numeric(0)
  )
  result <- filter_consistent_genes(dt, p_threshold = 0.05, consistency_threshold = 0.90)
  
  if (nrow(result) == 0) {
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
test_result <- test_consistency_filter()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
