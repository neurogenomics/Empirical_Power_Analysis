# ==============================================================================
# UNIT TEST: Chromosome Annotation
# ==============================================================================
# Tests the annotate_chromosomes() function
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

test_chromosome_annotation <- function() {
  cat("\n========================================\n")
  cat("UNIT TEST: Chromosome Annotation\n")
  cat("========================================\n\n")
  
  all_passed <- TRUE
  
  # Test 1: Known genes on autosomes
  cat("[TEST 1] Known genes on autosomes\n")
  # BRCA1 (chr17), TP53 (chr17), APOE (chr19)
  gene_ids <- c("ENSG00000012048", "ENSG00000141510", "ENSG00000130203")
  result <- annotate_chromosomes(gene_ids)
  
  # Check structure
  if (all(c("gene", "chromosome", "is_sex_chr") %in% names(result)) &&
      nrow(result) == 3 &&
      all(!result$is_sex_chr)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 3 autosomal genes, got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 2: Known genes on sex chromosomes
  cat("\n[TEST 2] Known genes on sex chromosomes\n")
  # XIST (chrX), SRY (chrY)
  gene_ids <- c("ENSG00000229807", "ENSG00000184895")
  result <- annotate_chromosomes(gene_ids)
  
  # Check that both are marked as sex chromosomes
  if (nrow(result) == 2 && all(result$is_sex_chr)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 2 sex chromosome genes, got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 3: Mixed chromosomes
  cat("\n[TEST 3] Mixed chromosomes (autosome + sex)\n")
  # BRCA1 (chr17), XIST (chrX)
  gene_ids <- c("ENSG00000012048", "ENSG00000229807")
  result <- annotate_chromosomes(gene_ids)
  
  if (nrow(result) == 2 &&
      sum(result$is_sex_chr) == 1 &&
      sum(!result$is_sex_chr) == 1) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected 1 autosome + 1 sex chr, got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 4: Invalid/unmapped gene ID
  cat("\n[TEST 4] Invalid gene ID (should return NA)\n")
  gene_ids <- c("ENSG00000000000")  # Fake ID
  result <- annotate_chromosomes(gene_ids)
  
  if (nrow(result) == 1 && is.na(result$chromosome[1])) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected NA chromosome, got:\n")
    print(result)
    all_passed <- FALSE
  }
  
  # Test 5: Empty input
  cat("\n[TEST 5] Empty input\n")
  gene_ids <- character(0)
  result <- annotate_chromosomes(gene_ids)
  
  if (nrow(result) == 0) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    all_passed <- FALSE
  }
  
  # Test 6: Output structure
  cat("\n[TEST 6] Output structure validation\n")
  gene_ids <- c("ENSG00000012048")
  result <- annotate_chromosomes(gene_ids)
  
  required_cols <- c("gene", "chromosome", "is_sex_chr")
  if (all(required_cols %in% names(result)) &&
      is.character(result$gene) &&
      is.character(result$chromosome) &&
      is.logical(result$is_sex_chr)) {
    cat("  ✓ PASS\n")
  } else {
    cat("  ✗ FAIL\n")
    cat("    Expected columns:", paste(required_cols, collapse = ", "), "\n")
    cat("    Got:", paste(names(result), collapse = ", "), "\n")
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
test_result <- test_chromosome_annotation()

# Exit with appropriate code
if (test_result) {
  quit(status = 0)
} else {
  quit(status = 1)
}
