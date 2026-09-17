#!/usr/bin/env Rscript
# ==============================================================================
# TEST: DESeq2 Worker Function
# PURPOSE: Validate run_deseq_worker() with real data
# STAGE: 1 (Unit Test)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 1 UNIT TEST: DESeq2 Worker Function\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(Matrix)
})

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",  # From project root
    "../config.R",                         # From tests/ directory
    "config.R"                             # From script directory
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
    stop("Could not find config.R. Ensure working directory is project root or tests/ directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

# Test counter
tests_passed <- 0
tests_failed <- 0

run_test <- function(name, expr) {
  cat(sprintf("[TEST] %s... ", name))
  result <- tryCatch({
    if (expr) {
      cat("PASSED\n")
      tests_passed <<- tests_passed + 1
      TRUE
    } else {
      cat("FAILED\n")
      tests_failed <<- tests_failed + 1
      FALSE
    }
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    tests_failed <<- tests_failed + 1
    FALSE
  })
  return(result)
}

# ==============================================================================
# 2. PREPARE TEST DATA
# ==============================================================================

cat("--- Preparing test data ---\n\n")

# Use Micro cell type
test_cell_type <- "Micro"
paths <- get_paths(test_cell_type)

if (!file.exists(paths$manifest_path)) {
  cat("\nERROR: Manifest not found. Cannot run DESeq2 tests.\n")
  quit(status = 1)
}

# Load donor metadata
donor_meta <- load_donor_metadata(paths$manifest_path, paths$shard_dir)

if (is.null(donor_meta) || nrow(donor_meta) < 6) {
  cat("\nERROR: Need at least 6 donors for DESeq2 tests.\n")
  quit(status = 1)
}

# Get 3 males and 3 females
males <- donor_meta$donor_id[donor_meta$sex == MALE_VAL]
females <- donor_meta$donor_id[donor_meta$sex == FEMALE_VAL]

if (length(males) < 3 || length(females) < 3) {
  cat("\nERROR: Need at least 3 donors of each sex.\n")
  quit(status = 1)
}

test_donors <- c(head(males, 3), head(females, 3))
cat(sprintf("Test donors: %d males, %d females\n\n", 3, 3))

# Create temp directory for test outputs
test_dir <- tempdir()
test_results_path <- file.path(test_dir, "test_deseq_results.qs")
test_log_path <- file.path(test_dir, "test_deseq_log.qs")

# ==============================================================================
# 3. AGGREGATE TEST DATA
# ==============================================================================

cat("--- Aggregating test data ---\n\n")

agg_result <- aggregate_from_shards(
  donor_ids = test_donors,
  shard_dir = paths$shard_dir,
  n_cells_to_sample = 100,  # Use 100 cells per donor for speed
  seed = 42,
  verbose = TRUE
)

if (is.null(agg_result)) {
  cat("\nERROR: Failed to aggregate test data.\n")
  quit(status = 1)
}

cat(sprintf("Aggregated SCE: %d genes x %d samples\n\n", 
            nrow(agg_result$sce), ncol(agg_result$sce)))

# ==============================================================================
# 4. TEST run_deseq_worker()
# ==============================================================================

cat("--- Testing run_deseq_worker() ---\n\n")

# Run DESeq2
run_test("DESeq2 worker completes successfully", {
  status <- run_deseq_worker(
    agg_result = agg_result,
    sex_col = SEX_COL,
    age_col = AGE_COL,
    cat_vars = CATEGORICAL_VARS,
    save_path_results = test_results_path,
    save_path_log = test_log_path,
    min_gene_count = MIN_GENE_COUNT
  )
  status == "SUCCESS"
})

run_test("Results file was created", {
  file.exists(test_results_path)
})

run_test("Log file was created", {
  file.exists(test_log_path)
})

# ==============================================================================
# 5. VALIDATE OUTPUT STRUCTURE
# ==============================================================================

cat("\n--- Validating output structure ---\n\n")

# Load results
deseq_results <- qs::qread(test_results_path)

run_test("Results is a DESeqResults object", {
  inherits(deseq_results, "DESeqResults")
})

run_test("Results has log2FoldChange column", {
  "log2FoldChange" %in% colnames(deseq_results)
})
run_test("Results has padj column", {
  "padj" %in% colnames(deseq_results)
})

run_test("Results has baseMean column", {
  "baseMean" %in% colnames(deseq_results)
})

run_test("Results has gene names as rownames", {
  length(rownames(deseq_results)) > 0 && all(nchar(rownames(deseq_results)) > 0)
})

# ==============================================================================
# 6. VALIDATE RESULT VALUES
# ==============================================================================

cat("\n--- Validating result values ---\n\n")

run_test("log2FoldChange values are finite (no Inf)", {
  lfc <- deseq_results$log2FoldChange
  all(is.na(lfc) | is.finite(lfc))
})

run_test("padj values are between 0 and 1 (or NA)", {
  padj <- deseq_results$padj
  all(is.na(padj) | (padj >= 0 & padj <= 1))
})

run_test("baseMean values are non-negative", {
  bm <- deseq_results$baseMean
  all(is.na(bm) | bm >= 0)
})

run_test("At least some genes have non-NA padj", {
  sum(!is.na(deseq_results$padj)) > 0
})

# Check for reasonable number of significant genes
n_sig <- sum(deseq_results$padj < 0.05, na.rm = TRUE)
cat(sprintf("  [INFO] Significant genes (padj < 0.05): %d\n", n_sig))

# ==============================================================================
# 7. VALIDATE SAMPLING LOG
# ==============================================================================

cat("\n--- Validating sampling log ---\n\n")

sampling_log <- qs::qread(test_log_path)

run_test("Sampling log is a list", {
  is.list(sampling_log)
})

run_test("Sampling log has entry for each aggregated donor", {
  # Note: Some donors may be skipped if they don't have enough cells
  # The log should match the number of samples in the SCE
  length(sampling_log) == ncol(agg_result$sce)
})

run_test("Each log entry contains cell indices", {
  all(sapply(sampling_log, function(x) is.numeric(x) && length(x) > 0))
})

run_test("Cell indices are positive integers", {
  all(sapply(sampling_log, function(x) all(x > 0) && all(x == floor(x))))
})

# ==============================================================================
# 8. TEST EDGE CASES
# ==============================================================================

cat("\n--- Testing edge cases ---\n\n")

# Test with minimal data (should still work)
run_test("Works with minimum viable sample size", {
  # Use just 2 males and 2 females
  minimal_donors <- c(head(males, 2), head(females, 2))
  minimal_agg <- aggregate_from_shards(minimal_donors, paths$shard_dir, 50, 42, FALSE)
  
  if (is.null(minimal_agg)) return(FALSE)
  
  minimal_results_path <- file.path(test_dir, "minimal_results.qs")
  minimal_log_path <- file.path(test_dir, "minimal_log.qs")
  
  status <- tryCatch({
    run_deseq_worker(minimal_agg, SEX_COL, AGE_COL, CATEGORICAL_VARS,
                     minimal_results_path, minimal_log_path, MIN_GENE_COUNT)
  }, error = function(e) e$message)
  
  status == "SUCCESS"
})

# ==============================================================================
# 9. CLEANUP
# ==============================================================================

# Remove test files
unlink(test_results_path)
unlink(test_log_path)
unlink(file.path(test_dir, "minimal_results.qs"))
unlink(file.path(test_dir, "minimal_log.qs"))

# ==============================================================================
# 10. SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed\n", tests_passed, tests_failed))
cat("================================================================================\n")
cat("\n")

if (tests_failed > 0) {
  cat("  STATUS: FAILED\n")
  quit(status = 1)
} else {
  cat("  STATUS: PASSED\n")
  quit(status = 0)
}
