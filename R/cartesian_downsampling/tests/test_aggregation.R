#!/usr/bin/env Rscript
# ==============================================================================
# TEST: Shard Aggregation Functions
# PURPOSE: Validate aggregate_from_shards() with real data
# STAGE: 1 (Unit Test)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 1 UNIT TEST: Shard Aggregation Functions\n")
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
# 2. FIND TEST DATA
# ==============================================================================

cat("--- Locating test data ---\n\n")

# Use Micro cell type (typically smallest)
test_cell_type <- "Micro"
paths <- get_paths(test_cell_type)

cat(sprintf("Shard directory: %s\n", paths$shard_dir))
cat(sprintf("Manifest path: %s\n", paths$manifest_path))

if (!file.exists(paths$manifest_path)) {
  cat("\nERROR: Manifest not found. Cannot run aggregation tests.\n")
  cat("Expected path: ", paths$manifest_path, "\n")
  quit(status = 1)
}

# Load manifest and get first 3 donors
manifest <- qs::qread(paths$manifest_path)$manifest
test_donors <- head(manifest$donor_id, 3)

cat(sprintf("Test donors: %s\n\n", paste(test_donors, collapse = ", ")))

# ==============================================================================
# 3. TEST aggregate_from_shards() - BASIC FUNCTIONALITY
# ==============================================================================

cat("--- Testing aggregate_from_shards() ---\n\n")

# Test with all cells (NA)
run_test("Aggregates 3 donors with all cells", {
  result <- aggregate_from_shards(
    donor_ids = test_donors,
    shard_dir = paths$shard_dir,
    n_cells_to_sample = NA,
    seed = 42,
    verbose = FALSE
  )
  !is.null(result) && !is.null(result$sce) && !is.null(result$sampling_log)
})

# Test output structure
run_test("Returns SingleCellExperiment object", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  inherits(result$sce, "SingleCellExperiment")
})

run_test("SCE has correct number of samples (donors)", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  ncol(result$sce) == length(test_donors)
})

run_test("SCE has counts assay", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  "counts" %in% assayNames(result$sce)
})

run_test("Counts are non-negative integers", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  counts <- assay(result$sce, "counts")
  all(counts >= 0) && all(counts == floor(counts))
})

run_test("Sampling log has entry for each donor", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  length(result$sampling_log) == length(test_donors)
})

# ==============================================================================
# 4. TEST CELL SAMPLING
# ==============================================================================

cat("\n--- Testing cell sampling ---\n\n")

# Get cell count for first donor
first_shard <- file.path(paths$shard_dir, paste0("donor_", gsub("[^[:alnum:]]", "_", test_donors[1]), ".qs"))
first_sce <- qs::qread(first_shard)
available_cells <- ncol(first_sce)

cat(sprintf("First donor has %d cells\n\n", available_cells))

# Test with specific cell count (if donor has enough)
if (available_cells >= 50) {
  run_test("Samples exactly 50 cells when requested", {
    result <- aggregate_from_shards(
      donor_ids = test_donors[1],
      shard_dir = paths$shard_dir,
      n_cells_to_sample = 50,
      seed = 42,
      verbose = FALSE
    )
    length(result$sampling_log[[1]]) == 50
  })
  
  run_test("Different seeds produce different samples", {
    result1 <- aggregate_from_shards(test_donors[1], paths$shard_dir, 50, 42, FALSE)
    result2 <- aggregate_from_shards(test_donors[1], paths$shard_dir, 50, 123, FALSE)
    !identical(result1$sampling_log[[1]], result2$sampling_log[[1]])
  })
  
  run_test("Same seed produces identical samples", {
    result1 <- aggregate_from_shards(test_donors[1], paths$shard_dir, 50, 42, FALSE)
    result2 <- aggregate_from_shards(test_donors[1], paths$shard_dir, 50, 42, FALSE)
    identical(result1$sampling_log[[1]], result2$sampling_log[[1]])
  })
}

# Test with more cells than available
run_test("Returns NULL when requesting more cells than available", {
  result <- aggregate_from_shards(
    donor_ids = test_donors[1],
    shard_dir = paths$shard_dir,
    n_cells_to_sample = available_cells + 1000,
    seed = 42,
    verbose = FALSE
  )
  is.null(result)
})

# ==============================================================================
# 5. TEST PSEUDOBULK AGGREGATION
# ==============================================================================

cat("\n--- Testing pseudobulk aggregation ---\n\n")

run_test("Pseudobulk counts are row sums of cell counts", {
  # Load raw shard
  raw_sce <- qs::qread(first_shard)
  expected_sums <- Matrix::rowSums(assay(raw_sce, "counts"))
  
  # Get aggregated result
  result <- aggregate_from_shards(test_donors[1], paths$shard_dir, NA, 42, FALSE)
  actual_sums <- assay(result$sce, "counts")[, 1]
  
  # Compare (allowing for gene name ordering)
  common_genes <- intersect(names(expected_sums), names(actual_sums))
  all(expected_sums[common_genes] == actual_sums[common_genes])
})

run_test("Gene names are preserved", {
  result <- aggregate_from_shards(test_donors[1], paths$shard_dir, NA, 42, FALSE)
  length(rownames(result$sce)) > 0 && all(nchar(rownames(result$sce)) > 0)
})

run_test("Donor IDs are in column names", {
  result <- aggregate_from_shards(test_donors, paths$shard_dir, NA, 42, FALSE)
  all(test_donors %in% colnames(result$sce))
})

# ==============================================================================
# 6. TEST ERROR HANDLING
# ==============================================================================

cat("\n--- Testing error handling ---\n\n")

run_test("Returns NULL for non-existent donor", {
  result <- aggregate_from_shards(
    donor_ids = c("FAKE_DONOR_12345"),
    shard_dir = paths$shard_dir,
    n_cells_to_sample = NA,
    seed = 42,
    verbose = FALSE
  )
  is.null(result)
})

run_test("Handles mix of valid and invalid donors gracefully", {
  result <- aggregate_from_shards(
    donor_ids = c(test_donors[1], "FAKE_DONOR_12345"),
    shard_dir = paths$shard_dir,
    n_cells_to_sample = NA,
    seed = 42,
    verbose = FALSE
  )
  # Should return result with just the valid donor
  !is.null(result) && ncol(result$sce) == 1
})

# ==============================================================================
# 7. SUMMARY
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
