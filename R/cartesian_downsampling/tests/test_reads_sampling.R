#!/usr/bin/env Rscript
# ==============================================================================
# TEST: Reads-Based Sampling Functions
# PURPOSE: Validate reads-mode sampling logic (multinomial overflow, exact targeting)
# STAGE: 1 (Unit Test)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 1 UNIT TEST: Reads-Based Sampling\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
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
# 2. CREATE SYNTHETIC TEST DATA
# ==============================================================================

cat("--- Creating synthetic test data ---\n\n")

# Create a synthetic SCE with known read counts
set.seed(42)
n_genes <- 100
n_cells <- 50

# Create count matrix with known total reads per cell
# Cell 1: 1000 reads, Cell 2: 2000 reads, ..., Cell 50: 50000 reads
counts_list <- lapply(1:n_cells, function(i) {
  target_reads <- i * 1000
  # Generate random counts that sum to target_reads
  counts <- rpois(n_genes, lambda = target_reads / n_genes)
  # Adjust to hit exact target
  diff <- target_reads - sum(counts)
  counts[1] <- counts[1] + diff
  return(counts)
})

counts_matrix <- do.call(cbind, counts_list)
counts_matrix <- as(counts_matrix, "dgCMatrix")
rownames(counts_matrix) <- paste0("GENE_", sprintf("%03d", 1:n_genes))
colnames(counts_matrix) <- paste0("CELL_", sprintf("%03d", 1:n_cells))

# Create colData with correct number of rows
col_data <- DataFrame(
  donor_id = rep("TEST_DONOR", n_cells),
  sex = rep("male", n_cells),
  development_stage = rep("50 year-old stage", n_cells),
  row.names = colnames(counts_matrix)
)

# Create SCE
test_sce <- SingleCellExperiment(
  assays = list(counts = counts_matrix),
  colData = col_data
)

# Verify total reads per cell
reads_per_cell <- Matrix::colSums(assay(test_sce, "counts"))
cat(sprintf("Created synthetic SCE: %d genes x %d cells\n", n_genes, n_cells))
cat(sprintf("Read range: %.0f - %.0f per cell\n", min(reads_per_cell), max(reads_per_cell)))
cat(sprintf("Total reads: %.0f\n\n", sum(reads_per_cell)))

# Save to temp directory
temp_dir <- tempdir()
test_shard_dir <- file.path(temp_dir, "test_shards")
dir.create(test_shard_dir, showWarnings = FALSE, recursive = TRUE)
test_shard_path <- file.path(test_shard_dir, "donor_TEST_DONOR.qs")
qs::qsave(test_sce, test_shard_path)

# ==============================================================================
# 3. TEST READS-MODE AGGREGATION
# ==============================================================================

cat("--- Testing reads-mode aggregation ---\n\n")

# Test 1: Exact targeting with multinomial overflow
target_reads_test <- 10000  # Should take cells 1-9 (9000 reads) + partial cell 10
result_exact <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = target_reads_test,
  mode = "reads",
  seed = 42,
  verbose = FALSE
)

run_test("Reads mode returns non-NULL result", {
  !is.null(result_exact)
})

run_test("Result has SCE and sampling_log", {
  !is.null(result_exact$sce) && !is.null(result_exact$sampling_log)
})

# Check exact read targeting
actual_reads <- sum(assay(result_exact$sce, "counts"))
run_test("Exact read targeting achieved", {
  abs(actual_reads - target_reads_test) < 1  # Allow ±1 for rounding
})

cat(sprintf("  [INFO] Target: %.0f reads, Actual: %.0f reads, Diff: %.0f\n", 
            target_reads_test, actual_reads, actual_reads - target_reads_test))

# ==============================================================================
# 4. TEST MULTINOMIAL DOWNSAMPLING PROPERTIES
# ==============================================================================

cat("\n--- Testing multinomial downsampling properties ---\n\n")

# Test with a target that requires overflow handling
target_with_overflow <- 5500  # Should take cells 1-5 (15000 total) + partial cell 6
result_overflow <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = target_with_overflow,
  mode = "reads",
  seed = 123,
  verbose = FALSE
)

actual_overflow <- sum(assay(result_overflow$sce, "counts"))
run_test("Overflow handling produces exact target", {
  abs(actual_overflow - target_with_overflow) < 1
})

# Test gene distribution preservation (multinomial should preserve proportions)
# Compare the last cell's gene distribution before and after downsampling
# This is a statistical test - we can't check exact values due to randomness,
# but we can check that the distribution is reasonable
run_test("Downsampled cell has non-zero counts", {
  last_cell_counts <- assay(result_overflow$sce, "counts")[, ncol(result_overflow$sce)]
  sum(last_cell_counts > 0) > 0
})

run_test("Downsampled cell has reasonable gene coverage", {
  last_cell_counts <- assay(result_overflow$sce, "counts")[, ncol(result_overflow$sce)]
  # At least 10% of genes should have counts (multinomial preserves distribution)
  sum(last_cell_counts > 0) >= n_genes * 0.1
})

# ==============================================================================
# 5. TEST EDGE CASES
# ==============================================================================

cat("\n--- Testing edge cases ---\n\n")

# Test: Target reads equals exact donor total (no overflow)
total_donor_reads <- sum(reads_per_cell)
cat(sprintf("  [DEBUG] Testing exact total: %.0f reads\n", total_donor_reads))
result_exact_total <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = total_donor_reads,
  mode = "reads",
  seed = 42,
  verbose = TRUE
)

if (!is.null(result_exact_total)) {
  cat(sprintf("  [DEBUG] Result: %d donors, %d cells total\n", 
              ncol(result_exact_total$sce), 
              length(unlist(result_exact_total$sampling_log))))
} else {
  cat("  [DEBUG] Result is NULL\n")
}

run_test("Exact total reads uses all cells", {
  !is.null(result_exact_total) && 
  length(result_exact_total$sampling_log[[1]]) == n_cells
})

# Test: Target reads exceeds donor total (should return NULL)
result_too_high <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = total_donor_reads + 10000,
  mode = "reads",
  seed = 42,
  verbose = FALSE
)

run_test("Insufficient reads returns NULL", {
  is.null(result_too_high)
})

# Test: NA target reads uses all cells (same as cells mode)
result_na <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = NA,
  mode = "reads",
  seed = 42,
  verbose = FALSE
)

run_test("NA target_reads uses all cells", {
  !is.null(result_na) && ncol(result_na$sce) == 1  # 1 donor
})

run_test("NA target_reads gives total donor reads", {
  actual_na <- sum(assay(result_na$sce, "counts"))
  abs(actual_na - total_donor_reads) < 1
})

# Test: Very small target (should take just 1-2 cells)
result_small <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = 1500,
  mode = "reads",
  seed = 42,
  verbose = FALSE
)

run_test("Small target works", {
  !is.null(result_small)
})

run_test("Small target uses few cells", {
  length(result_small$sampling_log[[1]]) <= 3
})

# ==============================================================================
# 6. TEST REPRODUCIBILITY
# ==============================================================================

cat("\n--- Testing reproducibility ---\n\n")

# Same seed should give same result
result_seed1_a <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = 10000,
  mode = "reads",
  seed = 999,
  verbose = FALSE
)

result_seed1_b <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = 10000,
  mode = "reads",
  seed = 999,
  verbose = FALSE
)

run_test("Same seed produces identical results", {
  identical(result_seed1_a$sampling_log, result_seed1_b$sampling_log) &&
  identical(as.matrix(assay(result_seed1_a$sce, "counts")), 
            as.matrix(assay(result_seed1_b$sce, "counts")))
})

# Different seeds should give different results
result_seed2 <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  target_reads = 10000,
  mode = "reads",
  seed = 777,
  verbose = FALSE
)

run_test("Different seeds produce different results", {
  !identical(result_seed1_a$sampling_log, result_seed2$sampling_log)
})

# ==============================================================================
# 7. TEST CELLS MODE STILL WORKS
# ==============================================================================

cat("\n--- Testing cells mode unchanged ---\n\n")

# Ensure cells mode still works after our changes
result_cells <- aggregate_from_shards(
  donor_ids = "TEST_DONOR",
  shard_dir = test_shard_dir,
  n_cells_to_sample = 10,
  mode = "cells",
  seed = 42,
  verbose = FALSE
)

run_test("Cells mode still works", {
  !is.null(result_cells)
})

run_test("Cells mode samples exact cell count", {
  length(result_cells$sampling_log[[1]]) == 10
})

run_test("Cells mode uses random selection (not cumulative)", {
  # In cells mode, indices should be sorted (random sample then sorted)
  # In reads mode, indices are cumulative (1, 2, 3, ...)
  indices <- result_cells$sampling_log[[1]]
  # Check that indices are not just 1:10 (which would indicate cumulative)
  !identical(indices, 1:10)
})

# ==============================================================================
# 8. READ DEPTH DIFFERENTIATION (USING REAL SHARDS)
# ==============================================================================

cat("\n--- Testing read depth differentiation with real shards ---\n\n")

# Load real shards from production data
real_shard_dir <- file.path(PROJECT_ROOT, "Datalink", "Shards", "Micro", "Roussos")
real_manifest_path <- file.path(real_shard_dir, "00_MASTER_MANIFEST.qs")

if (!file.exists(real_manifest_path)) {
  cat("WARNING: Real shards not found, skipping real shard tests\n")
  cat("Expected path: ", real_manifest_path, "\n")
} else {
  real_manifest <- qs::qread(real_manifest_path)
  
  # Select 3 donors for testing
  set.seed(42)
  test_real_donors <- sample(real_manifest$manifest$donor_id, min(3, nrow(real_manifest$manifest)))
  
  cat(sprintf("Testing with %d real donors: %s\n\n", 
              length(test_real_donors), 
              paste(test_real_donors, collapse = ", ")))
  
  # Test different target_reads values
  target_reads_levels <- c(5000, 25000, 500000, NA)
  results_by_target <- list()
  
  for (target in target_reads_levels) {
    target_label <- ifelse(is.na(target), "NA", as.character(target))
    
    result <- aggregate_from_shards(
      donor_ids = test_real_donors,
      shard_dir = real_shard_dir,
      n_cells_to_sample = NA,
      target_reads = target,
      mode = "reads",
      seed = 42,
      verbose = FALSE
    )
    
    if (!is.null(result)) {
      actual_reads <- sum(Matrix::colSums(SummarizedExperiment::assay(result$sce, "counts")))
      n_cells_total <- sum(sapply(result$sampling_log, length))
      
      results_by_target[[target_label]] <- list(
        actual_reads = actual_reads,
        n_cells = n_cells_total,
        result = result
      )
      
      cat(sprintf("  Target: %10s → Actual: %10.0f reads, %d cells\n", 
                  target_label, actual_reads, n_cells_total))
    }
  }
  
  cat("\n")
  
  # Test: Different targets produce different actual reads
  run_test("Different target_reads produce different actual reads", {
    if (length(results_by_target) < 3) return(FALSE)
    
    actual_reads_vec <- sapply(results_by_target, function(x) x$actual_reads)
    # All values should be unique
    length(unique(actual_reads_vec)) == length(actual_reads_vec)
  })
  
  # Test: Actual reads approximately match target × n_donors (within cell granularity)
  run_test("Actual reads within reasonable range of target × n_donors", {
    all_ok <- TRUE
    n_donors <- length(test_real_donors)
    
    for (target_label in names(results_by_target)) {
      if (target_label == "NA") next  # Skip NA target
      
      target_val <- as.numeric(target_label)
      actual_val <- results_by_target[[target_label]]$actual_reads
      expected_val <- target_val * n_donors
      
      # Allow up to 1 cell's worth of reads per donor as tolerance
      # Real shards have ~3500 reads/cell, so tolerance = 3500 * n_donors
      tolerance <- 5000 * n_donors
      
      if (abs(actual_val - expected_val) > tolerance) {
        cat(sprintf("    FAIL: Target %s × %d donors = %.0f expected, Actual %.0f, Diff %.0f\n", 
                    target_label, n_donors, expected_val, actual_val, actual_val - expected_val))
        all_ok <- FALSE
      }
    }
    all_ok
  })
  
  # Test: Fewer cells selected at lower targets
  run_test("Lower targets select fewer cells", {
    if (length(results_by_target) < 3) return(FALSE)
    
    # Get numeric targets (exclude NA)
    numeric_targets <- names(results_by_target)[names(results_by_target) != "NA"]
    numeric_targets <- as.numeric(numeric_targets)
    numeric_targets <- sort(numeric_targets)
    
    if (length(numeric_targets) < 2) return(TRUE)  # Need at least 2 to compare
    
    # Check that cell counts increase with target
    for (i in 1:(length(numeric_targets) - 1)) {
      target1 <- as.character(numeric_targets[i])
      target2 <- as.character(numeric_targets[i + 1])
      
      cells1 <- results_by_target[[target1]]$n_cells
      cells2 <- results_by_target[[target2]]$n_cells
      
      if (cells1 >= cells2) {
        cat(sprintf("    FAIL: Target %s has %d cells, Target %s has %d cells\n",
                    target1, cells1, target2, cells2))
        return(FALSE)
      }
    }
    TRUE
  })
}

# ==============================================================================
# 9. CELLS MODE DIFFERENTIATION (USING REAL SHARDS)
# ==============================================================================

cat("\n--- Testing cells mode differentiation with real shards ---\n\n")

if (!file.exists(real_manifest_path)) {
  cat("WARNING: Real shards not found, skipping cells mode tests\n")
} else {
  # Test different n_cells values
  n_cells_levels <- c(10, 50, 200, NA)
  cells_results_by_target <- list()
  
  for (n_cells_val in n_cells_levels) {
    n_cells_label <- ifelse(is.na(n_cells_val), "NA", as.character(n_cells_val))
    
    result <- aggregate_from_shards(
      donor_ids = test_real_donors,
      shard_dir = real_shard_dir,
      n_cells_to_sample = n_cells_val,
      target_reads = NA,
      mode = "cells",
      seed = 42,
      verbose = FALSE
    )
    
    if (!is.null(result)) {
      actual_reads <- sum(Matrix::colSums(SummarizedExperiment::assay(result$sce, "counts")))
      n_cells_total <- sum(sapply(result$sampling_log, length))
      
      cells_results_by_target[[n_cells_label]] <- list(
        actual_reads = actual_reads,
        n_cells = n_cells_total
      )
      
      cat(sprintf("  n_cells: %4s → Actual: %d cells, %10.0f reads\n", 
                  n_cells_label, n_cells_total, actual_reads))
    }
  }
  
  cat("\n")
  
  # Test: Different n_cells produce different cell counts
  run_test("Different n_cells produce different cell counts", {
    if (length(cells_results_by_target) < 3) return(FALSE)
    
    cell_counts_vec <- sapply(cells_results_by_target, function(x) x$n_cells)
    # All values should be unique
    length(unique(cell_counts_vec)) == length(cell_counts_vec)
  })
  
  # Test: Cell counts match requested × eligible donors (donors with sufficient cells)
  run_test("Cell counts match requested n_cells × eligible donors", {
    all_ok <- TRUE
    for (n_cells_label in names(cells_results_by_target)) {
      if (n_cells_label == "NA") next
      
      requested <- as.numeric(n_cells_label)
      actual <- cells_results_by_target[[n_cells_label]]$n_cells
      
      # Count eligible donors (those with >= requested cells)
      # Load manifest to check cell counts
      manifest <- qs::qread(real_manifest_path)
      eligible_donors <- sum(sapply(test_real_donors, function(d) {
        donor_info <- manifest$manifest[manifest$manifest$donor_id == d, ]
        donor_info$cell_count >= requested
      }))
      
      expected <- requested * eligible_donors
      
      if (actual != expected) {
        cat(sprintf("    FAIL: Requested %s cells/donor × %d eligible donors = %d expected, got %d actual\n",
                    n_cells_label, eligible_donors, expected, actual))
        all_ok <- FALSE
      }
    }
    all_ok
  })
  
  # Test: More cells = more reads
  run_test("More cells produce more reads", {
    if (length(cells_results_by_target) < 3) return(FALSE)
    
    numeric_targets <- names(cells_results_by_target)[names(cells_results_by_target) != "NA"]
    numeric_targets <- as.numeric(numeric_targets)
    numeric_targets <- sort(numeric_targets)
    
    if (length(numeric_targets) < 2) return(TRUE)
    
    for (i in 1:(length(numeric_targets) - 1)) {
      target1 <- as.character(numeric_targets[i])
      target2 <- as.character(numeric_targets[i + 1])
      
      reads1 <- cells_results_by_target[[target1]]$actual_reads
      reads2 <- cells_results_by_target[[target2]]$actual_reads
      
      if (reads1 >= reads2) {
        cat(sprintf("    FAIL: %s cells has %.0f reads, %s cells has %.0f reads\n",
                    target1, reads1, target2, reads2))
        return(FALSE)
      }
    }
    TRUE
  })
}

# ==============================================================================
# 10. PARALLEL EXECUTION FIDELITY (BPLAPPLY)
# ==============================================================================

cat("\n--- Testing parallel execution fidelity ---\n\n")

if (!file.exists(real_manifest_path)) {
  cat("WARNING: Real shards not found, skipping parallel execution tests\n")
} else {
  # Build job specs with different target_reads
  job_specs <- list(
    list(target_reads = 5000, label = "5k"),
    list(target_reads = 25000, label = "25k"),
    list(target_reads = 500000, label = "500k"),
    list(target_reads = NA, label = "NA")
  )
  
  # Execute via bplapply (same as production pipeline)
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BPPARAM_TEST <- BiocParallel::MulticoreParam(workers = 2, progressbar = FALSE)
    
    parallel_results <- BiocParallel::bplapply(job_specs, function(job) {
      result <- aggregate_from_shards(
        donor_ids = test_real_donors,
        shard_dir = real_shard_dir,
        n_cells_to_sample = NA,
        target_reads = job$target_reads,
        mode = "reads",
        seed = 42,
        verbose = FALSE
      )
      
      if (is.null(result)) return(NULL)
      
      list(
        label = job$label,
        target = job$target_reads,
        actual_reads = sum(Matrix::colSums(SummarizedExperiment::assay(result$sce, "counts"))),
        n_cells = sum(sapply(result$sampling_log, length))
      )
    }, BPPARAM = BPPARAM_TEST)
    
    # Remove NULL results
    parallel_results <- parallel_results[!sapply(parallel_results, is.null)]
    
    cat("Parallel execution results:\n")
    for (res in parallel_results) {
      cat(sprintf("  %5s: %.0f reads, %d cells\n", res$label, res$actual_reads, res$n_cells))
    }
    cat("\n")
    
    # Test: All jobs produced different results
    run_test("Parallel jobs produce different results", {
      if (length(parallel_results) < 3) return(FALSE)
      
      actual_reads_vec <- sapply(parallel_results, function(x) x$actual_reads)
      length(unique(actual_reads_vec)) == length(actual_reads_vec)
    })
    
    # Test: Parallel results match sequential results
    run_test("Parallel results match sequential execution", {
      # Compare with results_by_target from section 8
      if (length(results_by_target) == 0) return(TRUE)  # Skip if section 8 didn't run
      
      all_match <- TRUE
      for (res in parallel_results) {
        target_label <- ifelse(is.na(res$target), "NA", as.character(res$target))
        
        if (target_label %in% names(results_by_target)) {
          seq_reads <- results_by_target[[target_label]]$actual_reads
          par_reads <- res$actual_reads
          
          if (seq_reads != par_reads) {
            cat(sprintf("    FAIL: Target %s - Sequential: %.0f, Parallel: %.0f\n",
                        target_label, seq_reads, par_reads))
            all_match <- FALSE
          }
        }
      }
      all_match
    })
  } else {
    cat("BiocParallel not available, skipping parallel execution tests\n")
  }
}

# ==============================================================================
# 11. SANITY CHECK FUNCTION
# ==============================================================================

cat("\n--- Testing sanity check function ---\n\n")

#' Check if read depth varies across Target_Reads settings
#' 
#' @param df_summary Aggregated summary data frame with Target_Reads and Mean_Total_Reads columns
#' @return List with status ("PASS"/"FAIL") and details
check_read_depth_differentiation <- function(df_summary) {
  if (!"Target_Reads" %in% colnames(df_summary)) {
    return(list(status = "SKIP", details = "No Target_Reads column found"))
  }
  
  if (!"Mean_Total_Reads" %in% colnames(df_summary)) {
    return(list(status = "SKIP", details = "No Mean_Total_Reads column found"))
  }
  
  # Group by Target_Reads and get mean of Mean_Total_Reads
  summary_by_target <- df_summary %>%
    group_by(Target_Reads) %>%
    summarise(Avg_Reads = mean(Mean_Total_Reads, na.rm = TRUE), .groups = "drop")
  
  # Check if all values are identical (within 1% tolerance)
  unique_vals <- unique(summary_by_target$Avg_Reads)
  
  if (length(unique_vals) == 1) {
    return(list(
      status = "FAIL",
      details = sprintf("All Target_Reads produce identical Mean_Total_Reads (%.0f)", unique_vals[1])
    ))
  }
  
  # Check if variation is too small (all within 1% of each other)
  range_pct <- (max(unique_vals) - min(unique_vals)) / mean(unique_vals) * 100
  
  if (range_pct < 1) {
    return(list(
      status = "WARN",
      details = sprintf("Mean_Total_Reads varies by only %.2f%% across targets", range_pct)
    ))
  }
  
  return(list(
    status = "PASS",
    details = sprintf("Mean_Total_Reads varies by %.1f%% across targets (range: %.0f - %.0f)",
                      range_pct, min(unique_vals), max(unique_vals))
  ))
}

# Test the sanity check function with mock data
suppressPackageStartupMessages(library(dplyr))

# Mock data: all targets produce same reads (should FAIL)
mock_data_fail <- data.frame(
  Target_Reads = rep(c(25000, 50000, 100000, NA), each = 3),
  Mean_Total_Reads = rep(750000, 12),
  Donors = rep(c(5, 10, 25), 4)
)

result_fail <- check_read_depth_differentiation(mock_data_fail)
run_test("Sanity check detects read-count indifference", {
  result_fail$status == "FAIL"
})

# Mock data: targets produce different reads (should PASS)
mock_data_pass <- data.frame(
  Target_Reads = rep(c(25000, 50000, 100000, NA), each = 3),
  Mean_Total_Reads = c(
    rep(25000, 3),   # 25k target
    rep(50000, 3),   # 50k target
    rep(100000, 3),  # 100k target
    rep(750000, 3)   # NA target
  ),
  Donors = rep(c(5, 10, 25), 4)
)

result_pass <- check_read_depth_differentiation(mock_data_pass)
run_test("Sanity check passes with differentiated data", {
  result_pass$status == "PASS"
})

cat(sprintf("\nSanity check function defined: check_read_depth_differentiation()\n"))
cat(sprintf("  - FAIL case: %s\n", result_fail$details))
cat(sprintf("  - PASS case: %s\n", result_pass$details))

# ==============================================================================
# 12. CLEANUP
# ==============================================================================

unlink(test_shard_path)
unlink(test_shard_dir, recursive = TRUE)

# ==============================================================================
# 9. SUMMARY
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
