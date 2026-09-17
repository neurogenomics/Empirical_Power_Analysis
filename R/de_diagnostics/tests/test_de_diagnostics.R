# ==============================================================================
# TESTS: DE Diagnostics Pipeline - Unit Tests
# ==============================================================================
# Tests utility functions without requiring full data.
# Run from project root: Rscript R/de_diagnostics/tests/test_de_diagnostics.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Source config and utils
source("R/de_diagnostics/config.R")
source("R/de_diagnostics/utils.R")

# Track results
tests_passed <- 0
tests_failed <- 0

test <- function(name, expr) {
  result <- tryCatch({
    eval(expr)
    TRUE
  }, error = function(e) {
    message(sprintf("  [FAIL] %s: %s", name, e$message))
    FALSE
  })

  if (result) {
    message(sprintf("  [PASS] %s", name))
    tests_passed <<- tests_passed + 1
  } else {
    tests_failed <<- tests_failed + 1
  }
}

assert <- function(condition, msg = "Assertion failed") {
  if (!condition) stop(msg)
}

message("\n=== DE Diagnostics Unit Tests ===\n")

# ==============================================================================
# 1. GENE ID CLEANING
# ==============================================================================

message("--- Gene ID Cleaning ---")

test("clean_gene_ids strips version suffix", {
  result <- clean_gene_ids(c("ENSG00000229807.11", "ENSG00000000003.15"))
  assert(all(result == c("ENSG00000229807", "ENSG00000000003")),
         "Version suffixes not stripped correctly")
})

test("clean_gene_ids handles no-version IDs", {
  result <- clean_gene_ids(c("ENSG00000229807", "ENSG00000000003"))
  assert(all(result == c("ENSG00000229807", "ENSG00000000003")),
         "IDs without versions should be unchanged")
})

test("clean_gene_ids handles empty input", {
  result <- clean_gene_ids(character(0))
  assert(length(result) == 0, "Empty input should return empty output")
})

# ==============================================================================
# 2. AGE PARSING
# ==============================================================================

message("\n--- Age Parsing ---")

test("parse_age_to_years handles year format", {
  result <- parse_age_to_years("80 year-old stage")
  assert(result == 80, sprintf("Expected 80, got %s", result))
})

test("parse_age_to_years handles month format", {
  result <- parse_age_to_years("6 month-old stage")
  assert(abs(result - 0.5) < 0.01, sprintf("Expected ~0.5, got %s", result))
})

test("parse_age_to_years handles unknown", {
  result <- parse_age_to_years("unknown")
  assert(is.na(result), "Expected NA for 'unknown'")
})

test("parse_age_to_years handles vector input", {
  result <- parse_age_to_years(c("80 year-old stage", "unknown", "65 year-old stage"))
  assert(result[1] == 80, "First element should be 80")
  assert(is.na(result[2]), "Second element should be NA")
  assert(result[3] == 65, "Third element should be 65")
})

# ==============================================================================
# 3. STATS LABEL
# ==============================================================================

message("\n--- Stats Label ---")

test("get_stats_label computes medians", {
  df <- data.frame(
    val = c(1, 2, 3, 4, 5, 6),
    group = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  result <- get_stats_label(df, "val", "group")
  assert(grepl("A: 2.000", result), "Should contain median for A")
  assert(grepl("B: 5.000", result), "Should contain median for B")
})

# ==============================================================================
# 4. DENSITY PLOT FUNCTIONS
# ==============================================================================

message("\n--- Density Plot Functions ---")

test("plot_density_aggregated returns ggplot with explicit limits", {
  df <- data.frame(
    Abs_LFC = rnorm(100, mean = 3, sd = 0.5),
    Cell_Type = rep(c("Astro", "Oligo"), each = 50),
    stringsAsFactors = FALSE
  )
  p <- plot_density_aggregated(df, "Abs_LFC", "LFC", "Abs(LFC)", c(1, 5))
  assert(inherits(p, "ggplot"), "Should return a ggplot object")
})

test("plot_density_aggregated returns ggplot with NULL limits (auto-scale)", {
  df <- data.frame(
    Abs_LFC = rnorm(100, mean = 3, sd = 0.5),
    Cell_Type = rep(c("Astro", "Oligo"), each = 50),
    stringsAsFactors = FALSE
  )
  p <- plot_density_aggregated(df, "Abs_LFC", "LFC", "Abs(LFC)", NULL)
  assert(inherits(p, "ggplot"), "Should return a ggplot object with NULL limits")
})

test("plot_density_aggregated handles all-NA column", {
  df <- data.frame(
    Abs_LFC = rep(NA_real_, 10),
    Cell_Type = rep("Astro", 10),
    stringsAsFactors = FALSE
  )
  p <- plot_density_aggregated(df, "Abs_LFC", "LFC", "Abs(LFC)", c(1, 5))
  assert(is.null(p), "Should return NULL for all-NA data")
})

test("plot_density_stratified returns ggplot with explicit limits", {
  df <- data.frame(
    Abs_LFC = rnorm(100, mean = 3, sd = 0.5),
    Cell_Type = rep(c("Astro", "Oligo"), each = 50),
    Biotype_Class = rep(c("Coding", "Non-Coding"), 50),
    stringsAsFactors = FALSE
  )
  p <- plot_density_stratified(df, "Abs_LFC", "LFC", "Abs(LFC)", c(1, 5))
  assert(inherits(p, "ggplot"), "Should return a ggplot object")
})

test("plot_density_stratified returns ggplot with NULL limits (auto-scale)", {
  df <- data.frame(
    Abs_LFC = rnorm(100, mean = 3, sd = 0.5),
    Cell_Type = rep(c("Astro", "Oligo"), each = 50),
    Biotype_Class = rep(c("Coding", "Non-Coding"), 50),
    stringsAsFactors = FALSE
  )
  p <- plot_density_stratified(df, "Abs_LFC", "LFC", "Abs(LFC)", NULL)
  assert(inherits(p, "ggplot"), "Should return a ggplot object with NULL limits")
})

# ==============================================================================
# 4b. AUTO-LIMIT COMPUTATION
# ==============================================================================

message("\n--- Auto-Limit Computation ---")

test("compute_auto_limits returns padded range", {
  df <- data.frame(val = c(1, 2, 3, 4, 5))
  result <- compute_auto_limits(df, "val", padding = 0.1)
  # range = [1, 5], span = 4, pad = 0.4
  assert(abs(result[1] - 0.6) < 1e-6, sprintf("Lower should be 0.6, got %.3f", result[1]))
  assert(abs(result[2] - 5.4) < 1e-6, sprintf("Upper should be 5.4, got %.3f", result[2]))
})

test("compute_auto_limits respects lim_lower override", {
  df <- data.frame(val = c(1, 2, 3, 4, 5))
  result <- compute_auto_limits(df, "val", padding = 0.1, lim_lower = 0)
  assert(result[1] == 0, sprintf("Lower should be 0, got %.3f", result[1]))
  assert(abs(result[2] - 5.4) < 1e-6, sprintf("Upper should be 5.4, got %.3f", result[2]))
})

test("compute_auto_limits respects lim_upper override", {
  df <- data.frame(val = c(1, 2, 3, 4, 5))
  result <- compute_auto_limits(df, "val", padding = 0.1, lim_upper = 10)
  assert(abs(result[1] - 0.6) < 1e-6, sprintf("Lower should be 0.6, got %.3f", result[1]))
  assert(result[2] == 10, sprintf("Upper should be 10, got %.3f", result[2]))
})

test("compute_auto_limits returns NULL for all-NA data", {
  df <- data.frame(val = rep(NA_real_, 5))
  result <- compute_auto_limits(df, "val")
  assert(is.null(result), "Should return NULL for all-NA data")
})

test("compute_auto_limits handles Inf values", {
  df <- data.frame(val = c(1, 2, Inf, 3, -Inf))
  result <- compute_auto_limits(df, "val", padding = 0.1)
  # Should exclude Inf, range = [1, 3], span = 2, pad = 0.2
  assert(abs(result[1] - 0.8) < 1e-6, sprintf("Lower should be 0.8, got %.3f", result[1]))
  assert(abs(result[2] - 3.2) < 1e-6, sprintf("Upper should be 3.2, got %.3f", result[2]))
})

# ==============================================================================
# 5. CONFIG VALIDATION
# ==============================================================================

message("\n--- Config Validation ---")

test("CELL_TYPES is defined", {
  assert(exists("CELL_TYPES") && length(CELL_TYPES) == 3,
         "CELL_TYPES should have 3 elements")
})

test("DENSITY_METRICS has 5 entries", {
  assert(length(DENSITY_METRICS) == 5,
         sprintf("Expected 5 metrics, got %d", length(DENSITY_METRICS)))
})

test("Each metric has required fields", {
  for (m in DENSITY_METRICS) {
    assert(all(c("col", "suffix", "lab", "lim") %in% names(m)),
           sprintf("Metric %s missing required fields", m$suffix))
    # lim can be NULL (auto-scale) or a length-2 numeric vector
    assert(is.null(m$lim) || length(m$lim) == 2,
           sprintf("Metric %s lim should be NULL or have 2 elements", m$suffix))
  }
})

test("Path helper functions work", {
  de_path <- get_de_results_path("Astro")
  assert(grepl("Roussos_DE_results.csv", de_path), "DE path should end with CSV")

  dds_path <- get_dds_checkpoint_path("Astro")
  assert(grepl("dds_checkpoint.qs", dds_path), "DDS path should end with .qs")

  gt_path <- get_gt_metadata_path("Astro", "cells", 1)
  assert(grepl("Run_1", gt_path), "GT path should contain Run_1")

  density_dir <- get_density_output_dir()
  assert(grepl("Cross_CellType", density_dir), "Density dir should contain Cross_CellType")

  pca_dir <- get_pca_output_dir("cells", "Astro")
  assert(grepl("cells/Astro", pca_dir), "PCA dir should contain cells/Astro")
})

test("CARTESIAN_MODES contains cells and reads", {
  assert("cells" %in% CARTESIAN_MODES, "Should contain 'cells'")
  assert("reads" %in% CARTESIAN_MODES, "Should contain 'reads'")
})

# ==============================================================================
# 6. DATA AVAILABILITY CHECKS
# ==============================================================================

message("\n--- Data Availability ---")

test("Full GT DE results exist for all cell types", {
  for (ct in CELL_TYPES) {
    path <- get_de_results_path(ct)
    assert(file.exists(path),
           sprintf("DE results missing for %s: %s", ct, path))
  }
})

test("DDS checkpoints exist for all cell types", {
  for (ct in CELL_TYPES) {
    path <- get_dds_checkpoint_path(ct)
    assert(file.exists(path),
           sprintf("DDS checkpoint missing for %s: %s", ct, path))
  }
})

test("Cartesian GT metadata exists for at least Run 1", {
  for (ct in CELL_TYPES) {
    for (mode in CARTESIAN_MODES) {
      path <- get_gt_metadata_path(ct, mode, 1)
      assert(file.exists(path),
             sprintf("GT metadata missing: %s/%s/Run_1 for %s", mode, ct, path))
    }
  }
})

test("Shard directories exist for all cell types", {
  for (ct in CELL_TYPES) {
    shard_dir <- get_shard_dir(ct)
    assert(dir.exists(shard_dir),
           sprintf("Shard dir missing for %s: %s", ct, shard_dir))
  }
})

test("Ground_Truth.qs files exist for all cell types x modes x Run 1", {
  for (ct in CELL_TYPES) {
    for (mode in CARTESIAN_MODES) {
      path <- get_cartesian_gt_path(ct, mode, 1)
      assert(file.exists(path),
             sprintf("Ground_Truth.qs missing: %s/%s/Run_1 at %s", ct, mode, path))
    }
  }
})

# ==============================================================================
# 7. STEP 3 UTILITY FUNCTIONS
# ==============================================================================

message("\n--- Step 3: Overlap & Classification ---")

test("classify_chromosome handles autosomal", {
  result <- classify_chromosome(c("1", "22", "X", "Y", "MT"))
  assert(result[1] == "Autosomal", "Chr 1 should be Autosomal")
  assert(result[2] == "Autosomal", "Chr 22 should be Autosomal")
  assert(result[3] == "SexChr", "Chr X should be SexChr")
  assert(result[4] == "SexChr", "Chr Y should be SexChr")
  assert(result[5] == "Other", "Chr MT should be Other")
})

test("classify_chromosome handles NA", {
  result <- classify_chromosome(c(NA_character_, "1"))
  assert(is.na(result[1]) || result[1] == "Other",
         "NA chromosome should be NA or Other")
  assert(result[2] == "Autosomal", "Chr 1 should be Autosomal")
})

test("assign_lfc_bin assigns correct bins", {
  result <- assign_lfc_bin(c(0.5, 1.5, 2.5, 3.5))
  assert(is.na(result[1]), "0.5 should be NA (below threshold)")
  assert(result[2] == "Small (1-2)", "1.5 should be Small (1-2)")
  assert(result[3] == "Medium (2-3)", "2.5 should be Medium (2-3)")
  assert(result[4] == "Large (3+)", "3.5 should be Large (3+)")
})

test("assign_lfc_bin handles boundary values", {
  result <- assign_lfc_bin(c(1.0, 2.0, 3.0))
  assert(result[1] == "Small (1-2)", "1.0 should be Small (1-2)")
  assert(result[2] == "Medium (2-3)", "2.0 should be Medium (2-3)")
  assert(result[3] == "Large (3+)", "3.0 should be Large (3+)")
})

test("compute_overlap_stats returns correct structure", {
  result <- compute_overlap_stats(c("A", "B", "C"), c("B", "C", "D"))
  assert(result$N_FullGT == 3, "N_FullGT should be 3")
  assert(result$N_Split == 3, "N_Split should be 3")
  assert(result$N_Overlap == 2, "N_Overlap should be 2")
  assert(result$N_FullOnly == 1, "N_FullOnly should be 1")
  assert(result$N_SplitOnly == 1, "N_SplitOnly should be 1")
  assert(abs(result$Recall - 2/3) < 1e-6, "Recall should be 2/3")
  assert(abs(result$Precision - 2/3) < 1e-6, "Precision should be 2/3")
})

test("compute_overlap_stats handles empty sets", {
  result <- compute_overlap_stats(character(0), c("A", "B"))
  assert(result$N_FullGT == 0, "N_FullGT should be 0")
  assert(is.na(result$Recall), "Recall should be NA for empty full GT")
  assert(result$Precision == 0, "Precision should be 0 (0 overlap / 2 split)")
})

test("compute_overlap_stats handles perfect overlap", {
  result <- compute_overlap_stats(c("A", "B"), c("A", "B"))
  assert(result$N_Overlap == 2, "N_Overlap should be 2")
  assert(result$Recall == 1.0, "Recall should be 1.0")
  assert(result$Precision == 1.0, "Precision should be 1.0")
})

test("load_split_gt_results returns NULL for missing file", {
  result <- load_split_gt_results("/nonexistent/path.qs")
  assert(is.null(result), "Should return NULL for missing file")
})

message("\n--- Step 3: Config Validation ---")

test("Recovery thresholds are defined", {
  assert(exists("RECOVERY_FDR_THRESHOLD"), "RECOVERY_FDR_THRESHOLD should exist")
  assert(exists("RECOVERY_LFC_THRESHOLD"), "RECOVERY_LFC_THRESHOLD should exist")
  assert(RECOVERY_FDR_THRESHOLD == 0.05, "RECOVERY_FDR_THRESHOLD should be 0.05")
  assert(RECOVERY_LFC_THRESHOLD == 1.0, "RECOVERY_LFC_THRESHOLD should be 1.0")
})

test("Recovery LFC bins are defined correctly", {
  assert(length(RECOVERY_LFC_BINS) == 3, "Should have 3 LFC bins")
  assert("Small (1-2)" %in% names(RECOVERY_LFC_BINS), "Should have Small bin")
  assert("Medium (2-3)" %in% names(RECOVERY_LFC_BINS), "Should have Medium bin")
  assert("Large (3+)" %in% names(RECOVERY_LFC_BINS), "Should have Large bin")
})

test("Split recovery path helpers work", {
  top_dir <- get_split_recovery_output_dir()
  assert(grepl("Split_Recovery$", top_dir), "Top dir should end with Split_Recovery")

  ct_dir <- get_split_recovery_output_dir("Astro")
  assert(grepl("Split_Recovery/Astro$", ct_dir), "CT dir should end with Split_Recovery/Astro")

  gt_path <- get_cartesian_gt_path("Astro", "cells", 1)
  assert(grepl("Ground_Truth.qs$", gt_path), "GT path should end with Ground_Truth.qs")
  assert(grepl("Run_1", gt_path), "GT path should contain Run_1")
})

# ==============================================================================
# 8. STEP 4 CONFIG VALIDATION
# ==============================================================================

message("\n--- Step 4: Chromosome Grid Config ---")

test("Chromosome grid config constants exist", {
  assert(exists("CHR_GRID_WIDTH"), "CHR_GRID_WIDTH should exist")
  assert(exists("CHR_GRID_HEIGHT"), "CHR_GRID_HEIGHT should exist")
  assert(exists("CHR_GRID_COMBINED_WIDTH"), "CHR_GRID_COMBINED_WIDTH should exist")
  assert(exists("CHR_GRID_TEXT_SIZE"), "CHR_GRID_TEXT_SIZE should exist")
  assert(exists("CHR_GRID_FILL_LOW"), "CHR_GRID_FILL_LOW should exist")
  assert(exists("CHR_GRID_FILL_HIGH"), "CHR_GRID_FILL_HIGH should exist")
})

test("Chromosome grid path helper works", {
  grid_dir <- get_chr_grid_output_dir()
  assert(grepl("Chromosome_Grid$", grid_dir),
         "Grid dir should end with Chromosome_Grid")
})

# ==============================================================================
# SUMMARY
# ==============================================================================

message(sprintf("\n=== Test Summary: %d passed, %d failed ===\n",
                tests_passed, tests_failed))

if (tests_failed > 0) {
  quit(status = 1)
} else {
  message("All tests passed!")
  quit(status = 0)
}
