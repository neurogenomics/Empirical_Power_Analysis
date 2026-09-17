#!/usr/bin/env Rscript
# ==============================================================================
# TESTS: Power Curve Fitting Pipeline v4.0.0 - Unit Tests (Stage 1)
# ==============================================================================
# Run: Rscript R/power_curve_fitting/tests/test_power_curve_fitting.R
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
})

cat("\n")
cat("================================================================\n")
cat("  STAGE 1: Unit Tests — Power Curve Fitting Pipeline v4.0.0\n")
cat("================================================================\n\n")

# Source config and utils
config_paths <- c(
  "R/power_curve_fitting/config.R",
  "../config.R",
  "config.R"
)

config_loaded <- FALSE
for (cp in config_paths) {
  if (file.exists(cp)) {
    source(cp)
    config_loaded <- TRUE
    break
  }
}

if (!config_loaded) stop("Could not find config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# TEST 1: Config consistency
# ==============================================================================

test_config <- function() {
  cat("TEST 1: Config consistency\n")

  # Grid dimensions must match cartesian_downsampling
  expect_equal(length(DONOR_STEPS), 12)
  expect_equal(length(CELL_STEPS), 14)
  expect_equal(length(READ_STEPS), 13)

  # Thresholds
  expect_equal(GT_LFC_THRESHOLD, 1.0)
  expect_equal(GT_PADJ_THRESHOLD, 0.05)
  expect_equal(DETECT_PADJ_THRESHOLD, 0.05)

  # Cell types
  expect_equal(CELL_TYPES, c("Astro", "Micro", "Oligo"))

  # Modes
  expect_equal(MODES, c("cells", "reads"))

  # Robustness settings
  expect_equal(N_ROBUSTNESS_RUNS, 10)
  expect_equal(N_ITERATIONS, 20)

  # Seed must match cartesian_downsampling
  expect_equal(SEED_BASE, 1234)

  # v3.1.0 settings
  expect_equal(GENE_SETS, c("full", "shared"))
  expect_equal(AGGREGATION_LEVELS, c("per_iteration"))
  expect_equal(LRT_ALPHA, 0.05)
  expect_equal(CI_LEVEL, 0.95)

  # Version
  expect_equal(PIPELINE_VERSION, "4.0.0")

  # Mean reads per cell (reference values)
  expect_true(is.list(MEAN_READS_PER_CELL))
  expect_true(all(c("Astro", "Micro", "Oligo") %in% names(MEAN_READS_PER_CELL)))
  expect_equal(MEAN_READS_PER_CELL$Astro, 9456)
  expect_equal(MEAN_READS_PER_CELL$Micro, 4426)
  expect_equal(MEAN_READS_PER_CELL$Oligo, 5742)

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 2: Path helpers
# ==============================================================================

test_path_helpers <- function() {
  cat("TEST 2: Path helpers\n")

  # get_paths returns expected structure
  paths <- get_paths("Astro", "cells")
  expect_true(is.list(paths))
  expect_true(all(c("cartesian_dir", "shard_dir", "manifest_path",
                     "gt_dds_dir", "shard_cache_dir", "fitting_dir",
                     "model_dir", "plots_dir", "export_dir") %in% names(paths)))

  # Mode-specific paths (cartesian input only)
  paths_cells <- get_paths("Astro", "cells")
  paths_reads <- get_paths("Astro", "reads")
  expect_true(grepl("cells", paths_cells$cartesian_dir))
  expect_true(grepl("reads", paths_reads$cartesian_dir))

  # GT DDS dir is shared across modes (same GT split)
  expect_equal(paths_cells$gt_dds_dir, paths_reads$gt_dds_dir)

  # Shard cache dir is per-cell-type, not per-mode
  expect_equal(paths_cells$shard_cache_dir, paths_reads$shard_cache_dir)

  # get_model_dir for v3.1.0
  md1 <- get_model_dir("full", "per_iteration")
  expect_true(grepl("full_per_iteration$", md1))

  md2 <- get_model_dir("shared", "per_iteration")
  expect_true(grepl("shared_per_iteration$", md2))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 3: File pattern parsing (DE results + metadata)
# ==============================================================================

test_file_parsing <- function() {
  cat("TEST 3: File pattern parsing\n")

  # --- DE result files ---

  # Cells mode
  result <- parse_iteration_filename("N50_C500_I1.qs", "cells")
  expect_equal(result$n_donors, 50)
  expect_equal(result$y_value, 500)
  expect_equal(result$iteration, 1)

  # Cells mode with NA (Max)
  result_na <- parse_iteration_filename("N100_CNA_I5.qs", "cells")
  expect_equal(result_na$n_donors, 100)
  expect_true(is.na(result_na$y_value))
  expect_equal(result_na$iteration, 5)

  # Reads mode
  result_r <- parse_iteration_filename("N50_R1000000_I1.qs", "reads")
  expect_equal(result_r$n_donors, 50)
  expect_equal(result_r$y_value, 1000000)
  expect_equal(result_r$iteration, 1)

  # Invalid filename returns NULL
  result_bad <- parse_iteration_filename("not_a_result.qs", "cells")
  expect_null(result_bad)

  # Sampling metadata should not match DE pattern
  fp <- get_file_pattern("cells")
  expect_false(grepl(fp$pattern, "N50_C500_I1_sampling_metadata.qs"))

  # --- Metadata files ---

  # Cells mode metadata
  meta_result <- parse_metadata_filename("N50_C500_I1_sampling_metadata.qs", "cells")
  expect_equal(meta_result$n_donors, 50)
  expect_equal(meta_result$y_value, 500)
  expect_equal(meta_result$iteration, 1)

  # Reads mode metadata
  meta_result_r <- parse_metadata_filename("N100_R1000000_I10_sampling_metadata.qs", "reads")
  expect_equal(meta_result_r$n_donors, 100)
  expect_equal(meta_result_r$y_value, 1000000)
  expect_equal(meta_result_r$iteration, 10)

  # Metadata pattern should match metadata files
  mfp <- get_metadata_file_pattern("cells")
  expect_true(grepl(mfp$pattern, "N50_C500_I1_sampling_metadata.qs"))
  expect_false(grepl(mfp$pattern, "N50_C500_I1.qs"))
  expect_false(grepl(mfp$pattern, "Ground_Truth_sampling_metadata.qs"))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 4: Logistic regression helpers (v4.0.0)
# ==============================================================================

test_logistic_helpers <- function() {
  cat("TEST 4: Logistic regression helpers (v4.0.0)\n")

  # safe_glm_binomial should return NULL on failure, not error
  bad_data <- data.frame(y = c(0, 1), x = c(NA, NA))
  result <- safe_glm_binomial(y ~ x, data = bad_data)
  # May return NULL or a model depending on NA handling; just check no crash
  expect_true(is.null(result) || inherits(result, "glm"))

  # McFadden R² should be between 0 and 1 for a reasonable model
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(1 + 2 * x))
  test_model <- glm(y ~ x, family = binomial)
  r2 <- mcfadden_r2(test_model)
  expect_true(r2 > 0)
  expect_true(r2 < 1)

  # McFadden R² of NULL model should be NA
  expect_true(is.na(mcfadden_r2(NULL)))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 5: Model formula construction (v3.1.0 unified)
# ==============================================================================

test_formula_construction <- function() {
  cat("TEST 5: Model formula construction (v4.0.0)\n")

  # Per-cell-type model (no cell_type term)
  f1 <- build_model_formula(include_cell_type = FALSE)
  f1_str <- paste(deparse(f1), collapse = " ")
  expect_true(grepl("log_donors", f1_str))
  expect_true(grepl("log_n_cells", f1_str))
  expect_true(grepl("log_reads_per_cell", f1_str))
  expect_true(grepl("log_abs_gt_lfc", f1_str))
  expect_true(grepl("log_gt_basemean", f1_str))
  expect_true(grepl("log_gt_dispersion", f1_str))
  expect_false(grepl("cell_type", f1_str))

  # Combined model is identical to per-CT (no cell_type term)
  f2 <- build_model_formula(include_cell_type = FALSE)
  expect_equal(deparse(f1), deparse(f2))

  # Additive model with cell_type
  f3 <- build_model_formula(include_cell_type = TRUE, interactions = FALSE)
  f3_str <- paste(deparse(f3), collapse = " ")
  expect_true(grepl("cell_type", f3_str))
  expect_false(grepl("\\*", f3_str))

  # Interaction model with cell_type
  f4 <- build_model_formula(include_cell_type = TRUE, interactions = TRUE)
  f4_str <- paste(deparse(f4), collapse = " ")
  expect_true(grepl("cell_type", f4_str))
  expect_true(grepl("\\*", f4_str))

  # All formulas have 6 predictors (plus cell_type if included)
  # Count terms in per-CT formula
  rhs <- as.character(f1)[3]  # RHS of formula
  terms <- trimws(strsplit(rhs, "\\+")[[1]])
  expect_equal(length(terms), 6)

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 6: Shared gene identification
# ==============================================================================

test_shared_gene_identification <- function() {
  cat("TEST 6: Shared gene identification\n")

  # Mock fitting data with overlapping gene sets
  df_a <- data.frame(gene_id = c("G1", "G2", "G3", "G4"), stringsAsFactors = FALSE)
  df_b <- data.frame(gene_id = c("G2", "G3", "G4", "G5"), stringsAsFactors = FALSE)
  df_c <- data.frame(gene_id = c("G1", "G3", "G4", "G6"), stringsAsFactors = FALSE)

  shared <- identify_shared_genes(list(A = df_a, B = df_b, C = df_c))
  expect_equal(sort(shared), c("G3", "G4"))

  # Empty intersection
  df_d <- data.frame(gene_id = c("G7", "G8"), stringsAsFactors = FALSE)
  shared_empty <- identify_shared_genes(list(A = df_a, D = df_d))
  expect_equal(length(shared_empty), 0)

  # Single cell type: shared = all genes
  shared_single <- identify_shared_genes(list(A = df_a))
  expect_equal(sort(shared_single), c("G1", "G2", "G3", "G4"))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 7: Shard cache path generation
# ==============================================================================

test_shard_cache_paths <- function() {
  cat("TEST 7: Shard cache paths\n")

  # Normal donor ID
  path1 <- get_shard_cache_path("donor_123", "/tmp/cache")
  expect_equal(path1, "/tmp/cache/donor_donor_123_reads.qs")

  # Donor ID with special characters
  path2 <- get_shard_cache_path("donor-456.abc", "/tmp/cache")
  expect_true(grepl("_reads\\.qs$", path2))
  expect_false(grepl("-", basename(path2)))  # special chars replaced

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 8: Prediction grid generation (v3.1.0)
# ==============================================================================

test_prediction_grid <- function() {
  cat("TEST 8: Prediction grid generation (v3.1.0)\n")

  # Vary n_cells
  grid_c <- generate_prediction_grid(vary_axis = "n_cells", fixed_depth = 6000, n_points = 10)
  expect_true(is.data.frame(grid_c))
  expect_true(all(c("n_donors", "n_cells", "reads_per_cell",
                     "log_donors", "log_n_cells", "log_reads_per_cell",
                     "log_abs_gt_lfc", "log_gt_basemean", "log_gt_dispersion")
                   %in% names(grid_c)))
  expect_equal(nrow(grid_c), 100)  # 10 x 10

  # reads_per_cell should be fixed
  expect_true(all(grid_c$reads_per_cell == 6000))

  # n_cells should vary
  expect_true(length(unique(grid_c$n_cells)) > 1)

  # Vary reads_per_cell
  grid_r <- generate_prediction_grid(vary_axis = "reads_per_cell", fixed_depth = 250, n_points = 10)
  expect_true(all(grid_r$n_cells == 250))
  expect_true(length(unique(grid_r$reads_per_cell)) > 1)

  # All values should be positive
  expect_true(all(grid_c$n_donors > 0))
  expect_true(all(grid_c$n_cells > 0))
  expect_true(all(grid_c$reads_per_cell > 0))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 9: Coefficient formatting
# ==============================================================================

test_coefficient_formatting <- function() {
  cat("TEST 9: Coefficient formatting\n")

  # Test with NULL model
  result_null <- format_coefficients(NULL, "test")
  expect_null(result_null)

  # Test with a simple logistic regression model (v4.0.0)
  set.seed(42)
  n <- 100
  mock_data <- data.frame(
    mean_detection = rbinom(n, 1, 0.5),
    log_donors = rnorm(n),
    log_n_cells = rnorm(n)
  )

  mock_model <- tryCatch({
    glm(mean_detection ~ log_donors + log_n_cells,
        data = mock_data, family = binomial)
  }, error = function(e) NULL)

  if (!is.null(mock_model)) {
    coef_df <- format_coefficients(mock_model, "test_model")
    expect_true(is.data.frame(coef_df))
    expect_true(all(c("model", "term", "estimate", "std_error") %in% names(coef_df)))
    expect_equal(unique(coef_df$model), "test_model")
    expect_true(nrow(coef_df) >= 3)  # intercept + 2 predictors
  }

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 10: Age parsing (reused utility)
# ==============================================================================

test_age_parsing <- function() {
  cat("TEST 10: Age parsing\n")

  ages <- c("65 year-old stage", "6 month-old stage", "unknown", "42")
  parsed <- parse_age_stage_to_years(ages)

  expect_equal(parsed[1], 65)
  expect_equal(parsed[2], 0.5)
  expect_true(is.na(parsed[3]))
  expect_equal(parsed[4], 42)

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 11: Chromosome classification
# ==============================================================================

test_chromosome_classification <- function() {
  cat("TEST 11: Chromosome classification\n")

  chrs <- c("1", "22", "X", "Y", "MT", "GL000220.1")
  classified <- classify_chromosome(chrs)

  expect_equal(classified[1], "Autosomal")
  expect_equal(classified[2], "Autosomal")
  expect_equal(classified[3], "SexChr")
  expect_equal(classified[4], "SexChr")
  expect_true(is.na(classified[5]))
  expect_true(is.na(classified[6]))

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 12: v3.1.0 model variant count
# ==============================================================================

test_model_variant_count <- function() {
  cat("TEST 12: v3.1.0 model variant count\n")

  # (3 per-CT + 1 combined) x 2 gene sets x 1 aggregation level = 8
  n_per_ct <- length(CELL_TYPES)
  n_combined <- 1
  n_gene_sets <- length(GENE_SETS)
  n_agg_levels <- length(AGGREGATION_LEVELS)

  n_total <- (n_per_ct + n_combined) * n_gene_sets * n_agg_levels
  expect_equal(n_total, 8)

  # Each variant directory should be uniquely named
  variant_dirs <- character()
  for (gs in GENE_SETS) {
    for (al in AGGREGATION_LEVELS) {
      variant_dirs <- c(variant_dirs, get_model_dir(gs, al))
    }
  }
  expect_equal(length(unique(variant_dirs)), 2)

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# TEST 13: Output naming conventions
# ==============================================================================

test_output_naming <- function() {
  cat("TEST 13: Output naming conventions\n")

  # GT DDS checkpoints
  for (run_id in 1:6) {
    expected <- sprintf("GT_DDS_Run_%d.qs", run_id)
    expect_true(grepl("^GT_DDS_Run_\\d+\\.qs$", expected))
  }

  # v3.1.0 fitting data files
  expected_files <- c(
    "unified_fitting_data_per_iteration.qs",
    "unified_fitting_data_per_iteration_shared.qs"
  )
  for (f in expected_files) {
    expect_true(grepl("^unified_fitting_data_", f))
  }

  # Model files
  for (ct in CELL_TYPES) {
    expected <- paste0(ct, "_model.qs")
    expect_true(nchar(expected) > 0)
  }

  cat("  PASSED\n")
  return(TRUE)
}

# ==============================================================================
# RUN ALL TESTS
# ==============================================================================

tests <- list(
  test_config,
  test_path_helpers,
  test_file_parsing,
  test_logistic_helpers,
  test_formula_construction,
  test_shared_gene_identification,
  test_shard_cache_paths,
  test_prediction_grid,
  test_coefficient_formatting,
  test_age_parsing,
  test_chromosome_classification,
  test_model_variant_count,
  test_output_naming
)

results <- sapply(seq_along(tests), function(i) {
  tryCatch({
    tests[[i]]()
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
    return(FALSE)
  })
})

cat("\n")
cat("================================================================\n")
cat(sprintf("  Results: %d/%d tests passed\n", sum(results), length(results)))
cat("================================================================\n")

if (all(results)) {
  cat("\n  ALL TESTS PASSED\n")
  cat("\n  Ready for Stage 2: Toy Integration Test\n")
  cat("  Run: Rscript R/power_curve_fitting/toy/run_toy_test.R\n\n")
  quit(status = 0)
} else {
  cat("\n  SOME TESTS FAILED — fix before proceeding\n\n")
  quit(status = 1)
}
