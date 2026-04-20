# ==============================================================================
# TOY TEST: DE Diagnostics Pipeline - Integration Test with Synthetic Data
# ==============================================================================
# Creates minimal synthetic data and runs both pipeline steps to verify
# end-to-end functionality without requiring real data.
#
# Run from project root: Rscript R/de_diagnostics/toy/run_toy_test.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

message("\n=== DE Diagnostics Toy Integration Test ===\n")

# Source config and utils
source("R/de_diagnostics/config.R")
source("R/de_diagnostics/utils.R")

# Create temporary output directory
toy_output <- file.path(tempdir(), "de_diagnostics_toy")
if (dir.exists(toy_output)) unlink(toy_output, recursive = TRUE)
dir.create(toy_output, recursive = TRUE)

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

# ==============================================================================
# 1. TEST DENSITY PLOTS WITH SYNTHETIC DATA
# ==============================================================================

message("--- Step 1: Density Plots (Synthetic Data) ---")

set.seed(42)

# Create synthetic DE results mimicking real structure
n_genes <- 200
synthetic_de <- data.frame(
  gene_id = paste0("ENSG", sprintf("%011d", seq_len(n_genes))),
  log2FoldChange = c(rnorm(n_genes/2, mean = 3, sd = 1),
                     rnorm(n_genes/2, mean = -3, sd = 1)),
  lfcSE = runif(n_genes, 0.05, 0.4),
  baseMean = rlnorm(n_genes, meanlog = 5, sdlog = 2),
  padj = runif(n_genes, 0, 0.04),
  Cell_Type = rep(c("Astro", "Micro", "Oligo"), length.out = n_genes),
  dispersion = rlnorm(n_genes, meanlog = -2, sdlog = 1),
  Biotype_Class = rep(c("Coding", "Non-Coding"), length.out = n_genes),
  stringsAsFactors = FALSE
)

# Compute derived metrics
synthetic_de <- synthetic_de %>%
  mutate(
    Abs_LFC = abs(log2FoldChange),
    SNR = Abs_LFC / lfcSE,
    Log10_SNR = log10(SNR + 1e-6),
    Log10_BaseMean = log10(baseMean + 1),
    Log10_Dispersion = log10(dispersion + 1e-6)
  )

test("Synthetic data has correct dimensions", {
  assert(nrow(synthetic_de) == n_genes)
  assert(all(c("Abs_LFC", "lfcSE", "Log10_SNR", "Log10_BaseMean",
               "Log10_Dispersion") %in% colnames(synthetic_de)))
})

# Test each metric plot — resolve auto-limits the same way the real pipeline does
for (m in DENSITY_METRICS) {
  # Resolve limits: explicit > auto-scale
  if (!is.null(m$lim)) {
    xlim_val <- m$lim
  } else {
    xlim_val <- compute_auto_limits(
      synthetic_de, m$col,
      padding   = DENSITY_PADDING,
      lim_lower = m$lim_lower,
      lim_upper = m$lim_upper
    )
  }

  test(sprintf("Aggregated density plot for %s", m$suffix), {
    p <- plot_density_aggregated(synthetic_de, m$col, m$suffix, m$lab, xlim_val)
    assert(!is.null(p), "Plot should not be NULL")
    assert(inherits(p, "ggplot"), "Should be a ggplot")

    # Save to verify
    out_file <- file.path(toy_output, paste0("Density_", m$suffix, ".png"))
    ggsave(out_file, p, width = 8, height = 5, dpi = 72)
    assert(file.exists(out_file), "File should be saved")
    assert(file.info(out_file)$size > 0, "File should not be empty")
  })

  test(sprintf("Stratified density plot for %s", m$suffix), {
    p <- plot_density_stratified(synthetic_de, m$col, m$suffix, m$lab, xlim_val)
    assert(!is.null(p), "Plot should not be NULL")
    assert(inherits(p, "ggplot"), "Should be a ggplot")

    out_file <- file.path(toy_output, paste0("Density_", m$suffix, "_Stratified.png"))
    ggsave(out_file, p, width = 8, height = 7, dpi = 72)
    assert(file.exists(out_file), "File should be saved")
  })
}

# Test summary table
test("Summary table generation", {
  summary_df <- synthetic_de %>%
    group_by(Cell_Type) %>%
    summarise(
      N_Genes = n(),
      Median_AbsLFC = median(Abs_LFC, na.rm = TRUE),
      Median_SE = median(lfcSE, na.rm = TRUE),
      .groups = "drop"
    )
  assert(nrow(summary_df) == 3, "Should have 3 cell types")
  assert(all(summary_df$N_Genes > 0), "All cell types should have genes")

  csv_file <- file.path(toy_output, "Density_Summary.csv")
  write.csv(summary_df, csv_file, row.names = FALSE)
  assert(file.exists(csv_file), "CSV should be saved")
})

# ==============================================================================
# 2. TEST PCA HELPER FUNCTIONS
# ==============================================================================

message("\n--- Step 2: PCA Helpers (Synthetic Data) ---")

test("extract_gt_donor_ids handles missing file", {
  result <- extract_gt_donor_ids("/nonexistent/path.qs")
  assert(is.null(result), "Should return NULL for missing file")
})

test("get_all_donor_ids handles missing manifest", {
  result <- get_all_donor_ids("/nonexistent/shard_dir")
  assert(is.null(result), "Should return NULL for missing manifest")
})

test("Age parsing works for PCA metadata", {
  ages <- c("80 year-old stage", "65 year-old stage", "unknown",
            "72 year-old stage", "55 year-old stage")
  parsed <- parse_age_to_years(ages)
  assert(parsed[1] == 80, "First age should be 80")
  assert(is.na(parsed[3]), "Unknown should be NA")
  assert(parsed[4] == 72, "Fourth age should be 72")
})

# ==============================================================================
# 3. VERIFY REAL DATA PATHS (non-destructive)
# ==============================================================================

message("\n--- Data Path Verification ---")

test("DE results CSVs are readable", {
  for (ct in CELL_TYPES) {
    path <- get_de_results_path(ct)
    if (file.exists(path)) {
      df <- read.csv(path, nrows = 5)
      assert(nrow(df) > 0, sprintf("CSV for %s should have rows", ct))
      assert("padj" %in% colnames(df), sprintf("CSV for %s should have padj column", ct))
    }
  }
})

test("GT metadata files are readable", {
  for (ct in CELL_TYPES[1]) {  # Just test first cell type
    for (mode in CARTESIAN_MODES) {
      path <- get_gt_metadata_path(ct, mode, 1)
      if (file.exists(path)) {
        gt_ids <- extract_gt_donor_ids(path)
        assert(!is.null(gt_ids), sprintf("Should extract GT IDs for %s/%s", ct, mode))
        assert(length(gt_ids) > 0, sprintf("Should have >0 GT donors for %s/%s", ct, mode))
      }
    }
  }
})

# ==============================================================================
# CLEANUP & SUMMARY
# ==============================================================================

# Count output files
toy_files <- list.files(toy_output, recursive = TRUE)
message(sprintf("\nToy output files generated: %d", length(toy_files)))

# Clean up
unlink(toy_output, recursive = TRUE)

message(sprintf("\n=== Toy Test Summary: %d passed, %d failed ===\n",
                tests_passed, tests_failed))

if (tests_failed > 0) {
  quit(status = 1)
} else {
  message("All toy tests passed!")
  quit(status = 0)
}
