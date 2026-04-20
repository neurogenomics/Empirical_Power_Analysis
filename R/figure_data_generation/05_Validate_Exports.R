# ==============================================================================
# STEP 05: Validate All Exports
# ==============================================================================
# Checks that all expected CSV files exist in Paper_Source_Data/,
# validates row counts, and checks for unexpected NAs.
# ==============================================================================

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 05: VALIDATE EXPORTS ###\n")

# ==============================================================================
# Expected files and minimum row counts
# ==============================================================================

expected <- list(
  list(file = "dataset_statistics.csv", min_rows = 15, key_cols = c("Dataset", "N_Donors", "N_Cells")),
  list(file = "truth_gene_characteristics.csv", min_rows = 3, key_cols = c("Cell_Type", "N_Genes", "Median_LFC")),
  list(file = "model_coefficients.csv", min_rows = 6, key_cols = c("cell_type", "intercept", "pseudo_r2")),
  list(file = "model_summary.csv", min_rows = 6, key_cols = c("Model", "R2", "Req_Donors_80pct_moderate")),
  list(file = "fdp_summary.csv", min_rows = 39, key_cols = c("Cell_Type", "Source", "Group"))
)

n_pass <- 0
n_fail <- 0

for (e in expected) {
  fpath <- file.path(OUTPUT_BASE, e$file)

  # Check existence
  if (!file.exists(fpath)) {
    message(sprintf("  [FAIL] Missing: %s", e$file))
    n_fail <- n_fail + 1
    next
  }

  df <- read.csv(fpath)

  # Check row count
  if (nrow(df) < e$min_rows) {
    message(sprintf("  [FAIL] %s: %d rows (expected >= %d)", e$file, nrow(df), e$min_rows))
    n_fail <- n_fail + 1
    next
  }

  # Check key columns exist
  missing_cols <- setdiff(e$key_cols, colnames(df))
  if (length(missing_cols) > 0) {
    message(sprintf("  [FAIL] %s: missing columns: %s", e$file, paste(missing_cols, collapse = ", ")))
    n_fail <- n_fail + 1
    next
  }

  # Check for NAs in key columns (except FDP which legitimately has NAs for null tests)
  if (e$file != "fdp_summary.csv") {
    for (col in e$key_cols) {
      na_count <- sum(is.na(df[[col]]))
      if (na_count > 0) {
        message(sprintf("  [WARN] %s: %d NAs in column '%s'", e$file, na_count, col))
      }
    }
  }

  message(sprintf("  [PASS] %s: %d rows, %d columns", e$file, nrow(df), ncol(df)))
  n_pass <- n_pass + 1
}

message(sprintf("\n  Validation: %d passed, %d failed out of %d checks",
                n_pass, n_fail, length(expected)))

if (n_fail > 0) {
  message("\n  [WARN] Some validations failed. Run missing steps before using data in paper.")
} else {
  message("\n  All source data files validated. Paper can read from Paper_Source_Data/.")
}
