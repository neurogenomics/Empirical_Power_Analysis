# ==============================================================================
# SCRIPT: Export Power Equations (v4.0.0)
# PURPOSE: Extract fitted coefficients, generate lookup tables, and create
#          a standalone R prediction function for end users.
#          v4.0.0 uses logistic regression with the unified formula.
# INPUT: Fitted models from Step 2
# OUTPUT: Coefficient CSVs, lookup tables, predict_power.R function
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
    "config.R",
    "../config.R"
  )
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 5: Export Power Equations (v4.0.0)")

export_dir <- file.path(OUTPUT_BASE, "Exports")
safe_mkdir(export_dir)

# ==============================================================================
# 2. EXPORT COEFFICIENTS (baseline: full_per_iteration)
# ==============================================================================

log_info("Exporting coefficients from baseline variant (full_per_iteration)...")

checkpoint_path <- file.path(get_model_dir("full", "per_iteration"), "all_models_checkpoint.qs")
if (!file.exists(checkpoint_path)) {
  log_error("Model checkpoint not found: %s", checkpoint_path)
  stop("Missing model checkpoint for export")
}

checkpoint <- qs::qread(checkpoint_path)

all_equations <- list()

# Per-cell-type models
for (ct in names(checkpoint$per_cell_type)) {
  model <- checkpoint$per_cell_type[[ct]]
  if (is.null(model)) next

  coefs <- coef(model)
  s <- summary(model)

  eq_info <- list(
    model_type = "per_cell_type",
    cell_type = ct,
    intercept = coefs["(Intercept)"],
    beta_log_donors = coefs["log_donors"],
    beta_log_n_cells = coefs["log_n_cells"],
    beta_log_reads_per_cell = coefs["log_reads_per_cell"],
    beta_log_lfc = coefs["log_abs_gt_lfc"],
    beta_log_basemean = coefs["log_gt_basemean"],
    beta_log_dispersion = coefs["log_gt_dispersion"],
    pseudo_r2 = mcfadden_r2(model),
    n_obs = nobs(model),
    aic = AIC(model)
  )

  all_equations[[ct]] <- eq_info
}

# Combined model (no cell_type term)
if (!is.null(checkpoint$combined)) {
  model <- checkpoint$combined
  coefs <- coef(model)
  s <- summary(model)

  all_equations[["combined"]] <- list(
    model_type = "combined",
    cell_type = "all",
    intercept = coefs["(Intercept)"],
    beta_log_donors = coefs["log_donors"],
    beta_log_n_cells = coefs["log_n_cells"],
    beta_log_reads_per_cell = coefs["log_reads_per_cell"],
    beta_log_lfc = coefs["log_abs_gt_lfc"],
    beta_log_basemean = coefs["log_gt_basemean"],
    beta_log_dispersion = coefs["log_gt_dispersion"],
    pseudo_r2 = mcfadden_r2(model),
    n_obs = nobs(model),
    aic = AIC(model)
  )
}

# Save equations as CSV
if (length(all_equations) > 0) {
  eq_df <- bind_rows(lapply(all_equations, as.data.frame))
  write.csv(eq_df, file.path(export_dir, "power_equation_coefficients.csv"),
            row.names = FALSE)
  log_info("Saved equation coefficients: %d models", nrow(eq_df))
}

# ==============================================================================
# 3. GENERATE LOOKUP TABLE
# ==============================================================================

log_info("Generating lookup table...")

lookup_grid <- expand.grid(
  n_donors = LOOKUP_DONOR_STEPS,
  n_cells = LOOKUP_CELL_STEPS,
  reads_per_cell = LOOKUP_RPC_STEPS,
  abs_gt_lfc = LOOKUP_LFC_VALUES,
  gt_basemean = LOOKUP_BASEMEAN_VALUES,
  gt_dispersion = LOOKUP_DISPERSION_VALUES
)

lookup_grid$log_donors <- log(lookup_grid$n_donors)
lookup_grid$log_n_cells <- log(lookup_grid$n_cells)
lookup_grid$log_reads_per_cell <- log(lookup_grid$reads_per_cell)
lookup_grid$log_abs_gt_lfc <- log(lookup_grid$abs_gt_lfc)
lookup_grid$log_gt_basemean <- log(lookup_grid$gt_basemean)
lookup_grid$log_gt_dispersion <- log(lookup_grid$gt_dispersion)

# Predict for each per-CT model
for (ct in names(checkpoint$per_cell_type)) {
  model <- checkpoint$per_cell_type[[ct]]
  if (is.null(model)) next

  col_name <- paste0("power_", ct)
  lookup_grid[[col_name]] <- predict(model, newdata = lookup_grid, type = "response")
}

# Predict for combined model
if (!is.null(checkpoint$combined)) {
  lookup_grid$power_combined <- predict(checkpoint$combined, newdata = lookup_grid,
                                         type = "response")
}

# Clean up helper columns
lookup_clean <- lookup_grid %>%
  select(-starts_with("log_"))

write.csv(lookup_clean, file.path(export_dir, "power_lookup_table.csv"),
          row.names = FALSE)
log_info("Saved lookup table: %d rows", nrow(lookup_clean))

# ==============================================================================
# 4. GENERATE STANDALONE R PREDICTION FUNCTION
# ==============================================================================

log_info("Generating standalone prediction function...")

func_lines <- c(
  '# ==============================================================================',
  '# predict_power.R — Standalone Power Prediction Function (v4.0.0)',
  '# ==============================================================================',
  '# Generated by the power_curve_fitting pipeline',
  sprintf('# Date: %s', Sys.Date()),
  sprintf('# Pipeline version: %s', PIPELINE_VERSION),
  '#',
  '# Usage:',
  '#   source("predict_power.R")',
  '#   predict_power(n_donors = 100, n_cells = 500, reads_per_cell = 6000,',
  '#                 abs_lfc = 3.0, basemean = 100, dispersion = 0.5,',
  '#                 cell_type = "Astro")',
  '#',
  '# Arguments:',
  '#   n_donors        - Number of donors (integer, >= 5)',
  '#   n_cells         - Cells per donor',
  '#   reads_per_cell  - Mean reads per cell',
  sprintf('#   abs_lfc         - Absolute log2 fold change of the target gene (> %.0f)', GT_LFC_THRESHOLD),
  '#   basemean        - Mean expression level (DESeq2 baseMean)',
  '#   dispersion      - Gene-level dispersion (DESeq2 dispersion)',
  '#   cell_type       - One of "Astro", "Micro", "Oligo", or "combined"',
  '#',
  '# Returns:',
  '#   Predicted detection rate (power) in [0, 1]',
  '# ==============================================================================',
  '',
  '# Fitted coefficients (from logistic regression with logit link)',
  '# Model: logit(power) = intercept + b1*log(donors) + b2*log(n_cells)',
  '#                       + b3*log(reads_per_cell) + b4*log(|LFC|)',
  '#                       + b5*log(baseMean) + b6*log(dispersion)',
  ''
)

# Add coefficient tables
func_lines <- c(func_lines, '.power_coefficients <- list(')

for (eq_name in names(all_equations)) {
  eq <- all_equations[[eq_name]]

  func_lines <- c(func_lines, sprintf('  %s = c(', eq_name))
  func_lines <- c(func_lines, sprintf('    intercept = %.10f,', eq$intercept))
  func_lines <- c(func_lines, sprintf('    log_donors = %.10f,', eq$beta_log_donors))
  func_lines <- c(func_lines, sprintf('    log_n_cells = %.10f,', eq$beta_log_n_cells))
  func_lines <- c(func_lines, sprintf('    log_reads_per_cell = %.10f,', eq$beta_log_reads_per_cell))
  func_lines <- c(func_lines, sprintf('    log_abs_lfc = %.10f,', eq$beta_log_lfc))
  func_lines <- c(func_lines, sprintf('    log_basemean = %.10f,', eq$beta_log_basemean))
  func_lines <- c(func_lines, sprintf('    log_dispersion = %.10f', eq$beta_log_dispersion))
  func_lines <- c(func_lines, '  ),')
}

func_lines <- c(func_lines, '  .dummy = NULL', ')')
func_lines <- c(func_lines, '')

# Add the prediction function
func_lines <- c(func_lines,
  'predict_power <- function(n_donors,',
  '                          n_cells,',
  '                          reads_per_cell,',
  '                          abs_lfc,',
  '                          basemean,',
  '                          dispersion,',
  '                          cell_type = "combined") {',
  '',
  '  # Validate inputs',
  '  stopifnot(n_donors >= 1)',
  '  stopifnot(n_cells >= 1)',
  '  stopifnot(reads_per_cell > 0)',
  '  stopifnot(abs_lfc > 0)',
  '  stopifnot(basemean > 0)',
  '  stopifnot(dispersion > 0)',
  '  stopifnot(cell_type %in% c("Astro", "Micro", "Oligo",
                           "Astro_Shallow", "Astro_Shallow_Fewer",
                           "combined"))',
  '',
  '  # Look up coefficients',
  '  coefs <- .power_coefficients[[cell_type]]',
  '  if (is.null(coefs)) {',
  '    stop(sprintf("No coefficients found for cell_type=%s", cell_type))',
  '  }',
  '',
  '  # Compute linear predictor (logit scale)',
  '  eta <- coefs["intercept"] +',
  '         coefs["log_donors"] * log(n_donors) +',
  '         coefs["log_n_cells"] * log(n_cells) +',
  '         coefs["log_reads_per_cell"] * log(reads_per_cell) +',
  '         coefs["log_abs_lfc"] * log(abs_lfc) +',
  '         coefs["log_basemean"] * log(basemean) +',
  '         coefs["log_dispersion"] * log(dispersion)',
  '',
  '  # Inverse logit to get probability',
  '  power <- 1 / (1 + exp(-eta))',
  '',
  '  return(as.numeric(power))',
  '}',
  '',
  '# Vectorised version for data frames',
  'predict_power_vec <- function(df) {',
  '  mapply(predict_power,',
  '         n_donors = df$n_donors,',
  '         n_cells = df$n_cells,',
  '         reads_per_cell = df$reads_per_cell,',
  '         abs_lfc = df$abs_lfc,',
  '         basemean = df$basemean,',
  '         dispersion = df$dispersion,',
  '         cell_type = if ("cell_type" %in% names(df)) df$cell_type else "combined")',
  '}',
  ''
)

writeLines(func_lines, file.path(export_dir, "predict_power.R"))
log_info("Saved standalone prediction function: predict_power.R")

# ==============================================================================
# 5. GENERATE SUMMARY
# ==============================================================================

summary_dir <- file.path(OUTPUT_BASE, "Summary")
safe_mkdir(summary_dir)

summary_info <- list(
  pipeline_version = PIPELINE_VERSION,
  date = as.character(Sys.Date()),
  gt_lfc_threshold = GT_LFC_THRESHOLD,
  gt_padj_threshold = GT_PADJ_THRESHOLD,
  detect_padj_threshold = DETECT_PADJ_THRESHOLD,
  n_robustness_runs = N_ROBUSTNESS_RUNS,
  n_iterations = N_ITERATIONS,
  cell_types = CELL_TYPES,
  modes = MODES,
  gene_sets = GENE_SETS,
  aggregation_levels = AGGREGATION_LEVELS,
  mean_reads_per_cell = MEAN_READS_PER_CELL,
  n_models_fitted = length(all_equations),
  equations = all_equations
)

qs::qsave(summary_info, file.path(summary_dir, "pipeline_summary.qs"))

# Also as CSV
summary_flat <- data.frame(
  parameter = c("pipeline_version", "date", "gt_lfc_threshold", "gt_padj_threshold",
                 "detect_padj_threshold", "n_robustness_runs", "n_iterations",
                 "cell_types", "modes", "gene_sets", "aggregation_levels",
                 "n_models_exported",
                 "reads_per_cell_Astro", "reads_per_cell_Micro", "reads_per_cell_Oligo"),
  value = c(PIPELINE_VERSION, as.character(Sys.Date()),
            GT_LFC_THRESHOLD, GT_PADJ_THRESHOLD, DETECT_PADJ_THRESHOLD,
            N_ROBUSTNESS_RUNS, N_ITERATIONS,
            paste(CELL_TYPES, collapse = ","), paste(MODES, collapse = ","),
            paste(GENE_SETS, collapse = ","), paste(AGGREGATION_LEVELS, collapse = ","),
            length(all_equations),
            MEAN_READS_PER_CELL$Astro, MEAN_READS_PER_CELL$Micro, MEAN_READS_PER_CELL$Oligo),
  stringsAsFactors = FALSE
)

write.csv(summary_flat, file.path(summary_dir, "pipeline_summary.csv"), row.names = FALSE)

log_section("STEP 5 COMPLETE")
log_info("Power equations exported: %d models", length(all_equations))
log_info("Outputs in: %s", export_dir)
