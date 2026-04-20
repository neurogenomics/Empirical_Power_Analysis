#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Power Curve Fitting Pipeline v4.1.0 - Master Runner
# PURPOSE: Execute all pipeline steps in sequence
# USAGE: Rscript run_pipeline.R [steps] [n_cores]
#        Rscript run_pipeline.R           # Run all steps with default cores
#        Rscript run_pipeline.R 1,2       # Run only steps 1 and 2
#        Rscript run_pipeline.R all 48    # Run all steps with 48 cores
#        Rscript run_pipeline.R 0 24      # Run only step 0 (GT DDS generation)
#
# Step IDs:
#   0   = Generate GT DDS Checkpoints
#   0.5 = Cache Shard Read Counts (Step 00b)
#   1   = Build Unified Fitting Dataset
#   2   = Fit Models (8 variants)
#   3   = Visualise
#   4   = Cross-Validate Modes
#   5   = Export Equations
#   6   = Compare All Models (including Astro_Shallow)
# ==============================================================================

# ==============================================================================
# USER CONFIGURATION (Edit these values)
# ==============================================================================

N_CORES <- 40

# Steps to run (can be overridden by command line)
STEPS_TO_RUN <- "3,4,5,6"

# ==============================================================================
# ==============================================================================

# ==============================================================================
# COMMAND LINE ARGUMENT PARSING
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  STEPS_TO_RUN <- args[1]
}

if (length(args) >= 2) {
  N_CORES <- as.integer(args[2])
}

# ==============================================================================
# SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
    "config.R"
  )

  tryCatch({
    ofile_path <- sys.frame(1)$ofile
    if (!is.null(ofile_path) && nzchar(ofile_path)) {
      possible_config_paths <- c(possible_config_paths,
                                 file.path(dirname(ofile_path), "config.R"))
    }
  }, error = function(e) NULL)

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

# ==============================================================================
# STEP DEFINITIONS
# ==============================================================================

steps <- list(
  list(
    id = 0,
    id_str = "0",
    name = "Generate GT DDS Checkpoints",
    script = "00_Generate_GT_DDS_Checkpoints.R",
    description = "Re-run DESeq2 on GT halves to extract dispersion estimates"
  ),
  list(
    id = 0.5,
    id_str = "0.5",
    name = "Cache Shard Read Counts",
    script = "00b_Cache_Shard_Read_Counts.R",
    description = "Pre-compute per-cell read counts for all donor shards"
  ),
  list(
    id = 1,
    id_str = "1",
    name = "Build Unified Fitting Dataset",
    script = "01_Build_Fitting_Dataset.R",
    description = "Parse metadata, merge cells+reads modes, build unified datasets with actual n_cells and rpc"
  ),
  list(
    id = 2,
    id_str = "2",
    name = "Fit Models",
    script = "02_Fit_Models.R",
    description = "Fit per-CT + combined models across all datasets and gene sets"
  ),
  list(
    id = 3,
    id_str = "3",
    name = "Visualise",
    script = "03_Visualise.R",
    description = "Power curves, surfaces, cell type comparisons, residual diagnostics"
  ),
  list(
    id = 4,
    id_str = "4",
    name = "Cross-Validate Modes",
    script = "04_Cross_Validate_Modes.R",
    description = "Cross-validate cells and reads mode predictions"
  ),
  list(
    id = 5,
    id_str = "5",
    name = "Export Equations",
    script = "05_Export_Equations.R",
    description = "Export coefficients, lookup tables, and prediction function"
  ),
  list(
    id = 6,
    id_str = "6",
    name = "Compare All Models",
    script = "06_Compare_All_Models.R",
    description = "Generate comparison plots across all models including Astro_Shallow"
  )
)

# ==============================================================================
# DETERMINE WHICH STEPS TO RUN
# ==============================================================================

all_step_ids <- sapply(steps, `[[`, "id")

if (tolower(STEPS_TO_RUN) == "all") {
  steps_to_run <- all_step_ids
} else {
  # Parse comma-separated step IDs (supports 0.5 for Step 00b)
  steps_to_run <- as.numeric(unlist(strsplit(STEPS_TO_RUN, ",")))
  steps_to_run <- steps_to_run[steps_to_run %in% all_step_ids]
}

# ==============================================================================
# PRINT BANNER
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  POWER CURVE FITTING PIPELINE v", PIPELINE_VERSION, "\n", sep = "")
cat("================================================================================\n")
cat("\n")
  cat("  Configuration:\n")
  cat("    - Datasets:      ", length(DATASET_CONFIGS), "\n")
  for (ds_cfg in DATASET_CONFIGS) {
    ds_model_names <- sapply(ds_cfg$cell_types, function(ct) get_dataset_model_name(ct, ds_cfg))
    cat(sprintf("      %s: %s (combined=%s)\n",
                ds_cfg$dataset,
                paste(sprintf("%s->%s", ds_cfg$cell_types, ds_model_names), collapse = ", "),
                ds_cfg$fit_combined))
  }
  cat("    - Cell types:    ", paste(CELL_TYPES, collapse = ", "), "\n")
  cat("    - Modes:         ", paste(MODES, collapse = ", "), " (merged in v3.1.0)\n")
  cat("    - Cores:         ", N_CORES, "\n")
  cat("    - GT |LFC|:      >", GT_LFC_THRESHOLD, "\n")
  cat("    - GT padj:       <", GT_PADJ_THRESHOLD, "\n")
  cat("    - Detect padj:   <", DETECT_PADJ_THRESHOLD, "\n")
  cat("    - Runs:          ", N_ROBUSTNESS_RUNS, "\n")
  cat("    - Iterations:    ", N_ITERATIONS, "\n")
  cat("    - Force Recalc:  ", FORCE_RECALC, "\n")
  cat("    - Gene sets:     ", paste(GENE_SETS, collapse = ", "), "\n")
  cat("    - Agg levels:    ", paste(AGGREGATION_LEVELS, collapse = ", "), "\n")
cat("\n")
cat("  Steps to run: ", paste(steps_to_run, collapse = ", "), "\n")
cat("\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# EXECUTE STEPS
# ==============================================================================

start_time <- Sys.time()
step_times <- list()
step_status <- list()

for (step_id in steps_to_run) {

  # Find step definition
  step_idx <- which(sapply(steps, `[[`, "id") == step_id)
  if (length(step_idx) == 0) {
    log_warning("Step %s not found, skipping", step_id)
    next
  }
  step <- steps[[step_idx]]

  step_start <- Sys.time()

  cat("\n")
  log_section(sprintf("STEP %s: %s", step$id_str, step$name))
  log_info("Script: %s", step$script)
  log_info("Description: %s", step$description)

  script_path <- file.path(SCRIPT_DIR, step$script)

  if (!file.exists(script_path)) {
    log_error("Script not found: %s", script_path)
    step_status[[step$id_str]] <- "FAILED - Script not found"
    next
  }

  result <- tryCatch({
    source(script_path)
    "SUCCESS"
  }, error = function(e) {
    log_error("Step %s failed: %s", step$id_str, e$message)
    sprintf("FAILED - %s", e$message)
  })

  step_end <- Sys.time()
  step_duration <- difftime(step_end, step_start, units = "mins")

  step_times[[step$id_str]] <- step_duration
  step_status[[step$id_str]] <- result

  log_info("Step %s completed in %.1f minutes", step$id_str, as.numeric(step_duration))
}

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("================================================================================\n")
cat("  PIPELINE COMPLETION REPORT\n")
cat("================================================================================\n")
cat("\n")

n_success <- sum(sapply(step_status, function(x) x == "SUCCESS"))
n_failed <- length(step_status) - n_success

cat("  Summary:\n")
cat("    - Steps run:    ", length(steps_to_run), "\n")
cat("    - Successful:   ", n_success, "\n")
cat("    - Failed:       ", n_failed, "\n")
cat("    - Total time:   ", sprintf("%.1f minutes", as.numeric(total_duration)), "\n")
cat("\n")

cat("  Step Details:\n")
for (step_id in steps_to_run) {
  step_idx <- which(sapply(steps, `[[`, "id") == step_id)
  step <- steps[[step_idx]]
  status <- step_status[[step$id_str]]
  duration <- step_times[[step$id_str]]

  status_icon <- if (!is.null(status) && status == "SUCCESS") "[OK]" else "[FAIL]"

  cat(sprintf("    %s Step %s: %s (%.1f min)\n",
              status_icon, step$id_str, step$name,
              if (!is.null(duration)) as.numeric(duration) else 0))

  if (!is.null(status) && status != "SUCCESS") {
    cat(sprintf("        Error: %s\n", status))
  }
}

cat("\n")

if (n_failed == 0) {
  cat("  Status: SUCCESS\n")
  cat("  All steps completed successfully.\n")
  final_status <- "SUCCESS"
} else if (n_success > 0) {
  cat("  Status: PARTIAL SUCCESS\n")
  cat("  Some steps failed. Review errors above.\n")
  final_status <- "PARTIAL"
} else {
  cat("  Status: FAILED\n")
  cat("  All steps failed. Review errors above.\n")
  final_status <- "FAILED"
}

cat("\n")
cat("================================================================================\n")
cat("\n")

# Finalize logging
if (exists("finalize_logging") && is_logging_active()) {
  finalize_logging(
    status = final_status,
    steps_ok = sprintf("%d/%d", n_success, length(steps_to_run))
  )
}
