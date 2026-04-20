#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Full Ground Truth Pipeline - Master Runner
# PURPOSE: Execute all pipeline steps in sequence
# USAGE: Rscript run_pipeline.R [steps] [n_cores_agg] [n_cores_de]
#        Rscript run_pipeline.R           # Run all steps with default cores
#        Rscript run_pipeline.R 1,2       # Run only steps 1 and 2
#        Rscript run_pipeline.R all 12 12 # Run all steps with 12 cores each
#        Rscript run_pipeline.R 2 8 8     # Run step 2 with 8 cores
# ==============================================================================

# ==============================================================================
# USER CONFIGURATION (Edit these values)
# ==============================================================================

# Steps to run (can be overridden by command line)
# Options: "all", or comma-separated step numbers like "1,2,3"
STEPS_TO_RUN <- "1,2,3,4,5"

# ==============================================================================
# COMMAND LINE ARGUMENT PARSING
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  STEPS_TO_RUN <- args[1]
}

# ==============================================================================
# SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/full_ground_truth/config.R",  # From project root
    "config.R"                        # From script directory
  )

  # Try sys.frame fallback (may not work in all contexts)
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

# Override cores from command line if provided
if (length(args) >= 2) {
  N_CORES_AGG <- as.integer(args[2])
}
if (length(args) >= 3) {
  N_CORES_DE <- as.integer(args[3])
}

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

# Use shared logging infrastructure if available
if (file.exists(LOGGING_SCRIPT)) {
  init_logging(
    pipeline_name = "full_ground_truth",
    log_dir = project_path("Datalink", "Run_Output_Logs"),
    verbose = TRUE
  )
} else {
  # Fallback: manual log file
  LOG_DIR <- project_path("Datalink", "Run_Output_Logs")
  if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)

  log_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG_FILE <- file.path(LOG_DIR, paste0("full_ground_truth_", log_timestamp, ".log"))

  log_con <- file(LOG_FILE, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message", append = TRUE)

  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
    message(sprintf("\n[LOG] Output saved to: %s", LOG_FILE))
  }, add = TRUE)

  cat(sprintf("=== Log started: %s ===\n", Sys.time()))
  cat(sprintf("=== Log file: %s ===\n\n", LOG_FILE))
}

# ==============================================================================
# STEP DEFINITIONS
# ==============================================================================

steps <- list(
  list(
    id = 1,
    name = "Run DE All Datasets",
    script = "01_Run_DE_All_Datasets.R",
    description = "Pseudobulk aggregation + DESeq2 for all 7 datasets x 3 cell types"
  ),
  list(
    id = 2,
    name = "Cross-Dataset Diagnostics",
    script = "02_Cross_Dataset_Diagnostics.R",
    description = "UpSet plots, correlation matrices, concordance, forest plots, expanded LFC bounds, FDP"
  ),
  list(
    id = 3,
    name = "Visualise Summary",
    script = "03_Visualise_Summary.R",
    description = "Cross-cell-type summary tables and overview figures"
  ),
  list(
    id = 4,
    name = "Gene Audit Matrix",
    script = "04_Gene_Audit_Matrix.R",
    description = "Audit gene composition across datasets by chromosome, biotype, and PAR status"
  ),
  list(
    id = 5,
    name = "Misc Analyses",
    script = "05_Misc_Analyses.R",
    description = "Ad-hoc analyses: shared-gene pairwise LFC correlations across datasets"
  )
)

# ==============================================================================
# DETERMINE WHICH STEPS TO RUN
# ==============================================================================

if (tolower(STEPS_TO_RUN) == "all") {
  steps_to_run <- seq_along(steps)
} else {
  steps_to_run <- as.integer(unlist(strsplit(STEPS_TO_RUN, ",")))
  steps_to_run <- steps_to_run[steps_to_run >= 1 & steps_to_run <= length(steps)]
}

# ==============================================================================
# PRINT BANNER
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  FULL GROUND TRUTH PIPELINE v", PIPELINE_VERSION, "\n", sep = "")
cat("================================================================================\n")
cat("\n")
cat("  Configuration:\n")
cat("    - Datasets:    ", paste(RUN_DATASETS, collapse = ", "), "\n")
cat("    - Cell Types:  ", paste(CELL_TYPES, collapse = ", "), "\n")
cat("    - Cores (AGG): ", N_CORES_AGG, "\n")
cat("    - Cores (DE):  ", N_CORES_DE, "\n")
cat("    - Force Recalc:", FORCE_RECALC, "\n")
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

  step <- steps[[step_id]]
  step_start <- Sys.time()

  cat("\n")
  log_section(sprintf("STEP %d: %s", step$id, step$name))
  log_info("Script: %s", step$script)
  log_info("Description: %s", step$description)

  script_path <- file.path(SCRIPT_DIR, step$script)

  if (!file.exists(script_path)) {
    log_error("Script not found: %s", script_path)
    step_status[[as.character(step_id)]] <- "FAILED - Script not found"
    next
  }

  # Execute step
  result <- tryCatch({
    source(script_path)
    "SUCCESS"
  }, error = function(e) {
    log_error("Step %d failed: %s", step_id, e$message)
    sprintf("FAILED - %s", e$message)
  })

  step_end <- Sys.time()
  step_duration <- difftime(step_end, step_start, units = "mins")

  step_times[[as.character(step_id)]] <- step_duration
  step_status[[as.character(step_id)]] <- result

  log_info("Step %d completed in %.1f minutes", step_id, as.numeric(step_duration))
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
  step <- steps[[step_id]]
  status <- step_status[[as.character(step_id)]]
  duration <- step_times[[as.character(step_id)]]

  status_icon <- if (status == "SUCCESS") "[OK]" else "[FAIL]"

  cat(sprintf("    %s Step %d: %s (%.1f min)\n",
              status_icon, step_id, step$name, as.numeric(duration)))

  if (status != "SUCCESS") {
    cat(sprintf("        Error: %s\n", status))
  }
}

cat("\n")

if (n_failed == 0) {
  cat("  Status: SUCCESS\n")
  cat("  All steps completed successfully.\n")
} else if (n_success > 0) {
  cat("  Status: PARTIAL SUCCESS\n")
  cat("  Some steps failed. Review errors above.\n")
} else {
  cat("  Status: FAILED\n")
  cat("  All steps failed. Review errors above.\n")
}

cat("\n")
cat("================================================================================\n")
cat("\n")

# Finalize shared logging if active
if (exists("finalize_logging") && is.function(finalize_logging)) {
  finalize_logging(
    status = if (n_failed == 0) "SUCCESS" else "PARTIAL",
    steps_ok = n_success
  )
}
