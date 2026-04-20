# ==============================================================================
# MASTER RUNNER: Bulk DEG Processing Pipeline
# ==============================================================================
# This script runs the complete bulk DEG processing pipeline.
# 
# Usage (Command Line):
#   Rscript run_pipeline.R [steps]
#
# Usage (RStudio Background Job):
#   1. Open this file in RStudio
#   2. Click "Source" dropdown -> "Source as Background Job..."
#   3. Settings below will be used automatically
#
# Arguments (command line only):
#   steps:   Comma-separated list of steps to run (default: all)
#            1 = Extract bulk DEGs
#            2 = Collect consistent DEGs
#            3 = Annotate DEGs
#            4 = Sex chromosome analysis
#            5 = Visualise bulk vs SC concordance
#            6 = Recovery by expression level
#
# Examples:
#   Rscript run_pipeline.R                # Run all steps
#   Rscript run_pipeline.R 1,2            # Run steps 1-2 only
#   Rscript run_pipeline.R 3,4,5          # Run steps 3-5 only
# ==============================================================================

# ==============================================================================
# USER SETTINGS - EDIT THIS SECTION
# ==============================================================================
# These settings are used for both RStudio background jobs and command line
# (Command line arguments override these if provided)

STEPS_TO_RUN <- 1:6     # Which steps to run (1:6 = all, or e.g., c(1,2))

# ==============================================================================

# ==============================================================================
# ==============================================================================

# ==============================================================================
# DETECT EXECUTION CONTEXT AND RESOLVE PATHS
# ==============================================================================

# Detect if running from command line or interactively (RStudio)
is_interactive <- interactive()
args <- commandArgs(trailingOnly = TRUE)

# Override settings from command line arguments if provided
if (length(args) >= 1 && args[1] != "all") {
  STEPS_TO_RUN <- as.integer(unlist(strsplit(args[1], ",")))
}

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  # Search for config in common locations
  possible_config_paths <- c(
    "R/bulk_deg_processing/config.R",  # From project root
    "config.R"                          # From script directory
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

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# Define scripts
scripts <- c(
  "01_Extract_Bulk_DEGs.R",
  "02_Collect_Consistent_DEGs.R",
  "03_Annotate_DEGs.R",
  "04_Sex_Chr_Analysis.R",
  "05_Visualise_Bulk_vs_SC.R",
  "06_Recovery_By_Expression.R"
)

# Print banner
log_section("BULK DEG PROCESSING PIPELINE - MASTER RUNNER")
log_info("Execution mode: %s", if(is_interactive) "RStudio/Interactive" else "Command Line")
log_info("Script directory: %s", SCRIPT_DIR)
log_info("Steps to run: %s", paste(STEPS_TO_RUN, collapse = ", "))

# Print configuration
print_config()

# Verify script directory is valid
if (!file.exists(file.path(SCRIPT_DIR, "config.R"))) {
  log_error("Cannot find config.R in script directory: %s", SCRIPT_DIR)
  finalize_logging(status = "FAILED")
  stop("Please run from R/bulk_deg_processing/ or set working directory correctly.")
}

# Run selected scripts
n_success <- 0
n_failed <- 0
failed_steps <- c()

for (step in STEPS_TO_RUN) {
  if (step < 1 || step > length(scripts)) {
    log_warn("Invalid step %d. Skipping.", step)
    next
  }
  
  script_path <- file.path(SCRIPT_DIR, scripts[step])
  
  if (!file.exists(script_path)) {
    log_error("Script not found: %s", script_path)
    n_failed <- n_failed + 1
    failed_steps <- c(failed_steps, step)
    next
  }
  
  log_section(sprintf("RUNNING STEP %d: %s", step, scripts[step]))
  
  tryCatch({
    source(script_path)
    log_info("STEP %d COMPLETED SUCCESSFULLY", step)
    n_success <- n_success + 1
  }, error = function(e) {
    log_error("STEP %d FAILED: %s", step, e$message)
    n_failed <- n_failed + 1
    failed_steps <- c(failed_steps, step)
  })
}

# Final summary
log_section("PIPELINE EXECUTION COMPLETE")
log_info("Total steps attempted: %d", length(STEPS_TO_RUN))
log_info("Successful: %d", n_success)
log_info("Failed: %d", n_failed)

if (n_failed > 0) {
  log_error("Failed steps: %s", paste(failed_steps, collapse = ", "))
  finalize_logging(status = "PARTIAL", steps_ok = sprintf("%d/%d", n_success, length(STEPS_TO_RUN)))
} else {
  log_info("All steps completed successfully!")
  finalize_logging(status = "SUCCESS")
}
