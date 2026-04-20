#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Cartesian Downsampling Pipeline - Master Runner
# PURPOSE: Execute all pipeline steps in sequence
# USAGE: Rscript run_pipeline.R [steps] [n_cores]
#        Rscript run_pipeline.R           # Run all steps with default cores
#        Rscript run_pipeline.R 1,2       # Run only steps 1 and 2
#        Rscript run_pipeline.R all 48    # Run all steps with 48 cores
#        Rscript run_pipeline.R 3 16      # Run step 3 with 16 cores
# ==============================================================================

# ==============================================================================
# USER CONFIGURATION (Edit these values)
# ==============================================================================

# Number of cores for parallel processing
# This can be overridden by command line argument
N_CORES <- 18

# Steps to run (can be overridden by command line)
# Options: "all", or comma-separated step numbers like "1,2,3"
STEPS_TO_RUN <- "3,7"

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

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  # Search for config in common locations
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",  # From project root
    "config.R"                             # From script directory
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
# STEP DEFINITIONS
# ==============================================================================

steps <- list(
  list(
    id = 1,
    name = "Generate Ground Truth",
    script = "01_Generate_Ground_Truth.R",
    description = "Create ground truth DE results using 50% of donors"
  ),
  list(
    id = 2,
    name = "Run Downsampling Grid",
    script = "02_Run_Downsampling_Grid.R",
    description = "Execute Cartesian grid of donor x cell experiments"
  ),
  list(
    id = 3,
    name = "Aggregate Results",
    script = "03_Aggregate_Results.R",
    description = "Combine results and compute precision/recall metrics"
  ),
  list(
    id = 4,
    name = "Visualize Power",
    script = "04_Visualise_Power.R",
    description = "Generate heatmaps and summary plots"
  ),
  list(
    id = 5,
    name = "Cross-Cell-Type Aggregation",
    script = "05_Visualise_Cross_CellType.R",
    description = "Pool results across cell types for aggregated heatmaps"
  ),
  list(
    id = 6,
    name = "Astro Variants Visualization",
    script = "06_Visualise_Astro_Variants.R",
    description = "Generate comparison plots for Astro depth-equalized and fewer-cells variants"
  ),
  list(
    id = 7,
    name = "Continuous Grid Visualization",
    script = "07_Visualise_Continuous_Grid.R",
    description = "Generate interpolated continuous heatmaps with bilinear interpolation"
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
cat("  CARTESIAN DOWNSAMPLING PIPELINE v", PIPELINE_VERSION, "\n", sep = "")
cat("================================================================================\n")
cat("\n")
y_steps <- get_y_steps()
cat("  Configuration:\n")
cat("    - Dataset:     ", SHARD_DATASET, "\n")
cat("    - Mode:        ", DOWNSAMPLE_MODE, "\n")
cat("    - Cores:       ", N_CORES, "\n")
cat("    - Runs:        ", N_ROBUSTNESS_RUNS, "\n")
cat("    - Subsamples:  ", N_SUBSAMPLES, "\n")
cat("    - Grid:        ", length(DONOR_STEPS), " x ", length(y_steps), " = ", 
    length(DONOR_STEPS) * length(y_steps), " points\n", sep = "")
cat("    - Y-axis:      ", ifelse(DOWNSAMPLE_MODE == "reads", "Read steps", "Cell steps"), "\n")
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
  
  # Execute step (shared environment for RStudio sourcing compatibility)
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
