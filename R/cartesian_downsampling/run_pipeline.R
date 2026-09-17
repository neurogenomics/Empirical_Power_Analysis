#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Cartesian Downsampling Pipeline - Master Runner
# PURPOSE: Execute all pipeline steps in sequence
# USAGE: Rscript run_pipeline.R [steps] [n_cores]
#        Rscript run_pipeline.R           # Run all steps with default cores
#        Rscript run_pipeline.R 1,2       # Run only steps 1 and 2
#        Rscript run_pipeline.R all 48    # Run all steps with 48 cores
#        Rscript run_pipeline.R 3 16      # Run step 3 with 16 cores
#
# Step IDs:
#   1   = Generate Ground Truth
#   1.5 = Build HPC Job Manifest
#   2   = Run Downsampling Grid
#   3   = Aggregate Results
#   3.5 = Cross-Method Aggregation (pool DE methods -> all_methods)
#   4   = Visualize Power
#   5   = Continuous Grid Visualization
#   6   = Figure T / Figure S (small-effect figure of 4)
#
# To render the cross-method "all_methods" heatmaps in steps 4/5, set env
# CARTESIAN_INCLUDE_ALL_METHODS=true (step 3.5 must have run first).
# ==============================================================================

# ==============================================================================
# USER CONFIGURATION (Edit these values)
# ==============================================================================

# Number of cores for parallel processing
# This can be overridden by command line argument
N_CORES <- 18

# Steps to run (can be overridden by command line)
# Options: "all", or comma-separated step numbers like "1,2,3"
STEPS_TO_RUN <- "3,5"

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
    id_str = "1",
    name = "Generate Ground Truth",
    script = "01_Generate_Ground_Truth.R",
    description = "Create ground truth DE results using 50% of donors"
  ),
  list(
    id = 1.5,
    id_str = "1.5",
    name = "Build HPC Job Manifest",
    script = "01b_Build_HPC_Job_Manifest.R",
    description = "Create cost-balanced JOB_CHUNK manifest for PBS array runs"
  ),
  list(
    id = 2,
    id_str = "2",
    name = "Run Downsampling Grid",
    script = "02_Run_Downsampling_Grid.R",
    description = "Execute Cartesian grid of donor x cell experiments"
  ),
  list(
    id = 3,
    id_str = "3",
    name = "Aggregate Results",
    script = "03_Aggregate_Results.R",
    description = "Combine results and compute precision/recall metrics"
  ),
  list(
    id = 3.5,
    id_str = "3.5",
    name = "Cross-Method Aggregation",
    script = "08_Aggregate_CrossMethod.R",
    description = "Trial-weighted pool of DE methods -> all_methods checkpoint"
  ),
  list(
    id = 4,
    id_str = "4",
    name = "Visualize Power",
    script = "04_Visualise_Power.R",
    description = "Generate heatmaps and summary plots"
  ),
  list(
    id = 5,
    id_str = "5",
    name = "Continuous Grid Visualization",
    script = "07_Visualise_Continuous_Grid.R",
    description = "Generate interpolated continuous heatmaps with bilinear interpolation"
  ),
  list(
    id = 6,
    id_str = "6",
    name = "Figure T / Figure S",
    script = "09_Visualise_Figure_TS.R",
    description = "Small-effect figure of 4: LFC bin x 4 BaseMean categories"
  )
)

# ==============================================================================
# DETERMINE WHICH STEPS TO RUN
# ==============================================================================

all_step_ids <- sapply(steps, `[[`, "id")

if (tolower(STEPS_TO_RUN) == "all") {
  steps_to_run <- all_step_ids
} else {
  steps_to_run <- as.numeric(unlist(strsplit(STEPS_TO_RUN, ",")))
  steps_to_run <- steps_to_run[steps_to_run %in% all_step_ids]
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
  
  # Execute step (shared environment for RStudio sourcing compatibility)
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
  
  status_icon <- if (status == "SUCCESS") "[OK]" else "[FAIL]"
  
  cat(sprintf("    %s Step %s: %s (%.1f min)\n", 
              status_icon, step$id_str, step$name, as.numeric(duration)))
  
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

if (final_status != "SUCCESS") {
  quit(status = if (final_status == "PARTIAL") 2 else 1, save = "no")
}
