# ==============================================================================
# MASTER RUNNER: EDA Roussos Pipeline
# PURPOSE: Run all exploratory analysis steps in sequence with automatic logging
# USAGE: Rscript R/eda_roussos/run_pipeline.R [steps]
#        Optional: Specify steps as comma-separated numbers (e.g., "1,2,3")
#        Default: Uses STEPS_TO_RUN from config.R
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
})

# ==============================================================================
# 1. SETUP LOGGING
# ==============================================================================

# Load config - use tryCatch for Rscript compatibility
script_path <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # Fallback: try commandArgs --file, then hardcoded path
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    file.path(PROJECT_ROOT, "R", "eda_roussos")
  }
})

source(file.path(script_path, "config.R"))

message(sprintf("=== EDA ROUSSOS PIPELINE ==="))
message(sprintf("Started: %s", Sys.time()))

# ==============================================================================
# 2. PARSE COMMAND LINE ARGUMENTS
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  # Command-line args override config default
  steps_to_run <- as.integer(unlist(strsplit(args[1], "[, ]+")))
  message(sprintf("Running steps (from command line): %s",
                  paste(steps_to_run, collapse = ", ")))
} else {
  # Use config default
  steps_to_run <- STEPS_TO_RUN
  message(sprintf("Running steps (from config): %s",
                  paste(steps_to_run, collapse = ", ")))
}

# ==============================================================================
# 3. DEFINE PIPELINE STEPS
# ==============================================================================

pipeline_steps <- list(
  list(
    id = 1,
    name = "PCA Diagnostics",
    script = "01_PCA_Diagnostics.R",
    description = "Generate PCA plots colored by metadata variables"
  )
)

# ==============================================================================
# 4. RUN PIPELINE
# ==============================================================================

results <- list()

for (step in pipeline_steps) {
  
  if (!(step$id %in% steps_to_run)) {
    message(sprintf("\n[SKIP] Step %d: %s", step$id, step$name))
    next
  }
  
  message(sprintf("\n========================================"))
  message(sprintf("STEP %d: %s", step$id, step$name))
  message(sprintf("Description: %s", step$description))
  message(sprintf("Script: %s", step$script))
  message(sprintf("========================================\n"))
  
  step_start <- Sys.time()
  
  # Run step
  step_result <- tryCatch({
    script_path <- file.path(SCRIPT_DIR, step$script)
    
    if (!file.exists(script_path)) {
      stop(sprintf("Script not found: %s", script_path))
    }
    
    # Use parent = environment() so step scripts inherit SCRIPT_DIR from config
    step_env <- new.env(parent = environment())
    source(script_path, local = step_env)
    
    list(
      status = "SUCCESS",
      error = NULL,
      duration = as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    )
    
  }, error = function(e) {
    list(
      status = "FAILED",
      error = e$message,
      duration = as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    )
  })
  
  # Store result
  results[[step$id]] <- step_result
  
  # Report
  if (step_result$status == "SUCCESS") {
    message(sprintf("\n✓ Step %d completed successfully (%.1f seconds)", 
                    step$id, step_result$duration))
  } else {
    message(sprintf("\n✗ Step %d FAILED (%.1f seconds)", 
                    step$id, step_result$duration))
    message(sprintf("Error: %s", step_result$error))
  }
}

# ==============================================================================
# 5. SUMMARY
# ==============================================================================

message(sprintf("\n========================================"))
message(sprintf("PIPELINE SUMMARY"))
message(sprintf("========================================\n"))

n_success <- sum(sapply(results, function(x) x$status == "SUCCESS"))
n_failed <- sum(sapply(results, function(x) x$status == "FAILED"))
n_total <- length(results)

message(sprintf("Steps run: %d", n_total))
message(sprintf("Successful: %d", n_success))
message(sprintf("Failed: %d", n_failed))

if (n_failed > 0) {
  message("\nFailed steps:")
  for (i in seq_along(results)) {
    if (results[[i]]$status == "FAILED") {
      message(sprintf("  Step %d: %s", i, results[[i]]$error))
    }
  }
}

total_duration <- sum(sapply(results, function(x) x$duration))
message(sprintf("\nTotal duration: %.1f seconds (%.1f minutes)", 
                total_duration, total_duration / 60))

message(sprintf("\nFinished: %s", Sys.time()))
if (exists("log_file")) {
  message(sprintf("Log file: %s", log_file))
}

# ==============================================================================
# 6. CLEANUP
# ==============================================================================

# Determine final status
if (n_failed == 0) {
  final_status <- "SUCCESS"
} else if (n_success > 0) {
  final_status <- "PARTIAL"
} else {
  final_status <- "FAILED"
}

# Finalize logging (closes sinks and writes footer)
if (exists("finalize_logging") && is_logging_active()) {
  finalize_logging(
    status = final_status,
    steps_ok = sprintf("%d/%d", n_success, n_total)
  )
}

# Exit with appropriate code
if (n_failed > 0) {
  quit(status = 1)
} else {
  quit(status = 0)
}
