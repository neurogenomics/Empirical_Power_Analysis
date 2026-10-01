# ==============================================================================
# MASTER RUNNER: DE Diagnostics Pipeline
# PURPOSE: Run all diagnostic steps in sequence with automatic logging
# USAGE: Rscript R/de_diagnostics/run_pipeline.R [steps]
#        Optional: Specify steps as space-separated numbers (e.g., "4 5")
#        Default: Run all remaining steps (4, 5)
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
    file.path(PROJECT_ROOT, "R", "de_diagnostics")
  }
})

source(file.path(script_path, "config.R"))

message(sprintf("=== DE DIAGNOSTICS PIPELINE ==="))
message(sprintf("Started: %s", Sys.time()))

# ==============================================================================
# 2. PARSE COMMAND LINE ARGUMENTS
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  steps_to_run <- as.integer(strsplit(trimws(args[1]), "[, ]+")[[1]])
  message(sprintf("Running steps: %s", paste(steps_to_run, collapse = ", ")))
} else {
  steps_to_run <- c(4, 5)
  message("Running all steps (4, 5)")
}

# ==============================================================================
# 3. DEFINE PIPELINE STEPS
# ==============================================================================

pipeline_steps <- list(
  list(
    id = 4,
    name = "DEG Chromosome Grid",
    script = "04_DEG_Chromosome_Grid.R",
    description = "DEG counts by chromosome type (autosomal/sex) x LFC magnitude bin"
  ),
  list(
    id = 5,
    name = "Shared Gene Diagnostics",
    script = "05_Shared_Gene_Diagnostics.R",
    description = "Paired cross-CT comparisons of 23 shared truth genes (Wilcoxon, correlations, scatter plots)"
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
    step_script_path <- file.path(SCRIPT_DIR, step$script)

    if (!file.exists(step_script_path)) {
      stop(sprintf("Script not found: %s", step_script_path))
    }

    # Use parent = environment() so step scripts inherit config
    step_env <- new.env(parent = environment())
    source(step_script_path, local = step_env)

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
if (exists("finalize_logging") && exists("is_logging_active") && is_logging_active()) {
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
