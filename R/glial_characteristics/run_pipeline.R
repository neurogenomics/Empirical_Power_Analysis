#!/usr/bin/env Rscript
# ==============================================================================
# PIPELINE: Glial Characteristics Cross-Dataset Comparison
# PURPOSE: Compare cell-level and donor-level characteristics between Roussos
#          and DeJager datasets (both DLPFC) to test whether inter-cell-type
#          ratios are biological invariants or dataset-specific artifacts.
# ==============================================================================
# 
# USAGE:
#   Rscript run_pipeline.R           # Run all steps
#   Rscript run_pipeline.R 1         # Run only Step 1
#   Rscript run_pipeline.R 1 2       # Run Steps 1 and 2
#   Rscript run_pipeline.R 2 3       # Run Steps 2 and 3
#
# STEPS:
#   1. Extract Metrics    — Load shards, compute UMIs/cell, genes/cell, cells/donor
#   2. Cross-Dataset Comparison — Side-by-side histograms, violins, ratio table
#   3. Bootstrap Ratio Tests    — 95% CIs on ratios, forest plots, interpretation
#
# OUTPUTS:
#   Datalink/Glial_Characteristics/
#   ├── Metrics/   — Per-dataset/cell-type summary CSVs and distributions
#   ├── Plots/     — Histograms, violin plots, ratio comparison charts
#   ├── Tests/     — Bootstrap CIs, forest plots, ratio difference tests
#   └── Summary/   — Combined summary tables
#
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Get script directory (works when sourced or run via Rscript)
SCRIPT_DIR <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  file.path(PROJECT_ROOT, "R", "glial_characteristics")
})

# Load configuration
source(file.path(SCRIPT_DIR, "config.R"))

# Parse command-line arguments for step selection
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  STEPS_TO_RUN <- as.integer(args)
  message(sprintf("[PIPELINE] Running steps: %s", paste(STEPS_TO_RUN, collapse = ", ")))
} else {
  message(sprintf("[PIPELINE] Running all steps: %s", paste(STEPS_TO_RUN, collapse = ", ")))
}

# ==============================================================================
# 2. LOG SETUP
# ==============================================================================

# Create log directory if needed
if (!dir.exists(LOG_DIR)) {
  dir.create(LOG_DIR, recursive = TRUE)
}

# Create timestamped log file
log_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
log_file <- file.path(LOG_DIR, sprintf("%s_%s.log", PIPELINE_NAME, log_timestamp))

# Sink output to log file
sink(log_file, split = TRUE)

message("==============================================================================")
message(sprintf("PIPELINE: %s", PIPELINE_NAME))
message(sprintf("VERSION:  %s", CONFIG_VERSION))
message(sprintf("STARTED:  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
message("==============================================================================")

# ==============================================================================
# 3. RUN PIPELINE STEPS
# ==============================================================================

start_time <- Sys.time()

# Step 1: Extract Metrics
if (1 %in% STEPS_TO_RUN) {
  message("\n######################################################################")
  message("# STEP 1: EXTRACT METRICS")
  message("######################################################################\n")
  
  step1_start <- Sys.time()
  source(file.path(SCRIPT_DIR, "01_Extract_Metrics.R"))
  step1_end <- Sys.time()
  
  message(sprintf("\n[STEP 1] Completed in %.1f minutes", 
                  as.numeric(difftime(step1_end, step1_start, units = "mins"))))
}

# Step 2: Cross-Dataset Comparison
if (2 %in% STEPS_TO_RUN) {
  message("\n######################################################################")
  message("# STEP 2: CROSS-DATASET COMPARISON")
  message("######################################################################\n")
  
  step2_start <- Sys.time()
  source(file.path(SCRIPT_DIR, "02_Cross_Dataset_Comparison.R"))
  step2_end <- Sys.time()
  
  message(sprintf("\n[STEP 2] Completed in %.1f minutes",
                  as.numeric(difftime(step2_end, step2_start, units = "mins"))))
}

# Step 3: Bootstrap Ratio Tests
if (3 %in% STEPS_TO_RUN) {
  message("\n######################################################################")
  message("# STEP 3: BOOTSTRAP RATIO TESTS")
  message("######################################################################\n")
  
  step3_start <- Sys.time()
  source(file.path(SCRIPT_DIR, "03_Bootstrap_Ratio_Tests.R"))
  step3_end <- Sys.time()
  
  message(sprintf("\n[STEP 3] Completed in %.1f minutes",
                  as.numeric(difftime(step3_end, step3_start, units = "mins"))))
}

# ==============================================================================
# 4. PIPELINE SUMMARY
# ==============================================================================

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

message("\n==============================================================================")
message("PIPELINE COMPLETE")
message("==============================================================================")
message(sprintf("Total runtime: %.1f minutes", total_time))
message(sprintf("Log file: %s", log_file))
message(sprintf("Outputs: %s", OUTPUT_BASE))
message("==============================================================================")

# Close log sink
sink()

# Print final message to console
cat(sprintf("\n[PIPELINE] Completed in %.1f minutes. Log: %s\n", total_time, log_file))
