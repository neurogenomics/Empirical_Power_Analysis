#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Figure Data Generation Pipeline v1.0.0 - Master Runner
# PURPOSE: Generate all CSV source data for paper_final_draft.rmd
# USAGE: Rscript run_pipeline.R [steps] [n_cores]
#        Rscript run_pipeline.R           # Run all steps
#        Rscript run_pipeline.R 1         # Run only step 1 (Dataset Statistics)
#        Rscript run_pipeline.R 1,2       # Run steps 1 and 2
#        Rscript run_pipeline.R all 24    # Run all steps with 24 cores
#
# Step IDs:
#   1 = Dataset Statistics (from shards — ALL 5 datasets)
#   2 = Truth Gene Characteristics (from Full Ground Truth DE results)
#   3 = Model Coefficients & Summary (from Power Curve Fitting exports)
#   4 = FDP Summary (from Full Ground Truth diagnostics)
#   5 = Validate All Exports
#   6 = Autosomal vs Sex-Chromosome Truth Gene Properties
# ==============================================================================

N_CORES <- 12
STEPS_TO_RUN <- "all"

# ==============================================================================
# LOGGING
# ==============================================================================

# ==============================================================================
# COMMAND LINE ARGUMENTS
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) STEPS_TO_RUN <- args[1]
if (length(args) >= 2) N_CORES <- as.integer(args[2])

# ==============================================================================
# SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_paths <- c(
    "R/figure_data_generation/config.R",
    "config.R",
    file.path(dirname(tryCatch(sys.frame(1)$ofile, error = function(e) "")), "config.R")
  )
  for (p in possible_paths) {
    if (file.exists(p)) { source(p); break }
  }
}

if (!exists("PROJECT_ROOT")) stop("Could not load config.R")

# ==============================================================================
# STEP DEFINITIONS
# ==============================================================================

steps <- list(
  list(id = 1, name = "Dataset Statistics",
       script = file.path(SCRIPT_DIR, "01_Dataset_Statistics.R"),
       desc = "Compute N_Donors, N_Cells, UMIs/cell, genes/cell from shards (all 5 datasets)"),
  list(id = 2, name = "Truth Gene Characteristics",
       script = file.path(SCRIPT_DIR, "02_Truth_Gene_Characteristics.R"),
       desc = "Extract truth gene LFC, baseMean, dispersion from Roussos GT results"),
  list(id = 3, name = "Model Coefficients & Summary",
       script = file.path(SCRIPT_DIR, "03_Collect_Model_Data.R"),
       desc = "Read power_equation_coefficients.csv, compute model summary with predictions"),
  list(id = 4, name = "FDP Summary",
       script = file.path(SCRIPT_DIR, "04_Collect_FDP_Data.R"),
       desc = "Read FDP_Summary.csv, format for paper source data tables"),
  list(id = 5, name = "Validate Exports",
       script = file.path(SCRIPT_DIR, "05_Validate_Exports.R"),
       desc = "Check all expected CSVs exist, row counts correct, no NAs in key columns"),
  list(id = 6, name = "Autosomal vs Sex-Chr Summary",
       script = file.path(SCRIPT_DIR, "06_Autosomal_vs_SexChr_Summary.R"),
       desc = "Build _scratch_autosomal_vs_sexchr_summary.csv for paper Table (line 589)")
)

# ==============================================================================
# DETERMINE WHICH STEPS TO RUN
# ==============================================================================

if (STEPS_TO_RUN == "all") {
  run_ids <- sapply(steps, function(s) s$id)
} else {
  run_ids <- as.numeric(strsplit(STEPS_TO_RUN, ",")[[1]])
}

# ==============================================================================
# EXECUTE
# ==============================================================================

dir.create(OUTPUT_BASE, recursive = TRUE, showWarnings = FALSE)

results <- list()
for (step in steps) {
  if (!(step$id %in% run_ids)) next

  message(sprintf("\n%s\n", paste(rep("=", 70), collapse = "")))
  message(sprintf("STEP %d: %s", step$id, step$name))
  message(sprintf("  %s", step$desc))
  message(sprintf("%s\n", paste(rep("=", 70), collapse = "")))

  t0 <- Sys.time()
  status <- tryCatch({
    source(step$script, local = TRUE)
    "SUCCESS"
  }, error = function(e) {
    message(sprintf("[ERROR] Step %d failed: %s", step$id, conditionMessage(e)))
    paste("FAILED:", conditionMessage(e))
  })
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  results[[length(results) + 1]] <- list(
    id = step$id, name = step$name, status = status, elapsed = elapsed
  )
}

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("FIGURE DATA GENERATION PIPELINE — SUMMARY")
message(paste(rep("=", 70), collapse = ""), "\n")

for (r in results) {
  icon <- if (grepl("^SUCCESS", r$status)) "[OK]" else "[!!]"
  message(sprintf("  %s Step %d: %-35s %s (%.1fs)",
                  icon, r$id, r$name, r$status, r$elapsed))
}

n_ok <- sum(sapply(results, function(r) grepl("^SUCCESS", r$status)))
message(sprintf("\n  %d/%d steps completed successfully.", n_ok, length(results)))
message(sprintf("  Output: %s", OUTPUT_BASE))
message("\n[COMPLETE]")
