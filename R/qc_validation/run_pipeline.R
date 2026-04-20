# ==============================================================================
# SCRIPT: QC Validation Pipeline - Production Runner
# PURPOSE: Stage 3 — Full production run across all Roussos shards
# RUN: Rscript R/qc_validation/run_pipeline.R
# ==============================================================================

# ==============================================================================
# ==============================================================================

message("=" |> rep(70) |> paste0(collapse = ""))
message("QC VALIDATION PIPELINE — PRODUCTION RUN")
message(sprintf("Started: %s", Sys.time()))
message("=" |> rep(70) |> paste0(collapse = ""))

# ==============================================================================
# Source config
# ==============================================================================
source(file.path(PROJECT_ROOT, "R", "qc_validation", "config.R"))

# ==============================================================================
# Step 00: Inspect H5AD Layers
# ==============================================================================
log_section("Step 00: Inspect H5AD Layers")

h5ad_report <- file.path(OUTPUT_BASE, "h5ad_layer_report.txt")
if (file.exists(h5ad_report)) {
  log_info("H5AD layer report already exists. Skipping Step 00.")
  log_info("Delete %s to force re-inspection.", h5ad_report)
} else {
  tryCatch({
    source(file.path(SCRIPT_DIR, "00_Inspect_H5AD_Layers.R"))
    log_info("Step 00 complete.")
  }, error = function(e) {
    log_error("Step 00 failed: %s", e$message)
    log_warn("Continuing without H5AD inspection. Shard-level QC will still run.")
  })
}

# ==============================================================================
# Step 01: Compute QC Metrics
# ==============================================================================
log_section("Step 01: Compute QC Metrics")

tryCatch({
  source(file.path(SCRIPT_DIR, "01_Compute_QC_Metrics.R"))
  log_info("Step 01 complete.")
}, error = function(e) {
  log_error("Step 01 failed: %s", e$message)
  finalize_logging(status = "FAILED")
  stop(e)
})

# ==============================================================================
# Step 02: Generate QC Report
# ==============================================================================
log_section("Step 02: Generate QC Report")

tryCatch({
  source(file.path(SCRIPT_DIR, "02_Generate_QC_Report.R"))
  log_info("Step 02 complete.")
}, error = function(e) {
  log_error("Step 02 failed: %s", e$message)
  finalize_logging(status = "PARTIAL", steps_ok = "Step 01 OK, Step 02 FAILED")
  stop(e)
})

# ==============================================================================
# Final Summary
# ==============================================================================
log_section("Pipeline Complete")

message("\n=== QC VALIDATION PIPELINE COMPLETE ===")
message(sprintf("Output directory: %s", OUTPUT_BASE))
message("")
message("Key outputs:")
message(sprintf("  H5AD report: %s", h5ad_report))
for (ct in CELL_TYPES) {
  message(sprintf("  %s verdict: %s", ct, file.path(OUTPUT_BASE, ct, "QC_Verdict.txt")))
}
message(sprintf("  Global summary: %s", file.path(OUTPUT_BASE, "Global_QC_Summary.csv")))

finalize_logging(status = "SUCCESS")
