# ==============================================================================
# STEP 2: COLLECT CONSISTENT DEGs
# ==============================================================================
# PURPOSE: Filter genes significant in ≥90% of tissues, calculate averages
# INPUT: Datalink/Bulk_Analysis/Bulk_DE_genes_All.qs
# OUTPUT: Datalink/Bulk_Analysis/Filtered_Consistent_Genes.csv
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/bulk_deg_processing/config.R",  # From project root
    "config.R",                         # From script directory
    "../config.R"                       # From tests/ or toy/
  )
  
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

# Source utilities
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# MAIN LOGIC
# ==============================================================================

message("\n========================================")
message("STEP 2: COLLECT CONSISTENT DEGs")
message("========================================\n")

# Check if checkpoint exists and is valid
if (file.exists(FILTERED_FILE) && !FORCE_RECALC) {
  message("[CHECKPOINT] Loading cached data...")
  final_list <- fread(FILTERED_FILE)
  
  # Validate checkpoint
  is_valid <- validate_checkpoint(final_list, type = "filtered")
  
  if (!is_valid) {
    message("[CACHE INVALID] Regenerating...")
    final_list <- NULL
  } else {
    message("[CACHE VALID] Using cached data")
    message(sprintf("Loaded %d filtered genes", nrow(final_list)))
  }
}

# Generate fresh data if needed
if (!exists("final_list") || is.null(final_list)) {
  message("[COMPUTING] Loading bulk data from Step 1...")
  
  # Load data from Step 1
  if (!file.exists(BULK_ALL_FILE)) {
    stop(sprintf("Input file not found: %s\nRun Step 1 first.", BULK_ALL_FILE))
  }
  
  dt <- qread(BULK_ALL_FILE)
  setDT(dt)
  message(sprintf("Loaded %d rows", nrow(dt)))
  
  # Apply consistency filter
  message("\n[FILTERING] Applying consistency threshold...")
  final_list <- filter_consistent_genes(
    dt,
    p_threshold = P_VALUE_THRESHOLD,
    consistency_threshold = TISSUE_CONSISTENCY_THRESHOLD
  )
  
  # Apply LFC threshold if configured
  if (LFC_THRESHOLD > 0) {
    message(sprintf("\n[FILTERING] Applying LFC threshold (|LFC| > %.2f)...", LFC_THRESHOLD))
    n_before <- nrow(final_list)
    final_list <- final_list[abs(log2FoldChange) > LFC_THRESHOLD]
    n_after <- nrow(final_list)
    message(sprintf("Genes after LFC filter: %d (removed %d)", n_after, n_before - n_after))
  }
  
  # Sanity checks on filtered data
  filtered_checks <- list(
    "Non-empty" = function(data) {
      if (nrow(data) > 0) TRUE
      else "No genes passed filter"
    },
    "No NA in gene IDs" = function(data) {
      if (!any(is.na(data$gene))) TRUE
      else sprintf("%d genes have NA IDs", sum(is.na(data$gene)))
    },
    "Valid LFC range" = function(data) {
      max_lfc <- max(abs(data$log2FoldChange), na.rm = TRUE)
      if (max_lfc < 20) TRUE  # Sanity check: LFC shouldn't be extreme
      else sprintf("Max |LFC| = %.1f (unusually high)", max_lfc)
    },
    "Valid padj range" = function(data) {
      if (all(data$padj >= 0 & data$padj <= 1, na.rm = TRUE)) TRUE
      else "padj out of [0,1] range"
    }
  )
  
  sanity_passed <- sanity_check(final_list, "Filtered genes", filtered_checks)
  if (!sanity_passed) {
    warning("Filtered data failed sanity checks. Review results carefully.")
  }
  
  # Ensure column order
  setcolorder(final_list, c("gene", "log2FoldChange", "lfcSE", "padj"))
  
  # Save checkpoint
  message("\n[SAVING] Writing filtered gene list...")
  fwrite(final_list, FILTERED_FILE, quote = TRUE)
  message(sprintf("Saved to: %s", FILTERED_FILE))
}

# Final summary
message("\n========================================")
message("STEP 2 COMPLETE")
message("========================================")
message(sprintf("Filtered genes: %d", nrow(final_list)))
message(sprintf("Mean |LFC|: %.2f", mean(abs(final_list$log2FoldChange))))
message(sprintf("Median |LFC|: %.2f", median(abs(final_list$log2FoldChange))))
message(sprintf("Output: %s", FILTERED_FILE))
message("========================================\n")
