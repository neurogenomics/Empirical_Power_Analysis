# ==============================================================================
# STEP 1: EXTRACT BULK DEGs
# ==============================================================================
# PURPOSE: Load raw bulk DEG data and filter by tissue coverage threshold
# INPUT: Datalink/Bulk_Data/Bulk_Significant_DEGs.txt
# OUTPUT: Datalink/Bulk_Analysis/Bulk_DE_genes_All.qs
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
message("STEP 1: EXTRACT BULK DEGs")
message("========================================\n")

# Check if checkpoint exists and is valid
bulk_data <- NULL
if (file.exists(BULK_ALL_FILE) && !FORCE_RECALC) {
  message("[CHECKPOINT] Loading cached data...")
  bulk_data <- qread(BULK_ALL_FILE)
  
  # Validate checkpoint
  is_valid <- validate_checkpoint(bulk_data, type = "bulk_all")
  
  if (!is_valid) {
    message("[CACHE INVALID] Regenerating...")
    bulk_data <- NULL
  } else {
    message("[CACHE VALID] Using cached data")
    message(sprintf("Loaded %d rows for %d unique genes", 
                    nrow(bulk_data), uniqueN(bulk_data$gene)))
  }
}

# Generate fresh data if needed
if (is.null(bulk_data)) {
  message("[COMPUTING] Loading raw bulk data...")
  
  # Load raw data
  bulk_data <- fread(INPUT_FILE)
  message(sprintf("Loaded %d rows from raw file", nrow(bulk_data)))
  
  # Sanity checks on raw data
  raw_checks <- list(
    "Has required columns" = function(data) {
      required <- c("gene", "tissue", "effsize", "effsize_se", "lfsr")
      if (all(required %in% names(data))) TRUE
      else sprintf("Missing: %s", paste(setdiff(required, names(data)), collapse = ", "))
    },
    "Non-empty" = function(data) {
      if (nrow(data) > 0) TRUE
      else "Data has 0 rows"
    },
    "Valid lfsr range" = function(data) {
      if (all(data$lfsr >= 0 & data$lfsr <= 1, na.rm = TRUE)) TRUE
      else sprintf("lfsr out of [0,1] range")
    }
  )
  
  sanity_passed <- sanity_check(bulk_data, "Raw bulk data", raw_checks)
  if (!sanity_passed) {
    stop("Raw data failed sanity checks. Investigate input file.")
  }
  
  # Calculate tissue coverage
  total_tissues <- uniqueN(bulk_data$tissue)
  message(sprintf("Total unique tissues: %d", total_tissues))
  
  # Calculate cutoff (currently 0% - keep all genes)
  cutoff_percent <- 0  # Can be made configurable if needed
  cutoff <- cutoff_percent * total_tissues
  message(sprintf("Filtering for genes present in more than %.0f tissues", cutoff))
  
  # Count tissues per gene
  gene_counts <- count_tissues_per_gene(bulk_data)
  
  # Identify genes to keep
  keep_genes <- gene_counts[tissue_count > cutoff, gene]
  message(sprintf("Genes passing filter: %d", length(keep_genes)))
  
  # Filter data
  bulk_data <- bulk_data[gene %in% keep_genes]
  
  # Save checkpoint
  message("[SAVING] Writing checkpoint...")
  qsave(bulk_data, BULK_ALL_FILE)
  message(sprintf("Saved to: %s", BULK_ALL_FILE))
}

# Final summary
message("\n========================================")
message("STEP 1 COMPLETE")
message("========================================")
message(sprintf("Total rows: %d", nrow(bulk_data)))
message(sprintf("Unique genes: %d", uniqueN(bulk_data$gene)))
message(sprintf("Unique tissues: %d", uniqueN(bulk_data$tissue)))
message(sprintf("Output: %s", BULK_ALL_FILE))
message("========================================\n")
