# ==============================================================================
# STEP 3: ANNOTATE DEGs
# ==============================================================================
# PURPOSE: Add gene symbols and descriptions via BioMart
# INPUT: Datalink/Bulk_Analysis/Filtered_Consistent_Genes.csv
# OUTPUT: Datalink/Bulk_Analysis/Filtered_Consistent_Genes_with_Descriptions.csv
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
message("STEP 3: ANNOTATE DEGs")
message("========================================\n")

# Check if checkpoint exists and is valid
if (file.exists(ANNOTATED_FILE) && !FORCE_RECALC) {
  message("[CHECKPOINT] Loading cached data...")
  annotated_data <- fread(ANNOTATED_FILE)
  
  # Validate checkpoint
  is_valid <- validate_checkpoint(annotated_data, type = "annotated")
  
  if (!is_valid) {
    message("[CACHE INVALID] Regenerating...")
    annotated_data <- NULL
  } else {
    message("[CACHE VALID] Using cached data")
    message(sprintf("Loaded %d annotated genes", nrow(annotated_data)))
  }
}

# Generate fresh data if needed
if (!exists("annotated_data") || is.null(annotated_data)) {
  message("[COMPUTING] Loading filtered genes from Step 2...")
  
  # Load data from Step 2
  if (!file.exists(FILTERED_FILE)) {
    stop(sprintf("Input file not found: %s\nRun Step 2 first.", FILTERED_FILE))
  }
  
  bulk_data <- fread(FILTERED_FILE)
  message(sprintf("Loaded %d genes", nrow(bulk_data)))
  
  # Query BioMart for annotations
  message("\n[BIOMART] Fetching gene annotations...")
  gene_info <- query_biomart_safe(
    gene_ids = bulk_data$gene,
    dataset = BIOMART_DATASET,
    attributes = BIOMART_ATTRIBUTES
  )
  
  # Merge annotations with bulk data
  message("\n[MERGING] Combining annotations with bulk data...")
  annotated_data <- bulk_data %>%
    left_join(gene_info, by = c("gene" = "ensembl_gene_id"))
  
  # Convert back to data.table for consistent syntax
  annotated_data <- as.data.table(annotated_data)
  
  # Clean descriptions (remove BioMart source annotations)
  annotated_data$description <- stringr::str_replace(
    annotated_data$description, 
    " \\[Source:.*\\]", 
    ""
  )
  
  # Fill NAs for cleaner output
  annotated_data$description[is.na(annotated_data$description)] <- "No description available"
  annotated_data$external_gene_name[is.na(annotated_data$external_gene_name)] <- "Unknown"
  
  # Reorder columns (gene info first, then stats)
  cols <- colnames(annotated_data)
  new_order <- c(
    "gene", 
    "external_gene_name", 
    "description", 
    setdiff(cols, c("gene", "external_gene_name", "description"))
  )
  annotated_data <- annotated_data[, ..new_order]
  
  # Sanity checks on annotated data
  annotated_checks <- list(
    "Non-empty" = function(data) {
      if (nrow(data) > 0) TRUE
      else "No genes in annotated data"
    },
    "Most genes have symbols" = function(data) {
      unknown_pct <- mean(data$external_gene_name == "Unknown") * 100
      if (unknown_pct < 50) TRUE
      else sprintf("%.1f%% genes have no symbol (expect <50%%)", unknown_pct)
    },
    "Most genes have descriptions" = function(data) {
      no_desc_pct <- mean(data$description == "No description available") * 100
      if (no_desc_pct < 50) TRUE
      else sprintf("%.1f%% genes have no description (expect <50%%)", no_desc_pct)
    }
  )
  
  sanity_passed <- sanity_check(annotated_data, "Annotated genes", annotated_checks)
  if (!sanity_passed) {
    warning("Annotated data failed sanity checks. BioMart may have returned incomplete results.")
  }
  
  # Save checkpoint
  message("\n[SAVING] Writing annotated gene list...")
  fwrite(annotated_data, ANNOTATED_FILE, quote = TRUE)
  message(sprintf("Saved to: %s", ANNOTATED_FILE))
}

# Final summary
message("\n========================================")
message("STEP 3 COMPLETE")
message("========================================")
message(sprintf("Annotated genes: %d", nrow(annotated_data)))
message(sprintf("Genes with symbols: %d (%.1f%%)", 
                sum(annotated_data$external_gene_name != "Unknown"),
                mean(annotated_data$external_gene_name != "Unknown") * 100))
message(sprintf("Genes with descriptions: %d (%.1f%%)",
                sum(annotated_data$description != "No description available"),
                mean(annotated_data$description != "No description available") * 100))
message(sprintf("Output: %s", ANNOTATED_FILE))
message("========================================\n")

# Show first few results
message("First 3 annotated genes:")
print(head(annotated_data[, c("gene", "external_gene_name", "description", "log2FoldChange")], 3))
