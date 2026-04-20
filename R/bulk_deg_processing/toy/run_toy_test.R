# ==============================================================================
# TOY INTEGRATION TEST: Bulk DEG Processing Pipeline
# ==============================================================================
# Minimal end-to-end test using first 1000 genes from raw data
# ==============================================================================

# ==============================================================================
# ==============================================================================
# ==============================================================================
# SETUP
# ==============================================================================

log_section("TOY INTEGRATION TEST - SETUP")

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/bulk_deg_processing/config.R",  # From project root
    "../config.R",                      # From toy/ directory
    "config.R"                          # From script directory
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
    log_error("Could not find config.R")
    finalize_logging(status = "FAILED")
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

# Source utilities
source(file.path(SCRIPT_DIR, "utils.R"))

log_info("Config loaded successfully")
log_info("Project root: %s", PROJECT_ROOT)
log_info("Script directory: %s", SCRIPT_DIR)

# Override output paths for toy test
BULK_ALL_FILE <- file.path(TOY_OUTPUT_DIR, "Bulk_DE_genes_All.qs")
FILTERED_FILE <- file.path(TOY_OUTPUT_DIR, "Filtered_Consistent_Genes.csv")
ANNOTATED_FILE <- file.path(TOY_OUTPUT_DIR, "Filtered_Consistent_Genes_with_Descriptions.csv")
SEX_CHR_FILE <- file.path(TOY_OUTPUT_DIR, "Sex_Chr_Proportions_Intervals.csv")

log_info("Toy outputs will be saved to: %s", TOY_OUTPUT_DIR)

# ==============================================================================
# STEP 1: EXTRACT BULK DEGs (FIRST 1000 GENES)
# ==============================================================================

log_section("STEP 1: EXTRACT BULK DEGs (TOY DATA)")

log_info("Loading raw bulk data...")
dt_full <- fread(INPUT_FILE)
log_info("Loaded %d rows from raw file", nrow(dt_full))

# Get first 1000 unique genes
unique_genes <- unique(dt_full$gene)
toy_genes <- head(unique_genes, 1000)
log_info("Selected first %d genes for toy test", length(toy_genes))

# Filter to toy genes
dt <- dt_full[gene %in% toy_genes]
log_info("Toy dataset: %d rows, %d unique genes, %d unique tissues",
         nrow(dt), uniqueN(dt$gene), uniqueN(dt$tissue))

# Save toy checkpoint
qsave(dt, BULK_ALL_FILE)
log_info("Saved toy data to: %s", BULK_ALL_FILE)

# ==============================================================================
# STEP 2: COLLECT CONSISTENT DEGs
# ==============================================================================

log_section("STEP 2: COLLECT CONSISTENT DEGs")

log_info("Applying consistency filter...")
final_list <- filter_consistent_genes(
  dt,
  p_threshold = P_VALUE_THRESHOLD,
  consistency_threshold = TISSUE_CONSISTENCY_THRESHOLD
)

log_info("Filtered genes: %d", nrow(final_list))

if (nrow(final_list) == 0) {
  log_warn("No genes passed filter. Consider lowering thresholds for toy test.")
} else {
  log_info("Mean |LFC|: %.2f", mean(abs(final_list$log2FoldChange)))
  log_info("Median |LFC|: %.2f", median(abs(final_list$log2FoldChange)))
}

# Save filtered list
fwrite(final_list, FILTERED_FILE, quote = TRUE)
log_info("Saved filtered genes to: %s", FILTERED_FILE)

# ==============================================================================
# STEP 3: ANNOTATE DEGs
# ==============================================================================

log_section("STEP 3: ANNOTATE DEGs")

if (nrow(final_list) == 0) {
  log_warn("Skipping annotation (no genes to annotate)")
  annotated_data <- final_list
} else {
  # Limit to first 50 genes for toy test (BioMart can be slow)
  genes_to_annotate <- head(final_list$gene, 50)
  log_info("Annotating first %d genes (toy test limit)", length(genes_to_annotate))
  
  tryCatch({
    gene_info <- query_biomart_safe(genes_to_annotate)
    log_info("Retrieved annotations for %d genes", nrow(gene_info))
    
    # Merge with filtered list
    annotated_data <- final_list[gene %in% genes_to_annotate] %>%
      left_join(gene_info, by = c("gene" = "ensembl_gene_id"))
    
    # Clean descriptions
    annotated_data$description <- stringr::str_replace(
      annotated_data$description,
      " \\[Source:.*\\]",
      ""
    )
    
    # Fill NAs
    annotated_data$description[is.na(annotated_data$description)] <- "No description available"
    annotated_data$external_gene_name[is.na(annotated_data$external_gene_name)] <- "Unknown"
    
    # Reorder columns
    cols <- colnames(annotated_data)
    new_order <- c(
      "gene",
      "external_gene_name",
      "description",
      setdiff(cols, c("gene", "external_gene_name", "description"))
    )
    annotated_data <- annotated_data[, new_order]
    
    log_info("Annotated %d genes", nrow(annotated_data))
    
  }, error = function(e) {
    log_error("BioMart query failed: %s", e$message)
    log_warn("Continuing without annotations")
    annotated_data <- final_list[gene %in% genes_to_annotate]
    annotated_data$external_gene_name <- "Unknown"
    annotated_data$description <- "Annotation failed"
  })
}

# Save annotated list
fwrite(annotated_data, ANNOTATED_FILE, quote = TRUE)
log_info("Saved annotated genes to: %s", ANNOTATED_FILE)

# ==============================================================================
# STEP 4: SEX CHROMOSOME ANALYSIS
# ==============================================================================

log_section("STEP 4: SEX CHROMOSOME ANALYSIS")

# Clean gene IDs and aggregate
dt[, clean_gene := clean_gene_ids(gene)]
gene_summary <- dt[, .(mean_abs_lfc = mean(abs(effsize), na.rm = TRUE)), by = clean_gene]
log_info("Aggregated %d unique genes", nrow(gene_summary))

# Annotate chromosomes
chr_annotations <- annotate_chromosomes(gene_summary$clean_gene)
gene_summary <- merge(gene_summary, chr_annotations, by.x = "clean_gene", by.y = "gene")

# Assign intervals
gene_summary$interval <- bin_by_lfc_interval(gene_summary$mean_abs_lfc, LFC_INTERVALS)

# Calculate proportions
stats_list <- lapply(names(LFC_INTERVALS), function(label) {
  subset_genes <- gene_summary[interval == label]
  
  n_total <- nrow(subset_genes)
  n_sex <- sum(subset_genes$is_sex_chr, na.rm = TRUE)
  prop <- if (n_total > 0) n_sex / n_total else 0
  
  return(data.table(
    Interval = label,
    Total_Genes = n_total,
    Sex_Chr_Genes = n_sex,
    Proportion = round(prop, 4)
  ))
})

final_stats <- rbindlist(stats_list)

log_info("Sex chromosome proportions:")
print(final_stats)

# Save results
fwrite(final_stats, SEX_CHR_FILE)
log_info("Saved sex chromosome analysis to: %s", SEX_CHR_FILE)

# ==============================================================================
# VALIDATION AND SUMMARY
# ==============================================================================

log_section("TOY TEST VALIDATION")

# Check all outputs exist
outputs <- c(BULK_ALL_FILE, FILTERED_FILE, ANNOTATED_FILE, SEX_CHR_FILE)
all_exist <- all(file.exists(outputs))

if (all_exist) {
  log_info("✓ All output files generated successfully")
} else {
  log_error("✗ Some output files missing:")
  for (f in outputs) {
    if (!file.exists(f)) {
      log_error("  Missing: %s", f)
    }
  }
}

# Manual inspection prompts
log_section("MANUAL INSPECTION REQUIRED")
log_info("Please review the following outputs:")
log_info("1. %s", FILTERED_FILE)
log_info("   - Check: Reasonable number of genes, valid LFC values")
log_info("2. %s", ANNOTATED_FILE)
log_info("   - Check: Gene symbols present, descriptions meaningful")
log_info("3. %s", SEX_CHR_FILE)
log_info("   - Check: Proportions in reasonable range (0-50%%)")

# Final status
if (all_exist) {
  log_info("✓ TOY TEST COMPLETE - AWAITING USER APPROVAL")
  finalize_logging(status = "SUCCESS")
} else {
  log_error("✗ TOY TEST FAILED - MISSING OUTPUTS")
  finalize_logging(status = "FAILED")
}
