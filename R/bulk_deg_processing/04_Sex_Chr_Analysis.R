# ==============================================================================
# STEP 4: SEX CHROMOSOME ANALYSIS
# ==============================================================================
# PURPOSE: Analyze sex chromosome representation by LFC interval
# INPUT: Datalink/Bulk_Analysis/Bulk_DE_genes_All.qs
# OUTPUT: Datalink/Bulk_Analysis/Sex_Chr_Proportions_Intervals.csv
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
message("STEP 4: SEX CHROMOSOME ANALYSIS")
message("========================================\n")

# Check if output exists and skip if not forcing recalc
if (file.exists(SEX_CHR_FILE) && !FORCE_RECALC) {
  message("[CHECKPOINT] Output already exists. Use FORCE_RECALC=TRUE to regenerate.")
  final_stats <- fread(SEX_CHR_FILE)
  message(sprintf("Loaded existing results: %s", SEX_CHR_FILE))
} else {
  message("[COMPUTING] Loading bulk data from Step 1...")
  
  # Load data from Step 1
  if (!file.exists(BULK_ALL_FILE)) {
    stop(sprintf("Input file not found: %s\nRun Step 1 first.", BULK_ALL_FILE))
  }
  
  dt <- qread(BULK_ALL_FILE)
  setDT(dt)
  message(sprintf("Loaded %d rows", nrow(dt)))
  
  # Clean gene IDs and aggregate per gene
  message("\n[PROCESSING] Aggregating per gene (mean absolute LFC)...")
  dt[, clean_gene := clean_gene_ids(gene)]
  gene_summary <- dt[, .(mean_abs_lfc = mean(abs(effsize), na.rm = TRUE)), by = clean_gene]
  message(sprintf("Analyzed %d unique genes", nrow(gene_summary)))
  
  # Annotate chromosomes
  message("\n[ANNOTATING] Mapping genes to chromosomes...")
  chr_annotations <- annotate_chromosomes(gene_summary$clean_gene)
  gene_summary <- merge(gene_summary, chr_annotations, by.x = "clean_gene", by.y = "gene")
  
  # Assign LFC intervals
  message("\n[BINNING] Assigning genes to LFC intervals...")
  gene_summary$interval <- bin_by_lfc_interval(gene_summary$mean_abs_lfc, LFC_INTERVALS)
  
  # Calculate proportions by interval
  message("\n[CALCULATING] Computing sex chromosome proportions...")
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
  
  # Sanity checks on results
  sex_chr_checks <- list(
    "Non-empty" = function(data) {
      if (nrow(data) > 0) TRUE
      else "No intervals computed"
    },
    "All intervals present" = function(data) {
      if (nrow(data) == length(LFC_INTERVALS)) TRUE
      else sprintf("Expected %d intervals, got %d", length(LFC_INTERVALS), nrow(data))
    },
    "Valid proportions" = function(data) {
      if (all(data$Proportion >= 0 & data$Proportion <= 1)) TRUE
      else "Proportions out of [0,1] range"
    },
    "Reasonable sex chr proportion" = function(data) {
      # Human genome: ~5% of genes on sex chromosomes
      # Expect higher in sex-biased DEGs, but not >50%
      max_prop <- max(data$Proportion)
      if (max_prop <= 0.5) TRUE
      else sprintf("Max proportion = %.1f%% (unusually high)", max_prop * 100)
    }
  )
  
  sanity_passed <- sanity_check(final_stats, "Sex chromosome analysis", sex_chr_checks)
  if (!sanity_passed) {
    warning("Sex chromosome analysis failed sanity checks. Review results carefully.")
  }
  
  # Save results
  message("\n[SAVING] Writing sex chromosome analysis...")
  fwrite(final_stats, SEX_CHR_FILE)
  message(sprintf("Saved to: %s", SEX_CHR_FILE))
}

# Final summary
message("\n========================================")
message("STEP 4 COMPLETE")
message("========================================")
message("Sex chromosome proportions by LFC interval:")
print(final_stats)
message(sprintf("\nOutput: %s", SEX_CHR_FILE))
message("========================================\n")
