# ==============================================================================
# SCRIPT: DE Density Plots - Cross-Cell-Type Metric Distributions
# PURPOSE: Generate density plots of DE result metrics (|LFC|, SE, SNR,
#          baseMean, dispersion) overlaid across cell types, with optional
#          coding/non-coding biotype stratification.
# INPUT: Full ground truth DE results from Datalink/Full_Ground_Truth/
#        DDS checkpoints for dispersion extraction
#        BioMart for gene biotype annotation (cached)
# OUTPUT: 10 PNG density plots in Datalink/DE_Diagnostics/Density_Plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(qs)
  library(DESeq2)
  library(scales)
  library(RColorBrewer)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Load config - skip if already loaded (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR") || !grepl("de_diagnostics", SCRIPT_DIR)) {
  possible_config_paths <- c(
    "R/de_diagnostics/config.R",
    "config.R",
    "../config.R"
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
    stop("Could not find config.R. Run from project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 1: DE DENSITY PLOTS")
log_info("DE Dataset: %s", DE_DATASET)
log_info("Filters: padj < %g, |LFC| > %g", PADJ_CUTOFF, LFC_MIN_HIST)
log_info("Metrics: %s", paste(sapply(DENSITY_METRICS, `[[`, "suffix"), collapse = ", "))

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

output_dir <- get_density_output_dir()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  log_info("Created output directory: %s", output_dir)
}

# Check skip-if-exists
if (!FORCE_RECALC) {
  existing_plots <- list.files(output_dir, pattern = "\\.png$")
  expected_n <- length(DENSITY_METRICS) * 2  # aggregated + stratified per metric
  if (length(existing_plots) >= expected_n) {
    log_info("[SKIP] Density plots already exist (%d files). Set FORCE_RECALC=TRUE to regenerate.",
             length(existing_plots))
    # Still return success for pipeline runner
    invisible(NULL)
    # Use a flag to skip the rest
    .skip_step1 <- TRUE
  } else {
    .skip_step1 <- FALSE
  }
} else {
  .skip_step1 <- FALSE
}

if (!exists(".skip_step1") || !.skip_step1) {

# ==============================================================================
# 3. LOAD DE RESULTS FOR ALL CELL TYPES
# ==============================================================================

log_section("Loading DE results")

de_list <- lapply(CELL_TYPES, function(ct) {
  load_de_results(ct, padj_cutoff = PADJ_CUTOFF, lfc_min = LFC_MIN_HIST)
})
names(de_list) <- CELL_TYPES

# Remove NULLs
de_list <- de_list[!sapply(de_list, is.null)]

if (length(de_list) == 0) {
  stop("No DE results loaded for any cell type. Check paths and filters.")
}

combined_df <- do.call(rbind, de_list)
log_info("Combined: %d significant genes across %d cell types",
         nrow(combined_df), length(de_list))

# ==============================================================================
# 4. LOAD DISPERSIONS FROM DDS CHECKPOINTS
# ==============================================================================

log_section("Loading dispersions from DDS checkpoints")

for (ct in names(de_list)) {
  disp_vec <- load_dispersions(ct)
  if (!is.null(disp_vec)) {
    # Match dispersions to DE results
    ct_mask <- combined_df$Cell_Type == ct
    ct_genes <- combined_df$gene_id[ct_mask]
    common <- intersect(ct_genes, names(disp_vec))
    combined_df$dispersion[ct_mask] <- NA_real_
    if (length(common) > 0) {
      match_idx <- match(combined_df$gene_id[ct_mask], names(disp_vec))
      combined_df$dispersion[ct_mask] <- disp_vec[match_idx]
      log_info("  %s: matched dispersions for %d/%d genes",
               ct, sum(!is.na(disp_vec[match_idx])), sum(ct_mask))
    }
  } else {
    combined_df$dispersion[combined_df$Cell_Type == ct] <- NA_real_
  }
}

# ==============================================================================
# 5. BIOMART ANNOTATION
# ==============================================================================

log_section("Annotating gene biotypes")

cache_path <- file.path(OUTPUT_BASE, BIOMART_CACHE_FILE)

# Ensure parent directory exists for cache
if (!dir.exists(dirname(cache_path))) {
  dir.create(dirname(cache_path), recursive = TRUE)
}

all_gene_ids <- unique(combined_df$gene_id)
biotype_df <- get_biotype_annotations(all_gene_ids, cache_path)

if (!is.null(biotype_df) && nrow(biotype_df) > 0) {
  combined_df <- combined_df %>%
    dplyr::left_join(
      biotype_df[, c("ensembl_gene_id", "Biotype_Class")],
      by = c("gene_id" = "ensembl_gene_id")
    )

  # Fill NAs as "Unknown"
  combined_df$Biotype_Class[is.na(combined_df$Biotype_Class)] <- "Unknown"

  log_info("Biotype breakdown: %s",
           paste(sprintf("%s=%d", names(table(combined_df$Biotype_Class)),
                         as.integer(table(combined_df$Biotype_Class))),
                 collapse = ", "))
} else {
  log_warn("BioMart annotation failed — stratified plots will be skipped")
  combined_df$Biotype_Class <- NA_character_
}

# ==============================================================================
# 6. COMPUTE DERIVED METRICS
# ==============================================================================

log_section("Computing derived metrics")

combined_df <- combined_df %>%
  dplyr::mutate(
    Abs_LFC = abs(log2FoldChange),
    SNR = Abs_LFC / lfcSE,
    Log10_SNR = log10(SNR + 1e-6),
    Log10_BaseMean = log10(baseMean + 1),
    Log10_Dispersion = log10(dispersion + 1e-6)
  )

# Report NA counts for dispersion
na_disp <- sum(is.na(combined_df$Log10_Dispersion) |
                 is.infinite(combined_df$Log10_Dispersion))
if (na_disp > 0) {
  log_info("  %s genes have NA/Inf dispersion (will be excluded from dispersion plots)",
           format(na_disp, big.mark = ","))
}

# ==============================================================================
# 7. GENERATE DENSITY PLOTS
# ==============================================================================

log_section("Generating density plots")

for (m in DENSITY_METRICS) {
  log_info("  Plotting: %s", m$suffix)

  # --- Resolve x-axis limits ---
  # If lim is explicitly set, use it. Otherwise auto-scale from data.
  if (!is.null(m$lim)) {
    xlim_val <- m$lim
    log_info("    X-limits (manual): [%.2f, %.2f]", xlim_val[1], xlim_val[2])
  } else {
    xlim_val <- compute_auto_limits(
      combined_df, m$col,
      padding   = DENSITY_PADDING,
      lim_lower = m$lim_lower,  # NULL if not set (auto both sides)
      lim_upper = m$lim_upper   # NULL if not set (auto both sides)
    )
    if (!is.null(xlim_val)) {
      log_info("    X-limits (auto):   [%.2f, %.2f]", xlim_val[1], xlim_val[2])
    }
  }

  # --- Aggregated plot ---
  p_agg <- plot_density_aggregated(combined_df, m$col, m$suffix, m$lab, xlim_val)
  if (!is.null(p_agg)) {
    agg_file <- file.path(output_dir, paste0("Density_", m$suffix, ".png"))
    ggplot2::ggsave(agg_file, p_agg,
                    width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    log_info("    Saved: %s", basename(agg_file))
  } else {
    log_warn("    No data for aggregated %s plot", m$suffix)
  }

  # --- Stratified plot ---
  if (!all(is.na(combined_df$Biotype_Class))) {
    p_strat <- plot_density_stratified(combined_df, m$col, m$suffix, m$lab, xlim_val)
    if (!is.null(p_strat)) {
      strat_file <- file.path(output_dir, paste0("Density_", m$suffix, "_Stratified.png"))
      ggplot2::ggsave(strat_file, p_strat,
                      width = PLOT_WIDTH_STRATIFIED, height = PLOT_HEIGHT_STRATIFIED,
                      dpi = PLOT_DPI)
      log_info("    Saved: %s", basename(strat_file))
    } else {
      log_warn("    No data for stratified %s plot", m$suffix)
    }
  }
}

# ==============================================================================
# 8. SUMMARY TABLE
# ==============================================================================

log_section("Saving summary table")

summary_df <- combined_df %>%
  dplyr::group_by(Cell_Type) %>%
  dplyr::summarise(
    N_Genes = dplyr::n(),
    Median_AbsLFC = median(Abs_LFC, na.rm = TRUE),
    Median_SE = median(lfcSE, na.rm = TRUE),
    Median_SNR = median(SNR, na.rm = TRUE),
    Median_BaseMean = median(baseMean, na.rm = TRUE),
    Median_Dispersion = median(dispersion, na.rm = TRUE),
    Pct_Coding = 100 * mean(Biotype_Class == "Coding", na.rm = TRUE),
    .groups = "drop"
  )

summary_file <- file.path(output_dir, "Density_Summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
log_info("Saved: %s", basename(summary_file))

# Print summary
for (i in seq_len(nrow(summary_df))) {
  log_info("  %s: %d genes, median |LFC|=%.2f, median SE=%.3f, %.0f%% coding",
           summary_df$Cell_Type[i], summary_df$N_Genes[i],
           summary_df$Median_AbsLFC[i], summary_df$Median_SE[i],
           summary_df$Pct_Coding[i])
}

log_info("Step 1 complete: %d density plots generated",
         length(list.files(output_dir, pattern = "\\.png$")))

} # end skip check
