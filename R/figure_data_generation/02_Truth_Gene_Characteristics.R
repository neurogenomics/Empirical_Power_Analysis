# ==============================================================================
# STEP 02: Truth Gene Characteristics — From Full Ground Truth DE Results
# ==============================================================================
# Reads the Roussos DE results for each cell type, filters to significant genes
# (padj < 0.05, |LFC| > 1), and computes summary statistics.
#
# Output: Paper_Source_Data/truth_gene_characteristics.csv
#
# This replaces the hardcoded truth_chars_median and truth_chars_mean data frames
# in paper_final_draft.rmd.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 02: TRUTH GENE CHARACTERISTICS ###\n")

# ==============================================================================
# Read Roussos DE results for each cell type
# ==============================================================================

LFC_THRESHOLD <- 1.0
PADJ_THRESHOLD <- 0.05

all_rows <- list()

for (ct in CELL_TYPES) {
  de_path <- file.path(GT_BASE, ct, "Roussos_DE_results.csv")

  if (!file.exists(de_path)) {
    message(sprintf("  [SKIP] DE results not found: %s", de_path))
    next
  }

  message(sprintf("  Processing %s ...", ct))

  de <- read.csv(de_path)

  # Filter to significant truth genes
  sig <- de %>%
    filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(log2FoldChange) > LFC_THRESHOLD)

  n_genes <- nrow(sig)

  if (n_genes == 0) {
    message(sprintf("  [WARN] No significant genes for %s", ct))
    next
  }

  # Compute SNR = |LFC| / lfcSE (if lfcSE column exists)
  if ("lfcSE" %in% colnames(sig)) {
    snr <- abs(sig$log2FoldChange) / sig$lfcSE
  } else {
    snr <- rep(NA, n_genes)
  }

  row <- data.frame(
    Cell_Type = ct,
    N_Genes = n_genes,
    # Median stats
    Median_LFC = round(median(abs(sig$log2FoldChange)), 3),
    Mean_LFC = round(mean(abs(sig$log2FoldChange)), 3),
    Median_SNR = round(median(snr, na.rm = TRUE), 2),
    Mean_SNR = round(mean(snr, na.rm = TRUE), 2),
    Median_BaseMean = round(median(sig$baseMean), 2),
    Mean_BaseMean = round(mean(sig$baseMean), 2),
    Median_Dispersion = if ("dispersion" %in% colnames(sig)) round(median(sig$dispersion, na.rm = TRUE), 3) else NA,
    Mean_Dispersion = if ("dispersion" %in% colnames(sig)) round(mean(sig$dispersion, na.rm = TRUE), 3) else NA,
    stringsAsFactors = FALSE
  )

  all_rows[[length(all_rows) + 1]] <- row

  message(sprintf("  %s: %d genes, median |LFC| = %.3f, median baseMean = %.2f",
                  ct, n_genes, row$Median_LFC, row$Median_BaseMean))
}

# ==============================================================================
# Write
# ==============================================================================

truth_chars <- bind_rows(all_rows)

out_path <- file.path(OUTPUT_BASE, "truth_gene_characteristics.csv")
write.csv(truth_chars, out_path, row.names = FALSE)

message(sprintf("\n[DONE] Wrote %d rows to %s", nrow(truth_chars), out_path))
