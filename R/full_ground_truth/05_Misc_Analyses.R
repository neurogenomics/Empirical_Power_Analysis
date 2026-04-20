# ==============================================================================
# SCRIPT 05: Miscellaneous Analyses
# PURPOSE: Ad-hoc analyses that supplement the main cross-dataset diagnostics
# INPUT: DE results from Step 01 ({OUTPUT_BASE}/{CellType}/{Dataset}_DE_results.csv)
# OUTPUT: Summary tables and console output in {OUTPUT_BASE}/Misc/
# ==============================================================================
#
# Analysis A: Shared-Gene LFC Correlation
#   For each pair of datasets, restrict to genes significant in BOTH datasets
#   (padj < 0.05, |LFC| > 1), compute Pearson correlation of LFC values.
#   Tests whether Roussos and DeJager (both DLPFC) show the highest pairwise
#   correlation when gene universe differences are controlled for.
#
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/full_ground_truth/config.R",
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

# Source utils
utils_path <- file.path(SCRIPT_DIR, "utils.R")
if (file.exists(utils_path)) {
  source(utils_path)
} else {
  stop("utils.R not found at: ", utils_path)
}

suppressPackageStartupMessages({
  library(dplyr)
})

# Output directory
MISC_OUTPUT_DIR <- file.path(OUTPUT_BASE, "Misc")
if (!dir.exists(MISC_OUTPUT_DIR)) dir.create(MISC_OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# ANALYSIS A: Shared-Gene LFC Correlation
# ==============================================================================
# For each pair of datasets, restrict to genes significant in BOTH, then
# compute Pearson correlation of LFC. This controls for gene universe
# differences that confound raw concordance rates.
# ==============================================================================

log_section("ANALYSIS A: Shared-Gene Pairwise LFC Correlation")

all_pair_results <- list()

for (ct in CELL_TYPES) {

  log_info("Processing cell type: %s", ct)

  # --- Load DE results ---
  ct_output_dir <- file.path(OUTPUT_BASE, ct)
  all_results <- load_all_de_results(ct_output_dir)

  if (length(all_results) < 2) {
    log_warn("Fewer than 2 datasets for %s, skipping", ct)
    next
  }

  # --- Extract significant genes per dataset ---
  sig_list <- extract_significant_genes(
    all_results,
    padj_cutoff = FDR_THRESHOLD,
    lfc_cutoff  = LFC_THRESHOLD
  )

  dataset_names <- names(all_results)
  n_ds <- length(dataset_names)

  # --- Pairwise shared-gene correlations ---
  pair_rows <- list()

  for (i in seq_len(n_ds - 1)) {
    for (j in (i + 1):n_ds) {
      ds_a <- dataset_names[i]
      ds_b <- dataset_names[j]

      # Genes significant in BOTH datasets
      shared_sig <- intersect(sig_list[[ds_a]], sig_list[[ds_b]])
      n_shared <- length(shared_sig)

      if (n_shared < 3) {
        log_info("  %s vs %s: %d shared sig genes (too few for correlation)", ds_a, ds_b, n_shared)
        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          CellType  = ct,
          Dataset_A = ds_a,
          Dataset_B = ds_b,
          N_Shared_Sig = n_shared,
          Pearson_R = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      # Get LFC values for shared genes
      lfc_a <- all_results[[ds_a]][shared_sig, "log2FoldChange"]
      lfc_b <- all_results[[ds_b]][shared_sig, "log2FoldChange"]

      # Pearson correlation
      r <- cor(lfc_a, lfc_b, use = "complete.obs")

      log_info("  %s vs %s: %d shared sig genes, Pearson R = %.4f", ds_a, ds_b, n_shared, r)

      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        CellType  = ct,
        Dataset_A = ds_a,
        Dataset_B = ds_b,
        N_Shared_Sig = n_shared,
        Pearson_R = r,
        stringsAsFactors = FALSE
      )
    }
  }

  ct_pairs <- do.call(rbind, pair_rows)
  all_pair_results[[ct]] <- ct_pairs
}

# --- Combine across cell types ---
combined_pairs <- do.call(rbind, all_pair_results)

# --- Save ---
output_csv <- file.path(MISC_OUTPUT_DIR, "Shared_Gene_Pairwise_LFC_Correlations.csv")
write.csv(combined_pairs, output_csv, row.names = FALSE)
log_info("Saved pairwise correlations to: %s", output_csv)

# ==============================================================================
# SUMMARY: Is Roussos-DeJager the highest correlation?
# ==============================================================================

log_section("SUMMARY: Roussos-DeJager Shared-Gene Correlation Ranking")

for (ct in CELL_TYPES) {
  ct_data <- combined_pairs[combined_pairs$CellType == ct, ]
  ct_data <- ct_data[!is.na(ct_data$Pearson_R), ]

  if (nrow(ct_data) == 0) {
    log_info("[%s] No valid correlations", ct)
    next
  }

  # Sort by correlation descending
  ct_data <- ct_data[order(-ct_data$Pearson_R), ]

  # Find Roussos-DeJager pair
  rd_mask <- (ct_data$Dataset_A == "Roussos" & ct_data$Dataset_B == "DeJager") |
             (ct_data$Dataset_A == "DeJager" & ct_data$Dataset_B == "Roussos")
  rd_row <- ct_data[rd_mask, ]

  cat(sprintf("\n[%s] All pairwise shared-gene LFC correlations (ranked):\n", ct))
  for (row_i in seq_len(nrow(ct_data))) {
    row <- ct_data[row_i, ]
    marker <- if (any(rd_mask) && row_i == which(rd_mask[order(-ct_data$Pearson_R)])[1]) " <-- DLPFC pair" else ""
    # Simpler: just mark if this is the RD pair
    is_rd <- (row$Dataset_A == "Roussos" & row$Dataset_B == "DeJager") |
             (row$Dataset_A == "DeJager" & row$Dataset_B == "Roussos")
    marker <- if (is_rd) " <-- DLPFC pair" else ""
    cat(sprintf("  %d. %s vs %s: R = %.4f (n = %d shared sig genes)%s\n",
                row_i, row$Dataset_A, row$Dataset_B, row$Pearson_R, row$N_Shared_Sig, marker))
  }

  if (nrow(rd_row) > 0) {
    rank <- which(ct_data$Dataset_A == rd_row$Dataset_A[1] &
                  ct_data$Dataset_B == rd_row$Dataset_B[1])
    cat(sprintf("\n  Roussos-DeJager rank: %d / %d (R = %.4f, n = %d)\n",
                rank, nrow(ct_data), rd_row$Pearson_R[1], rd_row$N_Shared_Sig[1]))
  } else {
    cat("\n  Roussos-DeJager pair not found in results\n")
  }
}

cat("\n")
log_info("Analysis A complete.")

# ==============================================================================
# ANALYSIS B: Shared-Gene Concordance (controlling for gene universe)
# ==============================================================================
# The standard concordance metric (Diagnostic 5) uses:
#   concordance = |shared sig| / |Roussos sig|
# This penalises datasets with smaller gene universes because genes absent
# from their universe can never be detected as significant.
#
# Here we compute concordance restricted to genes present in BOTH datasets:
#   concordance_shared = |shared sig| / |Roussos sig genes in dataset universe|
#
# This controls for gene universe differences and tests whether DeJager
# (same brain region as Roussos: DLPFC) shows the highest concordance
# when the gene universe confound is removed.
# ==============================================================================

log_section("ANALYSIS B: Shared-Gene Concordance (Universe-Controlled)")

all_conc_results <- list()

for (ct in CELL_TYPES) {

  log_info("Processing cell type: %s", ct)

  ct_output_dir <- file.path(OUTPUT_BASE, ct)
  all_results <- load_all_de_results(ct_output_dir)

  if (length(all_results) < 2 || !REFERENCE_DATASET %in% names(all_results)) {
    log_warn("Missing reference or insufficient datasets for %s", ct)
    next
  }

  # Roussos significant genes
  ref_df <- all_results[[REFERENCE_DATASET]]
  ref_sig <- ref_df[!is.na(ref_df$padj) & ref_df$padj < FDR_THRESHOLD &
                      abs(ref_df$log2FoldChange) > LFC_THRESHOLD, ]
  ref_sig_genes <- rownames(ref_sig)

  conc_rows <- list()

  for (ds in names(all_results)) {
    if (ds == REFERENCE_DATASET) next

    ds_df <- all_results[[ds]]
    ds_sig <- ds_df[!is.na(ds_df$padj) & ds_df$padj < FDR_THRESHOLD &
                      abs(ds_df$log2FoldChange) > LFC_THRESHOLD, ]
    ds_sig_genes <- rownames(ds_sig)
    ds_universe <- rownames(ds_df)

    # Standard concordance: |shared sig| / |Roussos sig|
    shared_sig <- intersect(ref_sig_genes, ds_sig_genes)
    standard_conc <- length(shared_sig) / length(ref_sig_genes)

    # Universe-controlled: restrict Roussos sig to genes in dataset universe
    ref_sig_in_universe <- intersect(ref_sig_genes, ds_universe)
    shared_sig_in_universe <- intersect(ref_sig_in_universe, ds_sig_genes)

    if (length(ref_sig_in_universe) > 0) {
      controlled_conc <- length(shared_sig_in_universe) / length(ref_sig_in_universe)
    } else {
      controlled_conc <- NA_real_
    }

    # Also compute: what fraction of dataset's sig genes are Roussos sig?
    # (reverse concordance)
    if (length(ds_sig_genes) > 0) {
      reverse_conc <- length(shared_sig) / length(ds_sig_genes)
    } else {
      reverse_conc <- NA_real_
    }

    # Jaccard on shared universe
    union_sig_in_universe <- union(ref_sig_in_universe, ds_sig_genes)
    # Restrict union to genes in both universes for fair Jaccard
    ds_sig_in_ref_universe <- intersect(ds_sig_genes, rownames(ref_df))
    union_both <- union(ref_sig_in_universe, ds_sig_in_ref_universe)
    shared_both <- intersect(ref_sig_in_universe, ds_sig_in_ref_universe)
    jaccard <- if (length(union_both) > 0) length(shared_both) / length(union_both) else NA_real_

    log_info("  %s: standard=%.3f, controlled=%.3f, reverse=%.3f, jaccard=%.3f (ref_in_univ=%d, shared=%d)",
             ds, standard_conc, controlled_conc, reverse_conc, jaccard,
             length(ref_sig_in_universe), length(shared_sig_in_universe))

    conc_rows[[length(conc_rows) + 1]] <- data.frame(
      CellType = ct,
      Dataset = ds,
      N_Ref_Sig = length(ref_sig_genes),
      N_DS_Sig = length(ds_sig_genes),
      N_Ref_Sig_In_Universe = length(ref_sig_in_universe),
      N_Shared_Sig = length(shared_sig),
      Standard_Concordance = standard_conc,
      Controlled_Concordance = controlled_conc,
      Reverse_Concordance = reverse_conc,
      Jaccard = jaccard,
      stringsAsFactors = FALSE
    )
  }

  all_conc_results[[ct]] <- do.call(rbind, conc_rows)
}

combined_conc <- do.call(rbind, all_conc_results)

# Save
conc_csv <- file.path(MISC_OUTPUT_DIR, "Shared_Gene_Concordance_Universe_Controlled.csv")
write.csv(combined_conc, conc_csv, row.names = FALSE)
log_info("Saved universe-controlled concordance to: %s", conc_csv)

# ==============================================================================
# SUMMARY: Concordance rankings
# ==============================================================================

log_section("SUMMARY: Concordance Rankings (Standard vs Universe-Controlled)")

for (ct in CELL_TYPES) {
  ct_data <- combined_conc[combined_conc$CellType == ct, ]

  cat(sprintf("\n[%s] Concordance rankings:\n", ct))
  cat(sprintf("  %-12s  %8s  %10s  %10s  %8s  %s\n",
              "Dataset", "Standard", "Controlled", "Reverse", "Jaccard", "Ref_in_Univ"))

  # Sort by controlled concordance
  ct_data <- ct_data[order(-ct_data$Controlled_Concordance), ]

  for (row_i in seq_len(nrow(ct_data))) {
    row <- ct_data[row_i, ]
    marker <- if (row$Dataset == "DeJager") " <-- DLPFC" else ""
    cat(sprintf("  %-12s  %8.3f  %10.3f  %10.3f  %8.3f  %d/%d%s\n",
                row$Dataset, row$Standard_Concordance, row$Controlled_Concordance,
                row$Reverse_Concordance, row$Jaccard,
                row$N_Ref_Sig_In_Universe, row$N_Ref_Sig, marker))
  }
}

cat("\n")
log_info("Analysis B complete.")
