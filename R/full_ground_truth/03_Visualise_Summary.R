# ==============================================================================
# SCRIPT 03: Visualise Summary — Cross-Cell-Type Overview
# ==============================================================================
# Generates summary tables and overview figures spanning all cell types.
#
# Outputs:
#   1. Summary_All_CellTypes.csv           — sig gene counts, XIST stats
#   2. Summary_SigGene_Heatmap.png         — heatmap of sig genes (datasets x cell types)
#   3. Summary_XIST_Across_Datasets.png    — XIST LFC bar plot, faceted by cell type
#   4. Console completion report
#
# Prerequisites:
#   - Step 01: DE results in {OUTPUT_BASE}/{CellType}/{Dataset}_DE_results.csv
#   - Step 02: Diagnostic outputs (optional, not required by this script)
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
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

source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(viridis)
})

# Load shared utilities for prettified labels
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 3: CROSS-CELL-TYPE SUMMARY VISUALISATION")
log_info("Pipeline Version: %s", PIPELINE_VERSION)

# ==============================================================================
# 1. COLLECT SUMMARY DATA ACROSS ALL CELL TYPES
# ==============================================================================

log_section("Collecting DE results across cell types")

summary_rows <- list()

for (ct in CELL_TYPES) {
  paths <- get_paths(ct)

  if (!dir.exists(paths$output_dir)) {
    log_warn("Output directory missing for %s — skipping", ct)
    next
  }

  all_res <- load_all_de_results(paths$output_dir)

  for (ds in names(all_res)) {
    df <- all_res[[ds]]

    # Significant gene count
    sig <- df[!is.na(df$padj) & df$padj < FDR_THRESHOLD &
                abs(df$log2FoldChange) > LFC_THRESHOLD, ]
    n_sig <- nrow(sig)
    n_tested <- nrow(df)

    # XIST lookup
    df$gene_clean <- clean_gene_ids(df$gene)
    xist_row <- df[df$gene_clean == TARGET_GENE_ID, ]

    if (nrow(xist_row) > 0) {
      xist_row <- xist_row[1, ]
      xist_detected <- "yes"
      xist_padj <- xist_row$padj
      xist_lfc  <- xist_row$log2FoldChange
    } else {
      xist_detected <- "no"
      xist_padj <- NA_real_
      xist_lfc  <- NA_real_
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      CellType      = ct,
      Dataset       = ds,
      N_Significant = n_sig,
      N_Tested      = n_tested,
      XIST_Detected = xist_detected,
      XIST_padj     = xist_padj,
      XIST_LFC      = xist_lfc,
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_rows) == 0) {
  stop("No DE results found for any cell type. Run Step 01 first.")
}

summary_df <- do.call(rbind, summary_rows)
log_info("Collected results for %d dataset x cell-type combinations", nrow(summary_df))

# Save summary table
summary_csv <- file.path(OUTPUT_BASE, "Summary_All_CellTypes.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
log_info("Saved: %s", summary_csv)

# ==============================================================================
# 2. SIGNIFICANT GENE COUNT HEATMAP
# ==============================================================================

log_section("Generating significant gene count heatmap")

tryCatch({
  heat_df <- summary_df %>%
    dplyr::select(Dataset, CellType, N_Significant) %>%
    tidyr::pivot_wider(names_from = CellType, values_from = N_Significant) %>%
    as.data.frame()

  rownames(heat_df) <- prettify_label(heat_df$Dataset)
  heat_df$Dataset <- NULL
  colnames(heat_df) <- prettify_label(colnames(heat_df))
  heat_mat <- as.matrix(heat_df)

  heatmap_file <- file.path(OUTPUT_BASE, "Summary_SigGene_Heatmap.png")
  png(heatmap_file, width = 900, height = 700, res = 120)
  pheatmap::pheatmap(
    heat_mat,
    display_numbers = TRUE,
    number_format   = "%d",
    color           = viridis::viridis(50),
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    main            = sprintf("Significant Genes (padj < %g, |LFC| > %g)",
                              FDR_THRESHOLD, LFC_THRESHOLD),
    fontsize        = 14,
    fontsize_number = 12,
    fontsize_row    = 13,
    fontsize_col    = 13
  )
  dev.off()
  log_info("Saved: %s", heatmap_file)
}, error = function(e) {
  log_warn("Failed to generate heatmap: %s", e$message)
})

# ==============================================================================
# 3. XIST SUMMARY PLOT
# ==============================================================================

log_section("Generating XIST summary plot")

tryCatch({
  xist_df <- summary_df %>%
    dplyr::filter(XIST_Detected == "yes") %>%
    dplyr::mutate(
      Significant = ifelse(!is.na(XIST_padj) & XIST_padj < FDR_THRESHOLD,
                           "Significant", "Not Significant"),
      Dataset_Pretty = prettify_label(Dataset),
      CellType_Pretty = prettify_label(CellType)
    )
  xist_df$CellType_Pretty <- factor(xist_df$CellType_Pretty, 
                                     levels = prettify_label(CELL_TYPES))

  if (nrow(xist_df) == 0) {
    log_warn("XIST not detected in any dataset — skipping XIST plot")
  } else {
    p_xist <- ggplot(xist_df, aes(x = Dataset_Pretty, y = XIST_LFC, fill = Significant)) +
      geom_col(width = 0.7) +
      facet_wrap(~ CellType_Pretty, scales = "free_x") +
      scale_fill_manual(values = c("Significant" = "#2166AC",
                                   "Not Significant" = "#B2182B")) +
      labs(
        title    = sprintf("%s log2FoldChange Across Datasets", TARGET_GENE_NAME),
        x = NULL,
        y = "log2FoldChange"
      ) +
      theme_publication(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
      )

    xist_file <- file.path(OUTPUT_BASE, "Summary_XIST_Across_Datasets.png")
    ggsave(xist_file, p_xist, width = 12, height = 6, dpi = 150)
    log_info("Saved: %s", xist_file)
  }
}, error = function(e) {
  log_warn("Failed to generate XIST plot: %s", e$message)
})

# ==============================================================================
# 4. PIPELINE COMPLETION REPORT
# ==============================================================================

log_section("PIPELINE COMPLETION REPORT")

# Datasets per cell type
for (ct in CELL_TYPES) {
  ct_rows <- summary_df[summary_df$CellType == ct, ]
  log_info("  %s: %d datasets processed, %d total sig genes",
           ct, nrow(ct_rows), sum(ct_rows$N_Significant))
}

# Check for missing datasets
expected <- expand.grid(CellType = CELL_TYPES, Dataset = RUN_DATASETS,
                        stringsAsFactors = FALSE)
observed  <- summary_df[, c("CellType", "Dataset")]
missing   <- dplyr::anti_join(expected, observed, by = c("CellType", "Dataset"))

if (nrow(missing) > 0) {
  log_warn("Missing %d dataset x cell-type combinations:", nrow(missing))
  for (i in seq_len(nrow(missing))) {
    log_warn("  - %s / %s", missing$CellType[i], missing$Dataset[i])
  }
} else {
  log_info("All expected dataset x cell-type combinations present")
}

# Overall status
log_info("Total sig genes across pipeline: %d", sum(summary_df$N_Significant))
log_info("Summary CSV: %s", summary_csv)
log_info("Full Ground Truth Pipeline — Step 3 complete")
