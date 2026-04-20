# ==============================================================================
# SCRIPT: Shared Gene Diagnostics — Paired Cross-Cell-Type Comparisons
# PURPOSE: Formalise the paired gene-by-gene comparisons of the 23 shared truth
#          genes across cell types. Computes paired Wilcoxon tests, Pearson
#          correlations, and summary statistics for 5 metrics (dispersion, |LFC|,
#          lfcSE, SNR, baseMean). Also produces summary statistics for each CT's
#          full truth set. Outputs CSVs + scatter plots.
#
#          This provides the reproducible evidence chain for the Subtype
#          Heterogeneity Hypothesis: the same genes have similar intrinsic
#          characteristics across cell types, so power differences must come
#          from a cell-type-level mechanism (subtype composition).
#
# INPUT:  - Shared truth genes from Power Curve Fitting pipeline
#           (Datalink/Power_Curve_Fitting/Fitting_Data/shared_truth_genes.qs)
#         - Full Ground Truth DE results (Datalink/Full_Ground_Truth/{CT}/)
#         - Full Ground Truth DDS checkpoints (for dispersion)
# OUTPUT: - Shared_Gene_Paired_Tests.csv (test statistics for all pairs × metrics)
#         - Shared_Gene_Comparison_Matrix.csv (gene × CT × metric wide table)
#         - Full_Truth_Set_Summary.csv (per-CT medians for full significant genes)
#         - 15 scatter plots (5 metrics × 3 CT pairs)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(DESeq2)
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

log_section("STEP 5: SHARED GENE DIAGNOSTICS")

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

output_dir <- get_shared_gene_output_dir()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  log_info("Created output directory: %s", output_dir)
}

# ==============================================================================
# 3. LOAD SHARED TRUTH GENES
# ==============================================================================

log_info("Loading shared truth genes from: %s", basename(SHARED_GENES_PATH))

if (!file.exists(SHARED_GENES_PATH)) {
  stop(sprintf("Shared truth genes file not found: %s\n  Run power_curve_fitting Step 01 first.",
               SHARED_GENES_PATH))
}

shared_genes <- qs::qread(SHARED_GENES_PATH)
shared_genes <- clean_gene_ids(shared_genes)
n_shared <- length(shared_genes)
log_info("Loaded %d shared truth genes", n_shared)

# ==============================================================================
# 4. LOAD DE RESULTS AND DISPERSIONS FOR ALL CELL TYPES
# ==============================================================================

log_info("Loading DE results and dispersions for %d cell types...",
         length(CELL_TYPES))

# We load ALL tested genes (lfc_min = 0) so we can look up shared genes
# regardless of their significance in individual CTs
ct_data <- list()

for (ct in CELL_TYPES) {
  log_info("  Processing %s...", ct)

  # Load unfiltered DE results (all tested genes)
  de_all <- load_de_results(ct, padj_cutoff = 1.0, lfc_min = 0)
  if (is.null(de_all)) {
    stop(sprintf("Failed to load DE results for %s", ct))
  }

  # Load dispersions
  disp_vec <- load_dispersions(ct)
  if (is.null(disp_vec)) {
    stop(sprintf("Failed to load dispersions for %s", ct))
  }

  # Build per-gene table for shared genes
  shared_in_de <- intersect(shared_genes, de_all$gene_id)
  if (length(shared_in_de) < n_shared) {
    log_warn("  %s: only %d/%d shared genes found in DE results",
             ct, length(shared_in_de), n_shared)
  }

  de_shared <- de_all[de_all$gene_id %in% shared_genes, ]
  de_shared$dispersion <- disp_vec[de_shared$gene_id]
  de_shared$abs_lfc <- abs(de_shared$log2FoldChange)
  de_shared$snr <- de_shared$abs_lfc / de_shared$lfcSE

  ct_data[[ct]] <- list(
    de_shared = de_shared,
    de_all = de_all,
    disp_vec = disp_vec
  )

  log_info("  %s: %d shared genes loaded (dispersion available for %d)",
           ct, nrow(de_shared),
           sum(!is.na(de_shared$dispersion)))
}

# ==============================================================================
# 5. BUILD COMPARISON MATRIX (gene × CT × metric)
# ==============================================================================

log_info("Building comparison matrix...")

comparison_rows <- list()

for (ct in CELL_TYPES) {
  de <- ct_data[[ct]]$de_shared
  for (i in seq_len(nrow(de))) {
    comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
      gene_id     = de$gene_id[i],
      cell_type   = ct,
      baseMean    = de$baseMean[i],
      log2FC      = de$log2FoldChange[i],
      abs_lfc     = de$abs_lfc[i],
      lfcSE       = de$lfcSE[i],
      snr         = de$snr[i],
      dispersion  = de$dispersion[i],
      padj        = de$padj[i],
      stringsAsFactors = FALSE
    )
  }
}

comparison_long <- do.call(rbind, comparison_rows)
log_info("Comparison matrix: %d rows (%d genes × %d cell types)",
         nrow(comparison_long), n_shared, length(CELL_TYPES))

# Save long-format comparison matrix
comparison_path <- file.path(output_dir, "Shared_Gene_Comparison_Matrix.csv")
write.csv(comparison_long, comparison_path, row.names = FALSE)
log_info("Saved: %s", basename(comparison_path))

# Also create wide-format for easy inspection (one row per gene, columns per CT)
comparison_wide <- comparison_long %>%
  select(gene_id, cell_type, baseMean, abs_lfc, lfcSE, snr, dispersion) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = c(baseMean, abs_lfc, lfcSE, snr, dispersion),
    names_sep = "_"
  )

wide_path <- file.path(output_dir, "Shared_Gene_Comparison_Wide.csv")
write.csv(comparison_wide, wide_path, row.names = FALSE)
log_info("Saved: %s", basename(wide_path))

# ==============================================================================
# 6. PAIRED WILCOXON TESTS
# ==============================================================================

log_info("Running paired Wilcoxon signed-rank tests...")

test_results <- list()

for (pair in SHARED_GENE_CT_PAIRS) {
  ct_a <- pair[1]
  ct_b <- pair[2]
  pair_label <- paste(ct_a, "vs", ct_b)
  log_info("  Pair: %s", pair_label)

  # Get data for both CTs, aligned by gene_id
  de_a <- ct_data[[ct_a]]$de_shared
  de_b <- ct_data[[ct_b]]$de_shared

  # Ensure same genes in same order
  common_genes <- intersect(de_a$gene_id, de_b$gene_id)
  de_a <- de_a[match(common_genes, de_a$gene_id), ]
  de_b <- de_b[match(common_genes, de_b$gene_id), ]

  stopifnot(all(de_a$gene_id == de_b$gene_id))

  n_genes_pair <- length(common_genes)

  for (metric_info in SHARED_GENE_METRICS) {
    metric_col <- metric_info$col
    metric_lab <- metric_info$lab

    vals_a <- de_a[[metric_col]]
    vals_b <- de_b[[metric_col]]

    # Remove NAs (e.g., missing dispersion)
    valid <- !is.na(vals_a) & !is.na(vals_b)
    vals_a_clean <- vals_a[valid]
    vals_b_clean <- vals_b[valid]
    n_valid <- sum(valid)

    if (n_valid < 3) {
      log_warn("    %s: only %d valid pairs, skipping", metric_lab, n_valid)
      next
    }

    # Paired Wilcoxon signed-rank test
    wt <- tryCatch({
      wilcox.test(vals_a_clean, vals_b_clean, paired = TRUE,
                  exact = FALSE, conf.int = FALSE)
    }, error = function(e) {
      log_warn("    %s: Wilcoxon test failed: %s", metric_lab, e$message)
      list(statistic = NA, p.value = NA)
    })

    # Count how many genes have A > B
    n_a_greater <- sum(vals_a_clean > vals_b_clean)

    # Pearson correlation (for selected metrics)
    pearson_r <- NA
    pearson_p <- NA
    if (metric_col %in% SHARED_GENE_CORR_METRICS) {
      ct <- cor.test(vals_a_clean, vals_b_clean, method = "pearson")
      pearson_r <- ct$estimate
      pearson_p <- ct$p.value
    }

    test_results[[length(test_results) + 1]] <- data.frame(
      ct_a          = ct_a,
      ct_b          = ct_b,
      pair          = pair_label,
      metric        = metric_col,
      metric_label  = metric_lab,
      n_genes       = n_valid,
      n_a_greater   = n_a_greater,
      n_b_greater   = n_valid - n_a_greater,
      wilcox_V      = as.numeric(wt$statistic),
      wilcox_p      = wt$p.value,
      median_a      = median(vals_a_clean),
      median_b      = median(vals_b_clean),
      mean_a        = mean(vals_a_clean),
      mean_b        = mean(vals_b_clean),
      pearson_r     = pearson_r,
      pearson_p     = pearson_p,
      stringsAsFactors = FALSE
    )

    log_info("    %s: %s>%s in %d/%d genes, Wilcoxon p=%.2e%s",
             metric_lab, ct_a, ct_b, n_a_greater, n_valid, wt$p.value,
             ifelse(!is.na(pearson_r),
                    sprintf(", r=%.3f", pearson_r), ""))
  }
}

test_df <- do.call(rbind, test_results)

# Save paired test results
test_path <- file.path(output_dir, "Shared_Gene_Paired_Tests.csv")
write.csv(test_df, test_path, row.names = FALSE)
log_info("Saved: %s (%d rows)", basename(test_path), nrow(test_df))

# ==============================================================================
# 7. FULL TRUTH SET SUMMARY (per-CT, all significant genes)
# ==============================================================================

log_info("Computing full truth set summaries...")

full_summary_rows <- list()

for (ct in CELL_TYPES) {
  de_all <- ct_data[[ct]]$de_all
  disp_vec <- ct_data[[ct]]$disp_vec

  # Filter to significant genes (padj < 0.05, |LFC| > 2)
  de_sig <- de_all[!is.na(de_all$padj) &
                     de_all$padj < SHARED_GENE_PADJ_CUTOFF &
                     abs(de_all$log2FoldChange) > SHARED_GENE_LFC_MIN, ]

  de_sig$dispersion <- disp_vec[de_sig$gene_id]
  de_sig$abs_lfc <- abs(de_sig$log2FoldChange)
  de_sig$snr <- de_sig$abs_lfc / de_sig$lfcSE

  full_summary_rows[[ct]] <- data.frame(
    cell_type        = ct,
    n_genes          = nrow(de_sig),
    median_abs_lfc   = median(de_sig$abs_lfc, na.rm = TRUE),
    median_lfcSE     = median(de_sig$lfcSE, na.rm = TRUE),
    median_snr       = median(de_sig$snr, na.rm = TRUE),
    median_baseMean  = median(de_sig$baseMean, na.rm = TRUE),
    median_dispersion = median(de_sig$dispersion, na.rm = TRUE),
    mean_abs_lfc     = mean(de_sig$abs_lfc, na.rm = TRUE),
    mean_lfcSE       = mean(de_sig$lfcSE, na.rm = TRUE),
    mean_snr         = mean(de_sig$snr, na.rm = TRUE),
    mean_baseMean    = mean(de_sig$baseMean, na.rm = TRUE),
    mean_dispersion  = mean(de_sig$dispersion, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  log_info("  %s: %d significant genes (padj < %g, |LFC| > %g)",
           ct, nrow(de_sig), SHARED_GENE_PADJ_CUTOFF, SHARED_GENE_LFC_MIN)
}

full_summary_df <- do.call(rbind, full_summary_rows)

full_summary_path <- file.path(output_dir, "Full_Truth_Set_Summary.csv")
write.csv(full_summary_df, full_summary_path, row.names = FALSE)
log_info("Saved: %s", basename(full_summary_path))

# Also compute shared gene summary (same format, for the 23 genes only)
shared_summary_rows <- list()

for (ct in CELL_TYPES) {
  de <- ct_data[[ct]]$de_shared
  shared_summary_rows[[ct]] <- data.frame(
    cell_type        = ct,
    n_genes          = nrow(de),
    median_abs_lfc   = median(de$abs_lfc, na.rm = TRUE),
    median_lfcSE     = median(de$lfcSE, na.rm = TRUE),
    median_snr       = median(de$snr, na.rm = TRUE),
    median_baseMean  = median(de$baseMean, na.rm = TRUE),
    median_dispersion = median(de$dispersion, na.rm = TRUE),
    mean_abs_lfc     = mean(de$abs_lfc, na.rm = TRUE),
    mean_lfcSE       = mean(de$lfcSE, na.rm = TRUE),
    mean_snr         = mean(de$snr, na.rm = TRUE),
    mean_baseMean    = mean(de$baseMean, na.rm = TRUE),
    mean_dispersion  = mean(de$dispersion, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

shared_summary_df <- do.call(rbind, shared_summary_rows)

shared_summary_path <- file.path(output_dir, "Shared_Gene_Summary.csv")
write.csv(shared_summary_df, shared_summary_path, row.names = FALSE)
log_info("Saved: %s", basename(shared_summary_path))

# ==============================================================================
# 8. SCATTER PLOTS (5 metrics × 3 CT pairs = 15 plots)
# ==============================================================================

log_info("Generating scatter plots...")

plot_dir <- file.path(output_dir, "Scatter_Plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

n_plots <- 0

for (pair in SHARED_GENE_CT_PAIRS) {
  ct_a <- pair[1]
  ct_b <- pair[2]

  # Get aligned data
  de_a <- ct_data[[ct_a]]$de_shared
  de_b <- ct_data[[ct_b]]$de_shared
  common_genes <- intersect(de_a$gene_id, de_b$gene_id)
  de_a <- de_a[match(common_genes, de_a$gene_id), ]
  de_b <- de_b[match(common_genes, de_b$gene_id), ]

  for (metric_info in SHARED_GENE_METRICS) {
    metric_col <- metric_info$col
    metric_lab <- metric_info$lab
    use_log <- metric_info$log_scale

    vals_a <- de_a[[metric_col]]
    vals_b <- de_b[[metric_col]]

    # Skip if too many NAs
    valid <- !is.na(vals_a) & !is.na(vals_b)
    if (sum(valid) < 3) next

    plot_df <- data.frame(
      gene_id = common_genes[valid],
      x = vals_a[valid],
      y = vals_b[valid]
    )

    # Compute Pearson r for annotation
    r_val <- cor(plot_df$x, plot_df$y, method = "pearson")

    # Look up Wilcoxon p from test results
    test_row <- test_df[test_df$ct_a == ct_a &
                          test_df$ct_b == ct_b &
                          test_df$metric == metric_col, ]
    wilcox_p <- if (nrow(test_row) > 0) test_row$wilcox_p[1] else NA

    # Format annotation
    annot_text <- sprintf("r = %.3f", r_val)
    if (!is.na(wilcox_p)) {
      if (wilcox_p < 0.001) {
        annot_text <- paste0(annot_text, sprintf("\nWilcoxon p = %.1e", wilcox_p))
      } else {
        annot_text <- paste0(annot_text, sprintf("\nWilcoxon p = %.3f", wilcox_p))
      }
    }

    # Build plot
    p <- ggplot(plot_df, aes(x = x, y = y)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey50", linewidth = 0.5) +
      geom_point(size = SHARED_SCATTER_POINT_SIZE,
                 alpha = SHARED_SCATTER_ALPHA,
                 colour = "#2166AC") +
      annotate("text", x = -Inf, y = Inf, label = annot_text,
               hjust = -0.1, vjust = 1.3, size = 3.5, fontface = "italic") +
      labs(
        title = sprintf("%s: %s vs %s", metric_lab, ct_a, ct_b),
        subtitle = sprintf("23 shared truth genes — paired comparison"),
        x = ct_a,
        y = ct_b
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey40")
      )

    # Apply log scale if appropriate
    if (use_log) {
      p <- p +
        scale_x_log10() +
        scale_y_log10()
    }

    # Save
    fname <- sprintf("Scatter_%s_%s_vs_%s.png", metric_col, ct_a, ct_b)
    ggsave(file.path(plot_dir, fname), p,
           width = SHARED_SCATTER_WIDTH, height = SHARED_SCATTER_HEIGHT,
           dpi = SHARED_SCATTER_DPI)
    n_plots <- n_plots + 1
  }
}

log_info("Generated %d scatter plots in: %s", n_plots, basename(plot_dir))

# ==============================================================================
# 9. SUMMARY
# ==============================================================================

log_section("STEP 5 COMPLETE")
log_info("Shared genes analysed: %d", n_shared)
log_info("Paired tests computed: %d (across %d CT pairs × %d metrics)",
         nrow(test_df), length(SHARED_GENE_CT_PAIRS),
         length(SHARED_GENE_METRICS))
log_info("Scatter plots generated: %d", n_plots)
log_info("Output directory: %s", output_dir)

# Print key findings
log_info("")
log_info("=== KEY FINDINGS ===")

# Dispersion tests
disp_tests <- test_df[test_df$metric == "dispersion", ]
if (nrow(disp_tests) > 0) {
  log_info("Dispersion (paired Wilcoxon):")
  for (i in seq_len(nrow(disp_tests))) {
    row <- disp_tests[i, ]
    sig_marker <- ifelse(row$wilcox_p < 0.05, "*", "NS")
    log_info("  %s: p=%.3e [%s], r=%.3f",
             row$pair, row$wilcox_p, sig_marker,
             ifelse(!is.na(row$pearson_r), row$pearson_r, NA))
  }
}

# BaseMean correlations
bm_tests <- test_df[test_df$metric == "baseMean", ]
if (nrow(bm_tests) > 0) {
  log_info("BaseMean correlations:")
  for (i in seq_len(nrow(bm_tests))) {
    row <- bm_tests[i, ]
    log_info("  %s: r=%.3f", row$pair,
             ifelse(!is.na(row$pearson_r), row$pearson_r, NA))
  }
}

log_info("")
log_info("Files saved:")
log_info("  %s", basename(test_path))
log_info("  %s", basename(comparison_path))
log_info("  %s", basename(wide_path))
log_info("  %s", basename(full_summary_path))
log_info("  %s", basename(shared_summary_path))
log_info("  %d scatter plots in %s/", n_plots, basename(plot_dir))
