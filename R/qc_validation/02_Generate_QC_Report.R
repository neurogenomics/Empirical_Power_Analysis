# ==============================================================================
# SCRIPT: Generate QC Validation Report
# PURPOSE: Produce diagnostic plots and summary tables from computed QC metrics.
#          Overlays Fullard et al. 2025 thresholds on histograms and scatter plots.
# INPUT: QC metric files from Datalink/QC_Validation/{CellType}/
# OUTPUT: PNG plots and summary tables in Datalink/QC_Validation/{CellType}/
# ==============================================================================

if (!exists("PROJECT_ROOT")) source("config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

message("=" |> rep(70) |> paste0(collapse = ""))
message("STEP 02: Generate QC Validation Report")
message("=" |> rep(70) |> paste0(collapse = ""))

# ==============================================================================
# Helper: Plot histogram with threshold lines
# ==============================================================================

plot_qc_histogram <- function(values, metric_name, cell_type,
                              lower_thresh = NULL, upper_thresh = NULL,
                              x_label = NULL, log_scale = FALSE) {
  if (is.null(x_label)) x_label <- metric_name

  df <- data.frame(value = values)

  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 100, fill = "steelblue", colour = "grey30", alpha = 0.7) +
    labs(
      title = sprintf("%s Distribution — %s", metric_name, cell_type),
      subtitle = sprintf("n = %s cells", format(length(values), big.mark = ",")),
      x = x_label,
      y = "Count"
    ) +
    theme_minimal(base_size = 12)

  if (!is.null(lower_thresh)) {
    n_below <- sum(values < lower_thresh, na.rm = TRUE)
    pct_below <- round(n_below / length(values) * 100, 2)
    p <- p + geom_vline(xintercept = lower_thresh,
      colour = THRESHOLD_COLOUR, linetype = THRESHOLD_LINETYPE, linewidth = 1) +
      annotate("text", x = lower_thresh, y = Inf, vjust = 2, hjust = -0.1,
        label = sprintf("Min: %s\n(%d cells below, %.2f%%)",
          format(lower_thresh, big.mark = ","), n_below, pct_below),
        colour = THRESHOLD_COLOUR, size = 3.5)
  }

  if (!is.null(upper_thresh)) {
    n_above <- sum(values > upper_thresh, na.rm = TRUE)
    pct_above <- round(n_above / length(values) * 100, 2)
    p <- p + geom_vline(xintercept = upper_thresh,
      colour = THRESHOLD_COLOUR, linetype = THRESHOLD_LINETYPE, linewidth = 1) +
      annotate("text", x = upper_thresh, y = Inf, vjust = 2, hjust = 1.1,
        label = sprintf("Max: %s\n(%d cells above, %.2f%%)",
          format(upper_thresh, big.mark = ","), n_above, pct_above),
        colour = THRESHOLD_COLOUR, size = 3.5)
  }

  if (log_scale) {
    p <- p + scale_x_log10(labels = scales::comma)
  }

  return(p)
}

# ==============================================================================
# Helper: Scatter plot (genes vs UMIs, coloured by mito%)
# ==============================================================================

plot_genes_vs_umis <- function(qc_df, cell_type) {
  # Subsample for plotting if too many cells
  if (nrow(qc_df) > 100000) {
    set.seed(42)
    plot_df <- qc_df[sample(nrow(qc_df), 100000), ]
    subtitle <- sprintf("n = %s cells (100k subsample for plotting)",
      format(nrow(qc_df), big.mark = ","))
  } else {
    plot_df <- qc_df
    subtitle <- sprintf("n = %s cells", format(nrow(qc_df), big.mark = ","))
  }

  p <- ggplot(plot_df, aes(x = n_umis, y = n_genes, colour = pct_mito)) +
    geom_point(alpha = 0.1, size = 0.3) +
    scale_colour_viridis_c(name = "Mito %", limits = c(0, 10), oob = scales::squish) +
    # Threshold box
    geom_vline(xintercept = QC_MIN_UMIS, colour = THRESHOLD_COLOUR,
      linetype = THRESHOLD_LINETYPE) +
    geom_vline(xintercept = QC_MAX_UMIS, colour = THRESHOLD_COLOUR,
      linetype = THRESHOLD_LINETYPE) +
    geom_hline(yintercept = QC_MIN_GENES, colour = THRESHOLD_COLOUR,
      linetype = THRESHOLD_LINETYPE) +
    geom_hline(yintercept = QC_MAX_GENES, colour = THRESHOLD_COLOUR,
      linetype = THRESHOLD_LINETYPE) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = sprintf("Genes vs UMIs — %s", cell_type),
      subtitle = subtitle,
      x = "Total UMIs (log scale)",
      y = "Genes Detected (log scale)"
    ) +
    theme_minimal(base_size = 12) +
    # Annotate the threshold box
    annotate("rect",
      xmin = QC_MIN_UMIS, xmax = QC_MAX_UMIS,
      ymin = QC_MIN_GENES, ymax = QC_MAX_GENES,
      fill = NA, colour = THRESHOLD_COLOUR, linetype = "solid", linewidth = 0.5)

  return(p)
}

# ==============================================================================
# Process each cell type
# ==============================================================================

global_summary <- list()

for (cell_type in CELL_TYPES) {
  message(sprintf("\n%s", "=" |> rep(60) |> paste0(collapse = "")))
  message(sprintf("Generating report for: %s", cell_type))
  message("=" |> rep(60) |> paste0(collapse = ""))

  output_dir <- file.path(OUTPUT_BASE, cell_type)

  # Load QC data (check for partial or full)
  cells_file <- file.path(output_dir, "qc_all_cells.qs")
  partial_file <- file.path(output_dir, "qc_all_cells_PARTIAL.qs")

  if (file.exists(cells_file)) {
    qc_df <- qs::qread(cells_file)
    is_partial <- FALSE
  } else if (file.exists(partial_file)) {
    qc_df <- qs::qread(partial_file)
    is_partial <- TRUE
    message(sprintf("  [WARN] Using PARTIAL results for %s (fail-early triggered)", cell_type))
  } else {
    message(sprintf("  [SKIP] No QC data found for %s. Run Step 01 first.", cell_type))
    next
  }

  message(sprintf("  Loaded %s cells", format(nrow(qc_df), big.mark = ",")))

  # ------ Histograms ------

  # UMI distribution
  p <- plot_qc_histogram(qc_df$n_umis, "UMI Count", cell_type,
    lower_thresh = QC_MIN_UMIS, upper_thresh = QC_MAX_UMIS,
    x_label = "Total UMIs per Cell", log_scale = TRUE)
  ggsave(file.path(output_dir, "Histogram_UMIs.png"), p,
    width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

  # Gene count distribution
  p <- plot_qc_histogram(qc_df$n_genes, "Gene Count", cell_type,
    lower_thresh = QC_MIN_GENES, upper_thresh = QC_MAX_GENES,
    x_label = "Genes Detected per Cell")
  ggsave(file.path(output_dir, "Histogram_Genes.png"), p,
    width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

  # Mito % distribution
  if (!all(is.na(qc_df$pct_mito))) {
    p <- plot_qc_histogram(qc_df$pct_mito, "Mitochondrial %", cell_type,
      upper_thresh = QC_MAX_MITO_PCT,
      x_label = "Mitochondrial UMI %")
    ggsave(file.path(output_dir, "Histogram_Mito.png"), p,
      width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  }

  # ------ Scatter plot ------
  p <- plot_genes_vs_umis(qc_df, cell_type)
  ggsave(file.path(output_dir, "Scatter_Genes_vs_UMIs.png"), p,
    width = PLOT_WIDTH + 2, height = PLOT_HEIGHT, dpi = PLOT_DPI)

  # ------ Violation summary ------
  violations <- summarise_violations(qc_df)

  verdict_lines <- c(
    sprintf("QC VALIDATION VERDICT — %s", cell_type),
    sprintf("Generated: %s", Sys.time()),
    if (is_partial) "⚠️  PARTIAL RESULTS (fail-early triggered)" else "",
    "",
    sprintf("Total cells analysed: %s", format(violations$n_total, big.mark = ",")),
    "",
    "Threshold Violations:",
    sprintf("  UMIs < %s:     %s (%.4f%%)",
      format(QC_MIN_UMIS, big.mark = ","),
      format(violations$counts[["fail_low_umis"]], big.mark = ","),
      violations$percentages[["fail_low_umis"]]),
    sprintf("  UMIs > %s: %s (%.4f%%)",
      format(QC_MAX_UMIS, big.mark = ","),
      format(violations$counts[["fail_high_umis"]], big.mark = ","),
      violations$percentages[["fail_high_umis"]]),
    sprintf("  Genes < %s:    %s (%.4f%%)",
      format(QC_MIN_GENES, big.mark = ","),
      format(violations$counts[["fail_low_genes"]], big.mark = ","),
      violations$percentages[["fail_low_genes"]]),
    sprintf("  Genes > %s: %s (%.4f%%)",
      format(QC_MAX_GENES, big.mark = ","),
      format(violations$counts[["fail_high_genes"]], big.mark = ","),
      violations$percentages[["fail_high_genes"]]),
    sprintf("  Mito >= %.1f%%:    %s (%.4f%%)",
      QC_MAX_MITO_PCT,
      format(violations$counts[["fail_high_mito"]], big.mark = ","),
      violations$percentages[["fail_high_mito"]]),
    "",
    sprintf("  ANY violation:  %s (%.4f%%)",
      format(violations$counts[["fail_any"]], big.mark = ","),
      violations$percentages[["fail_any"]]),
    ""
  )

  # Determine verdict
  any_pct <- violations$percentages[["fail_any"]]
  if (any_pct == 0) {
    verdict <- "✅ PASS — All cells within Fullard et al. 2025 thresholds. Data is POST-QC."
  } else if (any_pct < 1) {
    verdict <- "✅ PASS (with minor edge cases) — <1% violations. Data is likely POST-QC with minor boundary cases."
  } else if (any_pct < 5) {
    verdict <- "⚠️  MARGINAL — 1-5% violations. Data may be partially QC'd or thresholds differ slightly."
  } else if (any_pct < 20) {
    verdict <- "🔴 CONCERN — 5-20% violations. Data may not have full QC applied."
  } else {
    verdict <- "🔴 FAIL — >20% violations. Data is likely RAW/PRE-QC."
  }

  verdict_lines <- c(verdict_lines, sprintf("VERDICT: %s", verdict))

  writeLines(verdict_lines, file.path(output_dir, "QC_Verdict.txt"))
  message(sprintf("\n  %s", verdict))

  # Store for global summary
  global_summary[[cell_type]] <- data.frame(
    cell_type = cell_type,
    n_cells = violations$n_total,
    n_fail_any = violations$counts[["fail_any"]],
    pct_fail_any = violations$percentages[["fail_any"]],
    n_fail_low_umis = violations$counts[["fail_low_umis"]],
    n_fail_high_umis = violations$counts[["fail_high_umis"]],
    n_fail_low_genes = violations$counts[["fail_low_genes"]],
    n_fail_high_genes = violations$counts[["fail_high_genes"]],
    n_fail_high_mito = violations$counts[["fail_high_mito"]],
    is_partial = is_partial,
    verdict = verdict,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Global summary
# ==============================================================================

if (length(global_summary) > 0) {
  global_df <- do.call(rbind, global_summary)
  global_file <- file.path(OUTPUT_BASE, "Global_QC_Summary.csv")
  write.csv(global_df, global_file, row.names = FALSE)
  message(sprintf("\n[DONE] Global summary written to: %s", global_file))

  # Print global summary
  message("\n=== GLOBAL QC VALIDATION SUMMARY ===")
  for (i in seq_len(nrow(global_df))) {
    row <- global_df[i, ]
    message(sprintf("  %s: %s cells, %.4f%% violations — %s",
      row$cell_type,
      format(row$n_cells, big.mark = ","),
      row$pct_fail_any,
      row$verdict))
  }
}

message("\n[DONE] Step 02 complete.")
