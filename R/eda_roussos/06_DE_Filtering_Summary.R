# ==============================================================================
# SCRIPT: DE Filtering Summary
# PURPOSE: Summarise how many genes were excluded by DESeq2 filtering and why,
#          for each cell type, binned by LFC magnitude
# INPUT:   Full Ground Truth DE results (CSV) and pseudobulk checkpoints (.qs)
#          from Datalink/Full_Ground_Truth/{CellType}/
# OUTPUT:  Per-cell-type and combined summary CSVs + stacked bar chart
#          in Datalink/Exploratory_Analysis/
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

# Load shared utilities for prettified labels
# Use hardcoded fallback path for Rscript compatibility
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Load config - skip if SCRIPT_DIR already set (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({ dirname(sys.frame(1)$ofile) }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "eda_roussos")
  })
  source(file.path(script_path, "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### STEP 06: DE FILTERING SUMMARY ###\n")

# ==============================================================================
# 2. CLASSIFY GENES BY FILTER CATEGORY
# ==============================================================================

#' Classify genes from DE results into filtering categories
#'
#' DESeq2 sets NA in different columns for different reasons:
#'   - baseMean == 0: gene has zero counts across all samples (all cols NA)
#'   - pvalue is NA (baseMean > 0): Cook's distance outlier
#'   - padj is NA but pvalue is not NA: independent filtering (low mean)
#'   - padj is not NA: fully tested
#'
#' Additionally, genes present in the pseudobulk but absent from DE results
#' were removed by the low-count pre-filter (rowSums < MIN_GENE_COUNT).
#'
#' @param de_df Data.frame of DE results (from CSV)
#' @param pb_genes Character vector of all gene IDs in the pseudobulk
#' @return Data.frame with columns: gene, filter_category, baseMean, log2FoldChange
classify_gene_filters <- function(de_df, pb_genes) {

  # --- Genes in DE results: classify by NA pattern ---
  de_classified <- de_df %>%
    mutate(
      filter_category = case_when(
        baseMean == 0                    ~ "Low-Count Pre-Filter",
        is.na(pvalue) & baseMean > 0     ~ "Cook's Distance Outlier",
        is.na(padj) & !is.na(pvalue)     ~ "Independent Filtering",
        !is.na(padj)                     ~ "Fully Tested",
        TRUE                             ~ "Unknown"
      )
    ) %>%
    select(gene, filter_category, baseMean, log2FoldChange)

  # --- Genes in pseudobulk but NOT in DE results: low-count pre-filter ---
  de_genes <- de_df$gene
  prefiltered_genes <- setdiff(pb_genes, de_genes)

  if (length(prefiltered_genes) > 0) {
    prefiltered_df <- data.frame(
      gene = prefiltered_genes,
      filter_category = "Low-Count Pre-Filter",
      baseMean = NA_real_,
      log2FoldChange = NA_real_,
      stringsAsFactors = FALSE
    )
    de_classified <- bind_rows(de_classified, prefiltered_df)
  }

  return(de_classified)
}

#' Bin genes by |LFC| magnitude
#'
#' @param classified_df Data.frame from classify_gene_filters()
#' @param breaks Numeric vector of bin boundaries (e.g., c(1, 2, 3, Inf))
#' @param labels Character vector of bin labels
#' @return Data.frame with added lfc_bin column (NA for genes without LFC)
bin_by_lfc <- function(classified_df, breaks, labels) {
  classified_df %>%
    mutate(
      abs_lfc = abs(log2FoldChange),
      lfc_bin = cut(abs_lfc, breaks = breaks, labels = labels,
                    include.lowest = TRUE, right = FALSE)
    ) %>%
    # Only keep genes with |LFC| >= min(breaks)
    # Genes below the first break or with NA LFC are excluded from binned summary
    select(-abs_lfc)
}

# ==============================================================================
# 3. PROCESS EACH CELL TYPE
# ==============================================================================

all_summaries <- list()
all_classified <- list()

for (cell_type in CELL_TYPES) {

  message(sprintf("\n=== Processing: %s ===", cell_type))

  # --- Load pseudobulk checkpoint (gene universe before DESeq2 filtering) ---
  pb_path <- file.path(FULL_GT_BASE, cell_type,
                        paste0(DE_DATASET, "_pb_checkpoint.qs"))

  if (!file.exists(pb_path)) {
    warning(sprintf("Pseudobulk checkpoint not found: %s. Skipping %s.",
                    pb_path, cell_type))
    next
  }

  pb_sce <- qs::qread(pb_path)
  pb_genes <- rownames(pb_sce)
  message(sprintf("  Pseudobulk gene universe: %d genes", length(pb_genes)))
  rm(pb_sce)
  gc(verbose = FALSE)

  # --- Load DE results CSV ---
  de_path <- file.path(FULL_GT_BASE, cell_type,
                        paste0(DE_DATASET, "_DE_results.csv"))

  if (!file.exists(de_path)) {
    warning(sprintf("DE results not found: %s. Skipping %s.", de_path, cell_type))
    next
  }

  de_df <- read.csv(de_path, stringsAsFactors = FALSE)
  message(sprintf("  DE results: %d genes", nrow(de_df)))

  # --- Classify genes ---
  classified <- classify_gene_filters(de_df, pb_genes)
  classified$cell_type <- cell_type

  # --- Summary: total counts per category ---
  total_summary <- classified %>%
    group_by(filter_category) %>%
    summarise(total_count = n(), .groups = "drop") %>%
    mutate(cell_type = cell_type)

  message(sprintf("  Filter category breakdown:"))
  for (i in seq_len(nrow(total_summary))) {
    message(sprintf("    %-30s %d", total_summary$filter_category[i],
                    total_summary$total_count[i]))
  }

  # --- Bin by LFC and summarise ---
  binned <- bin_by_lfc(classified, DE_FILTER_LFC_BINS, DE_FILTER_LFC_LABELS)

  # Only genes with |LFC| >= 1 (first bin boundary)
  binned_filtered <- binned %>%
    filter(!is.na(lfc_bin))

  lfc_summary <- binned_filtered %>%
    group_by(filter_category, lfc_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(cell_type = cell_type)

  message(sprintf("  Genes with |LFC| >= 1: %d (of %d total)",
                  nrow(binned_filtered), nrow(classified)))

  # --- Store results ---
  all_classified[[cell_type]] <- classified
  all_summaries[[cell_type]] <- list(
    total = total_summary,
    by_lfc = lfc_summary
  )

  # --- Save per-cell-type CSV ---
  output_paths <- get_output_paths(cell_type, dataset = DE_DATASET)
  create_output_dirs(output_paths)

  csv_path <- file.path(output_paths$tables, "DE_Filtering_Summary.csv")

  # Build wide-format table: rows = filter_category, columns = LFC bins + Total
  # Start from total_summary to ensure all categories appear (even those with
  # no LFC-binned genes, like Low-Count Pre-Filter)
  lfc_wide <- lfc_summary %>%
    select(-cell_type) %>%
    pivot_wider(names_from = lfc_bin, values_from = count, values_fill = 0)

  wide_summary <- total_summary %>%
    select(filter_category, total_count) %>%
    left_join(lfc_wide, by = "filter_category") %>%
    mutate(across(all_of(DE_FILTER_LFC_LABELS), ~ replace_na(.x, 0L))) %>%
    arrange(match(filter_category, DE_FILTER_CATEGORIES))

  write.csv(wide_summary, csv_path, row.names = FALSE)
  message(sprintf("  Saved: %s", csv_path))

  rm(de_df, classified, binned, binned_filtered)
  gc(verbose = FALSE)
}

# ==============================================================================
# 4. COMBINED SUMMARY CSV
# ==============================================================================

if (length(all_summaries) == 0) {
  warning("No cell types processed. Exiting Step 06.")
} else {

  # Combine all LFC summaries
  combined_lfc <- bind_rows(lapply(all_summaries, `[[`, "by_lfc"))
  combined_total <- bind_rows(lapply(all_summaries, `[[`, "total"))

  # Wide format: cell_type × filter_category × LFC bins
  # Start from total to ensure all categories appear (including Low-Count
  # Pre-Filter which has no LFC-binned genes)
  lfc_combined_wide <- combined_lfc %>%
    pivot_wider(
      id_cols = c(cell_type, filter_category),
      names_from = lfc_bin,
      values_from = count,
      values_fill = list(count = 0L)
    )

  combined_wide <- combined_total %>%
    select(cell_type, filter_category, total_count) %>%
    left_join(lfc_combined_wide, by = c("cell_type", "filter_category")) %>%
    mutate(across(any_of(DE_FILTER_LFC_LABELS), ~ replace_na(.x, 0L))) %>%
    arrange(cell_type, match(filter_category, DE_FILTER_CATEGORIES))

  combined_csv_path <- file.path(OUTPUT_BASE, "DE_Filtering_Summary_Combined.csv")
  if (!dir.exists(dirname(combined_csv_path))) {
    dir.create(dirname(combined_csv_path), recursive = TRUE)
  }
  write.csv(combined_wide, combined_csv_path, row.names = FALSE)
  message(sprintf("\nSaved combined summary: %s", combined_csv_path))

  # ==========================================================================
  # 5. STACKED BAR CHART (FACETED BY CELL TYPE)
  # ==========================================================================

  # Use the LFC-binned data (|LFC| >= 1 only)
  plot_data <- combined_lfc %>%
    mutate(
      filter_category = factor(filter_category, levels = DE_FILTER_CATEGORIES),
      lfc_bin = factor(lfc_bin, levels = DE_FILTER_LFC_LABELS),
      cell_type = factor(cell_type, levels = CELL_TYPES),
      cell_type_pretty = prettify_label(cell_type)
    )
  plot_data$cell_type_pretty <- factor(plot_data$cell_type_pretty, 
    levels = prettify_label(CELL_TYPES))

  # Remove categories with zero counts across all bins (e.g., Low-Count Pre-Filter
  # will have zero in LFC bins since those genes have no LFC estimate)
  plot_data <- plot_data %>%
    filter(count > 0)

  # --- Pre-compute stacked bar segment positions for hybrid labelling ---
  # Compute cumulative y-positions (ymin, ymax, ymid) per bar within each facet,
  # and flag segments as "small" if they occupy < 5% of the total bar height.
  SMALL_SEGMENT_THRESHOLD <- 0.05

  plot_data <- plot_data %>%
    arrange(cell_type, lfc_bin, filter_category) %>%
    group_by(cell_type, lfc_bin) %>%
    mutate(
      bar_total = sum(count),
      ymax = cumsum(count),
      ymin = ymax - count,
      ymid = (ymin + ymax) / 2,
      is_small = (count / bar_total) < SMALL_SEGMENT_THRESHOLD
    ) %>%
    ungroup()

  # Large segments: white bold text inside the bar (original style)
  large_data <- plot_data %>% filter(!is_small)
  # Small segments: repelled outside with leader lines
  small_data <- plot_data %>% filter(is_small)

  p <- ggplot(plot_data, aes(x = lfc_bin, y = count, fill = filter_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    # Large segments: labels inside
    geom_text(
      data = large_data,
      aes(x = lfc_bin, y = ymid, label = count),
      size = 4, colour = "white", fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Small segments: labels repelled outside with reduced nudge
    geom_text_repel(
      data = small_data,
      aes(x = lfc_bin, y = ymid, label = count),
      size = 3.5, colour = "grey30", fontface = "bold",
      direction = "y",
      nudge_x = 0.35,
      segment.size = 0.3,
      segment.colour = "grey50",
      min.segment.length = 0,
      box.padding = 0.2,
      point.padding = 0.1,
      force = 2,
      max.overlaps = Inf,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ cell_type_pretty, scales = "free_y", nrow = 1) +
    scale_fill_manual(
      values = DE_FILTER_COLOURS,
      name = "Filter Category",
      drop = TRUE
    ) +
    labs(
      title = "DESeq2 Gene Filtering by Category and LFC Magnitude",
      x = "|log2 Fold Change| Bin",
      y = "Number of Genes"
    ) +
    theme_publication(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  plot_path <- file.path(OUTPUT_BASE, "DE_Filtering_Summary_BarChart.png")
  ggsave(plot_path, p, width = 12, height = 5, dpi = 300)
  message(sprintf("Saved bar chart: %s", plot_path))

  # ==========================================================================
  # 6. TOTALS BAR CHART (ALL GENES, NOT BINNED BY LFC)
  # ==========================================================================

  # This shows the full picture including low-count pre-filtered genes
  totals_data <- combined_total %>%
    mutate(
      filter_category = factor(filter_category, levels = DE_FILTER_CATEGORIES),
      cell_type = factor(cell_type, levels = CELL_TYPES),
      cell_type_pretty = prettify_label(cell_type)
    )
  totals_data$cell_type_pretty <- factor(totals_data$cell_type_pretty,
    levels = prettify_label(CELL_TYPES))

  # --- Pre-compute stacked bar segment positions for hybrid labelling ---
  totals_data <- totals_data %>%
    arrange(cell_type, filter_category) %>%
    group_by(cell_type) %>%
    mutate(
      bar_total = sum(total_count),
      ymax = cumsum(total_count),
      ymin = ymax - total_count,
      ymid = (ymin + ymax) / 2,
      is_small = (total_count / bar_total) < SMALL_SEGMENT_THRESHOLD
    ) %>%
    ungroup()

  # Large segments: white bold text inside the bar (original style)
  totals_large <- totals_data %>% filter(!is_small)
  # Small segments: repelled outside with leader lines
  totals_small <- totals_data %>% filter(is_small)

  p_totals <- ggplot(totals_data,
                      aes(x = cell_type_pretty, y = total_count, fill = filter_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    # Large segments: labels inside
    geom_text(
      data = totals_large,
      aes(x = cell_type_pretty, y = ymid, label = format(total_count, big.mark = ",")),
      size = 4, colour = "white", fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Small segments: labels repelled outside with reduced nudge
    geom_text_repel(
      data = totals_small %>% mutate(cell_type_pretty = prettify_label(cell_type)),
      aes(x = cell_type_pretty, y = ymid, label = format(total_count, big.mark = ",")),
      size = 3.5, colour = "grey30", fontface = "bold",
      direction = "y",
      nudge_x = 0.35,
      segment.size = 0.3,
      segment.colour = "grey50",
      min.segment.length = 0,
      box.padding = 0.3,
      point.padding = 0.1,
      force = 3,
      max.overlaps = Inf,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = DE_FILTER_COLOURS,
      name = "Filter Category"
    ) +
    labs(
      title = "DESeq2 Gene Filtering — Total Genes by Category",
      x = "Cell Type",
      y = "Number of Genes"
    ) +
    theme_publication(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  totals_plot_path <- file.path(OUTPUT_BASE, "DE_Filtering_Summary_Totals.png")
  ggsave(totals_plot_path, p_totals, width = 8, height = 6, dpi = 300)
  message(sprintf("Saved totals bar chart: %s", totals_plot_path))

  # ==========================================================================
  # 7. CONSOLE SUMMARY
  # ==========================================================================

  message("\n========================================")
  message("DE FILTERING SUMMARY — ALL CELL TYPES")
  message("========================================\n")

  for (ct in CELL_TYPES) {
    if (!ct %in% names(all_summaries)) next

    total <- all_summaries[[ct]]$total
    n_total <- sum(total$total_count)

    message(sprintf("--- %s (%d total genes) ---", ct, n_total))
    for (i in seq_len(nrow(total))) {
      pct <- 100 * total$total_count[i] / n_total
      message(sprintf("  %-30s %6d (%5.1f%%)",
                      total$filter_category[i],
                      total$total_count[i], pct))
    }
    message("")
  }
}

message("\n### STEP 06 COMPLETE ###\n")
