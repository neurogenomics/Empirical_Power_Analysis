# ==============================================================================
# SCRIPT: Cross-Dataset Comparison - Visualise Metrics Across Datasets
# PURPOSE: Generate side-by-side histograms, violin plots, and ratio tables
# INPUT: Metrics from Step 01 (Datalink/Glial_Characteristics/Metrics/)
# OUTPUT: Comparison plots in Datalink/Glial_Characteristics/Plots/
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(gridExtra)
})

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Load config - skip if SCRIPT_DIR already set (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({ dirname(sys.frame(1)$ofile) }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "glial_characteristics")
  })
  source(file.path(script_path, "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### STEP 02: CROSS-DATASET COMPARISON ###\n")

# ==============================================================================
# 2. LOAD DATA FROM STEP 01
# ==============================================================================

message("Loading data from Step 01...")

# Load summary stats
all_stats <- read.csv(file.path(OUTPUT_SUMMARY, "All_Summary_Stats.csv"), 
                      stringsAsFactors = FALSE)
message(sprintf("  Loaded summary stats: %d rows", nrow(all_stats)))

# Load donor summaries
all_donors <- read.csv(file.path(OUTPUT_SUMMARY, "All_Donor_Summaries.csv"),
                       stringsAsFactors = FALSE)
message(sprintf("  Loaded donor summaries: %d donors", nrow(all_donors)))

# Load ratio comparison
ratio_comparison <- read.csv(file.path(OUTPUT_SUMMARY, "Ratio_Comparison.csv"),
                             stringsAsFactors = FALSE)

# Load cell-level distributions for plotting
cell_metrics_list <- list()
for (dataset in DATASET_ORDER) {
  for (cell_type in CELL_TYPES) {
    key <- paste(dataset, cell_type, sep = "_")
    cell_file <- file.path(OUTPUT_METRICS, 
                           sprintf("Cell_Metrics_%s_%s.qs", dataset, cell_type))
    if (file.exists(cell_file)) {
      cell_metrics_list[[key]] <- qs::qread(cell_file)
    }
  }
}

# Combine all cell metrics
all_cell_metrics <- bind_rows(cell_metrics_list)
message(sprintf("  Loaded cell metrics: %s cells", 
                format(nrow(all_cell_metrics), big.mark = ",")))

# Add prettified labels
all_cell_metrics$cell_type_pretty <- prettify_label(all_cell_metrics$cell_type)
all_cell_metrics$cell_type_pretty <- factor(all_cell_metrics$cell_type_pretty,
                                            levels = prettify_label(CELL_TYPES))
all_cell_metrics$dataset <- factor(all_cell_metrics$dataset, levels = DATASET_ORDER)

all_donors$cell_type_pretty <- prettify_label(all_donors$cell_type)
all_donors$cell_type_pretty <- factor(all_donors$cell_type_pretty,
                                      levels = prettify_label(CELL_TYPES))
all_donors$dataset <- factor(all_donors$dataset, levels = DATASET_ORDER)

# ==============================================================================
# 3. SIDE-BY-SIDE HISTOGRAMS
# ==============================================================================

message("\n=== Generating Side-by-Side Histograms ===")

# 3.1 UMIs per Cell - Faceted by Dataset
# --------------------------------------

p_umi_hist <- ggplot(all_cell_metrics, aes(x = umis, fill = cell_type_pretty)) +
  geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.5, position = "identity", color = NA) +
  facet_wrap(~dataset, ncol = 2, scales = "free_y") +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "UMIs per Cell — Roussos vs De Jager",
    subtitle = "Distribution comparison across DLPFC datasets",
    x = "UMIs per Cell",
    y = "Number of Cells",
    fill = "Cell Type"
  )

ggsave(file.path(OUTPUT_PLOTS, "Histogram_UMIs_Per_Cell_Comparison.png"),
       p_umi_hist, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Histogram_UMIs_Per_Cell_Comparison.png")

# Log10 version
p_umi_hist_log <- p_umi_hist + 
  scale_x_log10(labels = comma) +
  labs(
    title = "UMIs per Cell — Roussos vs De Jager (Log Scale)",
    x = "UMIs per Cell (log10)"
  )

ggsave(file.path(OUTPUT_PLOTS, "Histogram_UMIs_Per_Cell_Comparison_Log10.png"),
       p_umi_hist_log, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Histogram_UMIs_Per_Cell_Comparison_Log10.png")

# 3.2 Genes per Cell - Faceted by Dataset
# ----------------------------------------

p_genes_hist <- ggplot(all_cell_metrics, aes(x = genes, fill = cell_type_pretty)) +
  geom_histogram(bins = HISTOGRAM_BINS, alpha = 0.5, position = "identity", color = NA) +
  facet_wrap(~dataset, ncol = 2, scales = "free_y") +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Genes Detected per Cell — Roussos vs De Jager",
    subtitle = "Distribution comparison across DLPFC datasets",
    x = "Genes Detected per Cell",
    y = "Number of Cells",
    fill = "Cell Type"
  )

ggsave(file.path(OUTPUT_PLOTS, "Histogram_Genes_Per_Cell_Comparison.png"),
       p_genes_hist, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Histogram_Genes_Per_Cell_Comparison.png")

# 3.3 Cells per Donor - Faceted by Dataset
# -----------------------------------------

p_cells_hist <- ggplot(all_donors, aes(x = n_cells, fill = cell_type_pretty)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity", color = NA) +
  facet_wrap(~dataset, ncol = 2, scales = "free_y") +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Cells per Donor — Roussos vs De Jager",
    subtitle = "Distribution comparison across DLPFC datasets",
    x = "Cells per Donor",
    y = "Number of Donors",
    fill = "Cell Type"
  )

ggsave(file.path(OUTPUT_PLOTS, "Histogram_Cells_Per_Donor_Comparison.png"),
       p_cells_hist, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Histogram_Cells_Per_Donor_Comparison.png")

# Log10 version
p_cells_hist_log <- p_cells_hist + 
  scale_x_log10(labels = comma) +
  labs(
    title = "Cells per Donor — Roussos vs De Jager (Log Scale)",
    x = "Cells per Donor (log10)"
  )

ggsave(file.path(OUTPUT_PLOTS, "Histogram_Cells_Per_Donor_Comparison_Log10.png"),
       p_cells_hist_log, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Histogram_Cells_Per_Donor_Comparison_Log10.png")

# ==============================================================================
# 4. VIOLIN/BOX PLOTS
# ==============================================================================

message("\n=== Generating Violin/Box Plots ===")

# 4.1 UMIs per Cell - Split by Dataset
# -------------------------------------

# Subsample for violin plot (too many points otherwise)
# Subsample each group to at most 10,000 cells for violin plotting
set.seed(SUBSAMPLE_SEED)
cell_metrics_violin <- do.call(rbind, lapply(
  split(all_cell_metrics, list(all_cell_metrics$dataset, all_cell_metrics$cell_type)),
  function(df) {
    if (nrow(df) > 10000) df[sample.int(nrow(df), 10000), ] else df
  }
))

p_umi_violin <- ggplot(cell_metrics_violin, 
                       aes(x = cell_type_pretty, y = umis, fill = dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(0.8), scale = "width") +
  geom_boxplot(width = 0.15, position = position_dodge(0.8), 
               outlier.shape = NA, alpha = 0.8) +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = DATASET_COLOURS) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "UMIs per Cell — Dataset Comparison",
    subtitle = "Violin + boxplot overlay",
    x = "Cell Type",
    y = "UMIs per Cell",
    fill = "Dataset"
  )

ggsave(file.path(OUTPUT_PLOTS, "Violin_UMIs_Per_Cell_Comparison.png"),
       p_umi_violin, width = 10, height = 7, dpi = PLOT_DPI)
message("  Saved: Violin_UMIs_Per_Cell_Comparison.png")

# 4.2 Genes per Cell - Split by Dataset
# --------------------------------------

p_genes_violin <- ggplot(cell_metrics_violin,
                         aes(x = cell_type_pretty, y = genes, fill = dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(0.8), scale = "width") +
  geom_boxplot(width = 0.15, position = position_dodge(0.8),
               outlier.shape = NA, alpha = 0.8) +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = DATASET_COLOURS) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Genes Detected per Cell — Dataset Comparison",
    subtitle = "Violin + boxplot overlay",
    x = "Cell Type",
    y = "Genes Detected per Cell",
    fill = "Dataset"
  )

ggsave(file.path(OUTPUT_PLOTS, "Violin_Genes_Per_Cell_Comparison.png"),
       p_genes_violin, width = 10, height = 7, dpi = PLOT_DPI)
message("  Saved: Violin_Genes_Per_Cell_Comparison.png")

# 4.3 Cells per Donor - Split by Dataset
# ---------------------------------------

p_cells_violin <- ggplot(all_donors,
                         aes(x = cell_type_pretty, y = n_cells, fill = dataset)) +
  geom_violin(alpha = 0.6, position = position_dodge(0.8), scale = "width") +
  geom_boxplot(width = 0.15, position = position_dodge(0.8),
               outlier.shape = NA, alpha = 0.8) +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = DATASET_COLOURS) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Cells per Donor — Dataset Comparison",
    subtitle = "Violin + boxplot overlay",
    x = "Cell Type",
    y = "Cells per Donor",
    fill = "Dataset"
  )

ggsave(file.path(OUTPUT_PLOTS, "Violin_Cells_Per_Donor_Comparison.png"),
       p_cells_violin, width = 10, height = 7, dpi = PLOT_DPI)
message("  Saved: Violin_Cells_Per_Donor_Comparison.png")

# ==============================================================================
# 5. RATIO COMPARISON TABLE
# ==============================================================================

message("\n=== Generating Ratio Comparison Table ===")

# Create a formatted summary table
# Pivot all_stats to wide format: rows = cell_type × metric, cols = dataset

summary_wide <- all_stats %>%
  select(dataset, cell_type, metric, median) %>%
  pivot_wider(names_from = dataset, values_from = median, names_prefix = "median_") %>%
  mutate(
    diff_pct = 100 * (median_Roussos - median_DeJager) / median_DeJager,
    ratio_roussos_dejager = median_Roussos / median_DeJager
  )

summary_file <- file.path(OUTPUT_SUMMARY, "Median_Comparison_Wide.csv")
write.csv(summary_wide, summary_file, row.names = FALSE)
message(sprintf("  Saved: %s", basename(summary_file)))

# Create inter-cell-type ratio comparison plot
# This shows the ratio of ratios (e.g., Astro:Micro in Roussos vs DeJager)

ratio_plot_data <- ratio_comparison %>%
  mutate(
    metric_pretty = METRIC_LABELS[metric],
    ratio_label_pretty = gsub(":", " : ", ratio_label)
  )

# Pivot for plotting
ratio_long <- ratio_comparison %>%
  select(metric, ratio_label, Roussos, DeJager) %>%
  pivot_longer(cols = c(Roussos, DeJager), 
               names_to = "dataset", values_to = "ratio") %>%
  mutate(
    metric_pretty = METRIC_LABELS[metric],
    dataset = factor(dataset, levels = DATASET_ORDER)
  )

p_ratio_comparison <- ggplot(ratio_long, 
                             aes(x = ratio_label, y = ratio, fill = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~metric_pretty, ncol = 3, scales = "free_y") +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_fill_manual(values = DATASET_COLOURS) +
  labs(
    title = "Inter-Cell-Type Ratios — Roussos vs De Jager",
    subtitle = "Comparing median ratios relative to Microglia",
    x = "Cell Type Ratio",
    y = "Ratio (fold)",
    fill = "Dataset"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_PLOTS, "Ratio_Comparison_BarChart.png"),
       p_ratio_comparison, width = 12, height = 6, dpi = PLOT_DPI)
message("  Saved: Ratio_Comparison_BarChart.png")

# ==============================================================================
# 6. COMBINED SUMMARY FIGURE
# ==============================================================================

message("\n=== Generating Combined Summary Figure ===")

# Combine violin plots into one figure
combined_fig <- (p_umi_violin | p_genes_violin | p_cells_violin) +
  plot_annotation(
    title = "Glial Cell Characteristics — Cross-Dataset Comparison",
    subtitle = "Roussos (PsychAD, ~1,200 donors) vs De Jager (~450 donors), both DLPFC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_PLOTS, "Combined_Violin_Summary.png"),
       combined_fig, width = 18, height = 8, dpi = PLOT_DPI)
message("  Saved: Combined_Violin_Summary.png")

# ==============================================================================
# 7. PRINT SUMMARY
# ==============================================================================

message("\n=== Summary of Findings ===")

for (metric in unique(ratio_comparison$metric)) {
  metric_rows <- ratio_comparison %>% filter(metric == !!metric)
  message(sprintf("\n%s:", METRIC_LABELS[metric]))
  
  for (i in seq_len(nrow(metric_rows))) {
    r <- metric_rows[i, ]
    direction <- ifelse(r$Roussos > r$DeJager, "↑", "↓")
    diff_pct <- 100 * abs(r$Roussos / r$DeJager - 1)
    preserved <- diff_pct < 25  # <25% difference = "preserved"
    
    message(sprintf("  %s: Roussos=%.2f× vs DeJager=%.2f× %s (%.1f%% diff) %s",
                    r$ratio_label, r$Roussos, r$DeJager, direction, diff_pct,
                    ifelse(preserved, "[PRESERVED]", "[DIVERGENT]")))
  }
}

message("\n### STEP 02 COMPLETE ###")
