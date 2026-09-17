# ==============================================================================
# SCRIPT: DEG Chromosome Grid — Significant Gene Counts by Chr Type × LFC Bin
# PURPOSE: For each cell type, produce a 2×3 grid showing the number of
#          significant DEGs stratified by chromosome type (SexChr / Autosomal)
#          and LFC magnitude bin (Small 1-2, Medium 2-3, Large 3+).
#          Uses the full-cohort ground truth DE results from Roussos.
# INPUT: Full GT DE results from Datalink/Full_Ground_Truth/
#        BioMart cache for chromosome annotation
# OUTPUT: One grid plot per cell type + summary CSV
#         in Datalink/DE_Diagnostics/Chromosome_Grid/
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(scales)
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
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 4: DEG CHROMOSOME GRID")
log_info("Thresholds: padj < %g, |LFC| > %g",
         RECOVERY_FDR_THRESHOLD, RECOVERY_LFC_THRESHOLD)
log_info("LFC bins: %s", paste(names(RECOVERY_LFC_BINS), collapse = ", "))
log_info("Chr categories: Autosomal (1-22), SexChr (X, Y)")

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

output_dir <- get_chr_grid_output_dir()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  log_info("Created output directory: %s", output_dir)
}

# Check skip-if-exists
if (!FORCE_RECALC) {
  expected_plots <- paste0("DEG_Chr_Grid_", CELL_TYPES, ".png")
  existing <- file.exists(file.path(output_dir, expected_plots))
  if (all(existing) && file.exists(file.path(output_dir, "DEG_Chr_Grid_Summary.csv"))) {
    log_info("[SKIP] Chromosome grid outputs already exist. Set FORCE_RECALC=TRUE to regenerate.")
    .skip_step4 <- TRUE
  } else {
    .skip_step4 <- FALSE
  }
} else {
  .skip_step4 <- FALSE
}

if (!exists(".skip_step4") || !.skip_step4) {

# ==============================================================================
# 3. LOAD FULL GT DE RESULTS (ALL CELL TYPES)
# ==============================================================================

log_section("Loading Full Ground Truth DE results")

full_gt_list <- list()

for (ct in CELL_TYPES) {
  df <- load_de_results(ct,
                        padj_cutoff = RECOVERY_FDR_THRESHOLD,
                        lfc_min = RECOVERY_LFC_THRESHOLD)
  if (!is.null(df)) {
    full_gt_list[[ct]] <- df
  }
}

if (length(full_gt_list) == 0) {
  stop("No Full GT DE results loaded for any cell type. Check paths and filters.")
}

# ==============================================================================
# 4. GET CHROMOSOME ANNOTATIONS
# ==============================================================================

log_section("Loading chromosome annotations")

cache_path <- file.path(OUTPUT_BASE, BIOMART_CACHE_FILE)

# Ensure parent directory exists for cache
if (!dir.exists(dirname(cache_path))) {
  dir.create(dirname(cache_path), recursive = TRUE)
}

# Collect all significant gene IDs across cell types
all_sig_genes <- unique(unlist(lapply(full_gt_list, function(df) df$gene_id)))
log_info("Total unique significant genes across cell types: %d", length(all_sig_genes))

chr_annotations <- get_chromosome_annotations(all_sig_genes, cache_path)

if (is.null(chr_annotations) || nrow(chr_annotations) == 0) {
  stop("Chromosome annotation failed. Cannot generate grid plots.")
}

log_info("Chromosome annotations: %d genes annotated", nrow(chr_annotations))
chr_table <- table(chr_annotations$Chr_Type)
log_info("  Breakdown: %s",
         paste(sprintf("%s=%d", names(chr_table), as.integer(chr_table)),
               collapse = ", "))

# ==============================================================================
# 5. ANNOTATE AND CLASSIFY GENES
# ==============================================================================

log_section("Annotating genes with LFC bins and chromosome type")

grid_rows <- list()

for (ct in names(full_gt_list)) {
  df <- full_gt_list[[ct]]

  # Compute absolute LFC and assign bins
  df$Abs_LFC <- abs(df$log2FoldChange)
  df$LFC_Bin <- assign_lfc_bin(df$Abs_LFC)

  # Annotate with chromosome type
  chr_match <- chr_annotations[match(df$gene_id, chr_annotations$ensembl_gene_id), ]
  df$Chr_Type <- chr_match$Chr_Type

  # Keep only Autosomal and SexChr (drop Other/NA)
  df_clean <- df[!is.na(df$Chr_Type) & df$Chr_Type %in% c("Autosomal", "SexChr") &
                   !is.na(df$LFC_Bin), ]

  # Count genes per Chr_Type × LFC_Bin
  counts <- df_clean %>%
    dplyr::group_by(Chr_Type, LFC_Bin) %>%
    dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Cell_Type = ct)

  grid_rows[[ct]] <- counts

  log_info("  %s: %d genes classified (of %d significant)",
           ct, nrow(df_clean), nrow(df))
}

grid_df <- do.call(rbind, grid_rows)

# Ensure factor ordering for consistent plots
grid_df$LFC_Bin <- factor(grid_df$LFC_Bin,
                           levels = names(RECOVERY_LFC_BINS))
grid_df$Chr_Type <- factor(grid_df$Chr_Type,
                            levels = c("SexChr", "Autosomal"),
                            labels = prettify_label(c("SexChr", "Autosomal")))

# ==============================================================================
# 6. GENERATE GRID PLOTS (ONE PER CELL TYPE)
# ==============================================================================

log_section("Generating chromosome grid plots")

total_plots <- 0

for (ct in CELL_TYPES) {

  ct_data <- grid_df[grid_df$Cell_Type == ct, ]

  if (nrow(ct_data) == 0) {
    log_warn("  No data for %s — skipping", ct)
    next
  }

  # Ensure all 6 cells exist (fill missing combos with 0)
  chr_levels <- prettify_label(c("SexChr", "Autosomal"))
  full_grid <- expand.grid(
    Chr_Type = factor(chr_levels, levels = chr_levels),
    LFC_Bin = factor(names(RECOVERY_LFC_BINS), levels = names(RECOVERY_LFC_BINS)),
    stringsAsFactors = FALSE
  )
  ct_data <- merge(full_grid, ct_data, by = c("Chr_Type", "LFC_Bin"), all.x = TRUE)
  ct_data$N[is.na(ct_data$N)] <- 0
  ct_data$Cell_Type <- ct

  # Compute totals for subtitle
  n_auto <- sum(ct_data$N[ct_data$Chr_Type == prettify_label("Autosomal")])
  n_sex <- sum(ct_data$N[ct_data$Chr_Type == prettify_label("SexChr")])
  n_total <- n_auto + n_sex

  p <- ggplot(ct_data, aes(x = LFC_Bin, y = Chr_Type, fill = N)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = N), size = CHR_GRID_TEXT_SIZE, fontface = "bold") +
    scale_fill_gradient(low = CHR_GRID_FILL_LOW, high = CHR_GRID_FILL_HIGH,
                        name = "DEG Count") +
    labs(
      title = sprintf("%s: Significant DEGs by Chromosome Type and LFC Magnitude", prettify_label(ct)),
      subtitle = sprintf("Total: %d DEGs (Autosomal: %d, Sex Chr.: %d); padj < %g, |LFC| > %g",
                         n_total, n_auto, n_sex,
                         RECOVERY_FDR_THRESHOLD, RECOVERY_LFC_THRESHOLD),
      x = "LFC Magnitude Bin",
      y = "Chromosome Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )

  plot_file <- file.path(output_dir, sprintf("DEG_Chr_Grid_%s.png", ct))
  ggsave(plot_file, p,
         width = CHR_GRID_WIDTH, height = CHR_GRID_HEIGHT, dpi = PLOT_DPI)
  log_info("  Saved: %s", basename(plot_file))
  total_plots <- total_plots + 1
}

# ==============================================================================
# 7. GENERATE COMBINED PLOT (ALL CELL TYPES)
# ==============================================================================

log_info("Generating combined grid plot (all cell types)")

# Ensure all 18 cells exist (3 cell types × 2 chr × 3 bins)
full_combined <- expand.grid(
  Cell_Type = CELL_TYPES,
  Chr_Type = factor(prettify_label(c("SexChr", "Autosomal")), levels = prettify_label(c("SexChr", "Autosomal"))),
  LFC_Bin = factor(names(RECOVERY_LFC_BINS), levels = names(RECOVERY_LFC_BINS)),
  stringsAsFactors = FALSE
)
combined_data <- merge(full_combined, grid_df,
                       by = c("Cell_Type", "Chr_Type", "LFC_Bin"), all.x = TRUE)
combined_data$N[is.na(combined_data$N)] <- 0
combined_data$Cell_Type <- factor(combined_data$Cell_Type, levels = CELL_TYPES)

p_combined <- ggplot(combined_data, aes(x = LFC_Bin, y = Chr_Type, fill = N)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = N), size = CHR_GRID_TEXT_SIZE, fontface = "bold") +
  facet_wrap(~Cell_Type, nrow = 1, labeller = prettify_labeller()) +
  scale_fill_gradient(low = CHR_GRID_FILL_LOW, high = CHR_GRID_FILL_HIGH,
                      name = "DEG Count") +
  labs(
    title = "Significant DEGs by Chromosome Type × LFC Magnitude",
    subtitle = sprintf("padj < %g, |LFC| > %g; PsychAD full-cohort ground truth",
                       RECOVERY_FDR_THRESHOLD, RECOVERY_LFC_THRESHOLD),
    x = "LFC Magnitude Bin",
    y = "Chromosome Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

combined_file <- file.path(output_dir, "DEG_Chr_Grid_Combined.png")
ggsave(combined_file, p_combined,
       width = CHR_GRID_COMBINED_WIDTH, height = CHR_GRID_HEIGHT, dpi = PLOT_DPI)
log_info("  Saved: %s", basename(combined_file))
total_plots <- total_plots + 1

# ==============================================================================
# 8. SAVE SUMMARY CSV
# ==============================================================================

log_section("Saving summary table")

# Pivot wider for readable table: rows = Cell_Type × Chr_Type, cols = LFC bins
summary_wide <- grid_df %>%
  tidyr::pivot_wider(names_from = LFC_Bin, values_from = N, values_fill = 0) %>%
  dplyr::arrange(Cell_Type, desc(Chr_Type))

# Add row totals
lfc_cols <- names(RECOVERY_LFC_BINS)
summary_wide$Total <- rowSums(summary_wide[, lfc_cols, drop = FALSE], na.rm = TRUE)

summary_file <- file.path(output_dir, "DEG_Chr_Grid_Summary.csv")
write.csv(summary_wide, summary_file, row.names = FALSE)
log_info("Saved: %s", basename(summary_file))

# Print summary
for (ct in CELL_TYPES) {
  ct_rows <- summary_wide[summary_wide$Cell_Type == ct, ]
  for (i in seq_len(nrow(ct_rows))) {
    log_info("  %s / %s: %s | Total: %d",
             ct_rows$Cell_Type[i], ct_rows$Chr_Type[i],
             paste(sprintf("%s=%d", lfc_cols, as.integer(ct_rows[i, lfc_cols])),
                   collapse = ", "),
             ct_rows$Total[i])
  }
}

# ==============================================================================
# 9. COMPLETION REPORT
# ==============================================================================

log_section("STEP 4 COMPLETION REPORT")
log_info("Total plots generated: %d (3 per-cell-type + 1 combined)", total_plots)
log_info("Summary CSV: DEG_Chr_Grid_Summary.csv")
log_info("Step 4 complete")

} # end skip check
