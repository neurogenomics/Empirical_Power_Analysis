# =============================================================================
# SCRIPT: compose_recall_panels.R
# PURPOSE: Build a single multi-panel recall figure for the draft, stacking three
#          all_methods recall grids (Oligo, cells view by default):
#            A  Recall by padj quartile           (4 facets: Q1..Q4)
#            B  Recall across expression, Tiny bin (|LFC| 0.2-1.0, 4 BaseMean facets)
#            C  Recall across expression, Small bin (|LFC| 1.0-2.0, 4 BaseMean facets)
#          Panels A/B/C are tagged; a single shared legend is collected. Each grid
#          keeps a paper-presentable title (Chris reinstated titles 2026-08-17);
#          only subtitles and "(n = x genes)" facet annotations are dropped.
#
#          Facets are ordered so the highest category is TOP-LEFT and reads
#          row-major (top-left, top-right, bottom-left, bottom-right): panel A
#          Q1 first (most significant), panels B/C highest expression first.
#
# The three plots are rebuilt from the pipeline checkpoints (not the on-disk
# PNGs) so they are subtitle-free and carry NO "(n = x genes)" facet annotations,
# matching the panel conventions in compose_crossscore.R.
#
# SOURCES (paths derived from the cartesian pipeline layout):
#   A  per-method df_pval_q_summary, pooled across deseq2/edger/limma exactly as
#      04_Visualise_Power.R does for Mean_Recall_PValQuartile_AllMethods.png.
#   B/C  all_methods df_finecross_summary (Step 8 checkpoint), the Figure T/S source.
#
# Reads checkpoints only. Writes the composed figure to
#   R/figure_composition/output_recall_panels/<name>.{png,pdf}  (staging), and a
#   canonical copy to Datalink/Composed_Figures/figure_recall_panels_<cell>_<mode>.{png,pdf}
#   which the manuscript chunk fig:FigRecallPanels includes. Re-running refreshes
#   the draft figure. It does not edit the .rmd itself.
#
# Run (on a host where the .qs checkpoints and R live):
#   Rscript R/figure_composition/compose_recall_panels.R
#   Rscript R/figure_composition/compose_recall_panels.R --cell Oligo --mode cells
#   Rscript R/figure_composition/compose_recall_panels.R --datalink /path/to/Datalink
# =============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required. Install: install.packages('ggnewscale')")
  }
  library(ggnewscale)
})

# ---------------------------------------------------------------------------
# Paths / args
# ---------------------------------------------------------------------------

HERE <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE))
if (is.na(HERE) || HERE == "") HERE <- normalizePath(".", mustWork = FALSE)
PROJECT_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)

# Reuse the shared theme + pretty cell-type labels if available.
prettify_src <- file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R")
if (file.exists(prettify_src)) {
  source(prettify_src)
} else {
  prettify_label   <- function(x) x
  theme_publication <- function(base_size = 14, ...) theme_bw(base_size = base_size)
}

# Canonical grid config + axis cleaners (CELL_STEPS/DONOR_STEPS/READ_STEPS and
# clean_cells_factor / clean_donors_factor / clean_reads_factor). Using these
# rather than a homegrown numeric coercion is load-bearing: the top row of the
# grid is the "Max" (all-available cells) bucket, stored as NA in the checkpoint.
# The canonical cleaners relabel that NA as "Max"; a naive as.numeric() coercion
# turns it into a spurious "NA" row and drops it out of order.
Sys.setenv(POWER_ANALYSIS_ROOT = PROJECT_ROOT)
config_src <- file.path(PROJECT_ROOT, "R", "cartesian_downsampling", "config.R")
utils_src  <- file.path(PROJECT_ROOT, "R", "cartesian_downsampling", "utils.R")
source(config_src)
source(utils_src)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  pos <- which(args == flag)
  if (length(pos) && pos[1] < length(args)) return(args[pos[1] + 1])
  default
}

DATALINK <- normalizePath(get_arg("--datalink", file.path(PROJECT_ROOT, "Datalink")),
                          mustWork = FALSE)
CELL     <- get_arg("--cell", "Oligo")
MODE     <- get_arg("--mode", "cells")
# --dataset: shard dataset whose grid to draw (Roussos = sex contrast,
# Roussos_ADvsCtrl = AD vs control); --composed-dir: output dir under Datalink.
DATASET  <- get_arg("--dataset", "Roussos")
COMPOSED_SUBDIR <- get_arg("--composed-dir", "Composed_Figures")
if (!MODE %in% c("cells", "reads")) stop("--mode must be 'cells' or 'reads'")

Y_COL   <- if (MODE == "reads") "Target_Reads" else "Cells"
Y_LABEL <- if (MODE == "reads") "Target Counts per Donor" else "Cells per Donor"
clean_y_factor <- if (MODE == "reads") clean_reads_factor else clean_cells_factor

DE_METHODS         <- c("deseq2", "edger", "limma")
CROSS_METHOD_LABEL <- "all_methods"
BASEMEAN_ORDER     <- c("High (200+)", "Medium (50-200)", "Low (10-50)", "Very Low (<10)")
FINECROSS_LFC_BINS <- c("Tiny (0.2-1.0)", "Small (1.0-2.0)")

OUT_DIR <- file.path(HERE, "output_recall_panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

checkpoint_path <- function(cell, method) {
  file.path(DATALINK, "Cartesian_Downsampling", cell, DATASET, MODE, method,
            "Plots", "Aggregated_Summary_Data.qs")
}

msg <- function(...) cat(sprintf(...), "\n")

# ---------------------------------------------------------------------------
# Core recall-grid builder (subtitle-free, no gene-count facets)
#   df         : data with Donors, Y_COL, <facet_col>, Mean_Recall
#   facet_col  : facet variable name
#   facet_levs : ordered facet levels (top-left -> bottom-right, row-major ncol=2)
#   title_str  : paper-presentable per-grid title (kept; only subtitles/gene
#                counts were dropped by the general rule)
# ---------------------------------------------------------------------------

build_recall_grid <- function(df, facet_col, facet_levs, title_str) {
  df <- df %>%
    mutate(
      Donors = clean_donors_factor(Donors),
      !!Y_COL := clean_y_factor(.data[[Y_COL]]),
      .facet  = factor(.data[[facet_col]], levels = facet_levs)
    ) %>%
    filter(!is.na(.facet)) %>%
    # Complete the grid so incomplete regions render as grey (NA fill) tiles
    # rather than white gaps, matching every other grid in the pipeline.
    tidyr::complete(Donors, !!sym(Y_COL), .facet)

  # Grey out entire facets that carry essentially no ground-truth genes (mean
  # gene count below threshold): showing metrics for a near-empty bin (e.g. High
  # expression within a small-LFC bin) is misleading, so blank the whole facet.
  if ("Ref_GT_Count" %in% names(df)) {
    facet_means <- tapply(df$Ref_GT_Count, df$.facet, function(v) mean(v, na.rm = TRUE))
    empty_facets <- names(facet_means)[is.nan(facet_means) |
                                         facet_means < EMPTY_FACET_GENE_THRESHOLD]
    if (length(empty_facets) > 0) {
      df$Mean_Recall[as.character(df$.facet) %in% empty_facets] <- NA_real_
    }
  }

  p <- ggplot(df, aes(x = Donors, y = .data[[Y_COL]])) +
    geom_tile(aes(fill = Mean_Recall), color = "white") +
    facet_wrap(~.facet, ncol = 2) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50",
                         name = "Mean Recall", guide = guide_colorbar(order = 1)) +
    labs(title = title_str, x = "Donors", y = Y_LABEL) +
    theme_publication(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", size = 16),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      strip.text   = element_text(face = "bold", size = 12),
      panel.spacing = unit(0.1, "lines"),
      legend.position   = "right",
      legend.key.width  = unit(0.7, "cm"),
      legend.key.height = unit(1.3, "cm"),
      legend.text       = element_text(size = 16),
      legend.title      = element_text(size = 18)
    )

  # 0.8 threshold overlay (cyan outline on cells reaching 80% recall).
  df <- df %>%
    mutate(Threshold_Met = factor(
      ifelse(!is.na(Mean_Recall) & Mean_Recall >= 0.8, "\u2265 80% Recall", "< 80% Recall"),
      levels = c("\u2265 80% Recall", "< 80% Recall")
    ))
  any_met <- any(!is.na(df$Mean_Recall) & df$Mean_Recall >= 0.8)

  p <- p +
    new_scale_color() +
    geom_tile(data = df, aes(color = Threshold_Met), fill = NA, linewidth = 0.6,
              width = 0.85, height = 0.85) +
    scale_color_manual(
      name = "Threshold",
      values = setNames(c("#00FFFF", "grey70"), c("\u2265 80% Recall", "< 80% Recall")),
      limits = c("\u2265 80% Recall", "< 80% Recall"), drop = FALSE,
      guide = if (any_met) {
        guide_legend(order = 2, override.aes = list(fill = c(NA, "grey85"),
                     color = c("#00FFFF", "grey70"), linewidth = 1.2))
      } else "none"
    )

  # Adaptive value labels (dark on light tiles, white on dark).
  df$.tc <- ifelse(is.na(df$Mean_Recall) | df$Mean_Recall > 0.45, "black", "white")
  for (tc in unique(df$.tc)) {
    sub <- df[df$.tc == tc, ]
    if (nrow(sub) == 0) next
    p <- p + geom_text(data = sub,
                       aes(label = ifelse(is.na(Mean_Recall), "", round(Mean_Recall * 100))),
                       color = tc, size = 2.8, show.legend = FALSE)
  }
  p
}

# ---------------------------------------------------------------------------
# Panel A: recall by padj quartile, pooled across DE methods (mirrors 04)
# ---------------------------------------------------------------------------

build_panel_a <- function() {
  parts <- list()
  for (method in DE_METHODS) {
    cp_file <- checkpoint_path(CELL, method)
    if (!file.exists(cp_file)) { msg("[A] missing checkpoint: %s", cp_file); next }
    cp <- qs::qread(cp_file)
    dq <- cp$df_pval_q_summary
    if (is.null(dq) || nrow(dq) == 0) { msg("[A] no df_pval_q_summary for %s", method); next }
    parts[[method]] <- dq
  }
  if (length(parts) == 0) stop("Panel A: no per-method p-value-quartile summaries found.")

  df_all <- bind_rows(parts)

  # Trial-agnostic simple mean across methods, exactly as 04 does for the
  # aggregated all-methods pval-quartile plot.
  df_agg <- df_all %>%
    group_by(Donors, .data[[Y_COL]], PVal_Quartile) %>%
    summarise(
      Mean_Recall = mean(Mean_Recall, na.rm = TRUE),
      PQ_Lo       = if ("PQ_Lo" %in% names(df_all)) mean(PQ_Lo, na.rm = TRUE) else NA_real_,
      PQ_Hi       = if ("PQ_Hi" %in% names(df_all)) mean(PQ_Hi, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  # Relabel Q1-Q4 with mean padj bounds, inner boundaries only (matches 04).
  if (all(c("PQ_Lo", "PQ_Hi") %in% names(df_agg)) && !all(is.na(df_agg$PQ_Lo))) {
    bounds <- df_agg %>%
      group_by(PVal_Quartile) %>%
      summarise(PQ_Lo = mean(PQ_Lo, na.rm = TRUE),
                PQ_Hi = mean(PQ_Hi, na.rm = TRUE), .groups = "drop")
    lab <- setNames(
      mapply(function(nm, lo, hi) {
        if (grepl("^Q1", nm))      sprintf("%s (padj \u2264 %.1e)", nm, hi)
        else if (grepl("^Q4", nm)) sprintf("%s (padj \u2265 %.1e)", nm, lo)
        else                       sprintf("%s (%.1e\u2013%.1e)", nm, lo, hi)
      }, bounds$PVal_Quartile, bounds$PQ_Lo, bounds$PQ_Hi),
      bounds$PVal_Quartile
    )
    df_agg$PVal_Quartile <- lab[df_agg$PVal_Quartile]
  }

  # Q1 first (top-left, most significant / highest recall) -> Q4 last
  # (bottom-right), so the strongest panel reads top-left for consistency.
  facet_levs <- sort(unique(df_agg$PVal_Quartile), decreasing = FALSE)
  build_recall_grid(df_agg, "PVal_Quartile", facet_levs,
                    "Recall by adjusted p-value quartile")
}

# ---------------------------------------------------------------------------
# Panels B/C: recall across expression within a fixed LFC bin (all_methods)
# ---------------------------------------------------------------------------

build_finecross_panel <- function(lfc_bin) {
  cp_file <- checkpoint_path(CELL, CROSS_METHOD_LABEL)
  if (!file.exists(cp_file)) stop(sprintf("all_methods checkpoint missing: %s", cp_file))
  cp <- qs::qread(cp_file)
  df_fc <- cp$df_finecross_summary
  if (is.null(df_fc) || nrow(df_fc) == 0) {
    stop("df_finecross_summary absent from all_methods checkpoint (re-run Step 8).")
  }

  df_fc <- df_fc %>%
    filter(Gene_Category == "Total", LFC_Bin == lfc_bin) %>%
    mutate(
      Donors = clean_donors_factor(Donors),
      !!Y_COL := clean_y_factor(.data[[Y_COL]])
    ) %>%
    complete(Donors, !!sym(Y_COL), BaseMean_Bin)

  lfc_titles <- c(
    "Tiny (0.2-1.0)" = "Recall across gene counts, tiny effect sizes (|LFC| 0.2\u20131.0)",
    "Small (1.0-2.0)" = "Recall across gene counts, small effect sizes (|LFC| 1.0\u20132.0)"
  )
  build_recall_grid(df_fc, "BaseMean_Bin", BASEMEAN_ORDER, unname(lfc_titles[lfc_bin]))
}

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------

msg("Cell type: %s | mode: %s | datalink: %s", CELL, MODE, DATALINK)

pA <- build_panel_a()
pB <- build_finecross_panel(FINECROSS_LFC_BINS[1])  # Tiny (0.2-1.0)
pC <- build_finecross_panel(FINECROSS_LFC_BINS[2])  # Small (1.0-2.0)

combined <- (pA / pB / pC) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

base_name <- sprintf("recall_panels_%s_%s", tolower(CELL), MODE)
png_path  <- file.path(OUT_DIR, paste0(base_name, ".png"))
pdf_path  <- file.path(OUT_DIR, paste0(base_name, ".pdf"))

ggsave(png_path, combined, width = 13, height = 24, dpi = 150, limitsize = FALSE)
ggsave(pdf_path, combined, width = 13, height = 24, limitsize = FALSE)

msg("Saved: %s", png_path)
msg("Saved: %s", pdf_path)

# Canonical draft figure: emit into Datalink/Composed_Figures so the manuscript
# chunk (fig:FigRecallPanels) always reads a freshly-rendered file on re-render.
COMPOSED_DIR <- file.path(DATALINK, COMPOSED_SUBDIR)
dir.create(COMPOSED_DIR, recursive = TRUE, showWarnings = FALSE)
draft_png <- file.path(COMPOSED_DIR, sprintf("figure_recall_panels_%s_%s.png", tolower(CELL), MODE))
draft_pdf <- file.path(COMPOSED_DIR, sprintf("figure_recall_panels_%s_%s.pdf", tolower(CELL), MODE))
ggsave(draft_png, combined, width = 13, height = 24, dpi = 150, limitsize = FALSE)
ggsave(draft_pdf, combined, width = 13, height = 24, limitsize = FALSE)
msg("Draft figure: %s", draft_png)
msg("Draft figure: %s", draft_pdf)
msg("Done.")
