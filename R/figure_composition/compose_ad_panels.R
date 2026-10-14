# =============================================================================
# SCRIPT: compose_ad_panels.R
# PURPOSE: Composed figure for the AD vs control contrast (manuscript end
#          section), method-pooled (all_methods), cells view, one cell type:
#            A  Recall and precision of tiny-effect genes (|LFC| 0.2-1.0), side
#               by side. Nearly every AD vs control DEG with |LFC| > 0.2 sits in
#               this bin; the small/medium/large bins are near-empty.
#            B  Recall by ground-truth adjusted p-value quartile (as the sex
#               figure fig:FigRecallPanels panel A)
#            C  Recall across mean expression, tiny bin (as fig:FigRecallPanels
#               panel B)
#          Same grid styling as compose_recall_panels.R: shared 0-1 magma scale,
#          cyan borders at >= 80%, grey for unobserved design points, whole
#          facet greyed if its mean truth-gene count < EMPTY_FACET_GENE_THRESHOLD.
#
# SOURCES (cartesian checkpoints, Plots/Aggregated_Summary_Data.qs):
#   A  all_methods df_summary (LFC_Bin == Tiny, Gene_Category == Total)
#   B  per-method df_pval_q_summary, simple mean across deseq2/edger/limma
#   C  all_methods df_finecross_summary (LFC_Bin == Tiny, Gene_Category == Total)
#
# OUTPUT: Datalink/<composed-dir>/figure_ad_panels_<cell>_cells.{png,pdf}
#
# Run:
#   Rscript R/figure_composition/compose_ad_panels.R --cell Oligo
#   Rscript R/figure_composition/compose_ad_panels.R --cell Oligo \
#     --dataset Roussos_ADvsCtrl --composed-dir AD_vs_Control/Composed_Figures
# =============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(ggnewscale)
})

HERE <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE))
if (is.na(HERE) || HERE == "") HERE <- normalizePath(".", mustWork = FALSE)
PROJECT_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)

source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))
Sys.setenv(POWER_ANALYSIS_ROOT = PROJECT_ROOT)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  pos <- which(args == flag)
  if (length(pos) && pos[1] < length(args)) return(args[pos[1] + 1])
  default
}

DATALINK        <- normalizePath(get_arg("--datalink", file.path(PROJECT_ROOT, "Datalink")),
                                 mustWork = FALSE)
CELL            <- get_arg("--cell", "Oligo")
DATASET         <- get_arg("--dataset", "Roussos_ADvsCtrl")
COMPOSED_SUBDIR <- get_arg("--composed-dir", "AD_vs_Control/Composed_Figures")
MODE            <- "cells"

# config.R checks SHARD_DATASET against EPA_CONTRAST, so set both to match the
# grid being drawn before sourcing it (for the axis steps and cleaners).
Sys.setenv(SHARD_DATASET = DATASET,
           EPA_CONTRAST  = if (DATASET == "Roussos_ADvsCtrl") "ad_vs_control" else "sex")
source(file.path(PROJECT_ROOT, "R", "cartesian_downsampling", "config.R"))
source(file.path(PROJECT_ROOT, "R", "cartesian_downsampling", "utils.R"))

DE_METHODS     <- c("deseq2", "edger", "limma")
TINY_BIN       <- "Tiny (0.2-1.0)"
BASEMEAN_ORDER <- c("High (200+)", "Medium (50-200)", "Low (10-50)", "Very Low (<10)")

checkpoint <- function(method) {
  f <- file.path(DATALINK, "Cartesian_Downsampling", CELL, DATASET, MODE, method,
                 "Plots", "Aggregated_Summary_Data.qs")
  if (!file.exists(f)) stop("missing checkpoint: ", f)
  qs::qread(f)
}

# Grid of one metric (Mean_Recall or Mean_Precision) over donors x cells per
# donor, one facet per level of facet_col. The 80% legend is drawn only by grids
# with a tile at >= 80% (as compose_recall_panels.R), so collected guides merge.
build_grid <- function(df, metric, facet_col, facet_levs, title_str, ncol = 2) {
  df <- df %>%
    mutate(Donors = clean_donors_factor(Donors),
           Cells  = clean_cells_factor(Cells),
           .facet = factor(.data[[facet_col]], levels = facet_levs)) %>%
    filter(!is.na(.facet)) %>%
    tidyr::complete(Donors, Cells, .facet)
  df$.value <- df[[metric]]

  if ("Ref_GT_Count" %in% names(df)) {
    m <- tapply(df$Ref_GT_Count, df$.facet, function(v) mean(v, na.rm = TRUE))
    empty <- names(m)[is.nan(m) | m < EMPTY_FACET_GENE_THRESHOLD]
    df$.value[as.character(df$.facet) %in% empty] <- NA_real_
  }

  met_lab <- c("≥ 80%", "< 80%")
  df$Threshold_Met <- factor(ifelse(!is.na(df$.value) & df$.value >= 0.8, met_lab[1], met_lab[2]),
                             levels = met_lab)
  df$.tc <- ifelse(is.na(df$.value) | df$.value > 0.45, "black", "white")

  ggplot(df, aes(x = Donors, y = Cells)) +
    geom_tile(aes(fill = .value), color = "white") +
    facet_wrap(~.facet, ncol = ncol) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50",
                         name = "Recall or\nprecision", guide = guide_colorbar(order = 1)) +
    new_scale_color() +
    geom_tile(aes(color = Threshold_Met), fill = NA, linewidth = 0.6,
              width = 0.85, height = 0.85) +
    scale_color_manual(name = "Threshold", values = setNames(c("#00FFFF", "grey70"), met_lab),
                       limits = met_lab, drop = FALSE,
                       guide = if (any(df$.value >= 0.8, na.rm = TRUE)) guide_legend(order = 2, override.aes = list(
                         fill = c(NA, "grey85"), color = c("#00FFFF", "grey70"), linewidth = 1.2))
                               else "none") +
    new_scale_color() +
    geom_text(aes(label = ifelse(is.na(.value), "", round(.value * 100)), color = .tc),
              size = 2.8, show.legend = FALSE) +
    scale_color_identity() +
    labs(title = title_str, x = "Donors", y = "Cells per Donor") +
    theme_publication(base_size = 14) +
    theme(plot.title    = element_text(face = "bold", size = 16),
          axis.text.x   = element_text(angle = 45, hjust = 1),
          panel.grid    = element_blank(),
          strip.text    = element_text(face = "bold", size = 12),
          panel.spacing = unit(0.1, "lines"),
          legend.position   = "right",
          legend.key.width  = unit(0.7, "cm"),
          legend.key.height = unit(1.3, "cm"),
          legend.text       = element_text(size = 16),
          legend.title      = element_text(size = 18))
}

cat(sprintf("Cell type: %s | dataset: %s | datalink: %s\n", CELL, DATASET, DATALINK))
cp_all <- checkpoint("all_methods")

# ---- A: tiny-bin recall and precision -------------------------------------
df_tiny <- cp_all$df_summary %>%
  filter(Gene_Category == "Total", LFC_Bin == TINY_BIN) %>%
  mutate(.bin = "Tiny effect sizes (|LFC| 0.2–1.0)")
pA_rec  <- build_grid(df_tiny, "Mean_Recall",    ".bin", unique(df_tiny$.bin), "Recall", ncol = 1)
pA_prec <- build_grid(df_tiny, "Mean_Precision", ".bin", unique(df_tiny$.bin), "Precision", ncol = 1)

# ---- B: recall by adjusted p-value quartile (pooled as in 04) --------------
dq <- bind_rows(lapply(DE_METHODS, function(m) checkpoint(m)$df_pval_q_summary))
dq <- dq %>%
  group_by(Donors, Cells, PVal_Quartile) %>%
  summarise(Mean_Recall = mean(Mean_Recall, na.rm = TRUE),
            PQ_Lo = mean(PQ_Lo, na.rm = TRUE), PQ_Hi = mean(PQ_Hi, na.rm = TRUE),
            .groups = "drop")
bounds <- dq %>% group_by(PVal_Quartile) %>%
  summarise(PQ_Lo = mean(PQ_Lo, na.rm = TRUE), PQ_Hi = mean(PQ_Hi, na.rm = TRUE), .groups = "drop")
q_lab <- setNames(mapply(function(nm, lo, hi) {
  if (grepl("^Q1", nm))      sprintf("%s (padj ≤ %.1e)", nm, hi)
  else if (grepl("^Q4", nm)) sprintf("%s (padj ≥ %.1e)", nm, lo)
  else                       sprintf("%s (%.1e–%.1e)", nm, lo, hi)
}, bounds$PVal_Quartile, bounds$PQ_Lo, bounds$PQ_Hi), bounds$PVal_Quartile)
dq$PVal_Quartile <- q_lab[dq$PVal_Quartile]
pB <- build_grid(dq, "Mean_Recall", "PVal_Quartile", sort(unique(dq$PVal_Quartile)),
                 "Recall by adjusted p-value quartile")

# ---- C: tiny-bin recall across mean expression -----------------------------
df_fc <- cp_all$df_finecross_summary %>%
  filter(Gene_Category == "Total", LFC_Bin == TINY_BIN)
pC <- build_grid(df_fc, "Mean_Recall", "BaseMean_Bin", BASEMEAN_ORDER,
                 "Recall across expression, tiny effect sizes (|LFC| 0.2–1.0)")

combined <- ((pA_rec | pA_prec) / pB / pC) +
  plot_layout(guides = "collect", heights = c(1, 2, 2)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "C"))) &
  theme(plot.tag = element_text(face = "bold", size = 18))

out_dir <- file.path(DATALINK, COMPOSED_SUBDIR)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(out_dir, sprintf("figure_ad_panels_%s_%s", tolower(CELL), MODE))
ggsave(paste0(stem, ".png"), combined, width = 13, height = 22, dpi = 150, limitsize = FALSE)
ggsave(paste0(stem, ".pdf"), combined, width = 13, height = 22, limitsize = FALSE)
cat("Saved:", paste0(stem, c(".png", ".pdf")), sep = "\n")
