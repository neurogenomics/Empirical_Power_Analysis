# ==============================================================================
# SCRIPT: Figure 5 panels for the XGBoost power model
# PURPOSE: Draw the paper-styled panels from the CSV / OOF outputs of
#          03_Select_XGBoost.R, 08_Compare_Reweighted.R and
#          09_Trinary_Classification.R. No model refit. Reads:
#   trinary_classification/per_class_metrics.csv     -> precision / recall / F1 bars
#   shap_selection.csv + search/config_<k>/oof_predictions.qs -> calibration grid
#   comparison/shap_importance_across_models.csv     -> SHAP heatmap
#   model_selection.csv, model_selection_folds.csv   -> hyperparameter bracketing
# Writes into Models/xgboost/comparison/_panel_preview/ (no titles; the
# composite figure is assembled by R/figure_composition/compose_xgb_panels.R).
#
# Project root: $POWER_ANALYSIS_ROOT, else two levels up from the working
# directory (the script used to be run from R/power_curve_fitting/).
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT    <- Sys.getenv("POWER_ANALYSIS_ROOT", "../..")
DL      <- file.path(ROOT, "Datalink", "Power_Curve_Fitting", "Models", "xgboost")
out_dir <- file.path(DL, "comparison", "_panel_preview")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

shap_sel <- read.csv(file.path(DL, "shap_selection.csv"), stringsAsFactors = FALSE)

# Pretty labels for the 5 SHAP models, built from the selection itself so they
# stay right whichever configurations the search picks.
tag_pretty <- function(tag, alpha) {
  base <- c(median = "Median config", worst = "Worst config")
  ifelse(startsWith(tag, "winner_alpha_"),
         sprintf("Winner (α = %g)", alpha),
         sprintf("%s (α = %g)", base[tag], alpha))
}
shap_sel$pretty <- tag_pretty(shap_sel$tag, shap_sel$alpha)
tag_order <- c(sort(grep("^winner_alpha_", shap_sel$tag, value = TRUE)), "median", "worst")
shap_sel <- shap_sel[order(match(shap_sel$tag, tag_order)), ]

# ------------------------------------------------------------------------------
# PANEL C: precision / recall / F1 per class for the winning model, one facet
# per partition scheme.
# ------------------------------------------------------------------------------
pc <- read.csv(file.path(DL, "trinary_classification", "per_class_metrics.csv"),
               stringsAsFactors = FALSE)

# class label written by 09_Trinary_Classification.R -> paper label
class_label_map <- c(
  "Not Powered"         = "Unpowered",
  "Definitely Not"      = "Unpowered",
  "Underpowered"        = "Unpowered",
  "Likely Underpowered" = "Probably not",
  "Contested"           = "Contested",
  "Maybe"               = "Contested",
  "Likely Powered"      = "Probably yes",
  "Powered"             = "Powered",
  "Yes"                 = "Powered",
  "Well Powered"        = "Powered",
  "Q1: Lowest"          = "Lowest",
  "Q2: Low-Mid"         = "Low-mid",
  "Q3: High-Mid"        = "High-mid",
  "Q4: Highest"         = "Highest"
)
# global low -> high order that renders every facet correctly under free_x
class_levels <- c("Unpowered", "Probably not", "Contested", "Probably yes",
                  "Powered", "Lowest", "Low-mid", "High-mid", "Highest")

scheme_order <- c("binary", "trinary_01_09", "trinary_02_08",
                  "trinary_04_08", "quartile", "four_class")
scheme_label_map <- c(
  binary        = "Binary\n[0–0.5] [0.5–1]",
  trinary_01_09 = "Trinary strict\n[0–0.1] [0.1–0.9] [0.9–1]",
  trinary_02_08 = "Trinary standard\n[0–0.2] [0.2–0.8] [0.8–1]",
  trinary_04_08 = "Trinary skewed\n[0–0.4] [0.4–0.8] [0.8–1]",
  quartile      = "Quartiles\n[0–0.25] [0.25–0.5] [0.5–0.75] [0.75–1]",
  four_class    = "Four-class\n[0–0.2] [0.2–0.5] [0.5–0.8] [0.8–1]"
)

unmapped <- setdiff(unique(pc$class), names(class_label_map))
if (length(unmapped)) stop("Unmapped class labels in per_class_metrics.csv: ",
                           paste(unmapped, collapse = ", "))

prf_df <- pc %>%
  filter(role == "winner") %>%
  pivot_longer(c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  mutate(
    class_lbl  = factor(class_label_map[class], levels = class_levels),
    scheme_lbl = factor(scheme_label_map[scheme], levels = scheme_label_map[scheme_order]),
    metric     = factor(metric, levels = c("precision", "recall", "f1"),
                        labels = c("Precision", "Recall", "F1"))
  )

p_prf <- ggplot(prf_df, aes(x = class_lbl, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.8), vjust = -0.3,
            size = 2.6, colour = "black", na.rm = TRUE) +
  facet_wrap(~ scheme_lbl, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c(Precision = "#0072B2", Recall = "#E69F00", F1 = "#009E73")) +
  labs(x = NULL, y = "Score", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position  = "top",
        legend.text      = element_text(size = 13, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.y     = element_text(size = 14, colour = "black"),
        axis.text.y      = element_text(size = 12, colour = "black"),
        axis.text.x      = element_text(angle = 25, hjust = 1, size = 11, colour = "black"),
        strip.text       = element_text(size = 11, colour = "black", lineheight = 0.9))

ggsave(file.path(out_dir, "class_prf_grid.png"), p_prf, width = 8.5, height = 10, dpi = 150)
cat("Wrote class_prf_grid.png\n")

# ------------------------------------------------------------------------------
# PANEL A: OOF calibration curves for the 5 SHAP models. Each point (per
# predicted bin) carries a bar spanning the trial-weighted 10-90th percentile
# of the OBSERVED detection rate within that bin: wide in the middle (mid-power
# outcomes are intrinsically volatile), narrow at the extremes.
# ------------------------------------------------------------------------------
wtd_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w)
  p  <- (cw - 0.5 * w) / sum(w)
  stats::approx(p, x, xout = probs, rule = 2, ties = "ordered")$y
}

bin_oof <- function(oof) {
  oof %>%
    mutate(pred_bin = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
    group_by(pred_bin) %>%
    summarise(
      mean_pred = weighted.mean(pred, w_trials),
      mean_obs  = weighted.mean(obs,  w_trials),
      obs_p10   = wtd_quantile(obs, w_trials, 0.10),
      obs_p90   = wtd_quantile(obs, w_trials, 0.90),
      .groups   = "drop"
    ) %>%
    filter(is.finite(mean_pred), is.finite(mean_obs))
}

cal <- bind_rows(lapply(seq_len(nrow(shap_sel)), function(i) {
  f <- file.path(DL, "search", sprintf("config_%d", shap_sel$config_index[i]),
                 "oof_predictions.qs")
  if (!file.exists(f)) stop("OOF predictions missing: ", f)
  cc <- bin_oof(qs::qread(f)); cc$model <- shap_sel$pretty[i]; cc
}))
cal$model <- factor(cal$model, levels = shap_sel$pretty)

p_cal <- ggplot(cal, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = obs_p10, ymax = obs_p90),
                width = 0.02, colour = "steelblue", alpha = 0.55, linewidth = 0.4) +
  geom_line(colour = "steelblue", alpha = 0.6) +
  geom_point(size = 1.8, colour = "steelblue") +
  facet_wrap(~ model, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted detection rate (per bin)",
       y = "Observed detection rate (mean point, 10-90th percentile bar)") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text  = element_text(size = 12, colour = "black"),
        strip.text = element_text(size = 11.5, colour = "black"))

ggsave(file.path(out_dir, "calibration_curves_2x3.png"),
       p_cal, width = 7.2, height = 10, dpi = 150)
cat("Wrote calibration_curves_2x3.png\n")

# ------------------------------------------------------------------------------
# PANEL B: SHAP importance heatmap across the 5 models.
# ------------------------------------------------------------------------------
sha <- read.csv(file.path(DL, "comparison", "shap_importance_across_models.csv"),
                stringsAsFactors = FALSE)

feature_label_map <- c(
  abs_gt_lfc         = "|LFC|",
  n_donors           = "Donors",
  gt_basemean        = "Gene counts",
  gt_dispersion      = "Dispersion",
  n_cells            = "Cells per donor",
  reads_per_cell     = "UMIs per cell",
  `de_method=deseq2` = "DE method: DESeq2",
  `de_method=edger`  = "DE method: edgeR",
  `de_method=limma`  = "DE method: limma",
  `cell_type=Astro`  = "Cell type: Astro",
  `cell_type=Micro`  = "Cell type: Micro",
  `cell_type=Oligo`  = "Cell type: Oligo"
)
heat_tag_label <- setNames(sub("^(Median|Worst) config", "\\1", shap_sel$pretty), shap_sel$tag)

feat_order <- sha %>%
  group_by(feature) %>%
  summarise(avg_rank = mean(rank), .groups = "drop") %>%
  arrange(avg_rank) %>%
  pull(feature)

sha2 <- sha %>%
  mutate(
    feat_lbl = factor(feature_label_map[feature],
                      levels = rev(feature_label_map[feat_order])),
    tag_lbl  = factor(heat_tag_label[tag], levels = heat_tag_label[shap_sel$tag])
  )

p_shap <- ggplot(sha2, aes(x = tag_lbl, y = feat_lbl, fill = mean_abs_shap)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_abs_shap)),
            colour = "white", size = 3.6) +
  scale_fill_viridis_c(name = "mean|SHAP|") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 12,
                                   colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.position = "right")

ggsave(file.path(out_dir, "shap_importance_pretty.png"),
       p_shap, width = 8.2, height = 6.4, dpi = 150)
cat("Wrote shap_importance_pretty.png\n")

# ------------------------------------------------------------------------------
# PANEL B: hyperparameter bracketing. For each value of each searched
# hyperparameter, take the best configuration with that value (minimum
# bin-balanced MAE over every other axis) and plot how much worse it is than
# the selected configuration. Every configuration is scored on the same
# gene-grouped CV folds, so the difference and its standard error are
# fold-PAIRED (unpaired fold SEs are larger than the gaps between configs).
# The label above each point counts the folds on which it did worse than the
# selected configuration. Each hyperparameter gets its own y-scale so the large
# alpha effect does not flatten the others. If the selected value sits strictly
# inside the searched range, the grid brackets the optimum on that axis.
# Drawn narrow (half the figure width) to sit above panel C in the composite.
# ------------------------------------------------------------------------------
ms    <- read.csv(file.path(DL, "model_selection.csv"), stringsAsFactors = FALSE)
folds <- read.csv(file.path(DL, "model_selection_folds.csv"), stringsAsFactors = FALSE)
win   <- ms[which.min(ms$bin_balanced_MAE), ]
win_f <- folds[folds$config_index == win$config_index, c("fold", "bin_balanced_MAE")]
paired_diff <- function(cfg) {
  f <- merge(folds[folds$config_index == cfg, c("fold", "bin_balanced_MAE")], win_f,
             by = "fold", suffixes = c("", "_win"))
  d <- f$bin_balanced_MAE - f$bin_balanced_MAE_win
  c(diff = mean(d), se = if (cfg == win$config_index) 0 else sd(d) / sqrt(nrow(f)),
    worse = sum(d > 0), n = nrow(f))
}

hp_cols <- c(max_depth = "Max tree depth", eta = "Learning rate",
             min_child_weight = "Min child weight", alpha = "Reweighting \u03b1")
bracket <- bind_rows(lapply(names(hp_cols), function(h) {
  ms %>%
    group_by(value = .data[[h]]) %>%
    slice_min(bin_balanced_MAE, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(hp = hp_cols[[h]], value, config_index, selected = value == win[[h]])
}))
bracket <- bind_cols(bracket,
                     as.data.frame(t(vapply(bracket$config_index, paired_diff, numeric(4))))) %>%
  mutate(hp = factor(hp, levels = hp_cols),
         d1000 = diff * 1000, se1000 = se * 1000,
         lbl = ifelse(selected, "", sprintf("%d/%d", worse, n)))
bracket$value_lbl <- factor(format(bracket$value, drop0trailing = TRUE, trim = TRUE),
                            levels = unique(format(sort(unique(bracket$value)),
                                                   drop0trailing = TRUE, trim = TRUE)))
# Headroom per facet so the fold-count labels are never clipped.
headroom <- bracket %>%
  group_by(hp) %>%
  summarise(top = max(d1000 + se1000) * 1.3 + 1e-3, value_lbl = value_lbl[1])

p_br <- ggplot(bracket, aes(x = value_lbl, y = d1000, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_blank(data = headroom, aes(x = value_lbl, y = top), inherit.aes = FALSE) +
  geom_line(colour = "grey50") +
  geom_errorbar(data = function(d) d[d$se > 0, ],
                aes(ymin = d1000 - se1000, ymax = d1000 + se1000),
                width = 0.2, colour = "grey30") +
  geom_point(aes(fill = selected), shape = 21, size = 2.8, colour = "black",
             show.legend = FALSE) +
  geom_text(aes(y = d1000 + se1000, label = lbl), vjust = -0.6, size = 3.3,
            colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white")) +
  facet_wrap(~ hp, nrow = 1, scales = "free") +
  labs(x = NULL, y = "Extra MAE vs chosen\nmodel (\u00d71000)") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"),
        plot.margin = margin(4, 6, 2, 4))

ggsave(file.path(out_dir, "hyperparameter_bracketing.png"),
       p_br, width = 8.5, height = 3.0, dpi = 150)
cat("Wrote hyperparameter_bracketing.png\n")
cat("DONE\n")
