# ==============================================================================
# UTILITY (preview only): regenerate the two XGBoost paper panels from the
# already-committed CSVs -- no model refit, so the numbers stay identical to
# what the draft already cites. Reads:
#   trinary_classification/per_class_metrics.csv   -> F1 grid (reordered, relabelled)
#   comparison/unified_calibration_table.csv       -> calibration grid (2x3)
#   shap_selection.csv                             -> calibration model set + labels
# Writes previews into comparison/_panel_preview/. Titles/subtitles stripped.
# Once approved, this label/order logic gets folded into 08_ and 09_.
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
})

# Paths are relative to this script's directory (R/power_curve_fitting/).
DL      <- "../../Datalink/Power_Curve_Fitting/Models/xgboost/realistic"
per_class_path <- file.path(DL, "trinary_classification", "per_class_metrics.csv")
unified_path   <- file.path(DL, "comparison", "unified_calibration_table.csv")
shap_path      <- file.path(DL, "shap_selection.csv")
out_dir        <- file.path(DL, "comparison", "_panel_preview")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# PANEL: F1 grid (winner vs baseline), reordered + pretty labels, no underscores
# ------------------------------------------------------------------------------
pc <- read.csv(per_class_path, stringsAsFactors = FALSE)

# data class value -> pretty label (many-to-one; distinct across schemes)
class_label_map <- c(
  unpowered      = "Unpowered",
  definitely_not = "Unpowered",
  probably_not   = "Probably not",
  contested      = "Contested",
  maybe          = "Contested",
  probably_yes   = "Probably yes",
  powered        = "Powered",
  yes            = "Powered",
  q1_low         = "Lowest",
  q2_low_mid     = "Low-mid",
  q3_high_mid    = "High-mid",
  q4_high        = "Highest",
  # class labels written by the current 09_Trinary_Classification.R
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
  binary        = "Binary\n[0\u20130.5] [0.5\u20131]",
  trinary_01_09 = "Trinary strict\n[0\u20130.1] [0.1\u20130.9] [0.9\u20131]",
  trinary_02_08 = "Trinary standard\n[0\u20130.2] [0.2\u20130.8] [0.8\u20131]",
  trinary_04_08 = "Trinary skewed\n[0\u20130.4] [0.4\u20130.8] [0.8\u20131]",
  quartile      = "Quartiles\n[0\u20130.25] [0.25\u20130.5] [0.5\u20130.75] [0.75\u20131]",
  four_class    = "Four-class\n[0\u20130.2] [0.2\u20130.5] [0.5\u20130.8] [0.8\u20131]"
)

plot_df <- pc %>%
  filter(role %in% c("winner", "baseline")) %>%
  mutate(
    class_lbl   = factor(class_label_map[class], levels = class_levels),
    scheme      = factor(scheme, levels = scheme_order),
    scheme_lbl  = factor(scheme_label_map[as.character(scheme)],
                         levels = scheme_label_map[scheme_order]),
    model       = factor(model, levels = c("realistic_winner", "baseline"),
                         labels = c("Realistic winner", "Baseline"))
  )

p_f1 <- ggplot(plot_df, aes(x = class_lbl, y = f1, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.2f", f1)),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = 3.2, colour = "black", na.rm = TRUE) +
  facet_wrap(~ scheme_lbl, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c("Realistic winner" = "#F8766D",
                               "Baseline" = "#00BFC4")) +
  labs(x = NULL, y = "F1 score", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position  = "top",
        legend.text      = element_text(size = 13, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.y     = element_text(size = 14, colour = "black"),
        axis.text.y      = element_text(size = 12, colour = "black"),
        axis.text.x      = element_text(angle = 25, hjust = 1, size = 11,
                                        colour = "black"),
        strip.text       = element_text(size = 11, colour = "black",
                                        lineheight = 0.9))

ggsave(file.path(out_dir, "class_f1_grid_v2.png"),
       p_f1, width = 7.2, height = 10, dpi = 150)
cat("Wrote class_f1_grid_v2.png\n")

# ------------------------------------------------------------------------------
# PANEL: calibration curves, 2x3 (ncol=2), pretty facet labels, no title.
# Each point (per predicted bin) carries a vertical error bar spanning the
# trial-weighted 10-90th percentile of the OBSERVED detection rate within that
# bin -- wide in the middle (mid-power outcomes are intrinsically volatile),
# narrow at the extremes (near 0/1 the outcome is near-deterministic). This is
# the spread the mean-only point conceals. Bands are computed from the per-row
# OOF predictions (search/config_<k>/oof_predictions.qs for the 5 configs;
# Models/xgboost/oof_predictions.qs for the baseline, regenerated by
# 08b_Baseline_OOF.R). Weighting = n_trials, matching mean_obs.
# ------------------------------------------------------------------------------
shap_sel <- read.csv(shap_path, stringsAsFactors = FALSE)

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

chosen <- c("baseline",
            sprintf("realistic_c%d_a%g", shap_sel$config_index, shap_sel$alpha))

pretty_map <- c(baseline = "Baseline (full grid)")
tag_pretty <- c(winner_alpha_0 = "Winner (\u03b1 = 0)",
                winner_alpha_0.5 = "Winner (\u03b1 = 0.5)",
                winner_alpha_1 = "Winner (\u03b1 = 1)",
                median = "Median config (\u03b1 = 0.5)",
                worst  = "Worst config (\u03b1 = 1)")
oof_map <- c(baseline = file.path(DL, "..", "oof_predictions.qs"))
for (i in seq_len(nrow(shap_sel))) {
  key <- sprintf("realistic_c%d_a%g", shap_sel$config_index[i], shap_sel$alpha[i])
  pretty_map[key] <- tag_pretty[[shap_sel$tag[i]]]
  oof_map[key] <- file.path(DL, "search",
                            sprintf("config_%d", shap_sel$config_index[i]),
                            "oof_predictions.qs")
}

cal <- bind_rows(lapply(chosen, function(key) {
  f <- oof_map[key]
  if (!file.exists(f)) { cat(sprintf("WARN: OOF missing for %s (%s)\n", key, f)); return(NULL) }
  cc <- bin_oof(qs::qread(f)); cc$model <- key; cc
}))
present <- chosen[chosen %in% unique(cal$model)]
cal$model <- factor(cal$model, levels = present, labels = pretty_map[present])

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
# PANEL: SHAP importance heatmap across 5 models, pretty labels, no underscores,
# no title/subtitle. Reads comparison/shap_importance_across_models.csv.
# ------------------------------------------------------------------------------
shap_across_path <- file.path(DL, "comparison", "shap_importance_across_models.csv")
sha <- read.csv(shap_across_path, stringsAsFactors = FALSE)

feature_label_map <- c(
  abs_gt_lfc         = "|LFC|",
  n_donors           = "Donors",
  gt_basemean        = "Mean expression",
  gt_dispersion      = "Dispersion",
  n_cells            = "Cells per donor",
  reads_per_cell     = "Counts per cell",
  `de_method=deseq2` = "DE method: DESeq2",
  `de_method=edger`  = "DE method: edgeR",
  `de_method=limma`  = "DE method: limma",
  `cell_type=Astro`  = "Cell type: Astro",
  `cell_type=Micro`  = "Cell type: Micro",
  `cell_type=Oligo`  = "Cell type: Oligo"
)
tag_order <- c("winner_alpha_0", "winner_alpha_0.5", "winner_alpha_1",
               "median", "worst")
tag_label_map <- c(
  winner_alpha_0   = "Winner (\u03b1 = 0)",
  winner_alpha_0.5 = "Winner (\u03b1 = 0.5)",
  winner_alpha_1   = "Winner (\u03b1 = 1)",
  median           = "Median",
  worst            = "Worst"
)

# feature order: ascending average rank (most important at top)
feat_order <- sha %>%
  group_by(feature) %>%
  summarise(avg_rank = mean(rank), .groups = "drop") %>%
  arrange(avg_rank) %>%
  pull(feature)

sha2 <- sha %>%
  mutate(
    feat_lbl = factor(feature_label_map[feature],
                      levels = rev(feature_label_map[feat_order])),
    tag_lbl  = factor(tag_label_map[tag],
                      levels = tag_label_map[tag_order])
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
# PANEL: bin-balanced MAE across 12 configs x 3 alpha, no title/subtitle.
# Reads model_selection.csv. Config labels already underscore-free.
# ------------------------------------------------------------------------------
ms <- read.csv(file.path(DL, "model_selection.csv"), stringsAsFactors = FALSE)

ms2 <- ms %>%
  mutate(
    cfg_lbl = sprintf("depth %d\nlr %g \u00b7 mcw %d", max_depth, eta, min_child_weight),
    alpha_lbl = factor(alpha, levels = c(0, 0.5, 1),
                       labels = c("0  (no reweighting)",
                                  "0.5  (inverse \u221an)",
                                  "1  (inverse n)"))
  )
# order configs by the alpha=0 MAE (best first), matching the existing figure
cfg_levels <- ms2 %>%
  filter(alpha == 0) %>%
  arrange(bin_balanced_MAE) %>%
  pull(cfg_lbl)
ms2 <- ms2 %>% mutate(cfg_lbl = factor(cfg_lbl, levels = cfg_levels))

p_mae <- ggplot(ms2, aes(x = cfg_lbl, y = bin_balanced_MAE, fill = alpha_lbl)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_viridis_d(name = "bin-weight exponent  \u03b1", end = 0.9) +
  coord_cartesian(ylim = c(0.16, 0.20)) +
  labs(x = NULL, y = "Bin-balanced MAE") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_text(size = 13, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black",
                                   lineheight = 0.9))

ggsave(file.path(out_dir, "bin_balanced_MAE_pretty.png"),
       p_mae, width = 14, height = 5, dpi = 150)
cat("Wrote bin_balanced_MAE_pretty.png\n")
cat("DONE\n")
