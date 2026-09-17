# ==============================================================================
# SCRIPT: Compare reweighted XGBoost models vs baseline (M2 -- summary figures)
# PURPOSE: Consumes the 36-model realistic-grid sweep (03_Select_XGBoost.R in
#          XGB_MODEL_SUBDIR=realistic mode) alongside the baseline model, and
#          produces the paper-facing comparison artefacts:
#            1. Unified bin-wise calibration table (37 models x bins x metrics)
#            2. Bar chart of bin_balanced_MAE grouped by alpha
#            3. Side-by-side calibration curves for baseline + 5 SHAP models
#            4. SHAP importance heatmap across the 5 SHAP models
#
# INPUT:  Models/xgboost/calibration_residuals.csv                  (baseline)
#         Models/xgboost/realistic/search/config_*/calibration_residuals.csv
#         Models/xgboost/realistic/model_selection.csv               (36-row rank)
#         Models/xgboost/realistic/shap_selection.csv                (5 chosen)
#         Models/xgboost/realistic/shap_models/<tag>/shap_summary.csv
# OUTPUT: Models/xgboost/realistic/comparison/
#           unified_calibration_table.csv
#           bin_balanced_MAE_by_alpha.png
#           calibration_curves_grid.png
#           shap_importance_across_models.png
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
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
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 8: Compare reweighted models vs baseline")

XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "realistic")
baseline_base <- file.path(OUTPUT_BASE, "Models", "xgboost")
realistic_base <- file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
out_dir <- file.path(realistic_base, "comparison")
safe_mkdir(out_dir)

theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(plot.title       = element_text(face = "bold", size = rel(1.05)),
          plot.subtitle    = element_text(colour = "grey30", size = rel(0.82),
                                          margin = margin(b = 6)),
          axis.title       = element_text(face = "bold", size = rel(0.9)),
          axis.text        = element_text(colour = "grey20"),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          legend.title     = element_text(face = "bold", size = rel(0.8)),
          legend.text      = element_text(size = rel(0.75)))
}

# ==============================================================================
# 2. UNIFIED BIN-WISE CALIBRATION TABLE
# ==============================================================================

# Baseline: one calibration_residuals.csv at model_base root.
baseline_cr_path <- file.path(baseline_base, "calibration_residuals.csv")
if (!file.exists(baseline_cr_path)) {
  log_warning("Baseline calibration_residuals.csv missing: %s", baseline_cr_path)
  baseline_cr <- NULL
} else {
  baseline_cr <- read.csv(baseline_cr_path, stringsAsFactors = FALSE)
  baseline_cr$model <- "baseline"
  baseline_cr$alpha <- NA_real_
  baseline_cr$config_index <- NA_integer_
}

# Realistic: one calibration_residuals.csv per config_dir.
cfg_dirs <- list.files(file.path(realistic_base, "search"),
                       pattern = "^config_[0-9]+$", full.names = TRUE)
sel <- read.csv(file.path(realistic_base, "model_selection.csv"),
                stringsAsFactors = FALSE)
sel_lookup <- sel[, c("config_index", "alpha")]

realistic_cr_list <- lapply(cfg_dirs, function(d) {
  f <- file.path(d, "calibration_residuals.csv")
  if (!file.exists(f)) return(NULL)
  cr <- read.csv(f, stringsAsFactors = FALSE)
  ci <- as.integer(sub(".*config_", "", d))
  cr$config_index <- ci
  cr$alpha <- sel_lookup$alpha[match(ci, sel_lookup$config_index)]
  cr$model <- sprintf("realistic_c%d_a%g", ci, cr$alpha)
  cr
})
realistic_cr <- bind_rows(realistic_cr_list)

unified <- bind_rows(baseline_cr, realistic_cr)
write.csv(unified, file.path(out_dir, "unified_calibration_table.csv"), row.names = FALSE)
log_info("Unified calibration table: %d model x bin rows -> %s",
         nrow(unified), out_dir)

# ==============================================================================
# 3. GROUPED BAR CHART: bin_balanced_MAE, 12 HP configs x 3 alpha
# ==============================================================================
# bin_balanced_MAE was computed in 03_Select_XGBoost.R (realistic mode) and
# lives in model_selection.csv. One GROUP per hyperparameter config; within each
# group three dodged bars for the bin-weight exponent alpha. This collapses the
# old 3-panel facet into a single compact panel (one x-axis of 12 configs, one
# alpha legend). subsample and colsample_bytree are constant (0.8) across the
# whole sweep, so they are omitted from the labels.

alpha_levels <- sort(unique(sel$alpha))
alpha_pal    <- setNames(viridis(length(alpha_levels), option = "D", end = 0.88),
                         as.character(alpha_levels))
# Human-readable legend entries; fall back to the bare value for unknown alphas.
alpha_desc_map <- c("0"   = "0  (no reweighting)",
                    "0.5" = "0.5  (inverse \u221an)",
                    "1"   = "1  (inverse n)")
alpha_labels <- ifelse(as.character(alpha_levels) %in% names(alpha_desc_map),
                       alpha_desc_map[as.character(alpha_levels)],
                       sprintf("%g", alpha_levels))

bar_df <- sel %>%
  mutate(
    hp_key    = sprintf("depth %d\nlr %g \u00b7 mcw %g",
                        max_depth, eta, min_child_weight),
    alpha_lbl = factor(as.character(alpha), levels = as.character(alpha_levels))
  )

# Order the 12 groups left -> right by their best (usually alpha = 0) MAE, so the
# panel reads best-to-worst config.
grp_order <- bar_df %>%
  group_by(hp_key) %>%
  summarise(best_mae = min(bin_balanced_MAE), .groups = "drop") %>%
  arrange(best_mae) %>% pull(hp_key)
bar_df$hp_key <- factor(bar_df$hp_key, levels = grp_order)

# Zoom the y-axis to the data band (values sit in a narrow ~0.16-0.20 window),
# so the alpha effect within each group is legible.
y_lo <- floor(min(bar_df$bin_balanced_MAE) * 100) / 100
y_hi <- ceiling(max(bar_df$bin_balanced_MAE) * 100) / 100

p_bar <- ggplot(bar_df,
                aes(x = hp_key, y = bin_balanced_MAE, fill = alpha_lbl)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74,
           colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = alpha_pal, labels = alpha_labels,
                    name = expression(paste("bin-weight exponent  ", alpha))) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  labs(x = NULL, y = "Bin-balanced MAE",
       title = "Calibration error across 12 hyperparameter configs \u00d7 3 reweighting strengths") +
  theme_paper(11) +
  theme(legend.position = "top",
        legend.key.size  = unit(0.9, "lines"),
        axis.text.x      = element_text(size = rel(0.78), lineheight = 0.95),
        panel.grid.major.x = element_blank())
ggsave(file.path(out_dir, "bin_balanced_MAE_by_alpha.png"),
       p_bar, width = 13, height = 5.2, dpi = PLOT_DPI)
log_info("Wrote bin_balanced_MAE_by_alpha.png (grouped: 12 configs x %d alpha)",
         length(alpha_levels))

# ==============================================================================
# 4. CALIBRATION CURVES: baseline + 5 SHAP models (with observed 10-90 bands)
# ==============================================================================
# Each panel plots mean observed detection rate against mean predicted rate,
# per predicted-power bin, with y = x reference. Well-calibrated => hugs the
# diagonal.
#
# The point (mean_obs per bin) can look well-calibrated on average while hiding
# large within-bin spread of observed detection rate -- the mean is reassuring
# but individual outcomes in the mid range are volatile. To expose this, each
# point carries a vertical error bar spanning the trial-weighted 10th-90th
# PERCENTILE of the observed detection rate within that predicted bin. The band
# is wide in the middle (mid-power outcomes are intrinsically variable) and
# narrow at the extremes (near 0 and 1 the outcome is nearly deterministic),
# which is the shape the mean alone conceals.
#
# The band is computed from the per-row OOF predictions, not the binned summary:
#   * 5 realistic-grid SHAP configs -> search/config_<k>/oof_predictions.qs
#   * baseline (full grid)          -> Models/xgboost/oof_predictions.qs
#     (regenerated by 08b_Baseline_OOF.R; the full-search baseline persisted only
#      its binned summary, so 08b reproduces the per-row OOF identically).
# Weighting uses n_trials (w_trials), matching the mean_obs calibration weight
# in 02_Fit_XGBoost.R so the point and its band are on the same footing.

shap_sel <- read.csv(file.path(realistic_base, "shap_selection.csv"),
                     stringsAsFactors = FALSE)

#' Trial-weighted quantiles of a vector (continuity-corrected weighted ECDF).
#' Returns NA for empty input; clamps (rule = 2) for probs outside the support.
wtd_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w)
  p  <- (cw - 0.5 * w) / sum(w)
  stats::approx(p, x, xout = probs, rule = 2, ties = "ordered")$y
}

#' Per predicted-power bin: trial-weighted mean_pred, mean_obs, and 10-90 band
#' of observed detection rate. Bins match 02_Fit_XGBoost.R (seq(0, 1, 0.1)).
bin_oof <- function(oof) {
  oof %>%
    mutate(pred_bin = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
    group_by(pred_bin) %>%
    summarise(
      n         = dplyr::n(),
      mean_pred = weighted.mean(pred, w_trials),
      mean_obs  = weighted.mean(obs,  w_trials),
      obs_p10   = wtd_quantile(obs, w_trials, 0.10),
      obs_p90   = wtd_quantile(obs, w_trials, 0.90),
      .groups   = "drop"
    ) %>%
    filter(is.finite(mean_pred), is.finite(mean_obs))
}

# Model registry: key -> (pretty label, OOF path). Baseline first, then the 5
# chosen configs in shap_selection order.
model_keys <- c("baseline",
                sprintf("realistic_c%d_a%g", shap_sel$config_index, shap_sel$alpha))
model_pretty <- c(baseline = "Baseline (full grid)")
oof_paths <- c(baseline = file.path(baseline_base, "oof_predictions.qs"))
for (i in seq_len(nrow(shap_sel))) {
  key <- sprintf("realistic_c%d_a%g", shap_sel$config_index[i], shap_sel$alpha[i])
  model_pretty[key] <- sprintf("%s (config %d, alpha=%g)",
                               shap_sel$tag[i], shap_sel$config_index[i],
                               shap_sel$alpha[i])
  oof_paths[key] <- file.path(realistic_base, "search",
                              sprintf("config_%d", shap_sel$config_index[i]),
                              "oof_predictions.qs")
}

cal_list <- lapply(model_keys, function(key) {
  f <- oof_paths[key]
  if (!file.exists(f)) {
    log_warning("OOF file missing for %s (%s) -- panel will show points only, no band.",
                key, f)
    return(NULL)
  }
  oof <- qs::qread(f)
  cc <- bin_oof(oof)
  cc$model <- key
  cc
})
cal_curves <- bind_rows(cal_list)
if (nrow(cal_curves) == 0) stop("No OOF-derived calibration rows; check oof_predictions.qs paths.")

present_keys <- model_keys[model_keys %in% unique(cal_curves$model)]
cal_curves$model <- factor(cal_curves$model, levels = present_keys,
                           labels = model_pretty[present_keys])

p_cal <- ggplot(cal_curves, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = obs_p10, ymax = obs_p90),
                width = 0.02, colour = "steelblue", alpha = 0.55, linewidth = 0.4) +
  geom_line(colour = "steelblue", alpha = 0.6) +
  geom_point(size = 1.8, colour = "steelblue") +
  facet_wrap(~ model, ncol = 3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted detection rate (per bin)",
       y = "Observed detection rate (per bin: mean point, 10-90th percentile bar)",
       title = "OOF calibration: baseline vs 5 realistic-grid SHAP models",
       subtitle = "Error bars = trial-weighted 10-90th percentile of observed detection rate within each predicted bin") +
  theme_paper(11)
ggsave(file.path(out_dir, "calibration_curves_grid.png"),
       p_cal, width = 11, height = 7, dpi = PLOT_DPI)
write.csv(cal_curves, file.path(out_dir, "calibration_curves_bands.csv"), row.names = FALSE)
log_info("Wrote calibration_curves_grid.png (%d panels, with 10-90 observed bands)",
         length(present_keys))

# ==============================================================================
# 5. SHAP IMPORTANCE HEATMAP ACROSS 5 SHAP MODELS
# ==============================================================================
# Reads Models/xgboost/realistic/shap_models/<tag>/shap_summary.csv (written
# by 03_Select_XGBoost.R). Rank-normalise within each model so scales are
# comparable across differently-weighted training runs (design doc §10.8).

shap_imp_rows <- lapply(seq_len(nrow(shap_sel)), function(i) {
  tag <- shap_sel$tag[i]
  f <- file.path(realistic_base, "shap_models", tag, "shap_summary.csv")
  if (!file.exists(f)) { log_warning("Missing shap_summary for %s", tag); return(NULL) }
  ss <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(tag           = tag,
             feature       = ss$feature,
             mean_abs_shap = ss$mean_abs_shap,
             rank          = rank(-ss$mean_abs_shap, ties.method = "average"))
})
shap_imp <- bind_rows(shap_imp_rows)

if (nrow(shap_imp) > 0) {
  # Order features by average rank across models (best rank = 1 at top).
  feat_avg <- shap_imp %>%
    group_by(feature) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    arrange(desc(avg_rank))   # worst average rank at bottom => best at top after coord_flip
  shap_imp$feature <- factor(shap_imp$feature, levels = feat_avg$feature)
  shap_imp$tag <- factor(shap_imp$tag, levels = shap_sel$tag)

  p_heat <- ggplot(shap_imp, aes(x = tag, y = feature, fill = mean_abs_shap)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", mean_abs_shap)), size = 3, colour = "white") +
    scale_fill_viridis_c(name = "mean|SHAP|") +
    labs(title = "SHAP importance across 5 realistic-grid SHAP models",
         subtitle = "Colour = mean|SHAP| (log-odds units); features ordered by avg rank",
         x = NULL, y = NULL) +
    theme_paper(12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
  ggsave(file.path(out_dir, "shap_importance_across_models.png"),
         p_heat, width = 9, height = 6, dpi = PLOT_DPI)
  write.csv(shap_imp, file.path(out_dir, "shap_importance_across_models.csv"),
            row.names = FALSE)
  log_info("Wrote shap_importance_across_models.png")
}

log_section("STEP 8 COMPLETE")
log_info("Comparison artefacts written to %s", out_dir)
