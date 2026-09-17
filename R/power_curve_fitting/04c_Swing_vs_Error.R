# ==============================================================================
# SCRIPT: Cell-type swing vs model error floor (M2 -- follows 04b ICE)
# PURPOSE: 04b showed the cell-type ICE swing is large in the mid-power band
#          (median 0.39). The open question: is that swing a REAL cell-type
#          effect, or is it just the model's own resolution -- i.e. in the steep
#          mid-band EVERY input moves the prediction a lot, and cell type is no
#          more than the model's generic wobble there?
#
#          This script answers it by overlaying two per-real-grid-point
#          distributions and comparing them PAIRED, split by power region:
#            (E) OOF model error   e_i = |p_oof_i - p_obs_i|
#                (honest generalisation error: the winning realistic config's
#                 out-of-fold prediction vs the observed detection rate; the OOF
#                 preds were saved per-config by 02_Fit_XGBoost.R in ARRAY mode)
#            (S) cell-type swing   s_i = max_k p_ik - min_k p_ik  over
#                {Astro, Micro, Oligo}, computed on the SHIPPED exported booster
#                with everything except cell_type frozen at the row's native
#                design (identical recipe to 04b).
#
#          The decision-relevant statement is NOT "are the two marginals
#          similar" (they would look similar merely because both are large in
#          the mid-band and tiny on the 0/1 rails -- the confound). It is the
#          PAIRED comparison within each power region:
#            - if s_i ~ e_i at the same points (ratio ~1, large fraction with
#              s_i <= e_i), the cell-type swing sits inside the model's own
#              error floor: choosing a cell type is arguing below the model's
#              resolution.
#            - if s_i >> e_i, cell type moves the prediction by MORE than the
#              model's typical error, so it is a real, resolvable effect.
#
#          Region axis = the point's OOF predicted power (where the model thinks
#          it is), so error and swing are binned on the same anchor.
#
#          Caveat recorded in the summary: the OOF error e_i contains BOTH model
#          misspecification AND irreducible binomial label noise (obs is a
#          detection rate over n_trials). Comparing the swing against the FULL
#          OOF error (Chris's choice: model-error floor only) is therefore
#          slightly generous to the "within the floor" side, since label noise
#          inflates e_i. Stated so the reader can calibrate.
#
#          Runs off the exported bundle (swing) + the winner's saved OOF preds
#          (error). No HPC, no refit. base R + xgboost + ggplot2.
#
# INPUT:  Exports/xgboost/<subdir>/power_model_bundle.rds        (shipped model)
#         Models/xgboost/<subdir>/cv_metrics.csv                 (winner config_index)
#         Models/xgboost/<subdir>/search/config_<k>/oof_predictions.qs  (OOF error)
# OUTPUT: Models/xgboost/<subdir>/celltype_ice/swing_vs_error/
#           swing_vs_error_perpoint.csv     per-row e_i, s_i, ratio, region
#           swing_vs_error_summary.csv      paired stats overall + by region
#           swing_vs_error_ecdf.png         ECDF overlay (swing vs OOF error), by region
#           swing_vs_error_density.png      density overlay, by region
#           swing_vs_error_paired.png       paired hexbin s_i vs e_i with y=x
#           frac_swing_within_error.png     fraction s_i <= e_i by region (bar)
#           swing_vs_error_summary.md       paper-quotable interpretation
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(xgboost)
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

XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "realistic")
PLOT_SAMPLE      <- as.integer(Sys.getenv("SVE_PLOT_SAMPLE", "200000"))  # for hexbin only

set.seed(SEED_BASE)

log_section(sprintf("STEP 4c (XGBoost): Cell-type swing vs model-error floor [subdir=%s]",
                    if (nzchar(XGB_MODEL_SUBDIR)) XGB_MODEL_SUBDIR else "baseline"))

model_base <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Models", "xgboost")
}
export_dir <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Exports", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Exports", "xgboost")
}
out_dir <- file.path(model_base, "celltype_ice", "swing_vs_error")
safe_mkdir(out_dir)

# ==============================================================================
# 2. LOAD SHIPPED PREDICTOR (for the swing) -- same bundle 04b/04 use
# ==============================================================================

bundle_path <- file.path(export_dir, "power_model_bundle.rds")
if (!file.exists(bundle_path)) {
  stop(sprintf("Bundle not found: %s. Run 04_Export_XGBoost.R first.", bundle_path))
}
bundle  <- readRDS(bundle_path)
booster <- xgb.load.raw(bundle$booster_raw)

CONT_FEATURES <- bundle$cont_features
CAT_FEATURES  <- bundle$cat_features
cat_levels    <- bundle$cat_levels
feature_names <- bundle$feature_names
CELL_LEVELS   <- cat_levels$cell_type
METHOD_LEVELS <- cat_levels$de_method
DEFAULT_CT    <- bundle$defaults$cell_type

build_design <- function(df) {
  cont_mat <- as.matrix(df[, CONT_FEATURES, drop = FALSE])
  oh_blocks <- lapply(CAT_FEATURES, function(f) {
    lv <- cat_levels[[f]]
    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
    colnames(m) <- paste0(f, "=", lv)
    m
  })
  X <- do.call(cbind, c(list(cont_mat), oh_blocks))
  X[, feature_names, drop = FALSE]
}
predict_rows <- function(df) as.numeric(predict(booster, xgb.DMatrix(build_design(df))))

log_info("Loaded shipped bundle: %d rounds, cell types {%s}",
         bundle$metadata$best_nrounds, paste(CELL_LEVELS, collapse = ", "))

# ==============================================================================
# 3. LOAD THE WINNER'S OOF PREDICTIONS (for the error floor)
# ==============================================================================
# 03_Select_XGBoost.R promotes the winning config but does NOT copy its OOF
# preds to model_base; they stay in search/config_<k>/oof_predictions.qs. The
# winner's config_index is written to model_base/cv_metrics.csv. Recover k,
# then load that config's OOF table (obs, pred, native design + features).

cv_path <- file.path(model_base, "cv_metrics.csv")
if (!file.exists(cv_path)) {
  stop(sprintf("Winner cv_metrics.csv not found: %s. Run 03_Select_XGBoost.R first.", cv_path))
}
winner_cfg <- read.csv(cv_path, stringsAsFactors = FALSE)$config_index[1]
oof_path <- file.path(model_base, "search", sprintf("config_%d", winner_cfg),
                      "oof_predictions.qs")
if (!file.exists(oof_path)) {
  stop(sprintf("Winner OOF preds not found: %s.", oof_path))
}
oof <- qs::qread(oof_path)
log_info("Loaded winner (config %d) OOF preds: %d rows", winner_cfg, nrow(oof))

# Keep only finite, in-range rows with the design columns we need.
need_cols <- c("obs", "pred", "w_trials", "cell_type", "de_method",
               "n_donors", "n_cells", "reads_per_cell",
               "abs_gt_lfc", "gt_basemean", "gt_dispersion")
missing <- setdiff(need_cols, colnames(oof))
if (length(missing) > 0) {
  stop(sprintf("OOF table missing columns: %s", paste(missing, collapse = ", ")))
}
oof <- oof %>%
  filter(is.finite(obs), is.finite(pred), is.finite(w_trials), w_trials > 0,
         is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
         is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion),
         cell_type %in% CELL_LEVELS, de_method %in% METHOD_LEVELS)
n_rows <- nrow(oof)
log_info("Usable OOF grid points: %d", n_rows)

# ==============================================================================
# 4. THE TWO PAIRED QUANTITIES
# ==============================================================================
# (E) OOF model error at the point's NATIVE design.
oof_err <- abs(oof$pred - oof$obs)

# (S) cell-type swing on the shipped booster, native design frozen except cell_type.
base_df <- data.frame(
  n_donors       = oof$n_donors,
  n_cells        = oof$n_cells,
  reads_per_cell = oof$reads_per_cell,
  abs_gt_lfc     = oof$abs_gt_lfc,
  gt_basemean    = oof$gt_basemean,
  gt_dispersion  = oof$gt_dispersion,
  cell_type      = as.character(oof$cell_type),
  de_method      = as.character(oof$de_method),
  stringsAsFactors = FALSE
)
ct_preds <- matrix(NA_real_, nrow = n_rows, ncol = length(CELL_LEVELS),
                   dimnames = list(NULL, CELL_LEVELS))
for (ct in CELL_LEVELS) {
  df <- base_df
  df$cell_type <- ct
  ct_preds[, ct] <- predict_rows(df)
}
ct_swing <- apply(ct_preds, 1, max) - apply(ct_preds, 1, min)

# ==============================================================================
# 5. REGION AXIS + PAIRED RATIO
# ==============================================================================
# Bin on the OOF predicted power (where the model thinks the point sits): both
# error and swing are anchored to the same real row, so this is a single shared
# axis. Near 0/1 the surface saturates so both quantities vanish; the contest is
# in the mid band.
region_of <- function(p) {
  cut(p, breaks = c(-Inf, 0.1, 0.9, Inf),
      labels = c("low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]"))
}
region <- region_of(oof$pred)

eps_ratio  <- 1e-3   # floor to keep the ratio finite where error is ~0
swing_ratio <- ct_swing / pmax(oof_err, eps_ratio)
swing_within_error <- ct_swing <= oof_err   # swing sits inside the error floor

# ==============================================================================
# 6. PER-POINT TABLE
# ==============================================================================
perpoint <- data.frame(
  base_df,
  oof_pred            = oof$pred,
  oof_obs             = oof$obs,
  oof_error           = oof_err,
  celltype_swing      = ct_swing,
  swing_over_error    = swing_ratio,
  swing_within_error  = swing_within_error,
  region              = as.character(region),
  stringsAsFactors = FALSE
)
write.csv(perpoint, file.path(out_dir, "swing_vs_error_perpoint.csv"),
          row.names = FALSE)
log_info("Wrote per-point table -> swing_vs_error_perpoint.csv (%d rows)", n_rows)

# ==============================================================================
# 7. PAIRED SUMMARY (overall + by region)
# ==============================================================================
summarise_group <- function(idx, grp) {
  e <- oof_err[idx]; s <- ct_swing[idx]; r <- swing_ratio[idx]
  data.frame(
    region = grp, n = length(idx),
    err_median = median(e), err_p95 = as.numeric(quantile(e, 0.95)),
    swing_median = median(s), swing_p95 = as.numeric(quantile(s, 0.95)),
    ratio_median = median(r), ratio_p95 = as.numeric(quantile(r, 0.95)),
    frac_swing_le_error = mean(s <= e),
    frac_swing_2x_error = mean(s > 2 * e),
    stringsAsFactors = FALSE
  )
}
groups <- c(list(overall = seq_len(n_rows)),
            split(seq_len(n_rows), region))
summary_tbl <- bind_rows(lapply(names(groups), function(g) summarise_group(groups[[g]], g)))
# Order: overall, then regions in natural power order.
region_order <- c("overall", "low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]")
summary_tbl <- summary_tbl[match(region_order, summary_tbl$region), ]
summary_tbl <- summary_tbl[!is.na(summary_tbl$region), ]
write.csv(summary_tbl, file.path(out_dir, "swing_vs_error_summary.csv"), row.names = FALSE)
log_info("Wrote paired summary -> swing_vs_error_summary.csv")

mid_row <- summary_tbl[summary_tbl$region == "mid (0.1,0.9)", ]
all_row <- summary_tbl[summary_tbl$region == "overall", ]

# ==============================================================================
# 8. FIGURES
# ==============================================================================

theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(plot.title       = element_text(face = "bold", size = rel(1.05)),
          plot.subtitle    = element_text(colour = "grey30", size = rel(0.82),
                                          margin = margin(b = 6)),
          axis.title       = element_text(face = "bold", size = rel(0.9)),
          axis.text        = element_text(colour = "grey20"),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          strip.background = element_rect(fill = "grey92", colour = "grey40"),
          strip.text       = element_text(face = "bold", size = rel(0.8)),
          legend.title     = element_text(face = "bold", size = rel(0.8)),
          legend.text      = element_text(size = rel(0.75)))
}
QUANT_COLS <- c("Cell-type swing" = "#21908C", "OOF model error" = "#B8360E")

long_df <- bind_rows(
  data.frame(region = as.character(region), value = ct_swing,
             quantity = "Cell-type swing", stringsAsFactors = FALSE),
  data.frame(region = as.character(region), value = oof_err,
             quantity = "OOF model error", stringsAsFactors = FALSE)
)
long_df$region <- factor(long_df$region,
                         levels = c("low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]"))

# --- 8.1 ECDF overlay, faceted by region (headline) --------------------------
p_ecdf <- ggplot(long_df, aes(x = value, colour = quantity)) +
  stat_ecdf(linewidth = 0.7) +
  facet_wrap(~ region, nrow = 1) +
  scale_colour_manual(values = QUANT_COLS, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Value (power units)", y = "Cumulative fraction of grid points",
       title = "Is the cell-type swing bigger than the model's own error?",
       subtitle = sprintf("Per real grid point: cell-type swing vs OOF |pred - obs|. Mid band: median swing %.3f vs median error %.3f",
                          mid_row$swing_median, mid_row$err_median)) +
  theme_paper(12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "swing_vs_error_ecdf.png"), p_ecdf,
       width = 10, height = 4.2, dpi = PLOT_DPI)

# --- 8.2 density overlay, faceted by region ----------------------------------
p_den <- ggplot(long_df, aes(x = value, fill = quantity, colour = quantity)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  facet_wrap(~ region, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = QUANT_COLS, name = NULL) +
  scale_colour_manual(values = QUANT_COLS, name = NULL) +
  labs(x = "Value (power units)", y = "Density",
       title = "Distribution overlay: cell-type swing vs OOF model error",
       subtitle = "If the swing distribution sits inside the error distribution, cell type is below model resolution") +
  theme_paper(12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "swing_vs_error_density.png"), p_den,
       width = 10, height = 4.2, dpi = PLOT_DPI)

# --- 8.3 paired hexbin s_i vs e_i with y=x -----------------------------------
plot_idx <- if (n_rows > PLOT_SAMPLE) sort(sample.int(n_rows, PLOT_SAMPLE)) else seq_len(n_rows)
pair_df <- data.frame(error = oof_err[plot_idx], swing = ct_swing[plot_idx],
                      region = factor(as.character(region[plot_idx]),
                                      levels = c("low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]")))
p_pair <- ggplot(pair_df, aes(x = error, y = swing)) +
  geom_bin2d(bins = 60) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ region, nrow = 1) +
  scale_fill_viridis_c(name = "Grid\npoints", trans = "log10") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "OOF model error  |pred - obs|", y = "Cell-type swing  s",
       title = "Paired: cell-type swing vs OOF error at the same grid point",
       subtitle = "Above the dashed y=x line: swing exceeds the model's error (cell type is resolvable)") +
  theme_paper(12)
ggsave(file.path(out_dir, "swing_vs_error_paired.png"), p_pair,
       width = 10, height = 4.4, dpi = PLOT_DPI)

# --- 8.4 fraction swing <= error by region -----------------------------------
frac_df <- summary_tbl[summary_tbl$region != "overall", ]
frac_df$region <- factor(frac_df$region,
                         levels = c("low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]"))
p_frac <- ggplot(frac_df, aes(x = region, y = frac_swing_le_error, fill = region)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * frac_swing_le_error)),
            vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_fill_viridis_d(guide = "none", end = 0.85) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.05), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Fraction of grid points with swing <= OOF error",
       title = "Where does cell type sit inside the model's error floor?",
       subtitle = "High bar = the cell-type effect is within the noise floor there; low bar = it is resolvable") +
  theme_paper(12)
ggsave(file.path(out_dir, "frac_swing_within_error.png"), p_frac,
       width = 7, height = 5, dpi = PLOT_DPI)

log_info("Wrote 4 figures -> %s", out_dir)

# ==============================================================================
# 9. PAPER-QUOTABLE SUMMARY
# ==============================================================================
fmt_row <- function(r) sprintf("| %s | %d | %.4f | %.4f | %.4f | %.4f | %.3f | %.1f%% |",
                               r$region, r$n, r$err_median, r$err_p95,
                               r$swing_median, r$swing_p95, r$ratio_median,
                               100 * r$frac_swing_le_error)
md <- c(
  "# Cell-type swing vs OOF model-error floor",
  "",
  sprintf("Model: XGBoost `%s` (shipped bundle) for the swing; winner config %d OOF predictions for the error floor.",
          XGB_MODEL_SUBDIR, winner_cfg),
  sprintf("Real grid points: %d.", n_rows),
  "",
  "## Question",
  "04b showed the cell-type ICE swing is large in the mid-power band (median ~0.39).",
  "Is that a real cell-type effect, or just the model's own resolution -- in the",
  "steep mid-band every input moves the prediction, so is cell type any more than",
  "the model's generic wobble? Compare the swing against the model's OOF error at",
  "the SAME points, paired and split by power region.",
  "",
  "## Paired comparison",
  "",
  "| region | n | err median | err p95 | swing median | swing p95 | swing/err median | % swing <= err |",
  "|---|---|---|---|---|---|---|---|",
  fmt_row(all_row),
  paste(vapply(which(summary_tbl$region != "overall"),
               function(i) fmt_row(summary_tbl[i, ]), character(1)), collapse = "\n"),
  "",
  "## Reading",
  sprintf("- Mid band (the region you plan in): median swing %.3f vs median OOF error %.3f, ratio %.2f. %.1f%% of mid-band points have swing <= OOF error.",
          mid_row$swing_median, mid_row$err_median, mid_row$ratio_median,
          100 * mid_row$frac_swing_le_error),
  "- If swing/err >> 1 and few points have swing <= err, the cell-type effect is",
  "  RESOLVABLE: it moves the prediction by more than the model's typical error,",
  "  so it is not just mid-band steepness. Cell type should be kept, not dropped.",
  "- If swing/err ~ 1 with a large fraction swing <= err, the swing sits inside",
  "  the model's own error floor: committing to a specific cell type argues below",
  "  the model's resolution, and defaulting it (to Oligo, the conservative",
  "  reference) is the honest move.",
  "",
  "## Caveat (floor definition)",
  "The OOF error e_i mixes model misspecification with irreducible binomial label",
  "noise (obs is a detection rate over finite n_trials, max SE ~0.11 at 20 trials).",
  "Per Chris's choice the floor here is the FULL OOF error (model error only, no",
  "separate label-noise term), which is slightly generous to the 'within the",
  "floor' side because label noise inflates e_i. A model-independent floor (the",
  "detection-rate binomial SE) would be a stricter test if wanted."
)
writeLines(md, file.path(out_dir, "swing_vs_error_summary.md"))
log_info("Wrote paper-quotable summary -> swing_vs_error_summary.md")

log_section("STEP 4c (XGBoost) SWING-VS-ERROR COMPLETE")
log_info("Mid band: swing median=%.4f err median=%.4f ratio=%.2f | %% swing<=err mid=%.1f overall=%.1f",
         mid_row$swing_median, mid_row$err_median, mid_row$ratio_median,
         100 * mid_row$frac_swing_le_error, 100 * all_row$frac_swing_le_error)
log_info("Outputs in %s", out_dir)
