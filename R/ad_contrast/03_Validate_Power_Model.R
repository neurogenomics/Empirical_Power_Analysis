#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  03_Validate_Power_Model.R
# PURPOSE: Does the power model trained on the sex contrast generalise to a
#          disease contrast? Apply the selected XGBoost power model (trained
#          only on sex-biased detection outcomes) to the AD vs control grid,
#          used purely as an external test set, and compare its performance
#          with the model's own out-of-fold performance on the sex grid.
#
# TEST SET: AD vs control detection outcomes, aggregated exactly like the
#   training table (R/power_curve_fitting 00, 00b, 01, 01b, 01e with
#   EPA_CONTRAST=ad_vs_control): one row per (gene, cell type, DE method,
#   run, design point) with n_cells >= 100. Gene features come from the AD
#   ground-truth half of each run (|LFC|, mean expression and DESeq2
#   dispersion fitted on the AD truth-half donors), so the dispersion reflects
#   the AD design formula, not the sex one.
#
# METRICS (same definitions as the sex model's evaluation):
#   * trial-weighted MAE and bin-balanced MAE (03_Select_XGBoost.R)
#   * calibration per predicted-power bin with the 10-90th percentile of the
#     observed detection rate (08_Compare_Reweighted.R)
#   * per-class precision / recall / F1 for every partition scheme of
#     09_Trinary_Classification.R; the skewed trinary scheme
#     ([0-0.4] [0.4-0.8] [0.8-1]) is the headline
#   * share of test rows outside the training range of each feature
#
# INPUT:   Datalink/AD_vs_Control/Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data.qs
#          Datalink/Power_Curve_Fitting/Models/xgboost/{booster.ubj,feature_spec.qs}
#          Datalink/Power_Curve_Fitting/Models/xgboost/trinary_classification/per_class_metrics.csv
#          Datalink/Power_Curve_Fitting/Models/xgboost/comparison/calibration_curves_bands.csv
#          Datalink/Power_Curve_Fitting/Exports/xgboost/power_model_bundle.rds  (training ranges)
#          Datalink/Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data.qs (feature shift plot)
# OUTPUT:  Datalink/AD_vs_Control/Power_Model_Test/
#            test_metrics.csv              overall MAE, bin-balanced MAE, accuracy, macro-F1 per scheme
#            per_class_metrics.csv         precision / recall / F1, sex (OOF) and AD (test)
#            calibration_bands.csv
#            extrapolation.csv             share of AD rows outside each feature's training range
#            trinary_skewed_prf.{png,pdf}  sex OOF vs AD test, per class
#            calibration_sex_vs_ad.{png,pdf}
#            feature_shift.{png,pdf}       gene-feature distributions, sex training vs AD test
#            ad_predictions.qs             per-row predictions
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(xgboost)
})

ROOT    <- Sys.getenv("POWER_ANALYSIS_ROOT", ".")
DL      <- file.path(ROOT, "Datalink")
MODEL   <- file.path(DL, "Power_Curve_Fitting", "Models", "xgboost")
OUT_DIR <- file.path(DL, "AD_vs_Control", "Power_Model_Test")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
NTHREAD <- as.integer(Sys.getenv("XGB_NTHREAD", "8"))

# ---- model + test data --------------------------------------------------------------
spec    <- qs::qread(file.path(MODEL, "feature_spec.qs"))
booster <- xgboost::xgb.load(file.path(MODEL, "booster.ubj"))
xgboost::xgb.attr(booster, "nthread") <- as.character(NTHREAD)

clean <- function(d) {
  if ("gene_set" %in% names(d)) d <- d[d$gene_set == "full", ]
  d %>% filter(is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
               n_trials > 0, is.finite(n_donors), is.finite(n_cells),
               is.finite(reads_per_cell), is.finite(abs_gt_lfc),
               is.finite(gt_basemean), is.finite(gt_dispersion))
}
ad <- clean(qs::qread(file.path(DL, "AD_vs_Control", "Power_Curve_Fitting",
                                "Fitting_Data", "aggregated_fitting_data.qs")))
cat(sprintf("AD test rows: %d (%d genes)\n", nrow(ad), n_distinct(ad$gene_id)))

build_design <- function(df) {
  cont <- as.matrix(df[, spec$cont_features, drop = FALSE])
  oh <- lapply(spec$cat_features, function(f) {
    lv <- spec$cat_levels[[f]]
    unknown <- setdiff(unique(as.character(df[[f]])), lv)
    if (length(unknown)) stop("Unknown ", f, " levels in test data: ", paste(unknown, collapse = ", "))
    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
    colnames(m) <- paste0(f, "=", lv)
    m
  })
  do.call(cbind, c(list(cont), oh))[, spec$feature_names, drop = FALSE]
}
ad$pred <- predict(booster, xgboost::xgb.DMatrix(build_design(ad)))
qs::qsave(ad[, c("gene_id", "cell_type", "de_method", "run_id", "n_donors", "n_cells",
                 "reads_per_cell", "abs_gt_lfc", "gt_basemean", "gt_dispersion",
                 "n_trials", "detection_rate", "pred")],
          file.path(OUT_DIR, "ad_predictions.qs"))

# ---- continuous metrics ------------------------------------------------------------
bin_balanced_mae <- function(pred, obs, w) {
  bin <- cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  per_bin <- tapply(seq_along(pred), bin, function(ix) sum(w[ix] * abs(obs[ix] - pred[ix])) / sum(w[ix]))
  mean(per_bin, na.rm = TRUE)
}
mae_overall <- with(ad, sum(n_trials * abs(detection_rate - pred)) / sum(n_trials))
mae_bb      <- with(ad, bin_balanced_mae(pred, detection_rate, n_trials))

# ---- classification (schemes of 09_Trinary_Classification.R) -------------------------
SCHEMES <- list(
  binary        = list(breaks = 0.5,               labels = c("Not Powered", "Powered")),
  trinary_02_08 = list(breaks = c(0.2, 0.8),       labels = c("Not Powered", "Contested", "Powered")),
  trinary_04_08 = list(breaks = c(0.4, 0.8),       labels = c("Definitely Not", "Maybe", "Yes")),
  trinary_01_09 = list(breaks = c(0.1, 0.9),       labels = c("Not Powered", "Contested", "Powered")),
  four_class    = list(breaks = c(0.2, 0.5, 0.8),  labels = c("Underpowered", "Likely Underpowered",
                                                              "Likely Powered", "Well Powered")),
  quartile      = list(breaks = c(0.25, 0.5, 0.75), labels = c("Q1: Lowest", "Q2: Low-Mid",
                                                               "Q3: High-Mid", "Q4: Highest"))
)
classify <- function(x, breaks) findInterval(x, breaks, left.open = TRUE) + 1L

per_class <- list(); overall <- list()
for (nm in names(SCHEMES)) {
  s <- SCHEMES[[nm]]; K <- length(s$labels)
  truth <- classify(ad$detection_rate, s$breaks); pred <- classify(ad$pred, s$breaks)
  cm <- table(factor(truth, 1:K), factor(pred, 1:K))
  tp <- diag(cm); prec <- tp / colSums(cm); rec <- tp / rowSums(cm)
  f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  per_class[[nm]] <- data.frame(model = "winner", role = "ad_test", scheme = nm,
                                class = s$labels, support = as.integer(rowSums(cm)),
                                precision = as.numeric(prec), recall = as.numeric(rec),
                                f1 = as.numeric(f1))
  overall[[nm]] <- data.frame(scheme = nm, accuracy = sum(tp) / sum(cm),
                              macro_f1 = mean(f1, na.rm = TRUE))
  write.csv(as.data.frame(cm, responseName = "n") %>%
              rename(truth_class = Var1, pred_class = Var2) %>%
              mutate(truth_class = s$labels[truth_class], pred_class = s$labels[pred_class]),
            file.path(OUT_DIR, sprintf("confusion_ad_%s.csv", nm)), row.names = FALSE)
}
ad_pc <- bind_rows(per_class)
sex_pc <- read.csv(file.path(MODEL, "trinary_classification", "per_class_metrics.csv"),
                   stringsAsFactors = FALSE) %>%
  filter(role == "winner") %>% mutate(role = "sex_oof")
pc_all <- bind_rows(sex_pc, ad_pc)
write.csv(pc_all, file.path(OUT_DIR, "per_class_metrics.csv"), row.names = FALSE)

metrics <- bind_rows(overall) %>%
  mutate(test_set = "AD vs control", n_rows = nrow(ad), n_genes = n_distinct(ad$gene_id),
         mae = mae_overall, bin_balanced_mae = mae_bb) %>%
  relocate(test_set, n_rows, n_genes, mae, bin_balanced_mae)
write.csv(metrics, file.path(OUT_DIR, "test_metrics.csv"), row.names = FALSE)
print(metrics, row.names = FALSE)

# ---- extrapolation -------------------------------------------------------------------
bundle <- readRDS(file.path(DL, "Power_Curve_Fitting", "Exports", "xgboost", "power_model_bundle.rds"))
extrap <- bind_rows(lapply(names(bundle$training_ranges), function(f) {
  r <- bundle$training_ranges[[f]]
  data.frame(feature = f, train_min = r[["min"]], train_max = r[["max"]],
             pct_below = 100 * mean(ad[[f]] < r[["min"]]),
             pct_above = 100 * mean(ad[[f]] > r[["max"]]))
}))
write.csv(extrap, file.path(OUT_DIR, "extrapolation.csv"), row.names = FALSE)

# ---- figure: skewed trinary precision / recall / F1, sex OOF vs AD test --------------
class_map <- c("Definitely Not" = "Unpowered\n[0-0.4]", "Maybe" = "Contested\n[0.4-0.8]",
               "Yes" = "Powered\n[0.8-1]")
role_map <- c(sex_oof = "Sex (out-of-fold)", ad_test = "AD vs control (external test)")
tri <- pc_all %>%
  filter(scheme == "trinary_04_08") %>%
  group_by(role) %>% mutate(share = support / sum(support)) %>% ungroup() %>%
  pivot_longer(c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  mutate(class  = factor(class_map[class], levels = class_map),
         role   = factor(role_map[role], levels = role_map),
         metric = factor(metric, levels = c("precision", "recall", "f1"),
                         labels = c("Precision", "Recall", "F1")))
share_lbl <- tri %>% distinct(role, class, share) %>%
  mutate(lbl = sprintf("%.1f%% of rows", 100 * share))
p_tri <- ggplot(tri, aes(x = class, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.2f", value)), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.2, na.rm = TRUE) +
  geom_text(data = share_lbl, aes(x = class, y = 1.12, label = lbl), inherit.aes = FALSE,
            size = 3.3, colour = "grey30") +
  facet_wrap(~ role, nrow = 1) +
  scale_y_continuous(limits = c(0, 1.17), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c(Precision = "#0072B2", Recall = "#E69F00", F1 = "#009E73")) +
  labs(x = NULL, y = "Score", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"), strip.text = element_text(size = 12))
ggsave(file.path(OUT_DIR, "trinary_skewed_prf.png"), p_tri, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(OUT_DIR, "trinary_skewed_prf.pdf"), p_tri, width = 11, height = 5.5)

# ---- figure: calibration, sex OOF (winner) vs AD test ----------------------------------
wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  stats::approx((cumsum(w) - 0.5 * w) / sum(w), x, xout = p, rule = 2, ties = "ordered")$y
}
ad_cal <- ad %>%
  mutate(pred_bin = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(pred_bin) %>%
  summarise(n = n(), mean_pred = weighted.mean(pred, n_trials),
            mean_obs = weighted.mean(detection_rate, n_trials),
            obs_p10 = wq(detection_rate, n_trials, 0.1), obs_p90 = wq(detection_rate, n_trials, 0.9),
            .groups = "drop") %>%
  mutate(model = role_map[["ad_test"]])
sex_bands <- read.csv(file.path(MODEL, "comparison", "calibration_curves_bands.csv"),
                      stringsAsFactors = FALSE)
# 08_Compare_Reweighted.R labels each curve "<tag> (config <k>, alpha=<a>)";
# the promoted model is row 1 of model_selection.csv.
win_cfg   <- read.csv(file.path(MODEL, "model_selection.csv"))$config_index[1]
win_label <- unique(sex_bands$model[grepl(sprintf("(config %d,", win_cfg), sex_bands$model, fixed = TRUE)])
if (length(win_label) != 1) stop("Cannot find the winner's calibration curve in calibration_curves_bands.csv")
sex_cal <- sex_bands %>% filter(model == win_label) %>% mutate(model = role_map[["sex_oof"]])
cal <- bind_rows(sex_cal, ad_cal) %>% mutate(model = factor(model, levels = role_map))
write.csv(cal, file.path(OUT_DIR, "calibration_bands.csv"), row.names = FALSE)
p_cal <- ggplot(cal, aes(mean_pred, mean_obs, colour = model)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = obs_p10, ymax = obs_p90), width = 0.02, alpha = 0.6) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted detection rate (per bin)",
       y = "Observed detection rate (mean, 10-90th percentile bar)") +
  theme_bw(base_size = 13) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(OUT_DIR, "calibration_sex_vs_ad.png"), p_cal, width = 7, height = 7.5, dpi = 300)
ggsave(file.path(OUT_DIR, "calibration_sex_vs_ad.pdf"), p_cal, width = 7, height = 7.5)

# ---- figure: gene-feature shift, sex training vs AD test --------------------------------
gene_feats <- function(d, lbl) {
  d %>% distinct(gene_id, cell_type, de_method, run_id, abs_gt_lfc, gt_basemean, gt_dispersion) %>%
    mutate(set = lbl)
}
sex_train <- clean(qs::qread(file.path(DL, "Power_Curve_Fitting", "Fitting_Data",
                                       "aggregated_fitting_data.qs")))
feats <- bind_rows(gene_feats(sex_train, "Sex (training)"), gene_feats(ad, "AD vs control (test)")) %>%
  pivot_longer(c(abs_gt_lfc, gt_basemean, gt_dispersion), names_to = "feature", values_to = "value") %>%
  filter(value > 0) %>%
  mutate(feature = factor(feature, levels = c("abs_gt_lfc", "gt_basemean", "gt_dispersion"),
                          labels = c("Truth |LFC|", "Truth mean expression", "Truth DESeq2 dispersion")),
         set = factor(set, levels = c("Sex (training)", "AD vs control (test)")))
p_shift <- ggplot(feats, aes(value, fill = set)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, position = "identity", alpha = 0.55) +
  scale_x_log10() + facet_grid(cell_type ~ feature, scales = "free") +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Density (within set)") +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(OUT_DIR, "feature_shift.png"), p_shift, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "feature_shift.pdf"), p_shift, width = 12, height = 8)

cat("Done:", OUT_DIR, "\n")
