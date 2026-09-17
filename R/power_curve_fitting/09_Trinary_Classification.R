# ==============================================================================
# SCRIPT: N-ary classification metrics for the XGBoost power model
# PURPOSE: Reframe the continuous power model as a multi-class classifier
#          by sweeping several partition schemes and comparing precision /
#          recall / F1 / accuracy across classes, models, and schemes.
#
#          Motivation: mid-band MAE for the continuous model stays at ~0.22
#          across every configuration tried, but the extremes calibrate to
#          <2% MAE; treating the output as a decision (safe / borderline /
#          hopeless) rather than a decimal probability may be far more useful
#          for design planning. This script quantifies "how much more useful".
#
# SCHEMES EVALUATED (all at once, results keyed by scheme name):
#   binary          [0, 0.5] / (0.5, 1]
#                   Unpowered vs powered. Simplest split.
#
#   trinary_02_08   [0, 0.2] / (0.2, 0.8) / [0.8, 1]   ← original default
#                   Conservative contested band: only ~8% of rows contested.
#
#   trinary_04_08   [0, 0.4] / (0.4, 0.8) / [0.8, 1]   ← user-requested
#                   "Definitely not / maybe / yes" -- bottom 40% unpowered,
#                   middle 40% contested, top 20% powered.
#
#   trinary_01_09   [0, 0.1] / (0.1, 0.9) / [0.9, 1]
#                   Strict: only the near-certain extremes count as decided.
#
#   four_class      [0, 0.2] / (0.2, 0.5] / (0.5, 0.8) / [0.8, 1]
#                   Unpowered / probably-not / probably-yes / powered.
#                   Provides a fine-grained "borderline" split.
#
#   quartile        [0, 0.25] / (0.25, 0.5] / (0.5, 0.75] / (0.75, 1]
#                   Equal-width quartile partition.
#
# TWO MODES (env var XGB_TRINARY_MODE):
#   insample (default) -- apply saved booster.ubj to aggregated training data.
#                         Fast; overstates performance. Good for exploration.
#   oof                -- refit winner in 5-fold gene-grouped CV.
#                         ~1.5-2h for winner on iclneuro; honest for paper.
#
# WHICH MODELS (env var XGB_TRINARY_CONFIGS):
#   winner (default) -- realistic winner + baseline (2 models)
#   all              -- winner + baseline + all 36 config boosters
#
# INPUT:  Fitting_Data/aggregated_fitting_data_realistic.qs      (from 01e)
#         Models/xgboost/realistic/booster.ubj                     (winner)
#         Models/xgboost/booster.ubj                                (baseline, optional)
#         Models/xgboost/realistic/search/config_*/booster.ubj      (all 36, optional)
#
# OUTPUT: Models/xgboost/realistic/trinary_classification/
#           summary_by_scheme.csv       one row per (model, scheme): overall stats
#           per_class_metrics.csv       precision/recall/F1 per (model, scheme, class)
#           confusion_winner_<scheme>.csv  K×K confusion for the winner per scheme
#           class_f1_grid.png           faceted bar chart: scheme × class × model
#           README.md
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(xgboost)
  library(ggplot2)
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

log_section("STEP 9: N-ary classification metrics for XGBoost power model")

XGB_MODEL_SUBDIR    <- Sys.getenv("XGB_MODEL_SUBDIR",  "realistic")
XGB_FITTING_DATA    <- Sys.getenv("XGB_FITTING_DATA",  "aggregated_fitting_data_realistic.qs")
XGB_TRINARY_MODE    <- tolower(Sys.getenv("XGB_TRINARY_MODE",    "insample"))
XGB_TRINARY_CONFIGS <- tolower(Sys.getenv("XGB_TRINARY_CONFIGS", "winner"))
XGB_NFOLDS          <- as.integer(Sys.getenv("XGB_NFOLDS",   "5"))
XGB_NTHREAD         <- as.integer(Sys.getenv("XGB_NTHREAD",  as.character(min(8, N_CORES))))

stopifnot(XGB_TRINARY_MODE    %in% c("insample", "oof"))
stopifnot(XGB_TRINARY_CONFIGS %in% c("winner", "all"))

set.seed(SEED_BASE)

model_base_realistic <- file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
model_base_baseline  <- file.path(OUTPUT_BASE, "Models", "xgboost")
out_dir <- file.path(model_base_realistic, "trinary_classification")
safe_mkdir(out_dir)

log_info("Mode: %s | Configs: %s", XGB_TRINARY_MODE, XGB_TRINARY_CONFIGS)

# ==============================================================================
# 2. PARTITION SCHEMES
#    Each scheme is a list with:
#      name   -- short string key (used in filenames + CSVs)
#      breaks -- interior boundaries in (0, 1), strictly increasing
#      labels -- length(breaks)+1 class labels
# ==============================================================================

PARTITION_SCHEMES <- list(

  binary = list(
    name   = "binary",
    breaks = 0.5,
    labels = c("Not Powered", "Powered")
  ),

  trinary_02_08 = list(
    name   = "trinary_02_08",
    breaks = c(0.2, 0.8),
    labels = c("Not Powered", "Contested", "Powered")
  ),

  trinary_04_08 = list(
    name   = "trinary_04_08",
    breaks = c(0.4, 0.8),
    labels = c("Definitely Not", "Maybe", "Yes")
  ),

  trinary_01_09 = list(
    name   = "trinary_01_09",
    breaks = c(0.1, 0.9),
    labels = c("Not Powered", "Contested", "Powered")
  ),

  # four_class: full gradient — avoids the binary-sounding "powered/unpowered"
  # at the ends. Class 1 (detection <= 0.2) = Underpowered; class 4 (> 0.8) = Well Powered.
  four_class = list(
    name   = "four_class",
    breaks = c(0.2, 0.5, 0.8),
    labels = c("Underpowered", "Likely Underpowered", "Likely Powered", "Well Powered")
  ),

  quartile = list(
    name   = "quartile",
    breaks = c(0.25, 0.50, 0.75),
    labels = c("Q1: Lowest", "Q2: Low-Mid", "Q3: High-Mid", "Q4: Highest")
  )

)

log_info("Partition schemes: %s", paste(names(PARTITION_SCHEMES), collapse = ", "))

# ==============================================================================
# 3. LOAD FITTING DATA + REBUILD DESIGN MATRIX (mirrors 02_Fit_XGBoost.R)
# ==============================================================================

agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", XGB_FITTING_DATA)
if (!file.exists(agg_path)) {
  stop(sprintf("Fitting data not found: %s. Run 01e_Filter_Realistic.R first.", agg_path))
}
log_info("Loading %s", agg_path)
dat <- qs::qread(agg_path)

CONT_FEATURES <- c("n_donors", "n_cells", "reads_per_cell",
                   "abs_gt_lfc", "gt_basemean", "gt_dispersion")
CAT_FEATURES  <- c("cell_type", "de_method")

dat <- dat %>%
  filter(is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
         n_trials > 0,
         is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
         is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion))
log_info("Rows after finiteness filter: %d (genes: %d)",
         nrow(dat), dplyr::n_distinct(dat$gene_id))

cat_levels <- lapply(CAT_FEATURES, function(f) sort(unique(as.character(dat[[f]]))))
names(cat_levels) <- CAT_FEATURES

build_design <- function(df, cont_features, cat_features, cat_levels) {
  cont_mat <- as.matrix(df[, cont_features, drop = FALSE])
  oh_blocks <- lapply(cat_features, function(f) {
    lv <- cat_levels[[f]]
    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
    colnames(m) <- paste0(f, "=", lv)
    m
  })
  do.call(cbind, c(list(cont_mat), oh_blocks))
}

X <- build_design(dat, CONT_FEATURES, CAT_FEATURES, cat_levels)
y <- dat$detection_rate
log_info("Design matrix: %d rows x %d features", nrow(X), ncol(X))

# ==============================================================================
# 4. CLASSIFICATION HELPERS
# ==============================================================================

# Assign detection_rate values to integer classes 1..K using interior breaks.
# breaks must be strictly increasing values in (0, 1).
# Returns integer vector in 1..length(breaks)+1.
classify_nary <- function(p, breaks) {
  K <- length(breaks) + 1L
  cls <- rep(1L, length(p))
  for (i in seq_along(breaks)) {
    cls[p > breaks[i]] <- i + 1L
  }
  cls
}

# Given integer truth and prediction vectors in 1..K, compute:
#   confusion (K x K matrix, rows=truth, cols=pred)
#   per-class precision / recall / F1 / support
#   macro-F1, weighted-F1, accuracy
compute_metrics <- function(truth_cls, pred_cls, K, class_labels) {
  confusion <- table(factor(truth_cls, levels = seq_len(K)),
                     factor(pred_cls,  levels = seq_len(K)))
  storage.mode(confusion) <- "integer"
  dimnames(confusion) <- list(truth = class_labels, pred = class_labels)

  per_class <- lapply(seq_len(K), function(k) {
    tp <- confusion[k, k]
    fp <- sum(confusion[, k]) - tp
    fn <- sum(confusion[k, ]) - tp
    support   <- sum(confusion[k, ])
    precision <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)
    recall    <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
    f1        <- if (is.na(precision) || is.na(recall) || precision + recall == 0) NA_real_
                 else 2 * precision * recall / (precision + recall)
    data.frame(class = class_labels[k], support = as.integer(support),
               precision = precision, recall = recall, f1 = f1,
               stringsAsFactors = FALSE)
  })
  per_class_df <- do.call(rbind, per_class)

  total       <- sum(confusion)
  accuracy    <- sum(diag(confusion)) / total
  macro_f1    <- mean(per_class_df$f1, na.rm = TRUE)
  weights     <- per_class_df$support / total
  weighted_f1 <- sum(per_class_df$f1 * weights, na.rm = TRUE)

  list(confusion  = confusion,
       per_class  = per_class_df,
       summary    = data.frame(
         n           = as.integer(total),
         accuracy    = accuracy,
         macro_f1    = macro_f1,
         weighted_f1 = weighted_f1,
         stringsAsFactors = FALSE))
}

# ==============================================================================
# 5. COLLECT MODEL BOOSTERS
# ==============================================================================

model_specs <- list()

winner_booster <- file.path(model_base_realistic, "booster.ubj")
if (file.exists(winner_booster)) {
  model_specs[["realistic_winner"]] <- list(
    label = "realistic_winner", path = winner_booster, role = "winner")
} else {
  log_warning("Realistic winner booster not found: %s", winner_booster)
}

baseline_booster <- file.path(model_base_baseline, "booster.ubj")
if (file.exists(baseline_booster)) {
  model_specs[["baseline"]] <- list(
    label = "baseline", path = baseline_booster, role = "baseline")
} else {
  log_info("Baseline booster not found at %s -- skipping baseline", baseline_booster)
}

if (XGB_TRINARY_CONFIGS == "all") {
  cfg_dirs <- list.files(file.path(model_base_realistic, "search"),
                         pattern = "^config_[0-9]+$", full.names = TRUE)
  for (d in cfg_dirs) {
    b <- file.path(d, "booster.ubj")
    if (file.exists(b)) {
      cfg_idx <- as.integer(sub(".*config_", "", basename(d)))
      key <- paste0("config_", cfg_idx)
      model_specs[[key]] <- list(label = key, path = b, role = "grid")
    }
  }
}

if (length(model_specs) == 0) stop("No model boosters found. Unpack the model tars first.")
log_info("Evaluating %d models x %d schemes = %d classification runs",
         length(model_specs), length(PARTITION_SCHEMES),
         length(model_specs) * length(PARTITION_SCHEMES))

# ==============================================================================
# 6. PREDICT (insample) OR REFIT-CV (oof)
#    Predictions are model-level, not scheme-level — compute once per model,
#    then apply each scheme's break-points to discretise.
# ==============================================================================

get_predictions <- function(spec) {
  if (XGB_TRINARY_MODE == "insample") {
    booster <- xgboost::xgb.load(spec$path)
    xgboost::xgb.attr(booster, "nthread") <- as.character(XGB_NTHREAD)
    return(predict(booster, newdata = X))
  }
  # OOF mode: only supported for realistic winner.
  if (spec$label != "realistic_winner") {
    log_warning("OOF mode only for realistic_winner; skipping %s", spec$label)
    return(NULL)
  }
  sel_path <- file.path(model_base_realistic, "model_selection.csv")
  sel <- read.csv(sel_path, stringsAsFactors = FALSE)
  best <- sel[1, ]
  log_info("OOF refit: config %d alpha=%g depth=%d eta=%g nrounds=%d",
           best$config_index, best$alpha, best$max_depth, best$eta, best$best_nrounds)

  gene_ids <- unique(dat$gene_id)
  set.seed(SEED_BASE)
  gene_fold <- sample(rep_len(seq_len(XGB_NFOLDS), length(gene_ids)))
  names(gene_fold) <- gene_ids
  row_fold <- gene_fold[dat$gene_id]

  obs_bin <- cut(y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  bin_prop <- prop.table(table(obs_bin))
  bin_weight_base <- as.numeric(1 / bin_prop[as.character(obs_bin)])
  w_trials <- dat$n_trials
  alpha_val <- best$alpha
  w <- if (alpha_val == 0) w_trials else {
    bw <- bin_weight_base ^ alpha_val
    bw <- bw * length(bw) / sum(bw)
    w_trials * bw
  }

  monotone_map <- c(n_donors = 1L, n_cells = 1L, reads_per_cell = 1L,
                    abs_gt_lfc = 1L, gt_basemean = 1L, gt_dispersion = -1L)
  monotone_vec <- rep(0L, ncol(X)); names(monotone_vec) <- colnames(X)
  for (f in names(monotone_map)) monotone_vec[f] <- monotone_map[f]
  monotone_str <- paste0("(", paste(monotone_vec, collapse = ","), ")")

  params <- list(
    objective            = "binary:logistic",
    eval_metric          = "logloss",
    max_depth            = as.integer(best$max_depth),
    eta                  = as.numeric(best$eta),
    min_child_weight     = as.numeric(best$min_child_weight),
    subsample            = as.numeric(best$subsample),
    colsample_bytree     = as.numeric(best$colsample_bytree),
    tree_method          = "hist",
    monotone_constraints = monotone_str,
    nthread              = XGB_NTHREAD
  )

  oof <- numeric(nrow(X))
  for (k in seq_len(XGB_NFOLDS)) {
    tr <- which(row_fold != k); te <- which(row_fold == k)
    log_info("  fold %d/%d: train=%d test=%d", k, XGB_NFOLDS, length(tr), length(te))
    dtr <- xgboost::xgb.DMatrix(data = X[tr, , drop = FALSE], label = y[tr], weight = w[tr])
    m <- xgboost::xgb.train(params = params, data = dtr,
                            nrounds = as.integer(best$best_nrounds), verbose = 0)
    oof[te] <- predict(m, newdata = X[te, , drop = FALSE])
    rm(dtr, m); gc(verbose = FALSE)
  }
  oof
}

# Get predictions for every model first (each model predicted once, schemes applied later).
model_preds <- list()
for (nm in names(model_specs)) {
  spec <- model_specs[[nm]]
  log_info("Predicting: %s (%s)", spec$label, XGB_TRINARY_MODE)
  preds <- get_predictions(spec)
  if (!is.null(preds)) model_preds[[nm]] <- preds
}

# ==============================================================================
# 7. EVALUATE ALL (model × scheme) COMBINATIONS
# ==============================================================================

summary_rows   <- list()
per_class_rows <- list()
confusion_list <- list()   # keyed by scheme name, winner only

for (scheme_nm in names(PARTITION_SCHEMES)) {
  scheme <- PARTITION_SCHEMES[[scheme_nm]]
  K      <- length(scheme$breaks) + 1L

  truth_cls <- classify_nary(y, scheme$breaks)

  for (nm in names(model_preds)) {
    spec     <- model_specs[[nm]]
    preds    <- model_preds[[nm]]
    pred_cls <- classify_nary(preds, scheme$breaks)

    m <- compute_metrics(truth_cls, pred_cls, K, scheme$labels)

    s <- m$summary
    s$model       <- spec$label
    s$role        <- spec$role
    s$mode        <- XGB_TRINARY_MODE
    s$scheme      <- scheme_nm
    s$breaks      <- paste(scheme$breaks, collapse = "|")
    s$n_classes   <- K
    row_key <- paste0(scheme_nm, "__", nm)
    summary_rows[[row_key]] <- s[, c("model", "role", "mode", "scheme", "breaks",
                                     "n_classes", "n", "accuracy", "macro_f1", "weighted_f1")]

    pc <- m$per_class
    pc$model  <- spec$label
    pc$role   <- spec$role
    pc$scheme <- scheme_nm
    per_class_rows[[row_key]] <- pc[, c("model", "role", "scheme", "class",
                                        "support", "precision", "recall", "f1")]

    if (spec$role == "winner") {
      confusion_list[[scheme_nm]] <- m$confusion
    }

    log_info("  [%s] %s  acc=%.4f  macro_f1=%.4f  weighted_f1=%.4f",
             scheme_nm, spec$label, s$accuracy, s$macro_f1, s$weighted_f1)
  }
}

summary_df    <- do.call(rbind, summary_rows)
per_class_df  <- do.call(rbind, per_class_rows)

write.csv(summary_df,
          file.path(out_dir, "summary_by_scheme.csv"),   row.names = FALSE)
write.csv(per_class_df,
          file.path(out_dir, "per_class_metrics.csv"),   row.names = FALSE)

# Confusion matrices: one CSV per scheme for the winner
for (scheme_nm in names(confusion_list)) {
  conf_df <- as.data.frame.table(confusion_list[[scheme_nm]])
  colnames(conf_df) <- c("truth_class", "pred_class", "n")
  write.csv(conf_df,
            file.path(out_dir, sprintf("confusion_winner_%s.csv", scheme_nm)),
            row.names = FALSE)
}

# ==============================================================================
# 8. FIGURE: FACETED F1 BAR CHART (scheme × class, winner vs baseline)
# ==============================================================================

# Scheme display order (simpler → richer)
scheme_order <- c("binary", "trinary_01_09", "trinary_02_08",
                  "trinary_04_08", "quartile", "four_class")

plot_df <- per_class_df %>%
  filter(role %in% c("winner", "baseline")) %>%
  mutate(
    scheme = factor(scheme, levels = intersect(scheme_order, unique(scheme))),
    model  = factor(model, levels = c("realistic_winner", "baseline"),
                    labels = c("Realistic winner", "Baseline"))
  )

if (nrow(plot_df) > 0) {
  # Scheme label prettifier
  scheme_label_map <- c(
    binary         = "Binary (0.5)",
    trinary_01_09  = "Trinary strict\n[0–0.1] [0.1–0.9] [0.9–1]",
    trinary_02_08  = "Trinary standard\n[0–0.2] [0.2–0.8] [0.8–1]",
    trinary_04_08  = "Trinary skewed\n[0–0.4] [0.4–0.8] [0.8–1]",
    quartile       = "Quartiles\n[0–0.25] [0.25–0.5] [0.5–0.75] [0.75–1]",
    four_class     = "Four-class\n[0–0.2] [0.2–0.5] [0.5–0.8] [0.8–1]"
  )
  plot_df <- plot_df %>%
    mutate(scheme_label = factor(
      scheme_label_map[as.character(scheme)],
      levels = scheme_label_map[scheme_order[scheme_order %in% unique(as.character(scheme))]]))

  p <- ggplot(plot_df, aes(x = class, y = f1, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, na.rm = TRUE) +
    geom_text(aes(label = sprintf("%.2f", f1)),
              position = position_dodge(width = 0.7), vjust = -0.3,
              size = 2.5, na.rm = TRUE) +
    facet_wrap(~ scheme_label, scales = "free_x", ncol = 2) +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "F1 score", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "top",
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
          strip.text       = element_text(size = 8))

  ggsave(file.path(out_dir, "class_f1_grid.png"),
         p, width = 10, height = 11, dpi = 150)
  log_info("Saved class_f1_grid.png")
}

# Per-scheme individual bar plots (easier to inspect one at a time)
for (scheme_nm in names(PARTITION_SCHEMES)) {
  scheme_df <- plot_df %>% filter(as.character(scheme) == scheme_nm)
  if (nrow(scheme_df) == 0) next
  p_single <- ggplot(scheme_df, aes(x = class, y = f1, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, na.rm = TRUE) +
    geom_text(aes(label = sprintf("%.2f", f1)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3.5, na.rm = TRUE) +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "F1 score", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, sprintf("class_f1_%s.png", scheme_nm)),
         p_single, width = max(4, length(PARTITION_SCHEMES[[scheme_nm]]$labels) * 1.4),
         height = 4, dpi = 150)
}

# ==============================================================================
# 9. SUMMARY TABLE: BEST SCHEME BY WEIGHTED F1 FOR EACH ROLE
# ==============================================================================

best_per_role <- summary_df %>%
  group_by(role, model) %>%
  slice_max(weighted_f1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(model, role, scheme, n_classes, accuracy, macro_f1, weighted_f1)

write.csv(best_per_role,
          file.path(out_dir, "best_scheme_per_model.csv"), row.names = FALSE)

log_info("Best scheme per model (by weighted_f1):")
for (i in seq_len(nrow(best_per_role))) {
  r <- best_per_role[i, ]
  log_info("  %-20s  scheme=%-18s  wF1=%.4f  acc=%.4f",
           r$model, r$scheme, r$weighted_f1, r$accuracy)
}

# ==============================================================================
# 10. README
# ==============================================================================

scheme_summary <- paste(
  vapply(names(PARTITION_SCHEMES), function(nm) {
    s <- PARTITION_SCHEMES[[nm]]
    sprintf("  %-20s  breaks=[%s]  classes=%s",
            nm, paste(s$breaks, collapse = ", "), paste(s$labels, collapse = " / "))
  }, character(1)),
  collapse = "\n")

readme_text <- sprintf(
"# N-ary classification metrics

Generated %s by 09_Trinary_Classification.R.

## Partition schemes evaluated

%s

## Files

- `summary_by_scheme.csv`         one row per (model, scheme): accuracy, macro-F1, weighted-F1
- `per_class_metrics.csv`         precision / recall / F1 per (model, scheme, class)
- `confusion_winner_<scheme>.csv` confusion matrix for the realistic winner (K x K)
- `best_scheme_per_model.csv`     best-scoring scheme per model (by weighted-F1)
- `class_f1_grid.png`             faceted overview: all schemes x winner vs baseline
- `class_f1_<scheme>.png`         individual bar chart per scheme

## Mode

`XGB_TRINARY_MODE = %s`

%s

## Caveats

- The 'baseline' model was trained on the unfiltered 7.2M-row grid and is
  evaluated here on the 3.4M-row realistic subset. This is not a symmetric
  comparison: the baseline saw more training data of a different distribution.
  See xgboost_section_plan.md section 7 for a note on this asymmetry.
- 'Weighted-F1' weights each class by its support (row count). For schemes
  where unpowered rows dominate (>85%% of the grid), even a degenerate model
  that predicts 'unpowered' everywhere scores a high weighted-F1. Macro-F1
  is the fairer metric.
- Boundaries are set on the detection_rate (fraction of 20 iterations), not
  on the XGBoost predicted probability. Both are discretised using the same
  break-points for comparability.
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  scheme_summary,
  XGB_TRINARY_MODE,
  if (XGB_TRINARY_MODE == "insample")
    "**insample**: booster applied to its own training data. Overstates true generalisation."
  else
    "**oof**: 5-fold gene-grouped CV refit at the winner's best_nrounds. Honest estimate."
)

writeLines(readme_text, file.path(out_dir, "README.md"))

log_info("All N-ary classification outputs written to %s", out_dir)
log_section("STEP 9 complete")
