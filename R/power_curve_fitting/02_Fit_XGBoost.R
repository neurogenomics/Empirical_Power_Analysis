# ==============================================================================
# SCRIPT: Fit XGBoost Power Model (M2 -- replaces logistic regression)
# PURPOSE: Train a gradient-boosted regression-tree model of detection RATE over
#          the empirical design grid, replacing the 6-term additive logit
#          (02_Fit_Models.R). One model spans all cell types and DE methods;
#          tree structure captures the interactions the logit cannot.
#
# TWO RUN MODES (selected by env var XGB_CONFIG_INDEX):
#   * FULL SEARCH (XGB_CONFIG_INDEX unset / 0) -- local / single-node use.
#       Loops over the whole param_grid, picks the best config by grouped-CV
#       logloss, refits, computes OOF calibration, writes to Models/xgboost/.
#   * ARRAY MODE  (XGB_CONFIG_INDEX = k, 1..nrow(param_grid)) -- HPC Stage B.
#       Trains exactly ONE config (row k of param_grid): grouped CV to get the
#       score + early-stopped nrounds, then a full fit on all data at that
#       config. Writes to Models/xgboost/search/config_<k>/ so a PBS array can
#       fan the search out one config per task. Stage C (03_Select_XGBoost.R)
#       reads every config's cv_metrics.csv, picks the best, and runs SHAP +
#       calibration on the already-trained best booster.
#
# INPUT:  Fitting_Data/aggregated_fitting_data_all.qs  (from 01b)
# OUTPUT (full search): Models/xgboost/{booster.ubj,booster.qs,feature_spec.qs,
#           cv_metrics.csv,feature_importance.csv,calibration.csv}
# OUTPUT (array mode):  Models/xgboost/search/config_<k>/{booster.ubj,booster.qs,
#           feature_spec.qs,cv_metrics.csv,config.json}
#
# Modelling choices (locked 2026-07-15):
#   A run pooling -> handled in 01b (PCF_POOL_RUNS); this script is agnostic
#   B objective   -> weighted binary:logistic   [XGB_OBJECTIVE]
#   C monotone    -> ON                          [XGB_MONOTONE]
#   D de_method   -> single model, method as a feature (always)
#   full gene set only (shared dropped in 01b)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
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

# ------------------------------------------------------------------------------
# Config knobs (env-overridable; defaults = locked decisions)
# ------------------------------------------------------------------------------
XGB_OBJECTIVE    <- Sys.getenv("XGB_OBJECTIVE", "binary:logistic")   # decision B
XGB_MONOTONE     <- tolower(Sys.getenv("XGB_MONOTONE", "true")) %in% c("true", "1", "yes")  # decision C
XGB_NTHREAD      <- as.integer(Sys.getenv("XGB_NTHREAD", as.character(min(8, N_CORES))))
XGB_NFOLDS       <- as.integer(Sys.getenv("XGB_NFOLDS", "5"))
XGB_GENE_SET     <- Sys.getenv("XGB_GENE_SET", "full")   # only 'full' survives 01b
XGB_MAX_NROUNDS  <- as.integer(Sys.getenv("XGB_MAX_NROUNDS", "5000"))
XGB_EARLY_STOP   <- as.integer(Sys.getenv("XGB_EARLY_STOP", "50"))
# Array mode: 0 / unset => full local search; k>=1 => train only param_grid[k, ].
XGB_CONFIG_INDEX <- as.integer(Sys.getenv("XGB_CONFIG_INDEX", "0"))

# ---- realistic-grid reweighting knobs (xgboost_reweighting_design.md) --------
# XGB_FITTING_DATA: filename under Fitting_Data/ (default = baseline all-rows).
#   Set to "aggregated_fitting_data_realistic.qs" after running 01e_Filter_Realistic.R
#   to train on the n_cells >= 100 sub-grid.
# XGB_MODEL_SUBDIR: subdir under Models/xgboost/ ("" = baseline, "realistic" = the
#   reweighting sweep). Keeps the old baseline artefacts untouched for the
#   comparison figure in 08_Compare_Reweighted.R.
# XGB_ALPHA: per-config default exponent for the inverse-frequency bin weight.
#   In ARRAY mode, if the param_grid row carries an `alpha` column, THAT value
#   wins (the array indexes into the {config x alpha} joint grid); the env var
#   is only used when the grid has no alpha column (full-search legacy path).
XGB_FITTING_DATA <- Sys.getenv("XGB_FITTING_DATA", "aggregated_fitting_data_all.qs")
XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "")
XGB_ALPHA_ENV    <- as.numeric(Sys.getenv("XGB_ALPHA", "0"))

ARRAY_MODE <- !is.na(XGB_CONFIG_INDEX) && XGB_CONFIG_INDEX >= 1

log_section(sprintf("STEP 2 (XGBoost): Fit gradient-boosted power model [%s]",
                    if (ARRAY_MODE) sprintf("ARRAY config %d", XGB_CONFIG_INDEX) else "FULL SEARCH"))

set.seed(SEED_BASE)

log_info("Objective=%s | monotone=%s | nthread=%d | folds=%d | gene_set=%s | max_nrounds=%d | early_stop=%d",
         XGB_OBJECTIVE, XGB_MONOTONE, XGB_NTHREAD, XGB_NFOLDS, XGB_GENE_SET,
         XGB_MAX_NROUNDS, XGB_EARLY_STOP)

# ==============================================================================
# 2. LOAD AGGREGATED DATA
# ==============================================================================

agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", XGB_FITTING_DATA)
if (!file.exists(agg_path)) {
  stop(sprintf("Fitting data not found: %s. Run 01b_Aggregate_Fitting_Data.R (baseline) or 01e_Filter_Realistic.R (realistic sub-grid) first.", agg_path))
}
log_info("Reading fitting data: %s", agg_path)

dat <- qs::qread(agg_path)
log_info("Loaded aggregated table: %d rows", nrow(dat))

if ("gene_set" %in% colnames(dat)) {
  dat <- dat %>% filter(gene_set == XGB_GENE_SET)
  log_info("Filtered to gene_set='%s': %d rows", XGB_GENE_SET, nrow(dat))
}

# Continuous + categorical features.
# Raw (un-logged) predictors: XGBoost splits on rank order, so a monotone log
# transform is inert for the tree model (identical partitions, identical
# predictions). Kept raw so SHAP / interpretation axes are in native units
# (donors, |LFC|, baseMean, ...) rather than log space.
CONT_FEATURES <- c("n_donors", "n_cells", "reads_per_cell",
                   "abs_gt_lfc", "gt_basemean", "gt_dispersion")
CAT_FEATURES  <- c("cell_type", "de_method")

# Monotone direction per continuous feature (decision C): power is
# non-decreasing in donors/cells/reads/|LFC|/baseMean, non-increasing in dispersion.
CONT_MONOTONE <- c(n_donors = 1L, n_cells = 1L, reads_per_cell = 1L,
                   abs_gt_lfc = 1L, gt_basemean = 1L, gt_dispersion = -1L)

# Drop rows with non-finite predictors or target
dat <- dat %>%
  filter(
    is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
    n_trials > 0,
    is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
    is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion)
  )
log_info("After finiteness filter: %d rows, %d genes",
         nrow(dat), dplyr::n_distinct(dat$gene_id))

# ==============================================================================
# 3. BUILD DESIGN MATRIX (one-hot categoricals; monotone-aligned)
# ==============================================================================

# One-hot encode categoricals with a stored level map so predict_power can
# reproduce columns exactly. Monotone constraint is 0 for one-hot columns.
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
w_trials <- dat$n_trials    # binomial-trust weight per row; used inside a bin
feature_names <- colnames(X)

# ---- inverse-frequency bin-weight base (design doc §4) ----------------------
# Observed detection rate binned into 10 equal-width bins. bin_weight_base is
# 1 / (bin proportion) per row; the per-config alpha exponent is applied inside
# fit_one_config so the same base can be reused across alphas without reload.
obs_bin           <- cut(y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
bin_prop          <- prop.table(table(obs_bin))
bin_weight_base   <- as.numeric(1 / bin_prop[as.character(obs_bin)])
log_info("Observed-power bin proportions: %s",
         paste(sprintf("%s=%.3f", names(bin_prop), as.numeric(bin_prop)), collapse = " "))

#' Combine per-row trial weight with inverse-frequency bin weight (design doc §4).
#'
#' alpha=0 collapses to w = n_trials (baseline). alpha in (0,1] tempers the
#' inverse-frequency weight; alpha=1 = full class-balanced. The bin_weight is
#' normalised so its mean is 1, preserving the total effective sample size and
#' keeping XGBoost's min_child_weight / regularisation scaling comparable
#' across alphas.
make_weights <- function(alpha) {
  if (!is.finite(alpha) || alpha == 0) return(w_trials)
  bw <- bin_weight_base ^ alpha
  bw <- bw * length(bw) / sum(bw)   # normalise so sum(bw) == length(bw)
  w_trials * bw
}

# Monotone constraint vector aligned to feature_names
monotone_vec <- rep(0L, length(feature_names))
names(monotone_vec) <- feature_names
if (XGB_MONOTONE) {
  for (f in names(CONT_MONOTONE)) monotone_vec[f] <- CONT_MONOTONE[f]
}
monotone_str <- paste0("(", paste(monotone_vec, collapse = ","), ")")

log_info("Design matrix: %d rows x %d features", nrow(X), ncol(X))
log_info("Features: %s", paste(feature_names, collapse = ", "))

# ==============================================================================
# 4. GENE-GROUPED FOLDS (design doc §10.7, Option A)
# ==============================================================================
# Pure gene grouping (all rows of a gene stay in one fold) prevents gene
# leakage between train and test. No stratification: with millions of rows
# across 5 folds, random assignment gives a good enough bin balance for the
# purposes here (Chris's call).

gene_ids <- unique(dat$gene_id)
set.seed(SEED_BASE)
gene_fold <- sample(rep_len(seq_len(XGB_NFOLDS), length(gene_ids)))
names(gene_fold) <- gene_ids
row_fold <- gene_fold[dat$gene_id]
folds <- lapply(seq_len(XGB_NFOLDS), function(k) which(row_fold == k))
log_info("Gene-grouped %d-fold on %d genes (fold sizes: %s)",
         XGB_NFOLDS, length(gene_ids), paste(sapply(folds, length), collapse = ", "))

# ==============================================================================
# 5. HYPERPARAMETER GRID (single source of truth)
# ==============================================================================
# The HPC search array's task count MUST equal nrow(param_grid). Widen/narrow
# here and the array -J bound follows (see hpc/cx3-xgb-search.pbs).
#
# For the realistic-grid reweighting sweep the grid is expanded to include an
# `alpha` column ({0, 0.5, 1.0}), giving 12 hyperparameter configs x 3 alphas
# = 36 rows. Baseline runs (XGB_MODEL_SUBDIR="") keep the 12-row grid via a
# single alpha=0 value so XGB_CONFIG_INDEX mapping is unchanged from the
# original 12-config search.

XGB_REWEIGHT <- nzchar(XGB_MODEL_SUBDIR)
alpha_values <- if (XGB_REWEIGHT) c(0, 0.5, 1.0) else 0

param_grid <- expand.grid(
  alpha            = alpha_values,
  max_depth        = c(3L, 4L, 6L),
  eta              = c(0.03, 0.1),
  min_child_weight = c(1, 10),
  subsample        = 0.8,
  colsample_bytree = 0.8,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
# Put alpha last for legibility (matches the "hyperparameter config x alpha"
# framing in the design doc); order is stable so XGB_CONFIG_INDEX still maps
# deterministically to a param_grid row.
param_grid <- param_grid[, c("max_depth", "eta", "min_child_weight",
                              "subsample", "colsample_bytree", "alpha")]
N_CONFIGS <- nrow(param_grid)
log_info("Param grid: %d configs (%s alphas x %d hyperparameter configs)",
         N_CONFIGS,
         paste(alpha_values, collapse = ","),
         N_CONFIGS / length(alpha_values))

make_params <- function(p) {
  list(
    objective            = XGB_OBJECTIVE,
    eval_metric          = "logloss",
    max_depth            = p$max_depth,
    eta                  = p$eta,
    min_child_weight     = p$min_child_weight,
    subsample            = p$subsample,
    colsample_bytree     = p$colsample_bytree,
    tree_method          = "hist",
    monotone_constraints = monotone_str,
    nthread              = XGB_NTHREAD
  )
}

#' Grouped CV for one config -> (cv_logloss, best_nrounds, OOF preds), then full fit.
#'
#' Each config now carries its own alpha (bin-weight exponent), so the training
#' weight vector is rebuilt per call. OOF predictions are rebuilt by a manual
#' per-fold refit at best_iter rounds (Fix A, 2026-07-31): xgb.cv(prediction=TRUE)
#' returns cv$pred == NULL in xgboost 3.x, so we can't rely on it. The manual
#' refit is version-proof and roughly doubles per-config wall time (still under
#' the realistic-search 1h walltime for the sizes we see, ~16 min -> ~30 min).
fit_one_config <- function(p) {
  params    <- make_params(p)
  alpha_val <- if (!is.null(p$alpha)) as.numeric(p$alpha) else XGB_ALPHA_ENV
  w_cfg     <- make_weights(alpha_val)
  dtrain    <- xgb.DMatrix(data = X, label = y, weight = w_cfg)

  cv <- xgb.cv(
    params                = params,
    data                  = dtrain,
    nrounds               = XGB_MAX_NROUNDS,
    folds                 = folds,
    early_stopping_rounds = XGB_EARLY_STOP,
    verbose               = 0
  )
  metric_col <- grep("test_.*_mean", colnames(cv$evaluation_log), value = TRUE)[1]
  # Robustly resolve the best iteration across xgboost versions. In xgboost 3.x
  # xgb.cv() no longer returns a length-1 $best_iteration, so nrounds would be
  # length-zero and xgb.train() aborts ("if (nrounds == 0): argument is of
  # length zero"). Fall back to the attribute, then to which.min on the eval log
  # (guaranteed valid: metric_col is already used for best_score below).
  best_iter <- cv$best_iteration
  if (is.null(best_iter) || length(best_iter) != 1L || is.na(best_iter))
    best_iter <- attr(cv, "best_iteration")
  if (is.null(best_iter) || length(best_iter) != 1L || is.na(best_iter))
    best_iter <- which.min(cv$evaluation_log[[metric_col]])
  best_iter  <- as.integer(best_iter)
  best_score <- cv$evaluation_log[[metric_col]][best_iter]

  # ---- Manual OOF refit (xgboost-3.x-safe) --------------------------------
  # Refit each fold's TRAIN portion at best_iter rounds, predict its held-out
  # rows, assemble the length-N OOF vector keyed by original row_id. No
  # early-stopping here (best_iter is already the CV-selected round count).
  oof_pred <- numeric(nrow(X))
  n_all    <- nrow(X)
  for (k in seq_along(folds)) {
    test_idx  <- folds[[k]]
    train_idx <- setdiff(seq_len(n_all), test_idx)
    dtr_k <- xgb.DMatrix(data   = X[train_idx, , drop = FALSE],
                         label  = y[train_idx],
                         weight = w_cfg[train_idx])
    dte_k <- xgb.DMatrix(data   = X[test_idx,  , drop = FALSE])
    m_k <- xgb.train(params = params, data = dtr_k,
                     nrounds = best_iter, verbose = 0)
    oof_pred[test_idx] <- predict(m_k, dte_k)
    rm(dtr_k, dte_k, m_k); invisible(gc(verbose = FALSE))
  }
  if (any(!is.finite(oof_pred))) {
    stop(sprintf("Non-finite OOF prediction(s) after manual refit (n=%d)",
                 sum(!is.finite(oof_pred))))
  }

  model <- xgb.train(params = params, data = dtrain,
                     nrounds = best_iter, verbose = 0)
  list(model        = model,
       params       = params,
       cv_logloss   = best_score,
       best_nrounds = best_iter,
       alpha        = alpha_val,
       oof_pred     = oof_pred)
}

feature_spec_of <- function(best_params, best_nrounds, alpha_val = NA_real_) {
  list(
    cont_features = CONT_FEATURES,
    cat_features  = CAT_FEATURES,
    cat_levels    = cat_levels,
    feature_names = feature_names,
    monotone      = monotone_vec,
    objective     = XGB_OBJECTIVE,
    gene_set      = XGB_GENE_SET,
    best_params   = best_params,
    best_nrounds  = best_nrounds,
    alpha         = alpha_val,
    fitting_data  = XGB_FITTING_DATA,
    model_subdir  = XGB_MODEL_SUBDIR
  )
}

# Resolve model root once so both ARRAY_MODE and FULL_SEARCH write to the
# same subtree (Models/xgboost/ or Models/xgboost/<subdir>/).
xgb_model_base <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Models", "xgboost")
}

# ==============================================================================
# 6a. ARRAY MODE -- train exactly one config, write to search/config_<k>/
# ==============================================================================

if (ARRAY_MODE) {
  if (XGB_CONFIG_INDEX > N_CONFIGS) {
    stop(sprintf("XGB_CONFIG_INDEX=%d exceeds param_grid size (%d)", XGB_CONFIG_INDEX, N_CONFIGS))
  }
  p <- param_grid[XGB_CONFIG_INDEX, ]
  alpha_val <- if (!is.null(p$alpha)) as.numeric(p$alpha) else XGB_ALPHA_ENV
  log_section(sprintf("Array config %d/%d: depth=%d eta=%.2f mcw=%g subsample=%g colsample=%g alpha=%g",
                      XGB_CONFIG_INDEX, N_CONFIGS, p$max_depth, p$eta,
                      p$min_child_weight, p$subsample, p$colsample_bytree, alpha_val))

  res <- fit_one_config(p)
  log_info("cv_logloss=%.5f @ %d rounds  (alpha=%g)",
           res$cv_logloss, res$best_nrounds, res$alpha)

  out_dir <- file.path(xgb_model_base, "search", sprintf("config_%d", XGB_CONFIG_INDEX))
  safe_mkdir(out_dir)

  xgb.save(res$model, file.path(out_dir, "booster.ubj"))
  qs::qsave(res$model, file.path(out_dir, "booster.qs"))
  qs::qsave(feature_spec_of(res$params, res$best_nrounds, res$alpha),
            file.path(out_dir, "feature_spec.qs"))

  cv_row <- cbind(config_index = XGB_CONFIG_INDEX, p,
                  best_nrounds = res$best_nrounds, cv_logloss = res$cv_logloss)
  # Ensure alpha is always present in the CV row even in baseline (non-reweight)
  # mode so downstream selection (03) can uniformly key on it.
  if (!"alpha" %in% colnames(cv_row)) cv_row$alpha <- res$alpha
  write.csv(cv_row, file.path(out_dir, "cv_metrics.csv"), row.names = FALSE)

  cfg_json <- sprintf(
    '{"config_index":%d,"max_depth":%d,"eta":%g,"min_child_weight":%g,"subsample":%g,"colsample_bytree":%g,"alpha":%g,"best_nrounds":%d,"cv_logloss":%.6f,"fitting_data":"%s"}',
    XGB_CONFIG_INDEX, p$max_depth, p$eta, p$min_child_weight, p$subsample,
    p$colsample_bytree, res$alpha, res$best_nrounds, res$cv_logloss,
    XGB_FITTING_DATA)
  writeLines(cfg_json, file.path(out_dir, "config.json"))

  # OOF predictions keyed by original row index. Downstream (03_Select) reads
  # these to compute bin-balanced calibration metrics without refitting folds.
  # Trial weight (n_trials) is stored for calibration reporting; the training
  # weight used to fit is n_trials * bin_weight(alpha) and is recoverable from
  # bin_weight(y, alpha) if needed.
  meta_cols <- intersect(c("gene_id", "cell_type", "de_method",
                            "n_donors", "n_cells", "reads_per_cell",
                            "abs_gt_lfc", "gt_basemean", "gt_dispersion"),
                          colnames(dat))
  oof_df <- data.frame(
    row_id  = seq_len(nrow(X)),
    fold    = row_fold,
    obs     = y,
    pred    = res$oof_pred,
    w_trials = w_trials,
    alpha   = res$alpha,
    config_index = XGB_CONFIG_INDEX
  )
  oof_df <- cbind(oof_df, dat[, meta_cols, drop = FALSE])
  qs::qsave(oof_df, file.path(out_dir, "oof_predictions.qs"))
  log_info("Saved OOF predictions: %d rows -> %s", nrow(oof_df),
           file.path(out_dir, "oof_predictions.qs"))

  log_section("STEP 2 (XGBoost) ARRAY TASK COMPLETE")
  log_info("Config %d written to %s", XGB_CONFIG_INDEX, out_dir)
  quit(save = "no", status = 0)
}

# ==============================================================================
# 6b. FULL SEARCH -- loop all configs, pick best (local / single-node)
# ==============================================================================

log_section(sprintf("Hyperparameter search (%d configs, grouped CV)", N_CONFIGS))

cv_records <- list()
best <- list(score = Inf, params = NULL, nrounds = NA_integer_,
             model = NULL, alpha = NA_real_, oof_pred = NULL)

for (i in seq_len(N_CONFIGS)) {
  p <- param_grid[i, ]
  res <- fit_one_config(p)
  cv_records[[i]] <- cbind(config_index = i, p,
                           best_nrounds = res$best_nrounds, cv_logloss = res$cv_logloss)
  log_info("  cfg %2d/%2d depth=%d eta=%.2f mcw=%g alpha=%g -> logloss=%.5f @ %d rounds",
           i, N_CONFIGS, p$max_depth, p$eta, p$min_child_weight, res$alpha,
           res$cv_logloss, res$best_nrounds)
  if (res$cv_logloss < best$score) {
    best$score    <- res$cv_logloss
    best$params   <- res$params
    best$nrounds  <- res$best_nrounds
    best$grid     <- p
    best$model    <- res$model
    best$alpha    <- res$alpha
    best$oof_pred <- res$oof_pred
  }
}

cv_metrics <- bind_rows(cv_records)
log_info("Best config: depth=%d eta=%.2f mcw=%g alpha=%g | logloss=%.5f @ %d rounds",
         best$grid$max_depth, best$grid$eta, best$grid$min_child_weight,
         best$alpha, best$score, best$nrounds)

# ------------------------------------------------------------------------------
# OOF calibration for the winner (already produced by xgb.cv above -- no refit)
# ------------------------------------------------------------------------------
oof_pred    <- best$oof_pred
final_model <- best$model
final_w     <- make_weights(best$alpha)   # trial x bin weight used for training

wrmse <- sqrt(sum(final_w * (oof_pred - y)^2, na.rm = TRUE) / sum(final_w))
eps <- 1e-7
pc <- pmin(pmax(oof_pred, eps), 1 - eps)
wlogloss <- -sum(final_w * (y * log(pc) + (1 - y) * log(1 - pc)), na.rm = TRUE) / sum(final_w)
log_info("XGBoost OOF (weighted): RMSE=%.4f  logloss=%.4f", wrmse, wlogloss)

# Calibration reporting uses trial weights (not the training bin weights) so it
# stays comparable across alphas.
calib <- data.frame(pred = oof_pred, obs = y, w = w_trials) %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    n         = dplyr::n(),
    mean_pred = weighted.mean(pred, w),
    mean_obs  = weighted.mean(obs, w),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
model_out_dir <- xgb_model_base
safe_mkdir(model_out_dir)

xgb.save(final_model, file.path(model_out_dir, "booster.ubj"))
qs::qsave(final_model, file.path(model_out_dir, "booster.qs"))
qs::qsave(feature_spec_of(best$params, best$nrounds, best$alpha),
          file.path(model_out_dir, "feature_spec.qs"))

write.csv(cv_metrics, file.path(model_out_dir, "cv_metrics.csv"), row.names = FALSE)
write.csv(calib, file.path(model_out_dir, "calibration.csv"), row.names = FALSE)

imp <- xgb.importance(model = final_model)
write.csv(imp, file.path(model_out_dir, "feature_importance.csv"), row.names = FALSE)

log_section("STEP 2 (XGBoost) COMPLETE")
log_info("Model + diagnostics saved to: %s", model_out_dir)
log_info("OOF weighted RMSE=%.4f logloss=%.4f", wrmse, wlogloss)
