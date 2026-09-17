# ==============================================================================
# SCRIPT: Regenerate baseline per-row OOF predictions (M2 -- calibration bands)
# PURPOSE: The baseline full-grid XGBoost model was fitted in FULL SEARCH mode
#          (02_Fit_XGBoost.R section 6b), which computes an out-of-fold (OOF)
#          prediction vector in memory but only PERSISTS the binned calibration
#          summary (calibration.csv / calibration_residuals.csv) -- not the
#          per-row OOF. The 5 realistic-grid SHAP configs, fitted in ARRAY mode,
#          DO persist per-row OOF (oof_predictions.qs). To draw a consistent
#          10-90th percentile band of observed detection rate on ALL SIX
#          calibration panels in 08_Compare_Reweighted.R, the baseline needs the
#          same per-row OOF the other five already have.
#
#          This script regenerates ONLY the baseline OOF, cheaply and without
#          touching the promoted booster: it reloads the baseline winning config
#          from the saved feature_spec.qs and reproduces 02's exact
#          design-matrix + gene-grouped 5-fold + per-fold refit logic
#          (set.seed(SEED_BASE), deterministic), then writes oof_predictions.qs
#          alongside the baseline model with the same schema as the config OOFs.
#
#          CPU-only, no HPC. ~10-30 min (5 fold refits at best_nrounds on the
#          full grid). Does NOT refit the promoted model or alter any existing
#          artefact -- it only ADDS Models/xgboost/oof_predictions.qs.
#
# INPUT:  Models/xgboost/feature_spec.qs                 (baseline winning config)
#         Fitting_Data/aggregated_fitting_data_all.qs    (baseline full grid)
# OUTPUT: Models/xgboost/oof_predictions.qs
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

log_section("STEP 8b: Regenerate baseline per-row OOF for calibration bands")

# Mirror 02_Fit_XGBoost.R knobs that shape the design matrix / folds.
XGB_NFOLDS   <- as.integer(Sys.getenv("XGB_NFOLDS", "5"))
XGB_NTHREAD  <- as.integer(Sys.getenv("XGB_NTHREAD", as.character(min(8, N_CORES))))
XGB_GENE_SET <- Sys.getenv("XGB_GENE_SET", "full")

xgb_model_base <- file.path(OUTPUT_BASE, "Models", "xgboost")
fs_path <- file.path(xgb_model_base, "feature_spec.qs")
if (!file.exists(fs_path)) stop(sprintf("Baseline feature_spec.qs not found: %s", fs_path))
fs <- qs::qread(fs_path)

# Baseline fitting data: feature_spec may carry an empty fitting_data string for
# the original (pre-provenance-field) baseline run; default to the all-rows grid.
baseline_fitting <- if (!is.null(fs$fitting_data) && nzchar(fs$fitting_data)) {
  fs$fitting_data
} else {
  "aggregated_fitting_data_all.qs"
}
best_nrounds <- as.integer(fs$best_nrounds)
best_params  <- fs$best_params
# Baseline is the unreweighted model (alpha = 0 => training weight == n_trials).
alpha_val <- if (!is.null(fs$alpha) && length(fs$alpha) == 1L && is.finite(fs$alpha)) as.numeric(fs$alpha) else 0

log_info("Baseline config: depth=%d eta=%g mcw=%g | nrounds=%d | alpha=%g | data=%s",
         best_params$max_depth, best_params$eta, best_params$min_child_weight,
         best_nrounds, alpha_val, baseline_fitting)

# Ensure nthread on the reused param list matches this host.
best_params$nthread <- XGB_NTHREAD

# ==============================================================================
# 2. LOAD DATA + REPRODUCE 02's FILTER (identical order -> identical folds)
# ==============================================================================

agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", baseline_fitting)
if (!file.exists(agg_path)) stop(sprintf("Fitting data not found: %s", agg_path))
log_info("Reading fitting data: %s", agg_path)
dat <- qs::qread(agg_path)
log_info("Loaded aggregated table: %d rows", nrow(dat))

if ("gene_set" %in% colnames(dat)) {
  dat <- dat %>% filter(gene_set == XGB_GENE_SET)
  log_info("Filtered to gene_set='%s': %d rows", XGB_GENE_SET, nrow(dat))
}

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
# 3. DESIGN MATRIX (verbatim from 02_Fit_XGBoost.R section 3)
# ==============================================================================

CONT_FEATURES <- c("n_donors", "n_cells", "reads_per_cell",
                   "abs_gt_lfc", "gt_basemean", "gt_dispersion")
CAT_FEATURES  <- c("cell_type", "de_method")

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
w_trials <- dat$n_trials
# Baseline: alpha = 0 -> training weight is the raw per-row trial weight.
w_cfg <- w_trials

log_info("Design matrix: %d rows x %d features", nrow(X), ncol(X))
if (!identical(colnames(X), fs$feature_names)) {
  log_warning("Design columns differ from feature_spec$feature_names; reordering to match.")
  X <- X[, fs$feature_names, drop = FALSE]
}

# ==============================================================================
# 4. GENE-GROUPED FOLDS (verbatim from 02_Fit_XGBoost.R section 4)
# ==============================================================================

gene_ids <- unique(dat$gene_id)
set.seed(SEED_BASE)
gene_fold <- sample(rep_len(seq_len(XGB_NFOLDS), length(gene_ids)))
names(gene_fold) <- gene_ids
row_fold <- gene_fold[dat$gene_id]
folds <- lapply(seq_len(XGB_NFOLDS), function(k) which(row_fold == k))
log_info("Gene-grouped %d-fold on %d genes (fold sizes: %s)",
         XGB_NFOLDS, length(gene_ids), paste(sapply(folds, length), collapse = ", "))

# ==============================================================================
# 5. PER-FOLD REFIT -> OOF (verbatim from 02_Fit_XGBoost.R fit_one_config)
# ==============================================================================

oof_pred <- numeric(nrow(X))
n_all    <- nrow(X)
for (k in seq_along(folds)) {
  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(n_all), test_idx)
  dtr_k <- xgb.DMatrix(data   = X[train_idx, , drop = FALSE],
                       label  = y[train_idx],
                       weight = w_cfg[train_idx])
  dte_k <- xgb.DMatrix(data   = X[test_idx,  , drop = FALSE])
  m_k <- xgb.train(params = best_params, data = dtr_k,
                   nrounds = best_nrounds, verbose = 0)
  oof_pred[test_idx] <- predict(m_k, dte_k)
  log_info("  fold %d/%d: refit on %d rows, predicted %d held-out rows",
           k, XGB_NFOLDS, length(train_idx), length(test_idx))
  rm(dtr_k, dte_k, m_k); invisible(gc(verbose = FALSE))
}
if (any(!is.finite(oof_pred))) {
  stop(sprintf("Non-finite OOF prediction(s) after baseline refit (n=%d)",
               sum(!is.finite(oof_pred))))
}

# ==============================================================================
# 6. SAVE -- same schema as the config OOF files (config_index = NA for baseline)
# ==============================================================================

meta_cols <- intersect(c("gene_id", "cell_type", "de_method",
                          "n_donors", "n_cells", "reads_per_cell",
                          "abs_gt_lfc", "gt_basemean", "gt_dispersion"),
                        colnames(dat))
oof_df <- data.frame(
  row_id       = seq_len(nrow(X)),
  fold         = row_fold,
  obs          = y,
  pred         = oof_pred,
  w_trials     = w_trials,
  alpha        = alpha_val,
  config_index = NA_integer_
)
oof_df <- cbind(oof_df, dat[, meta_cols, drop = FALSE])
qs::qsave(oof_df, file.path(xgb_model_base, "oof_predictions.qs"))
log_info("Saved baseline OOF predictions: %d rows -> %s",
         nrow(oof_df), file.path(xgb_model_base, "oof_predictions.qs"))

log_section("STEP 8b COMPLETE")
