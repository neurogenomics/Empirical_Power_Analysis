#!/usr/bin/env Rscript
# =============================================================================
# SMOKE TEST: Fix A -- manual OOF refit path in 02_Fit_XGBoost.R
# =============================================================================
# Purpose: verify that fit_one_config's new per-fold refit produces a valid
#          OOF prediction vector under xgboost 3.x (where xgb.cv()$pred is NULL).
#          Runs on a subsampled slice of the realistic .qs so it finishes in a
#          couple of minutes on a login node.
#
# Run (inside the container, from EPA repo root):
#   apptainer exec --bind $PWD:/workspace ./epa-cx3-2026-07-19.sif \
#       Rscript /workspace/hpc/smoke-oof-refit.R \
#           /workspace/Datalink/Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_realistic.qs
#
# Success criteria (all must pass):
#   [OK] length(oof_pred) == nrow(X)
#   [OK] no non-finite values in oof_pred
#   [OK] all oof_pred in [0, 1] (binary:logistic output range)
#   [OK] Pearson correlation with y > 0.3 (sanity: model actually learned something)
#   [OK] every fold covered (no fold has zero-length prediction span)
# =============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(xgboost)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: smoke-oof-refit.R <path-to-aggregated_fitting_data_realistic.qs> [n_genes_subsample=300]")
}
DATA_PATH  <- args[[1]]
N_GENES    <- as.integer(if (length(args) >= 2L) args[[2]] else "300")
SEED       <- 42L
MAX_ROUNDS <- 100L
EARLY_STOP <- 10L
NFOLDS     <- 5L
NTHREAD    <- as.integer(Sys.getenv("XGB_NTHREAD", "4"))

cat(sprintf("[smoke] xgboost version: %s\n", as.character(packageVersion("xgboost"))))
cat(sprintf("[smoke] R version:       %s\n", R.version.string))
cat(sprintf("[smoke] reading:         %s\n", DATA_PATH))
if (!file.exists(DATA_PATH)) stop("data file not found: ", DATA_PATH)
dat <- qs::qread(DATA_PATH)
cat(sprintf("[smoke] loaded rows:     %d  genes: %d\n",
            nrow(dat), dplyr::n_distinct(dat$gene_id)))

if ("gene_set" %in% colnames(dat)) dat <- dplyr::filter(dat, gene_set == "full")

dat <- dat %>% dplyr::filter(
  is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
  n_trials > 0,
  is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
  is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion)
)

# Subsample by gene so fold grouping stays legitimate.
set.seed(SEED)
gene_pool <- unique(dat$gene_id)
if (length(gene_pool) > N_GENES) gene_pool <- sample(gene_pool, N_GENES)
dat <- dat[dat$gene_id %in% gene_pool, , drop = FALSE]
cat(sprintf("[smoke] subsampled to:   %d rows  across %d genes\n",
            nrow(dat), length(gene_pool)))

# Design matrix (mirrors 02_Fit_XGBoost.R build_design).
CONT_FEATURES <- c("n_donors", "n_cells", "reads_per_cell",
                   "abs_gt_lfc", "gt_basemean", "gt_dispersion")
CAT_FEATURES  <- c("cell_type", "de_method")
cat_levels    <- lapply(CAT_FEATURES, function(f) sort(unique(as.character(dat[[f]]))))
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
cat(sprintf("[smoke] design matrix:   %d x %d\n", nrow(X), ncol(X)))

# Gene-grouped folds (Option A, matches 02).
gene_ids   <- as.character(dat$gene_id)
uniq_genes <- unique(gene_ids)
gene_fold  <- sample(rep_len(seq_len(NFOLDS), length(uniq_genes)))
names(gene_fold) <- uniq_genes
row_fold   <- gene_fold[gene_ids]
folds <- lapply(seq_len(NFOLDS), function(k) which(row_fold == k))
cat(sprintf("[smoke] fold sizes:      %s\n",
            paste(vapply(folds, length, integer(1)), collapse = ", ")))

# Monotone constraints (only for continuous columns; 0 for one-hot).
CONT_MONOTONE <- c(n_donors = 1L, n_cells = 1L, reads_per_cell = 1L,
                   abs_gt_lfc = 1L, gt_basemean = 1L, gt_dispersion = -1L)
monotone_vec <- integer(ncol(X))
for (i in seq_along(colnames(X))) {
  cn <- colnames(X)[i]
  if (cn %in% names(CONT_MONOTONE)) monotone_vec[i] <- CONT_MONOTONE[[cn]]
}
monotone_str <- paste0("(", paste(monotone_vec, collapse = ","), ")")

# Trivial "config" to exercise the code path fast.
params <- list(
  objective            = "binary:logistic",
  eval_metric          = "logloss",
  max_depth            = 3,
  eta                  = 0.1,
  min_child_weight     = 5,
  subsample            = 0.8,
  colsample_bytree     = 0.8,
  tree_method          = "hist",
  monotone_constraints = monotone_str,
  nthread              = NTHREAD
)

# Unit weights (no alpha reweighting) -- the OOF refit path is alpha-agnostic.
w <- rep(1, nrow(X))
dtrain <- xgb.DMatrix(data = X, label = y, weight = w)

cat("[smoke] running xgb.cv (early-stopped, no prediction=)... \n")
t0 <- Sys.time()
cv <- xgb.cv(params = params, data = dtrain,
             nrounds = MAX_ROUNDS, folds = folds,
             early_stopping_rounds = EARLY_STOP, verbose = 0)
metric_col <- grep("test_.*_mean", colnames(cv$evaluation_log), value = TRUE)[1]
best_iter <- cv$best_iteration
if (is.null(best_iter) || length(best_iter) != 1L || is.na(best_iter))
  best_iter <- attr(cv, "best_iteration")
if (is.null(best_iter) || length(best_iter) != 1L || is.na(best_iter))
  best_iter <- which.min(cv$evaluation_log[[metric_col]])
best_iter <- as.integer(best_iter)
cat(sprintf("[smoke] xgb.cv done: %.1fs, best_iter=%d, best_%s=%.5f\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            best_iter, metric_col, cv$evaluation_log[[metric_col]][best_iter]))

# --- Confirm the bug we're fixing: cv$pred should be NULL under xgboost 3.x --
cat(sprintf("[smoke] xgb.cv $pred is NULL?  %s  (length: %s)\n",
            is.null(cv$pred), length(cv$pred)))

# --- Fix A: manual per-fold refit ------------------------------------------
cat("[smoke] running manual OOF refit (Fix A)...\n")
t0 <- Sys.time()
oof_pred <- numeric(nrow(X))
for (k in seq_along(folds)) {
  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  dtr_k <- xgb.DMatrix(data   = X[train_idx, , drop = FALSE],
                       label  = y[train_idx],
                       weight = w[train_idx])
  dte_k <- xgb.DMatrix(data   = X[test_idx,  , drop = FALSE])
  m_k <- xgb.train(params = params, data = dtr_k,
                   nrounds = best_iter, verbose = 0)
  oof_pred[test_idx] <- predict(m_k, dte_k)
  rm(dtr_k, dte_k, m_k); invisible(gc(verbose = FALSE))
}
cat(sprintf("[smoke] refit done: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- Assertions -----------------------------------------------------------
ok <- TRUE
report <- function(msg, cond) {
  cat(sprintf("[%s] %s\n", if (cond) "OK" else "FAIL", msg))
  if (!cond) ok <<- FALSE
}
report(sprintf("length(oof_pred) == nrow(X)  (%d == %d)", length(oof_pred), nrow(X)),
       length(oof_pred) == nrow(X))
report(sprintf("no NA / NaN / Inf  (found %d)", sum(!is.finite(oof_pred))),
       all(is.finite(oof_pred)))
report(sprintf("range in [0,1]  [%.4f, %.4f]", min(oof_pred), max(oof_pred)),
       min(oof_pred) >= 0 & max(oof_pred) <= 1)
cor_yhat <- suppressWarnings(cor(oof_pred, y))
report(sprintf("Pearson cor(oof_pred, y) > 0.30  (got %.3f)", cor_yhat),
       is.finite(cor_yhat) && cor_yhat > 0.30)
covered <- sapply(folds, function(idx) any(oof_pred[idx] != 0))
report(sprintf("every fold covered (%d/%d)", sum(covered), length(folds)),
       all(covered))

cat("\n[smoke] mean OOF pred by observed-power bin:\n")
bin <- cut(y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
print(round(tapply(oof_pred, bin, mean), 3))

if (!ok) {
  cat("\n[smoke] FAIL -- do NOT resubmit the full array.\n")
  quit(status = 1)
}
cat(sprintf("\n[smoke] PASS -- Fix A OOF refit produces a valid length-%d prediction vector.\n",
            length(oof_pred)))
cat("[smoke] Safe to resubmit the 36-config search array.\n")
quit(status = 0)
