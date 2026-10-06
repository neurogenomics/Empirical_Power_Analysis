# ==============================================================================
# SCRIPT: Select the XGBoost power model + SHAP refits (M2 -- HPC Stage C)
# PURPOSE: After the search array (02_Fit_XGBoost.R in ARRAY mode) has trained
#          one booster per {tree hyperparameters x alpha} configuration into
#          Models/xgboost/search/config_<k>/, this script:
#            1. scores every configuration by out-of-fold bin-balanced MAE
#               (overall and per CV fold), writes calibration_residuals.csv per
#               configuration, and ranks them (model_selection.csv)
#            2. picks 5 models for SHAP (best per alpha, overall median,
#               overall worst), refits each on the full training table and
#               computes TreeSHAP on a uniform row subsample
#            3. promotes the overall winner to Models/xgboost/ (the paper model)
#
# INPUT:  Fitting_Data/aggregated_fitting_data.qs                   (from 01e)
#         Models/xgboost/search/config_*/{cv_metrics.csv,oof_predictions.qs,feature_spec.qs}
# OUTPUT: Models/xgboost/model_selection.csv        all configs ranked (Supplementary Data 4)
#         Models/xgboost/model_selection_folds.csv  bin-balanced MAE per config x CV fold
#         Models/xgboost/shap_selection.csv         the 5 SHAP models
#         Models/xgboost/shap_models/<tag>/...      refit boosters + SHAP
#         Models/xgboost/{booster.ubj,booster.qs,feature_spec.qs,shap_values.qs,
#                         shap_summary.csv,feature_importance.csv,
#                         calibration_residuals.csv,cv_metrics.csv}   (winner)
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

log_section("STEP 3 (XGBoost): Select best config + SHAP/interpretation")

set.seed(SEED_BASE)

XGB_NFOLDS   <- as.integer(Sys.getenv("XGB_NFOLDS", "5"))
XGB_NTHREAD  <- as.integer(Sys.getenv("XGB_NTHREAD", as.character(min(8, N_CORES))))
XGB_GENE_SET <- Sys.getenv("XGB_GENE_SET", "full")

# Must match those used in 02_Fit_XGBoost.R.
XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "")
XGB_FITTING_DATA <- Sys.getenv("XGB_FITTING_DATA", "aggregated_fitting_data.qs")

model_base <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Models", "xgboost")
}
search_dir <- file.path(model_base, "search")
if (!dir.exists(search_dir)) {
  stop(sprintf("No search dir: %s. Run the Stage B array (02_Fit_XGBoost.R, ARRAY mode) first.", search_dir))
}
log_info("Model base: %s | fitting data: %s", model_base, XGB_FITTING_DATA)

# --------------------------------------------------------------------------
# Bin-balanced MAE helper (design doc §10.9). Mean over predicted-power bins
# of the trial-weighted MAE within each bin; bins with n=0 are excluded from
# the mean. Bins are 10 equal-width intervals over [0, 1].
# --------------------------------------------------------------------------
compute_bin_balanced <- function(pred, obs, w) {
  ok <- is.finite(pred) & is.finite(obs) & is.finite(w) & w > 0
  pred <- pred[ok]; obs <- obs[ok]; w <- w[ok]
  bin <- cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  err <- abs(obs - pred)
  per_bin <- tapply(seq_along(err), bin, function(ix) {
    if (length(ix) == 0) return(c(mae = NA_real_, p90 = NA_real_))
    ww <- w[ix]
    ee <- err[ix]
    mae <- sum(ww * ee) / sum(ww)
    o   <- order(ee); ee_o <- ee[o]; ww_o <- ww[o]
    cw  <- (cumsum(ww_o) - 0.5 * ww_o) / sum(ww_o)
    p90 <- approx(cw, ee_o, xout = 0.90, rule = 2)$y
    c(mae = mae, p90 = p90)
  })
  mat <- do.call(rbind, per_bin)
  list(
    bin_balanced_MAE      = mean(mat[, "mae"], na.rm = TRUE),
    bin_balanced_p90_err  = mean(mat[, "p90"], na.rm = TRUE),
    per_bin               = mat
  )
}

# --------------------------------------------------------------------------
# Weighted summary helpers (used to build one calibration_residuals.csv per
# configuration).
# --------------------------------------------------------------------------
wsd_h <- function(x, wt) {
  m <- weighted.mean(x, wt); sqrt(sum(wt * (x - m)^2) / sum(wt))
}
wquantile_h <- function(x, wt, probs) {
  o <- order(x); x <- x[o]; wt <- wt[o]
  cw <- (cumsum(wt) - 0.5 * wt) / sum(wt)
  approx(cw, x, xout = probs, rule = 2)$y
}
build_calibration_residuals <- function(pred, obs, w) {
  ok <- is.finite(pred) & is.finite(obs) & is.finite(w) & w > 0
  d <- data.frame(pred = pred[ok], obs = obs[ok], w = w[ok])
  d <- transform(
    d,
    err          = obs - pred,
    abs_err      = abs(obs - pred),
    pred_bin     = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE),
    obs_bin      = cut(obs,  breaks = seq(0, 1, 0.1), include.lowest = TRUE),
    pred_bin_idx = as.integer(cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)),
    obs_bin_idx  = as.integer(cut(obs,  breaks = seq(0, 1, 0.1), include.lowest = TRUE))
  )
  resid_row <- function(dd) {
    sw <- sum(dd$w)
    data.frame(
      n                   = nrow(dd),
      sum_w               = sw,
      mean_pred           = weighted.mean(dd$pred, dd$w),
      mean_obs            = weighted.mean(dd$obs, dd$w),
      bias                = weighted.mean(dd$obs, dd$w) - weighted.mean(dd$pred, dd$w),
      sd_error            = wsd_h(dd$err, dd$w),
      min_error           = min(dd$err),
      max_error           = max(dd$err),
      mae                 = weighted.mean(dd$abs_err, dd$w),
      p50_abs_error       = wquantile_h(dd$abs_err, dd$w, 0.50),
      p90_abs_error       = wquantile_h(dd$abs_err, dd$w, 0.90),
      p95_abs_error       = wquantile_h(dd$abs_err, dd$w, 0.95),
      pct_within_005      = sum(dd$w[dd$abs_err <= 0.05]) / sw,
      pct_within_010      = sum(dd$w[dd$abs_err <= 0.10]) / sw,
      pct_same_bin        = sum(dd$w[dd$obs_bin_idx == dd$pred_bin_idx]) / sw,
      pct_same_or_adj_bin = sum(dd$w[abs(dd$obs_bin_idx - dd$pred_bin_idx) <= 1]) / sw
    )
  }
  by_bin <- d %>%
    group_by(pred_bin) %>%
    group_modify(~ resid_row(.x)) %>%
    ungroup() %>%
    mutate(pred_bin = as.character(pred_bin))
  overall <- resid_row(d); overall$pred_bin <- "ALL"
  bind_rows(by_bin, overall) %>% select(pred_bin, everything())
}

# ==============================================================================
# 2. RANK, SELECT, REFIT, PROMOTE (design doc §6.3)
# ==============================================================================
# Rank by bin-balanced MAE computed from the per-config OOF predictions (saved
# in ARRAY mode by 02), select 5 models for SHAP (3 alpha-winners + median +
# worst by bin_balanced), refit each on the full training table and compute
# SHAP, then promote the overall winner to model_base.

{
  log_section("Rank by bin_balanced_MAE, 5-model SHAP refit")

  cfg_dirs <- list.files(search_dir, pattern = "^config_[0-9]+$", full.names = TRUE)
  cfg_dirs <- cfg_dirs[order(as.integer(sub(".*config_", "", cfg_dirs)))]
  if (length(cfg_dirs) == 0) stop(sprintf("No config_* dirs under %s", search_dir))

  # --- per-config metrics + per-config calibration_residuals.csv ---------------
  fold_rows <- list()
  per_cfg_rows <- lapply(cfg_dirs, function(d) {
    cv_f  <- file.path(d, "cv_metrics.csv")
    oof_f <- file.path(d, "oof_predictions.qs")
    if (!file.exists(cv_f) || !file.exists(oof_f)) {
      log_warning("Missing metrics/OOF in %s -- skipping", d); return(NULL)
    }
    cv  <- read.csv(cv_f, stringsAsFactors = FALSE)
    oof <- qs::qread(oof_f)
    bb  <- compute_bin_balanced(oof$pred, oof$obs, oof$w_trials)
    cr  <- build_calibration_residuals(oof$pred, oof$obs, oof$w_trials)
    write.csv(cr, file.path(d, "calibration_residuals.csv"), row.names = FALSE)
    # Per-fold bin-balanced MAE: the folds are gene-grouped and identical for
    # every configuration (same seed and gene list), so per-fold scores can be
    # compared pairwise to ask whether a gap between configurations is larger
    # than fold-to-fold noise.
    fold_mae <- vapply(sort(unique(oof$fold)), function(k) {
      ix <- oof$fold == k
      compute_bin_balanced(oof$pred[ix], oof$obs[ix], oof$w_trials[ix])$bin_balanced_MAE
    }, numeric(1))
    fold_rows[[length(fold_rows) + 1]] <<- data.frame(
      config_index = cv$config_index[1], fold = sort(unique(oof$fold)),
      bin_balanced_MAE = fold_mae)
    data.frame(
      config_index         = cv$config_index[1],
      config_dir           = d,
      alpha                = if ("alpha" %in% colnames(cv)) cv$alpha[1] else NA_real_,
      max_depth            = cv$max_depth[1],
      eta                  = cv$eta[1],
      min_child_weight     = cv$min_child_weight[1],
      subsample            = cv$subsample[1],
      colsample_bytree     = cv$colsample_bytree[1],
      best_nrounds         = cv$best_nrounds[1],
      cv_logloss           = cv$cv_logloss[1],
      bin_balanced_MAE     = bb$bin_balanced_MAE,
      bin_balanced_p90_err = bb$bin_balanced_p90_err,
      bin_balanced_MAE_fold_sd = sd(fold_mae),
      stringsAsFactors     = FALSE
    )
  })
  per_cfg <- bind_rows(per_cfg_rows)
  if (nrow(per_cfg) == 0) stop("No configs with OOF predictions -- did Stage B complete?")
  per_cfg <- per_cfg %>% arrange(bin_balanced_MAE, bin_balanced_p90_err)
  log_info("Ranked %d configs by bin_balanced_MAE (best=%.5f, worst=%.5f)",
           nrow(per_cfg), per_cfg$bin_balanced_MAE[1],
           per_cfg$bin_balanced_MAE[nrow(per_cfg)])
  write.csv(per_cfg %>% select(-config_dir),
            file.path(model_base, "model_selection.csv"), row.names = FALSE)
  write.csv(bind_rows(fold_rows) %>% arrange(config_index, fold),
            file.path(model_base, "model_selection_folds.csv"), row.names = FALSE)

  # --- select 5 SHAP models -----------------------------------------------------
  # 1) best per alpha (3 rows if 3 alphas)
  alpha_winners <- per_cfg %>%
    group_by(alpha) %>%
    slice_min(order_by = bin_balanced_MAE, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(selection_reason = sprintf("winner_alpha_%g", alpha),
           tag              = sprintf("winner_alpha_%g", alpha))
  chosen <- alpha_winners
  # 2) overall median (by bin_balanced_MAE); skip if it duplicates a winner
  ranked  <- per_cfg %>% arrange(bin_balanced_MAE)
  mid_ix  <- ceiling(nrow(ranked) / 2)
  median_row <- ranked[mid_ix, ] %>%
    mutate(selection_reason = "overall_median", tag = "median")
  if (!median_row$config_index %in% chosen$config_index) {
    chosen <- bind_rows(chosen, median_row)
  } else {
    # walk outward until we find a non-overlapping row
    for (off in seq_len(nrow(ranked))) {
      for (cand_ix in unique(c(mid_ix + off, mid_ix - off))) {
        if (cand_ix < 1 || cand_ix > nrow(ranked)) next
        cand <- ranked[cand_ix, ]
        if (!cand$config_index %in% chosen$config_index) {
          chosen <- bind_rows(chosen,
                              cand %>% mutate(selection_reason = "overall_median_shifted",
                                              tag = "median"))
          off <- Inf; break
        }
      }
      if (is.infinite(off)) break
    }
  }
  # 3) overall worst; skip if it duplicates a winner
  worst_row <- ranked[nrow(ranked), ] %>%
    mutate(selection_reason = "overall_worst", tag = "worst")
  if (!worst_row$config_index %in% chosen$config_index) {
    chosen <- bind_rows(chosen, worst_row)
  } else {
    for (i in seq(nrow(ranked) - 1, 1)) {
      cand <- ranked[i, ]
      if (!cand$config_index %in% chosen$config_index) {
        chosen <- bind_rows(chosen,
                            cand %>% mutate(selection_reason = "overall_worst_shifted",
                                            tag = "worst"))
        break
      }
    }
  }
  chosen <- chosen[!duplicated(chosen$config_index), ]
  write.csv(chosen %>% select(-config_dir),
            file.path(model_base, "shap_selection.csv"), row.names = FALSE)
  log_info("Selected %d SHAP models: %s",
           nrow(chosen), paste(chosen$tag, collapse = ", "))

  # --- rebuild design matrix on FILTERED data ----------------------------------
  agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", XGB_FITTING_DATA)
  if (!file.exists(agg_path)) stop(sprintf("Fitting data not found: %s", agg_path))
  dat <- qs::qread(agg_path)
  if ("gene_set" %in% colnames(dat)) dat <- dat %>% filter(gene_set == XGB_GENE_SET)
  dat <- dat %>%
    filter(is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
           n_trials > 0,
           is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
           is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion))
  log_info("Rebuilt filtered fitting set: %d rows, %d genes",
           nrow(dat), dplyr::n_distinct(dat$gene_id))

  # Feature spec pulled from any one config (they share the design).
  ref_spec <- qs::qread(file.path(chosen$config_dir[1], "feature_spec.qs"))
  CONT_FEATURES <- ref_spec$cont_features
  CAT_FEATURES  <- ref_spec$cat_features
  cat_levels    <- ref_spec$cat_levels
  feature_names <- ref_spec$feature_names

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
  X <- X[, feature_names, drop = FALSE]
  y <- dat$detection_rate
  w_trials <- dat$n_trials

  # Bin-weight base for reweighting during refit (see 02_Fit_XGBoost.R §4).
  obs_bin_r         <- cut(y, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  bin_prop_r        <- prop.table(table(obs_bin_r))
  bin_weight_base_r <- as.numeric(1 / bin_prop_r[as.character(obs_bin_r)])
  make_w <- function(alpha) {
    if (!is.finite(alpha) || alpha == 0) return(w_trials)
    bw <- bin_weight_base_r ^ alpha
    bw <- bw * length(bw) / sum(bw)
    w_trials * bw
  }

  # --- refit each chosen config on full filtered data + SHAP -------------------
  shap_models_dir <- file.path(model_base, "shap_models")
  safe_mkdir(shap_models_dir)
  SHAP_MAX_ROWS <- as.integer(Sys.getenv("XGB_SHAP_MAX_ROWS", "200000"))
  force_shap    <- tolower(Sys.getenv("XGB_FORCE_SHAP", "false")) %in% c("true", "1", "yes")

  for (i in seq_len(nrow(chosen))) {
    row <- chosen[i, ]
    tag_dir <- file.path(shap_models_dir, row$tag)
    safe_mkdir(tag_dir)

    out_ubj  <- file.path(tag_dir, "booster.ubj")
    out_shap <- file.path(tag_dir, "shap_values.qs")
    if (!force_shap && file.exists(out_ubj) && file.exists(out_shap)) {
      log_info("[%s] SHAP artefacts present, skipping refit (XGB_FORCE_SHAP=true to force)", row$tag)
      next
    }

    log_info("[%s] refit: config %d (depth=%d eta=%g mcw=%g alpha=%g), %d rounds",
             row$tag, row$config_index, row$max_depth, row$eta,
             row$min_child_weight, row$alpha, row$best_nrounds)

    monotone_str <- paste0("(", paste(ref_spec$monotone, collapse = ","), ")")
    params <- list(
      objective            = ref_spec$objective,
      eval_metric          = "logloss",
      max_depth            = as.integer(row$max_depth),
      eta                  = as.numeric(row$eta),
      min_child_weight     = as.numeric(row$min_child_weight),
      subsample            = as.numeric(row$subsample),
      colsample_bytree     = as.numeric(row$colsample_bytree),
      tree_method          = "hist",
      monotone_constraints = monotone_str,
      nthread              = XGB_NTHREAD
    )
    w_cfg  <- make_w(row$alpha)
    dtrain <- xgb.DMatrix(data = X, label = y, weight = w_cfg)
    bst    <- xgb.train(params = params, data = dtrain,
                        nrounds = as.integer(row$best_nrounds), verbose = 0)

    xgb.save(bst, out_ubj)
    qs::qsave(bst, file.path(tag_dir, "booster.qs"))
    spec_i <- ref_spec
    spec_i$best_params  <- params
    spec_i$best_nrounds <- as.integer(row$best_nrounds)
    spec_i$alpha        <- as.numeric(row$alpha)
    spec_i$fitting_data <- XGB_FITTING_DATA
    spec_i$model_subdir <- XGB_MODEL_SUBDIR
    spec_i$model_tag    <- row$tag
    qs::qsave(spec_i, file.path(tag_dir, "feature_spec.qs"))

    # SHAP on a uniform subsample so exact TreeSHAP stays tractable
    if (nrow(X) > SHAP_MAX_ROWS) {
      set.seed(SEED_BASE)
      shap_idx <- sort(sample.int(nrow(X), SHAP_MAX_ROWS))
    } else {
      shap_idx <- seq_len(nrow(X))
    }
    log_info("[%s] SHAP on %d rows (exact TreeSHAP)", row$tag, length(shap_idx))
    shap <- predict(bst, xgb.DMatrix(X[shap_idx, , drop = FALSE]), predcontrib = TRUE)
    colnames(shap) <- c(feature_names, "BIAS")
    qs::qsave(shap, out_shap)

    shap_feat <- shap[, feature_names, drop = FALSE]
    shap_summary <- data.frame(
      feature       = feature_names,
      mean_abs_shap = colMeans(abs(shap_feat)),
      mean_shap     = colMeans(shap_feat),
      row.names = NULL
    ) %>% arrange(desc(mean_abs_shap))
    write.csv(shap_summary, file.path(tag_dir, "shap_summary.csv"), row.names = FALSE)

    imp <- xgb.importance(feature_names = feature_names, model = bst)
    write.csv(imp, file.path(tag_dir, "feature_importance.csv"), row.names = FALSE)
  }

  # --- promote overall winner to canonical model_base --------------------------
  win <- per_cfg[1, ]   # already sorted ascending by bin_balanced_MAE
  log_info("Promoting overall winner (config %d, alpha=%g, bin_balanced_MAE=%.5f)",
           win$config_index, win$alpha, win$bin_balanced_MAE)

  win_tag <- (chosen %>% filter(config_index == win$config_index) %>% pull(tag))[1]
  if (is.na(win_tag)) {
    # Winner wasn't in the 5 SHAP models (shouldn't happen with our selection),
    # so materialise it now.
    win_tag <- "overall_winner"
    warning("Overall winner not among 5 SHAP models -- refitting to promote.")
  }
  win_dir <- file.path(shap_models_dir, win_tag)
  file.copy(file.path(win_dir, "booster.ubj"),
            file.path(model_base, "booster.ubj"), overwrite = TRUE)
  file.copy(file.path(win_dir, "booster.qs"),
            file.path(model_base, "booster.qs"), overwrite = TRUE)
  file.copy(file.path(win_dir, "feature_spec.qs"),
            file.path(model_base, "feature_spec.qs"), overwrite = TRUE)
  file.copy(file.path(win_dir, "shap_values.qs"),
            file.path(model_base, "shap_values.qs"), overwrite = TRUE)
  file.copy(file.path(win_dir, "shap_summary.csv"),
            file.path(model_base, "shap_summary.csv"), overwrite = TRUE)
  file.copy(file.path(win_dir, "feature_importance.csv"),
            file.path(model_base, "feature_importance.csv"), overwrite = TRUE)
  # Winner's calibration & residuals (from its own config_dir OOF preds).
  file.copy(file.path(win$config_dir, "calibration_residuals.csv"),
            file.path(model_base, "calibration_residuals.csv"), overwrite = TRUE)
  # cv_metrics for the winner + trimmed model_selection for the paper.
  write.csv(win %>% select(-config_dir),
            file.path(model_base, "cv_metrics.csv"), row.names = FALSE)

  log_section("STEP 3 (XGBoost) COMPLETE")
  log_info("Winner promoted to %s (bin_balanced_MAE=%.4f)",
           model_base, win$bin_balanced_MAE)
  quit(save = "no", status = 0)
}
