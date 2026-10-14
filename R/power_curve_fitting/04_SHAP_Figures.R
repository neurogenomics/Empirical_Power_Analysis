# ==============================================================================
# SCRIPT: SHAP figures for the XGBoost power model (M2 -- HPC Stage D)
# PURPOSE: Turn the Stage B (search) + Stage C (select) artifacts into paper
#          figures, no model re-fitting. For each of the 5 SHAP models chosen
#          by 03_Select_XGBoost.R (best per alpha, median, worst) it draws a
#          beeswarm, violin and per-feature dependence panel into
#          shap_figures/<tag>/, plus one cross-model mean|SHAP| heatmap.
#
# INPUT:  Fitting_Data/aggregated_fitting_data.qs                     (01e)
#         Models/xgboost/{booster.ubj,feature_spec.qs,shap_values.qs, (Stage C)
#                         shap_summary.csv,model_selection.csv}
#         Models/xgboost/{shap_selection.csv,shap_models/<tag>/{booster.ubj,shap_values.qs}}  (Stage C)
# OUTPUT: Models/xgboost/shap_figures/*.png
#
# ENV KNOBS (all optional):
#   XGB_GENE_SET       gene set filter               (default "full")
#   XGB_SHAP_MAX_ROWS  MUST match Stage C to reuse   (default 200000)
#                      the saved shap_values.qs subsample; if it doesn't line
#                      up the script recomputes TreeSHAP on the winner.
#   XGB_BEESWARM_ROWS  points drawn per beeswarm     (default 30000)
#   XGB_NTHREAD        predict threads               (default min(8,N_CORES))
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(xgboost)
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

log_section("STEP 4 (XGBoost): SHAP figures (winner beeswarm + per-config overlay)")

set.seed(SEED_BASE)

XGB_GENE_SET     <- Sys.getenv("XGB_GENE_SET", "full")
XGB_NTHREAD      <- as.integer(Sys.getenv("XGB_NTHREAD", as.character(min(8, N_CORES))))
SHAP_MAX_ROWS    <- as.integer(Sys.getenv("XGB_SHAP_MAX_ROWS", "200000"))
BEESWARM_ROWS    <- as.integer(Sys.getenv("XGB_BEESWARM_ROWS", "30000"))

# Must match 02/03.
XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "")
XGB_FITTING_DATA <- Sys.getenv("XGB_FITTING_DATA", "aggregated_fitting_data.qs")

model_base <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Models", "xgboost")
}
search_dir <- file.path(model_base, "search")
fig_dir    <- file.path(model_base, "shap_figures")
safe_mkdir(fig_dir)

if (!dir.exists(model_base))
  stop(sprintf("No model dir: %s. Run Stage B + Stage C first.", model_base))

# Human-readable feature labels (native units; no log space after the raw-feature
# fix). Includes the one-hot categorical columns so the beeswarm axis reads e.g.
# "DESeq2" / "Oligodendrocytes" instead of the raw "de_method=deseq2".
FEATURE_LABELS <- c(
  n_donors            = "Donors",
  n_cells             = "Cells / donor",
  reads_per_cell      = "Counts / cell",
  abs_gt_lfc          = "|LFC|",
  gt_basemean         = "Median expression",
  gt_dispersion       = "Dispersion",
  `de_method=deseq2`  = "DESeq2",
  `de_method=edger`   = "edgeR",
  `de_method=limma`   = "limma-voom",
  `cell_type=Astro`   = "Astrocytes",
  `cell_type=Micro`   = "Microglia",
  `cell_type=Oligo`   = "Oligodendrocytes"
)
pretty_feature <- function(f) ifelse(f %in% names(FEATURE_LABELS), FEATURE_LABELS[f], f)

# Shared paper-ready theme: clean serif-neutral look, no minor grid, framed
# panel, legend on the right. Reused across every figure in this script so the
# whole SHAP figure set is visually consistent for the manuscript.
theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = "grey30", size = rel(0.82),
                                      margin = margin(b = 6)),
      plot.title.position = "plot",
      axis.title       = element_text(face = "bold", size = rel(0.9)),
      axis.text        = element_text(colour = "grey20"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      legend.title     = element_text(face = "bold", size = rel(0.8)),
      legend.text      = element_text(size = rel(0.75)),
      legend.key.size  = unit(0.9, "lines"),
      plot.margin      = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# 2. REBUILD DESIGN MATRIX (identical recipe to 02/03 -- for feature VALUES)
#    shap_values.qs stores only contributions, so we rebuild X to recover the
#    feature values needed to colour the beeswarm and drive the dependence plots.
# ==============================================================================

winner_spec <- qs::qread(file.path(model_base, "feature_spec.qs"))
CONT_FEATURES <- winner_spec$cont_features
CAT_FEATURES  <- winner_spec$cat_features
cat_levels    <- winner_spec$cat_levels
feature_names <- winner_spec$feature_names

agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", XGB_FITTING_DATA)
if (!file.exists(agg_path)) stop(sprintf("Fitting data not found: %s", agg_path))
dat <- qs::qread(agg_path)
if ("gene_set" %in% colnames(dat)) dat <- dat %>% filter(gene_set == XGB_GENE_SET)

dat <- dat %>%
  filter(
    is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
    n_trials > 0,
    is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
    is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion)
  )

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
log_info("Design matrix rebuilt: %d rows x %d features", nrow(X), ncol(X))

# ==============================================================================
# 3. FIGURES (design doc §6.4) — loop over the 5 SHAP models
# ==============================================================================
# 03_Select_XGBoost.R has already refit 5 SHAP models (3 alpha-winners + median + worst) under
# model_base/shap_models/<tag>/ and produced shap_values.qs alongside each
# booster. This branch loops over shap_selection.csv, produces the beeswarm /
# violin / size x alpha / per-feature dependence figures for each tag into
# shap_figures/<tag>/, and one cross-model mean|SHAP| heatmap at the root.

{
  sel_path <- file.path(model_base, "shap_selection.csv")
  if (!file.exists(sel_path))
    stop(sprintf("shap_selection.csv missing: %s -- run 03_Select_XGBoost.R first", sel_path))
  selection <- read.csv(sel_path, stringsAsFactors = FALSE)
  log_info("Rendering figures for %d SHAP models (%s)",
           nrow(selection), paste(selection$tag, collapse = ", "))

  # --- reusable helpers ------------------------------------------------------
  scale01 <- function(v) {
    rng <- range(v, finite = TRUE)
    if (!is.finite(rng[1]) || diff(rng) == 0) return(rep(0.5, length(v)))
    (v - rng[1]) / diff(rng)
  }
  build_bee <- function(vtrans, feat_order, shap_feat, Xs, draw_idx) {
    frames <- lapply(feat_order, function(f) {
      data.frame(feature = f,
                 shap    = shap_feat[draw_idx, f],
                 fval    = scale01(vtrans(Xs[draw_idx, f])),
                 stringsAsFactors = FALSE)
    })
    out <- bind_rows(frames)
    out$feature <- factor(out$feature, levels = rev(feat_order),
                          labels = rev(pretty_feature(feat_order)))
    out
  }
  BEE_DOT_SIZE  <- 0.28; BEE_DOT_ALPHA <- 0.32
  make_beeswarm <- function(bee, colour_name,
                            dot_size = BEE_DOT_SIZE, dot_alpha = BEE_DOT_ALPHA) {
    ggplot(bee, aes(x = feature, y = shap, colour = fval)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_point(size = dot_size, alpha = dot_alpha, stroke = 0,
                 position = position_jitter(width = 0.4, height = 0, seed = SEED_BASE)) +
      coord_flip() +
      scale_colour_viridis_c(name = colour_name, limits = c(0, 1), breaks = c(0, 1),
                             labels = c("Low", "High"),
                             guide = guide_colourbar(barwidth = unit(0.5, "lines"),
                                                     barheight = unit(7, "lines"),
                                                     ticks = FALSE, title.position = "top")) +
      labs(x = NULL, y = "SHAP value (impact on log-odds of detection)") +
      theme_paper(base_size = 12) +
      theme(panel.grid.major.y = element_blank())
  }
  make_shap_violin <- function(bee, colour_name,
                               dot_size = 0.20, dot_alpha = 0.28) {
    ggplot(bee, aes(x = feature, y = shap)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
      geom_violin(aes(fill = feature), colour = "grey40", linewidth = 0.25,
                  alpha = 0.18, scale = "width", trim = TRUE, show.legend = FALSE) +
      scale_fill_manual(values = rep("grey70", length(levels(bee$feature)))) +
      geom_point(aes(colour = fval), size = dot_size, alpha = dot_alpha, stroke = 0,
                 position = position_jitter(width = 0.28, height = 0, seed = SEED_BASE)) +
      coord_flip() +
      scale_colour_viridis_c(name = colour_name, limits = c(0, 1), breaks = c(0, 1),
                             labels = c("Low", "High"),
                             guide = guide_colourbar(barwidth = unit(0.5, "lines"),
                                                     barheight = unit(7, "lines"),
                                                     ticks = FALSE, title.position = "top")) +
      labs(x = NULL, y = "SHAP value (impact on log-odds of detection)") +
      theme_paper(base_size = 12) +
      theme(panel.grid.major.y = element_blank())
  }
  N_BINS_DEP <- 40

  # --- per-tag figures ---------------------------------------------------------
  cross_imp <- list()
  for (i in seq_len(nrow(selection))) {
    row <- selection[i, ]
    tag <- row$tag
    tag_src_dir <- file.path(model_base, "shap_models", tag)
    tag_fig_dir <- file.path(fig_dir, tag)
    safe_mkdir(tag_fig_dir)

    # Self-skip guard scoped per tag (design doc §6.4).
    tag_marker <- file.path(tag_fig_dir, "shap_beeswarm.png")
    if (file.exists(tag_marker) && !identical(tolower(Sys.getenv("XGB_FORCE_FIG", "false")),
                                              "true")) {
      log_info("[%s] figures present, skipping (set XGB_FORCE_FIG=true to redo)", tag)
      # Still read SHAP for cross-model heatmap
      sh <- qs::qread(file.path(tag_src_dir, "shap_values.qs"))
      cross_imp[[tag]] <- colMeans(abs(sh[, feature_names, drop = FALSE]))
      next
    }

    booster_path <- file.path(tag_src_dir, "booster.ubj")
    shap_path    <- file.path(tag_src_dir, "shap_values.qs")
    if (!file.exists(booster_path) || !file.exists(shap_path)) {
      log_warning("[%s] missing booster or shap_values, skipping", tag); next
    }

    # Reproduce 03's deterministic subsample so shap_values line up with X.
    if (nrow(X) > SHAP_MAX_ROWS) {
      set.seed(SEED_BASE)
      shap_idx <- sort(sample.int(nrow(X), SHAP_MAX_ROWS))
    } else {
      shap_idx <- seq_len(nrow(X))
    }
    shap      <- qs::qread(shap_path)
    if (nrow(shap) != length(shap_idx)) {
      log_warning("[%s] shap_values rows (%d) != subsample (%d); recomputing",
                  tag, nrow(shap), length(shap_idx))
      bst  <- xgb.load(booster_path)
      shap <- predict(bst, xgb.DMatrix(X[shap_idx, , drop = FALSE]), predcontrib = TRUE)
      colnames(shap) <- c(feature_names, "BIAS")
    }
    shap_feat <- shap[, feature_names, drop = FALSE]
    Xs        <- X[shap_idx, feature_names, drop = FALSE]

    feat_order <- names(sort(colMeans(abs(shap_feat)), decreasing = TRUE))
    cross_imp[[tag]] <- colMeans(abs(shap_feat))

    n_draw <- min(BEESWARM_ROWS, nrow(shap_feat))
    set.seed(SEED_BASE)
    draw_idx <- sort(sample.int(nrow(shap_feat), n_draw))

    bee_raw <- build_bee(identity, feat_order, shap_feat, Xs, draw_idx)
    bee_log <- build_bee(log1p,    feat_order, shap_feat, Xs, draw_idx)

    ggsave(file.path(tag_fig_dir, "shap_beeswarm.png"),
           make_beeswarm(bee_raw, "Feature value\n(raw)"),
           width = 9, height = 5.5, dpi = PLOT_DPI)
    ggsave(file.path(tag_fig_dir, "shap_beeswarm_log.png"),
           make_beeswarm(bee_log, "Feature value\n(log scale)"),
           width = 9, height = 5.5, dpi = PLOT_DPI)
    ggsave(file.path(tag_fig_dir, "shap_violin.png"),
           make_shap_violin(bee_raw, "Feature value\n(raw)"),
           width = 9, height = 5.5, dpi = PLOT_DPI)
    ggsave(file.path(tag_fig_dir, "shap_violin_log.png"),
           make_shap_violin(bee_log, "Feature value\n(log scale)"),
           width = 9, height = 5.5, dpi = PLOT_DPI)

    grid_dir <- file.path(tag_fig_dir, "beeswarm_grid")
    safe_mkdir(grid_dir)
    for (sz in c(0.18, 0.25, 0.35, 0.50)) {
      for (al in c(0.20, 0.32, 0.45, 0.65)) {
        ggsave(file.path(grid_dir,
                         sprintf("shap_beeswarm_log_size%0.2f_alpha%0.2f.png", sz, al)),
               make_beeswarm(bee_log, "Feature value\n(log scale)",
                             dot_size = sz, dot_alpha = al),
               width = 9, height = 5.5, dpi = PLOT_DPI)
      }
    }

    # Per-feature SHAP dependence: bin feature value into quantiles, mean SHAP
    # per bin. One panel per continuous feature; single line (no cross-config
    # overlay in per-tag figures -- cross-model comparison lives in 08).
    cont_present <- intersect(CONT_FEATURES, feature_names)
    dep_panels <- lapply(cont_present, function(f) {
      v <- Xs[, f]
      qs_break <- unique(quantile(v, probs = seq(0, 1, length.out = N_BINS_DEP + 1),
                                  na.rm = TRUE, type = 8))
      if (length(qs_break) < 3) return(NULL)
      bin     <- cut(v, breaks = qs_break, include.lowest = TRUE)
      bin_mid <- tapply(v, bin, median)
      ms      <- tapply(shap_feat[, f], bin, mean)
      df <- data.frame(x = as.numeric(bin_mid[names(ms)]),
                        shap = as.numeric(ms))
      ggplot(df, aes(x = x, y = shap)) +
        geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
        geom_line(linewidth = 0.6, colour = "steelblue") +
        labs(x = pretty_feature(f), y = "Mean SHAP (log-odds)") +
        theme_paper(base_size = 11)
    })
    dep_panels <- Filter(Negate(is.null), dep_panels)
    if (length(dep_panels) > 0) {
      combined <- wrap_plots(dep_panels, ncol = 2)
      ggsave(file.path(tag_fig_dir, "shap_dependence.png"), combined,
             width = 10, height = 3.0 * ceiling(length(dep_panels) / 2),
             dpi = PLOT_DPI)
    }
    log_info("[%s] wrote figures to %s", tag, tag_fig_dir)
  }

  # --- cross-model importance heatmap -----------------------------------------
  if (length(cross_imp) > 0) {
    cm <- do.call(rbind, lapply(names(cross_imp), function(tg) {
      data.frame(tag = tg, feature = names(cross_imp[[tg]]),
                 mean_abs_shap = as.numeric(cross_imp[[tg]]),
                 stringsAsFactors = FALSE)
    }))
    feat_levels <- names(sort(tapply(cm$mean_abs_shap, cm$feature, mean),
                              decreasing = FALSE))
    cm$feature  <- factor(cm$feature, levels = feat_levels,
                          labels = pretty_feature(feat_levels))
    cm$tag      <- factor(cm$tag, levels = selection$tag)
    p_cross <- ggplot(cm, aes(x = tag, y = feature, fill = mean_abs_shap)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_viridis_c(name = "mean|SHAP|") +
      labs(title = "SHAP importance across 5 selected models (reweighted grid)",
           x = NULL, y = NULL) +
      theme_paper(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            panel.grid = element_blank())
    ggsave(file.path(fig_dir, "shap_importance_across_models.png"), p_cross,
           width = 9, height = 5.5, dpi = PLOT_DPI)
    write.csv(cm, file.path(fig_dir, "shap_importance_across_models.csv"),
              row.names = FALSE)
  }

  n_png <- length(list.files(fig_dir, pattern = "\\.png$", recursive = TRUE))
  log_section("STEP 4 (XGBoost) COMPLETE")
  log_info("%d PNGs across %d SHAP models -> %s", n_png, nrow(selection), fig_dir)
}
