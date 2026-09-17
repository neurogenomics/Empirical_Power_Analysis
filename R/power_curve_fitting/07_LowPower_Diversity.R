# ==============================================================================
# SCRIPT: Low-power mass diversity exploration (M2 -- diagnostic, no model)
# PURPOSE: The realistic-grid fitting set (n_cells >= 100) is still dominated
#          by rows with observed detection rate in [0, 0.1] (~unpowered). If
#          this "low-power mass" is HOMOGENEOUS (one archetype: e.g. tiny LFC
#          + adequate expression at any donor count) then class-reweighting
#          the XGBoost model is fighting an easy battle. If it is DIVERSE
#          (multiple archetypes: some are tiny LFC, some are low expression,
#          some are few donors, etc.) the reweighting has more room to help.
#          This diagnostic characterises which regime we're in via:
#            - feature marginals of low-power vs global (§5.1)
#            - 2D density heatmaps coloured by mean observed power (§5.2)
#            - k-means sub-clustering on standardised features (§5.3)
#
#          Runs on the FILTERED aggregated data (n_cells >= 100). No model
#          fitting; uses base R + ggplot2 only. UMAP intentionally SKIPPED
#          this pass (design doc §5.4, Chris opted out on 2026-07-30).
#
# INPUT:  Fitting_Data/aggregated_fitting_data_realistic.qs   (from 01e)
# OUTPUT: Models/xgboost/realistic/lowpower_diversity/
#           feature_marginals/{feature}_hist.png    (5 features)
#           density2d/{x}_vs_{y}_meanpower.png       (4 pairs)
#           kmeans_cluster_profile.csv               (cluster mean/median per feature)
#           kmeans_pca_scatter.png                    (points on first 2 PCs)
#           kmeans_summary.md                         (paper-quotable interpretation)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(viridis)
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

log_section("STEP 7: Low-power mass diversity exploration")

XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "realistic")
XGB_FITTING_DATA <- Sys.getenv("XGB_FITTING_DATA", "aggregated_fitting_data_realistic.qs")
LOW_POWER_MAX    <- as.numeric(Sys.getenv("XGB_LOW_POWER_MAX", "0.1"))
KMEANS_K         <- as.integer(Sys.getenv("XGB_KMEANS_K", "4"))

set.seed(SEED_BASE)

model_base <- file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
out_dir    <- file.path(model_base, "lowpower_diversity")
safe_mkdir(out_dir)
safe_mkdir(file.path(out_dir, "feature_marginals"))
safe_mkdir(file.path(out_dir, "density2d"))

# ==============================================================================
# 2. LOAD FILTERED DATA
# ==============================================================================

agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", XGB_FITTING_DATA)
if (!file.exists(agg_path)) {
  stop(sprintf("Filtered data not found: %s. Run 01e_Filter_Realistic.R first.", agg_path))
}
dat <- qs::qread(agg_path)
dat <- dat %>%
  filter(is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
         is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
         is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion))
log_info("Filtered rows: %d (from %s)", nrow(dat), XGB_FITTING_DATA)

low_mask <- dat$detection_rate <= LOW_POWER_MAX
n_low <- sum(low_mask)
log_info("Low-power (obs <= %.2f): %d rows (%.1f%% of filtered set)",
         LOW_POWER_MAX, n_low, 100 * n_low / nrow(dat))

# Paper-style theme.
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
# 3. FEATURE MARGINALS (§5.1)
# ==============================================================================

FEATURES <- c(
  abs_gt_lfc     = "|LFC|",
  gt_basemean    = "Mean expression (log axis)",
  n_donors       = "Donors",
  n_cells        = "Cells / donor",
  gt_dispersion  = "Dispersion"
)
LOG_FEATURES <- c("gt_basemean")   # apply log axis for readability

for (f in names(FEATURES)) {
  df <- data.frame(
    value  = c(dat[[f]], dat[[f]][low_mask]),
    subset = c(rep("all filtered rows", nrow(dat)),
               rep(sprintf("obs <= %.2f", LOW_POWER_MAX), n_low)),
    stringsAsFactors = FALSE
  )
  df$subset <- factor(df$subset, levels = c(sprintf("obs <= %.2f", LOW_POWER_MAX),
                                             "all filtered rows"))
  p <- ggplot(df, aes(x = value, fill = subset)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   position = "identity", alpha = 0.45, colour = NA) +
    scale_fill_manual(values = c("steelblue", "grey60"), name = NULL) +
    labs(x = FEATURES[[f]], y = "Density",
         title = sprintf("Marginal distribution: %s", FEATURES[[f]]),
         subtitle = sprintf("Low-power (obs <= %.2f) vs full filtered set",
                            LOW_POWER_MAX)) +
    theme_paper(12) +
    theme(legend.position = "top")
  if (f %in% LOG_FEATURES) p <- p + scale_x_log10()
  ggsave(file.path(out_dir, "feature_marginals",
                   sprintf("%s_hist.png", f)),
         p, width = 7, height = 4.5, dpi = PLOT_DPI)
}
log_info("Wrote %d marginal histograms -> feature_marginals/", length(FEATURES))

# ==============================================================================
# 4. 2D DENSITY HEATMAPS (§5.2)
# ==============================================================================

# Coloured by mean observed power (viridis). Uses stat_summary_2d for a
# regular grid; more legible than raw scatter with millions of points.
DENSITY_PAIRS <- list(
  c("abs_gt_lfc",  "n_donors"),
  c("abs_gt_lfc",  "n_cells"),
  c("gt_basemean", "n_donors"),
  c("gt_basemean", "n_cells")
)

for (pr in DENSITY_PAIRS) {
  xf <- pr[1]; yf <- pr[2]
  x_lab <- FEATURES[[xf]]; y_lab <- FEATURES[[yf]]

  df <- data.frame(x = dat[[xf]], y = dat[[yf]], power = dat$detection_rate)
  p <- ggplot(df, aes(x = x, y = y, z = power)) +
    stat_summary_2d(fun = mean, bins = 60) +
    scale_fill_viridis_c(name = "Mean observed\ndetection rate",
                         limits = c(0, 1), oob = scales::squish) +
    labs(x = x_lab, y = y_lab,
         title = sprintf("Mean observed power over %s x %s", x_lab, y_lab),
         subtitle = "Bottom-left corner homogeneous / spread diverse") +
    theme_paper(12)
  if (xf %in% LOG_FEATURES) p <- p + scale_x_log10()
  if (yf %in% LOG_FEATURES) p <- p + scale_y_log10()
  ggsave(file.path(out_dir, "density2d",
                   sprintf("%s_vs_%s_meanpower.png", xf, yf)),
         p, width = 7, height = 5, dpi = PLOT_DPI)
}
log_info("Wrote %d 2D density heatmaps -> density2d/", length(DENSITY_PAIRS))

# ==============================================================================
# 5. K-MEANS SUB-CLUSTERING OF THE LOW-POWER BIN (§5.3)
# ==============================================================================

lp <- dat[low_mask, ]
km_feat_raw <- data.frame(
  abs_gt_lfc     = lp$abs_gt_lfc,
  gt_basemean    = log1p(lp$gt_basemean),   # log for skew
  n_donors       = lp$n_donors,
  n_cells        = lp$n_cells,
  gt_dispersion  = lp$gt_dispersion
)
# Standardise (each feature mean 0, sd 1) so k-means isn't dominated by scale.
km_feat_std <- scale(km_feat_raw)

set.seed(42)
km <- kmeans(km_feat_std, centers = KMEANS_K, nstart = 25, iter.max = 50)
lp$cluster <- factor(km$cluster)

# --- cluster profile CSV ------------------------------------------------------
profile_rows <- lapply(levels(lp$cluster), function(k) {
  sub <- lp[lp$cluster == k, ]
  data.frame(
    cluster       = k,
    n             = nrow(sub),
    pct_low_power = 100 * nrow(sub) / n_low,
    mean_abs_lfc  = mean(sub$abs_gt_lfc,     na.rm = TRUE),
    median_abs_lfc= median(sub$abs_gt_lfc,   na.rm = TRUE),
    mean_basemean = mean(sub$gt_basemean,    na.rm = TRUE),
    median_basemean = median(sub$gt_basemean, na.rm = TRUE),
    mean_donors   = mean(sub$n_donors,       na.rm = TRUE),
    mean_cells    = mean(sub$n_cells,        na.rm = TRUE),
    mean_dispersion = mean(sub$gt_dispersion, na.rm = TRUE)
  )
})
profile <- bind_rows(profile_rows)
write.csv(profile, file.path(out_dir, "kmeans_cluster_profile.csv"), row.names = FALSE)
log_info("Cluster profile saved (K=%d)", KMEANS_K)

# --- PCA scatter --------------------------------------------------------------
pca <- prcomp(km_feat_std, center = FALSE, scale. = FALSE)
pc_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cluster = lp$cluster)
p_pca <- ggplot(pc_df, aes(x = PC1, y = PC2, colour = cluster)) +
  geom_point(size = 0.3, alpha = 0.35, stroke = 0) +
  scale_colour_brewer(palette = "Set1", name = "Cluster") +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * pca$sdev[1]^2 / sum(pca$sdev^2)),
       y = sprintf("PC2 (%.1f%%)", 100 * pca$sdev[2]^2 / sum(pca$sdev^2)),
       title = sprintf("Low-power mass k-means (K=%d)", KMEANS_K),
       subtitle = "Standardised {|LFC|, log(baseMean), donors, cells, dispersion}") +
  theme_paper(12) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "kmeans_pca_scatter.png"),
       p_pca, width = 7.5, height = 5.5, dpi = PLOT_DPI)

# --- paper-quotable interpretation lines --------------------------------------
# Heuristic archetype labels from the standardised centres.
centres <- km$centers   # rows = cluster; cols = features on std scale
label_cluster <- function(row) {
  parts <- c()
  if (row["abs_gt_lfc"]     < -0.4) parts <- c(parts, "low |LFC|")
  if (row["abs_gt_lfc"]     >  0.4) parts <- c(parts, "high |LFC|")
  if (row["gt_basemean"]    < -0.4) parts <- c(parts, "low expression")
  if (row["gt_basemean"]    >  0.4) parts <- c(parts, "high expression")
  if (row["n_donors"]       < -0.4) parts <- c(parts, "few donors")
  if (row["n_donors"]       >  0.4) parts <- c(parts, "many donors")
  if (row["n_cells"]        < -0.4) parts <- c(parts, "few cells")
  if (row["n_cells"]        >  0.4) parts <- c(parts, "many cells")
  if (row["gt_dispersion"]  >  0.4) parts <- c(parts, "high dispersion")
  if (length(parts) == 0) parts <- c("mixed / near-centroid")
  paste(parts, collapse = ", ")
}
lines <- vapply(seq_len(nrow(centres)), function(i) {
  sprintf("- Cluster %d (%d rows, %.1f%%): %s",
          i, profile$n[i], profile$pct_low_power[i],
          label_cluster(centres[i, ]))
}, character(1))

writeLines(c(
  sprintf("# Low-power mass k-means (K=%d)", KMEANS_K),
  "",
  sprintf("Low-power mass = rows with observed detection rate <= %.2f", LOW_POWER_MAX),
  sprintf("N low-power rows: %d (of %d filtered)", n_low, nrow(dat)),
  "",
  "## Archetypes (auto-labelled from standardised cluster centres)",
  lines,
  "",
  "Interpretation: if archetypes span several failure modes (small effect,",
  "low expression, few donors) the low-power mass is DIVERSE and reweighting",
  "has more room to help. If all archetypes concentrate on one axis it is",
  "HOMOGENEOUS and reweighting is fighting an easy battle."
), file.path(out_dir, "kmeans_summary.md"))
log_info("Wrote k-means summary + PCA scatter -> %s", out_dir)

log_section("STEP 7 COMPLETE")
log_info("Diagnostics written to %s", out_dir)
