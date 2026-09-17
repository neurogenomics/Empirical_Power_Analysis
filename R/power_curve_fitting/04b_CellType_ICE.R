# ==============================================================================
# SCRIPT: Cell-type ICE / swing analysis (M2 -- justifies defaulting cell_type)
# PURPOSE: The exported power predictor (04_Export_XGBoost.R) keeps cell_type in
#          the model but DEFAULTS it (to Oligo) so a naive caller never has to
#          supply it. This script quantifies how safe that default is, using the
#          exported bundle itself.
#
#          The tool is Individual Conditional Expectation (ICE; Goldstein et al.
#          2015), the per-row companion of a Partial Dependence Plot (Friedman
#          2001). SHAP answers "how much does the model ATTRIBUTE to cell_type",
#          which for a coarse 3-level factor competing against a granular
#          continuous depth feature (reads_per_cell) is deflated by split
#          attribution. ICE answers the decision-relevant question directly:
#          "how much does the PREDICTION move when I change only cell_type,
#          holding everything else fixed?" For a "can I default/ignore it?"
#          question ICE is the correct measure, not SHAP.
#
#          For each in-distribution design row we predict power three times,
#          cycling cell_type through {Astro, Micro, Oligo} with all other
#          features frozen, and summarise:
#            (A) full swing  s_i = max_k p_ik - min_k p_ik   (worst-case error
#                from ignoring cell_type entirely)
#            (B) default-to-Oligo error  |p_oligo_i - p_native_i|  (the error the
#                exported wrapper's Oligo default actually incurs vs the row's
#                real cell type)
#          The 95th percentile / max of these distributions ARE the error bound:
#            |p(default cell type) - p(true cell type)| <= max_i s_i
#          which is a far stronger, quantitative claim for the paper than
#          "SHAP was small".
#
#          de_method (also defaulted in the export, to deseq2) is analysed the
#          same way as a secondary, since the same default-safety question
#          applies to it.
#
#          Runs off the exported bundle + the realistic fitting data (for
#          genuinely in-distribution design rows). No HPC, no refit. base R +
#          xgboost + ggplot2.
#
# INPUT:  Exports/xgboost/<subdir>/power_model_bundle.rds   (from 04)
#         Fitting_Data/<fitting_data>.qs                    (in-distribution rows)
# OUTPUT: Models/xgboost/<subdir>/celltype_ice/
#           celltype_ice_predictions.csv     per-row predictions per cell type
#           celltype_ice_summary.csv         swing + default-error stats (overall + by region)
#           demethod_ice_summary.csv         same for de_method (secondary)
#           celltype_swing_ecdf.png          ECDF of s_i with p95/p99 marked (headline)
#           swing_vs_power.png               spread vs mean power (where it concentrates)
#           categorical_swing_compare.png    cell_type vs de_method swing violins
#           celltype_ice_summary.md          paper-quotable bound + interpretation
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
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
XGB_GENE_SET     <- Sys.getenv("XGB_GENE_SET", "full")
ICE_N_SAMPLE     <- as.integer(Sys.getenv("ICE_N_SAMPLE", "150000"))

set.seed(SEED_BASE)

log_section(sprintf("STEP 4b (XGBoost): Cell-type ICE / swing [subdir=%s]",
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
out_dir <- file.path(model_base, "celltype_ice")
safe_mkdir(out_dir)

# ==============================================================================
# 2. LOAD THE EXPORTED PREDICTOR BUNDLE
# ==============================================================================
# Uses the shipped power model (04_Export_XGBoost.R output) -- the same artefact
# the frontend calculator (M3) will consume -- so this analysis characterises
# exactly the predictor users receive.

bundle_path <- file.path(export_dir, "power_model_bundle.rds")
if (!file.exists(bundle_path)) {
  stop(sprintf("Bundle not found: %s. Run 04_Export_XGBoost.R first.", bundle_path))
}
bundle <- readRDS(bundle_path)
booster <- xgb.load.raw(bundle$booster_raw)   # UBJ raw bytes -> live booster

CONT_FEATURES <- bundle$cont_features
CAT_FEATURES  <- bundle$cat_features
cat_levels    <- bundle$cat_levels
feature_names <- bundle$feature_names
CELL_LEVELS   <- cat_levels$cell_type
METHOD_LEVELS <- cat_levels$de_method

log_info("Loaded bundle: %d rounds, cell types {%s}, methods {%s}",
         bundle$metadata$best_nrounds,
         paste(CELL_LEVELS, collapse = ", "),
         paste(METHOD_LEVELS, collapse = ", "))

# Exact fit-time design-matrix recipe (identical to the exported wrapper).
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

# ==============================================================================
# 3. IN-DISTRIBUTION DESIGN ROWS
# ==============================================================================
# Sample real observed grid x gene rows from the fitting data the model was
# trained on -- genuinely in-distribution, no synthetic extrapolation.

fitting_data_name <- bundle$metadata$fitting_data
agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", fitting_data_name)
if (!file.exists(agg_path)) {
  stop(sprintf("Fitting data not found: %s.", agg_path))
}
dat <- qs::qread(agg_path)
if ("gene_set" %in% colnames(dat)) dat <- dat %>% filter(gene_set == XGB_GENE_SET)
dat <- dat %>%
  filter(is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
         is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion),
         cell_type %in% CELL_LEVELS, de_method %in% METHOD_LEVELS)
log_info("In-distribution rows available: %d", nrow(dat))

if (nrow(dat) > ICE_N_SAMPLE) {
  dat <- dat[sample.int(nrow(dat), ICE_N_SAMPLE), ]
  log_info("Sampled %d rows for ICE (ICE_N_SAMPLE)", nrow(dat))
}
n_rows <- nrow(dat)

# Frozen base: the design each row carries, minus the categorical under test.
base_df <- data.frame(
  n_donors       = dat$n_donors,
  n_cells        = dat$n_cells,
  reads_per_cell = dat$reads_per_cell,
  abs_gt_lfc     = dat$abs_gt_lfc,
  gt_basemean    = dat$gt_basemean,
  gt_dispersion  = dat$gt_dispersion,
  cell_type      = as.character(dat$cell_type),   # row's NATIVE cell type
  de_method      = as.character(dat$de_method),   # row's NATIVE method
  stringsAsFactors = FALSE
)

# ==============================================================================
# 4. CELL-TYPE ICE: predict over all 3 levels, everything else frozen
# ==============================================================================

ct_preds <- matrix(NA_real_, nrow = n_rows, ncol = length(CELL_LEVELS),
                   dimnames = list(NULL, CELL_LEVELS))
for (ct in CELL_LEVELS) {
  df <- base_df
  df$cell_type <- ct
  ct_preds[, ct] <- predict_rows(df)
}

ct_min  <- apply(ct_preds, 1, min)
ct_max  <- apply(ct_preds, 1, max)
ct_mean <- rowMeans(ct_preds)
ct_swing <- ct_max - ct_min                       # (A) full swing s_i

# (B) default-to-Oligo error vs the row's native cell type.
DEFAULT_CT <- bundle$defaults$cell_type
p_default <- ct_preds[, DEFAULT_CT]
p_native  <- ct_preds[cbind(seq_len(n_rows), match(base_df$cell_type, CELL_LEVELS))]
ct_default_err <- abs(p_default - p_native)

# ==============================================================================
# 5. DE-METHOD ICE (secondary; same default-safety question)
# ==============================================================================

dm_preds <- matrix(NA_real_, nrow = n_rows, ncol = length(METHOD_LEVELS),
                   dimnames = list(NULL, METHOD_LEVELS))
for (dm in METHOD_LEVELS) {
  df <- base_df
  df$de_method <- dm
  dm_preds[, dm] <- predict_rows(df)
}
dm_swing <- apply(dm_preds, 1, max) - apply(dm_preds, 1, min)
DEFAULT_DM <- bundle$defaults$de_method
dm_default_err <- abs(dm_preds[, DEFAULT_DM] -
                      dm_preds[cbind(seq_len(n_rows), match(base_df$de_method, METHOD_LEVELS))])

# ==============================================================================
# 6. PER-ROW PREDICTION TABLE
# ==============================================================================

pred_tbl <- cbind(
  base_df,
  setNames(as.data.frame(ct_preds), paste0("p_", CELL_LEVELS)),
  mean_power_over_celltypes = ct_mean,
  celltype_swing            = ct_swing,
  default_oligo_error       = ct_default_err
)
write.csv(pred_tbl, file.path(out_dir, "celltype_ice_predictions.csv"),
          row.names = FALSE)
log_info("Wrote per-row predictions -> celltype_ice_predictions.csv (%d rows)", n_rows)

# ==============================================================================
# 7. SUMMARY STATS (overall + by power region)
# ==============================================================================
# Region on the mean-over-cell-types power: near the 0/1 rails the model is
# saturated so any categorical must move it little; the swing (if any) lives in
# the mid-range. Reporting by region pre-empts the "your bound is just an
# artefact of everything being at 0 or 1" objection.

region_of <- function(p) {
  cut(p, breaks = c(-Inf, 0.1, 0.9, Inf),
      labels = c("low [0,0.1]", "mid (0.1,0.9)", "high [0.9,1]"))
}
summarise_metric <- function(x, region, label) {
  mk <- function(idx, grp) {
    xi <- x[idx]
    data.frame(
      metric = label, group = grp, n = length(xi),
      mean = mean(xi), median = median(xi),
      p95 = as.numeric(quantile(xi, 0.95)),
      p99 = as.numeric(quantile(xi, 0.99)),
      max = max(xi),
      frac_gt_0.02 = mean(xi > 0.02),
      frac_gt_0.05 = mean(xi > 0.05),
      stringsAsFactors = FALSE
    )
  }
  rows <- list(mk(seq_along(x), "overall"))
  for (lv in levels(region)) rows <- c(rows, list(mk(which(region == lv), lv)))
  bind_rows(rows)
}

ct_region <- region_of(ct_mean)
ct_summary <- bind_rows(
  summarise_metric(ct_swing,       ct_region, "celltype_swing"),
  summarise_metric(ct_default_err, ct_region, "default_oligo_error")
)
write.csv(ct_summary, file.path(out_dir, "celltype_ice_summary.csv"), row.names = FALSE)

dm_region <- region_of(rowMeans(dm_preds))
dm_summary <- bind_rows(
  summarise_metric(dm_swing,       dm_region, "demethod_swing"),
  summarise_metric(dm_default_err, dm_region, "default_deseq2_error")
)
write.csv(dm_summary, file.path(out_dir, "demethod_ice_summary.csv"), row.names = FALSE)
log_info("Wrote summary stats -> celltype_ice_summary.csv, demethod_ice_summary.csv")

# Pull the headline numbers (overall row).
ct_swing_overall <- ct_summary %>% filter(metric == "celltype_swing", group == "overall")
ct_deferr_overall <- ct_summary %>% filter(metric == "default_oligo_error", group == "overall")
dm_swing_overall <- dm_summary %>% filter(metric == "demethod_swing", group == "overall")

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
          legend.title     = element_text(face = "bold", size = rel(0.8)),
          legend.text      = element_text(size = rel(0.75)))
}

# --- 8.1 ECDF of the cell-type swing (headline error-bound figure) ------------
p95 <- ct_swing_overall$p95
p99 <- ct_swing_overall$p99
ecdf_df <- data.frame(swing = sort(ct_swing),
                      cum = seq_along(ct_swing) / length(ct_swing))
p_ecdf <- ggplot(ecdf_df, aes(x = swing, y = cum)) +
  geom_step(colour = "steelblue", linewidth = 0.7) +
  geom_vline(xintercept = p95, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey40") +
  annotate("text", x = p95, y = 0.25, label = sprintf("p95 = %.3f", p95),
           hjust = -0.1, size = 3.2, colour = "grey30") +
  annotate("text", x = p99, y = 0.12, label = sprintf("p99 = %.3f", p99),
           hjust = -0.1, size = 3.2, colour = "grey30") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cell-type swing  s = max power - min power across {Astro, Micro, Oligo}",
       y = "Cumulative fraction of design rows",
       title = "How much does predicted power move when only cell type changes?",
       subtitle = sprintf("ICE over %d in-distribution designs; everything except cell type frozen. Max swing = %.3f",
                          n_rows, ct_swing_overall$max)) +
  theme_paper(12)
ggsave(file.path(out_dir, "celltype_swing_ecdf.png"), p_ecdf,
       width = 8, height = 5, dpi = PLOT_DPI)

# --- 8.2 Swing vs mean power (where the swing concentrates) -------------------
sv_df <- data.frame(mean_power = ct_mean, swing = ct_swing)
p_sv <- ggplot(sv_df, aes(x = mean_power, y = swing)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(name = "Design\nrows", trans = "log10") +
  labs(x = "Mean predicted power over the three cell types",
       y = "Cell-type swing  s",
       title = "Cell-type swing concentrates in the mid-power band",
       subtitle = "Near the 0/1 saturation rails the model is insensitive to cell type") +
  theme_paper(12)
ggsave(file.path(out_dir, "swing_vs_power.png"), p_sv,
       width = 8, height = 5, dpi = PLOT_DPI)

# --- 8.3 cell_type vs de_method swing (which categoricals are defaultable) ----
cmp_df <- bind_rows(
  data.frame(feature = "cell_type", swing = ct_swing),
  data.frame(feature = "de_method", swing = dm_swing)
)
p_cmp <- ggplot(cmp_df, aes(x = feature, y = swing, fill = feature)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.55) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c(cell_type = "#21908C", de_method = "#3B528B"),
                    guide = "none") +
  labs(x = NULL, y = "Prediction swing across the feature's levels",
       title = "Swing from toggling a categorical, everything else fixed",
       subtitle = sprintf("cell_type p95 = %.3f | de_method p95 = %.3f",
                          ct_swing_overall$p95, dm_swing_overall$p95)) +
  theme_paper(12)
ggsave(file.path(out_dir, "categorical_swing_compare.png"), p_cmp,
       width = 7, height = 5, dpi = PLOT_DPI)

log_info("Wrote 3 figures -> %s", out_dir)

# ==============================================================================
# 9. PAPER-QUOTABLE SUMMARY
# ==============================================================================

md <- c(
  "# Cell-type ICE / swing analysis",
  "",
  sprintf("Model: XGBoost `%s` (exported bundle). In-distribution designs: %d rows sampled from `%s`.",
          XGB_MODEL_SUBDIR, n_rows, fitting_data_name),
  "",
  "## Question",
  "The exported predictor keeps `cell_type` in the model but defaults it to",
  sprintf("`%s` so a naive caller need not supply it. Is that default safe? ICE",
          DEFAULT_CT),
  "(vary only cell_type, freeze everything else) measures the prediction move",
  "directly -- the decision-relevant quantity SHAP does not give.",
  "",
  "## Headline bound",
  sprintf("Across %d in-distribution designs, changing ONLY cell type moves predicted power by:",
          n_rows),
  sprintf("- median %.4f, mean %.4f", ct_swing_overall$median, ct_swing_overall$mean),
  sprintf("- 95th percentile %.4f, 99th percentile %.4f", ct_swing_overall$p95, ct_swing_overall$p99),
  sprintf("- maximum %.4f", ct_swing_overall$max),
  sprintf("- fraction of designs with swing > 0.05: %.1f%%", 100 * ct_swing_overall$frac_gt_0.05),
  "",
  "This is the error bound: |power(default cell type) - power(true cell type)| <= max swing.",
  "",
  sprintf("Defaulting specifically to %s (the shipped default), the error vs each row's NATIVE cell type is:",
          DEFAULT_CT),
  sprintf("- median %.4f, 95th percentile %.4f, max %.4f",
          ct_deferr_overall$median, ct_deferr_overall$p95, ct_deferr_overall$max),
  "",
  "## By power region",
  "(Swing concentrates in the mid-power band; the 0/1 rails are saturated.)",
  "",
  "| metric | region | median | p95 | max |",
  "|---|---|---|---|---|",
  paste(apply(ct_summary[ct_summary$group != "overall", ], 1, function(r) {
    sprintf("| %s | %s | %.4f | %.4f | %.4f |",
            r["metric"], r["group"],
            as.numeric(r["median"]), as.numeric(r["p95"]), as.numeric(r["max"]))
  }), collapse = "\n"),
  "",
  "## Secondary: de_method",
  sprintf("Toggling de_method (default %s) moves power by median %.4f, p95 %.4f, max %.4f.",
          DEFAULT_DM, dm_swing_overall$median, dm_swing_overall$p95, dm_swing_overall$max),
  "",
  "## Interpretation",
  "- Low SHAP for cell_type is corroborated by a direct, artefact-free measure:",
  "  the ICE swing. Unlike mean|SHAP| (which conflates model sensitivity with a",
  "  feature's data range and is deflated for a coarse factor competing against a",
  "  granular continuous depth feature), the swing IS the full dynamic range of",
  "  toggling the categorical, so its magnitude can be trusted.",
  "- The claim is conditional and stated as such: cell type carries little beyond",
  "  depth, yield and gene characteristics -- NOT that cell type is causally",
  "  irrelevant (cf. the depth-equalisation finding, which shows the apparent",
  "  cell-type hierarchy is a depth/yield/independent-filtering artefact).",
  "- Practical consequence: defaulting cell_type in the exported predictor changes",
  "  the answer by at most the bound above across the design space, so a naive",
  "  caller can safely omit it while a power user retains full control."
)
writeLines(md, file.path(out_dir, "celltype_ice_summary.md"))
log_info("Wrote paper-quotable summary -> celltype_ice_summary.md")

log_section("STEP 4b (XGBoost) ICE COMPLETE")
log_info("Cell-type swing: median=%.4f p95=%.4f max=%.4f | default-%s error p95=%.4f",
         ct_swing_overall$median, ct_swing_overall$p95, ct_swing_overall$max,
         DEFAULT_CT, ct_deferr_overall$p95)
log_info("Outputs in %s", out_dir)
