# =============================================================================
# _paper_values.R  --  Canonical numeric-reference binding (val() system)
# -----------------------------------------------------------------------------
# Single source of truth for data-derived STATISTICS quoted in the prose of
# BOTH `nature_comms_draft.rmd` (main text) and `nature_comms_supplementary.rmd`.
#
# WHY: hand-typed "magic numbers" in the prose (e.g. "CV logloss of 0.100",
# "1,303 donors") drift out of sync every time the pipeline / source CSVs are
# regenerated. This is the number-analogue of the SF() figure-referencing system
# in `_supp_fig_index.R`: author with a stable symbolic KEY and let val(key)
# resolve to the current value, read live from the source CSV, at knit time.
# No stale number can survive a rebuild.
#
# RELATIONSHIP TO SF():
#   SF()  only needs an ORDERING          -> resolver is match(key, keys).
#   val() must READ DATA and COMPUTE       -> resolver reads a source CSV via a
#                                             tiny procurement fn, then formats.
#
# USAGE in an .Rmd (identical wiring to the SF() source line):
#   ```{r} source("_paper_values.R") ```                    # in the setup chunk
#   ... reached `r val("xgb_cv_logloss_realistic", digits = 3)` cross-validated
#       logloss across `r val("n_donors_oligo")` oligodendrocyte donors ...
#
#   Inside a chunk `fig.cap=` you CANNOT use inline backtick-r; build the string
#   with paste0(..., val("key"), ...) exactly as the SF() caption refs do.
#
# TO ADD A NUMBER:
#   1. Find the source CSV and the exact filter that yields the single value.
#   2. register_val("<metric>_<celltype>_<condition>", function() { ... }).
#   3. Replace the literal in the draft with `r val("key", digits=, pct=)`.
#
# GUARDRAILS (mirrors SF()):
#   - val() errors loudly on an unknown key  -> a typo fails the knit.
#   - load_csv() stops with the path if the source CSV is absent -> a moved or
#     renamed source is caught at knit, never shipped blank.
#   - register_val() errors on a duplicate key -> no ambiguous binding.
#   - NEVER fabricate a binding: if a number cannot be traced to a source CSV,
#     leave it hardcoded in the prose and flag it. A fake data-binding is worse
#     than an honest hardcode.
#
# SCOPE NOTE (2026-08-17): this file currently registers the EASY tier only
#   (single-CSV lookups: PsychAD cohort sizes + realistic XGBoost diagnostics).
#   Donor-for-80%-power numbers are deliberately NOT bound: they were model
#   inversions off the now-archived logistic fit and are pending the XGBoost
#   Methods rewrite. The medium tier (specific recall-grid cells) is deferred
#   until a tidy long-format grid CSV is exported. See
#   `numeric_reference_binding_plan.md` for the full taxonomy.
# =============================================================================

suppressPackageStartupMessages(library(here))

# --- Datalink path helper ----------------------------------------------------
.dl <- function(...) here::here("Datalink", ...)

# --- Memoised CSV loader (each source file read at most once per knit) --------
.csv_cache <- new.env(parent = emptyenv())
load_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("val(): source CSV not found: %s", path), call. = FALSE)
  }
  k <- normalizePath(path, mustWork = FALSE)
  if (is.null(.csv_cache[[k]])) {
    .csv_cache[[k]] <- utils::read.csv(
      path, check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  .csv_cache[[k]]
}

# --- Registry: key -> function returning ONE raw (unformatted) numeric --------
.val_registry <- new.env(parent = emptyenv())
register_val <- function(key, fn) {
  if (!is.null(.val_registry[[key]])) {
    stop(sprintf("register_val(): duplicate key '%s' in _paper_values.R", key),
         call. = FALSE)
  }
  .val_registry[[key]] <- fn
  invisible()
}

# --- Resolver: read + optionally format. Errors loudly on an unknown key. -----
# digits : if non-NULL, formatC to fixed decimals (returns a string).
# pct    : multiply a proportion (<=1) by 100 before formatting; a value already
#          on a 0-100 scale is passed through unchanged.
# big    : format with thousands separators (returns a string; implies integer).
# Keys not registered below are looked up in Paper_Source_Data/paper_values.csv
# (written by R/figure_data_generation/09_Paper_Values.R, which documents each
# key's definition in the table itself).
.paper_values_table <- function(key) {
  pv <- load_csv(.dl("Paper_Source_Data", "paper_values.csv"))
  x <- pv$value[pv$key == key]
  if (length(x) != 1) {
    stop(sprintf("val(): unknown key '%s' (not in _paper_values.R or paper_values.csv)", key),
         call. = FALSE)
  }
  x
}
val <- function(key, digits = NULL, pct = FALSE, big = FALSE) {
  fn <- .val_registry[[key]]
  x <- if (is.null(fn)) .paper_values_table(key) else fn()
  if (pct) x <- x * ifelse(x <= 1, 100, 1)
  if (big) return(formatC(round(x), format = "d", big.mark = ","))
  if (!is.null(digits)) return(formatC(x, format = "f", digits = digits))
  x
}

# =============================================================================
# EASY TIER — pure single-CSV lookups. High value, low effort, zero risk.
# =============================================================================

# --- PsychAD filtered cohort sizes -------------------------------------------
# Source: Datalink/Paper_Source_Data/dataset_statistics.csv (multi-dataset;
# filter Dataset == "PsychAD" AND the cell type). Astro 1293 / Micro 1214 /
# Oligo 1303.
.psychad_n_donors <- function(cell_type) {
  d <- load_csv(.dl("Paper_Source_Data", "dataset_statistics.csv"))
  v <- d$N_Donors[d$Dataset == "PsychAD" & d$Cell_Type == cell_type]
  if (length(v) != 1L) {
    stop(sprintf("val(): expected 1 PsychAD %s cohort row, got %d",
                 cell_type, length(v)), call. = FALSE)
  }
  v
}
register_val("n_donors_astro", function() .psychad_n_donors("Astro"))
register_val("n_donors_micro", function() .psychad_n_donors("Micro"))
register_val("n_donors_oligo", function() .psychad_n_donors("Oligo"))

# --- Per-dataset donor-count ranges (across the three cell types) -------------
# Source: Datalink/Paper_Source_Data/dataset_statistics.csv. Each dataset has one
# row per cell type; the min and max N_Donors over those rows give the range
# quoted in the Supplementary dataset-summary table. Binding these prevents the
# stale hand-typed ranges that previously drifted from the CSV (e.g. SEA-AD and
# Rexach).
.dataset_donor_minmax <- function(dataset, which = c("min", "max")) {
  which <- match.arg(which)
  d <- load_csv(.dl("Paper_Source_Data", "dataset_statistics.csv"))
  v <- d$N_Donors[d$Dataset == dataset]
  if (length(v) == 0L) {
    stop(sprintf("val(): no dataset_statistics rows for Dataset == '%s'", dataset),
         call. = FALSE)
  }
  if (which == "min") min(v) else max(v)
}
for (.ds in list(c("seaad", "SEA-AD"), c("rosmap", "ROSMAP"),
                 c("johansen", "Johansen"), c("rexach", "Rexach"))) {
  local({
    key <- .ds[1]; name <- .ds[2]
    register_val(paste0("n_donors_min_", key), function() .dataset_donor_minmax(name, "min"))
    register_val(paste0("n_donors_max_", key), function() .dataset_donor_minmax(name, "max"))
  })
}

# --- Realistic XGBoost model diagnostics -------------------------------------
# Source: Datalink/Power_Curve_Fitting/Models/xgboost/cv_metrics.csv.
# (Keys keep their historical "_realistic" suffix because the manuscript cites
# them by name; there is now only one power model.)
# The realistic (production) model is the promoted winner. The file currently
# holds only that winner row; if a future rerun writes all configs, select the
# winner by the documented realistic ranking (minimum bin_balanced_MAE).
.xgb_realistic_winner <- function() {
  m <- load_csv(.dl("Power_Curve_Fitting", "Models", "xgboost",
                    "cv_metrics.csv"))
  m[which.min(m$bin_balanced_MAE), , drop = FALSE]
}
register_val("xgb_winner_config_realistic", function() .xgb_realistic_winner()$config_index)
register_val("xgb_cv_logloss_realistic",    function() .xgb_realistic_winner()$cv_logloss)
register_val("xgb_mae_realistic",           function() .xgb_realistic_winner()$bin_balanced_MAE)
register_val("xgb_nrounds_realistic",       function() .xgb_realistic_winner()$best_nrounds)

# --- Realistic XGBoost feature gain (training-loss attribution) --------------
# Source: .../Models/xgboost/feature_importance.csv. Gain per feature. NOTE gain is
# redundancy-sensitive; SHAP tells a different story for low-gain features (see
# project.md) -- bind gain only where the prose actually quotes gain.
.xgb_realistic_gain <- function(feature) {
  fi <- load_csv(.dl("Power_Curve_Fitting", "Models", "xgboost",
                     "feature_importance.csv"))
  v <- fi$Gain[fi$Feature == feature]
  if (length(v) != 1L) {
    stop(sprintf("val(): expected 1 gain row for feature '%s', got %d",
                 feature, length(v)), call. = FALSE)
  }
  v
}
register_val("xgb_gain_lfc",        function() .xgb_realistic_gain("abs_gt_lfc"))
register_val("xgb_gain_dispersion", function() .xgb_realistic_gain("gt_dispersion"))
register_val("xgb_gain_ndonors",    function() .xgb_realistic_gain("n_donors"))
register_val("xgb_gain_basemean",   function() .xgb_realistic_gain("gt_basemean"))

# =============================================================================
# SEX-BIASED DE TRUTH SET (padj < 0.05, |LFC| > 0.2) — Results intro paragraph
# =============================================================================
# Source A: Datalink/_scratch_autosomal_vs_sexchr_summary.csv
#   Long format, one row per (cell_type, bin, chr_group) with n + baseMean /
#   dispersion / absLFC quantiles. The truth set is padj < 0.05 & |LFC| > 0.2;
#   summing n across all effect-size bins for a cell type gives its total, and
#   summing within chr_group in {"Autosomal","Sex_chr"} gives the split. This is
#   the same CSV that backs Supplementary Table 2 (autosomal vs sex-chromosome).
# Source B: Datalink/Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv
#   Tested gene universe (PsychAD/Roussos): Total, X_chr, Y_chr, PAR, autosomal.
#   Total / X_chr / Y_chr / PAR are identical across the three cell-type rows for
#   the Roussos dataset, so any Roussos row is representative.

.sexbias_summary <- function() {
  load_csv(.dl("_scratch_autosomal_vs_sexchr_summary.csv"))
}
# n summed across effect-size bins, optionally within a single chr_group.
.sexbias_n <- function(cell_type, chr_group = NULL) {
  d <- .sexbias_summary()
  d <- d[d$cell_type == cell_type, , drop = FALSE]
  if (!is.null(chr_group)) d <- d[d$chr_group == chr_group, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop(sprintf("val(): no sexbias rows for cell_type='%s' chr_group='%s'",
                 cell_type, ifelse(is.null(chr_group), "<any>", chr_group)),
         call. = FALSE)
  }
  sum(d$n)
}
register_val("sexbias_n_total_astro", function() .sexbias_n("Astro"))
register_val("sexbias_n_total_micro", function() .sexbias_n("Micro"))
register_val("sexbias_n_total_oligo", function() .sexbias_n("Oligo"))
register_val("sexbias_n_auto_astro",  function() .sexbias_n("Astro", "Autosomal"))
register_val("sexbias_n_auto_micro",  function() .sexbias_n("Micro", "Autosomal"))
register_val("sexbias_n_auto_oligo",  function() .sexbias_n("Oligo", "Autosomal"))
register_val("sexbias_n_sex_astro",   function() .sexbias_n("Astro", "Sex_chr"))
register_val("sexbias_n_sex_micro",   function() .sexbias_n("Micro", "Sex_chr"))
register_val("sexbias_n_sex_oligo",   function() .sexbias_n("Oligo", "Sex_chr"))

# Percentage autosomal, on a 0-100 scale (quote with digits = 0).
.sexbias_pct_auto <- function(cell_type) {
  100 * .sexbias_n(cell_type, "Autosomal") / .sexbias_n(cell_type)
}
register_val("sexbias_pct_auto_astro", function() .sexbias_pct_auto("Astro"))
register_val("sexbias_pct_auto_micro", function() .sexbias_pct_auto("Micro"))
register_val("sexbias_pct_auto_oligo", function() .sexbias_pct_auto("Oligo"))

# --- Tested gene universe (PsychAD / Roussos) --------------------------------
.gene_universe <- function() {
  u <- load_csv(.dl("Full_Ground_Truth", "Gene_Audit_Full_Universe_Combined.csv"))
  u <- u[u$Dataset == "Roussos", , drop = FALSE]
  if (nrow(u) == 0L) {
    stop("val(): no Roussos row in Gene_Audit_Full_Universe_Combined.csv",
         call. = FALSE)
  }
  u[1L, , drop = FALSE]
}
register_val("gene_universe_total", function() .gene_universe()$Total)
# Allosomal (X + Y + PAR) share of the tested universe, 0-100 scale.
register_val("gene_universe_allosomal_pct", function() {
  u <- .gene_universe()
  100 * (u$X_chr + u$Y_chr + u$PAR) / u$Total
})

# =============================================================================
# MEDIUM TIER — specific recall/precision grid cells.
# =============================================================================
# Source: Datalink/Paper_Source_Data/recall_precision_grid.csv, the tidy
# long-format flatten of the Aggregated_Summary_Data.qs checkpoints
# (produced by figure_data_generation/08_Recall_Precision_Grid.R). One row per
#   (cell_type, mode, method, n_donors, y_value, lfc_bin, gene_category)
# with mean_recall / mean_precision (proportions, 0-1) and ref_gt_count.
#
# GRAIN for the paper's headline recall/precision claims: the pooled
# `all_methods` model over ALL ground-truth genes (`Total` gene_category), cells
# mode. lfc_bin is matched by prefix ("Tiny", "Small", ...) so the exact bin
# label (e.g. "Tiny (0.2-1.0)") can change without breaking the key.
.recall_grid <- function() {
  load_csv(.dl("Paper_Source_Data", "recall_precision_grid.csv"))
}
# Resolve ONE grid cell to a single metric value (a 0-1 proportion).
.grid_cell <- function(cell_type, n_donors, y_value, lfc_prefix,
                       metric = c("mean_recall", "mean_precision"),
                       gene_category = "Total", method = "all_methods",
                       mode = "cells") {
  metric <- match.arg(metric)
  d <- .recall_grid()
  d <- d[d$cell_type == cell_type & d$mode == mode & d$method == method &
           d$gene_category == gene_category & d$n_donors == n_donors &
           d$y_value == y_value & startsWith(d$lfc_bin, lfc_prefix), , drop = FALSE]
  if (nrow(d) != 1L) {
    stop(sprintf(
      "val(): expected 1 grid cell (%s %s d=%d y=%d %s* %s/%s) got %d",
      cell_type, metric, n_donors, y_value, lfc_prefix, method, gene_category,
      nrow(d)), call. = FALSE)
  }
  d[[metric]]
}
# Tiny-effect bin recall at 300 donors x 100 cells (Total genes, all_methods).
# Verified 2026-08-22: Micro 0.206, Astro 0.124, Oligo 0.122 (quote pct = TRUE).
register_val("recall_tiny_d300_c100_astro",
             function() .grid_cell("Astro", 300, 100, "Tiny", "mean_recall"))
register_val("recall_tiny_d300_c100_micro",
             function() .grid_cell("Micro", 300, 100, "Tiny", "mean_recall"))
register_val("recall_tiny_d300_c100_oligo",
             function() .grid_cell("Oligo", 300, 100, "Tiny", "mean_recall"))
