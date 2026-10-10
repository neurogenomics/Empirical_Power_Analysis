# ==============================================================================
# STEP 09: Paper Values -- every data-derived number quoted in the manuscript
# ==============================================================================
# Writes ONE table, Datalink/Paper_Source_Data/paper_values.csv, with columns
#   key, value, unit, definition, source
# The manuscript reads it through val("key") (_paper_values.R falls back to this
# table for any key it does not register itself). Each definition is written
# out in full so a reader can check exactly what a quoted number means.
#
# Units: "percent" values are on a 0-100 scale; "proportion" on 0-1; counts
# and model quantities are raw.
#
# Inputs (all produced by earlier stages):
#   Paper_Source_Data/dataset_statistics.csv, cohort_prefilter.csv   (eda)
#   Sex_Validation/sex_exclusion_combined.csv                          (eda)
#   _scratch_autosomal_vs_sexchr_summary.csv                           (tables)
#   Paper_Source_Data/recall_precision_grid.csv                        (tables)
#   Depth_Eq_DEG_Loss/deg_loss_table.csv                               (deg_loss)
#   Power_Curve_Fitting/Models/xgboost/{model_selection,calibration_residuals}.csv (xgb_select)
#   Power_Curve_Fitting/Models/xgboost/comparison/unified_calibration_table.csv   (xgb_post)
#   Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data.qs        (pcf_aggregate)
#   Gene_Classification/oligo_gene_group_power_summary.csv             (gene_class)
#   R/depth_equalization/config.R (the downsampling factor actually used)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 09: PAPER VALUES ###\n")

dl <- function(...) project_path("Datalink", ...)
rd <- function(...) utils::read.csv(dl(...), check.names = FALSE, stringsAsFactors = FALSE)

values <- list()
add <- function(key, value, unit, definition, source) {
  if (length(value) != 1 || !is.finite(value)) stop("paper value ", key, " is not a single finite number")
  if (key %in% names(values)) stop("duplicate paper value key: ", key)
  values[[key]] <<- data.frame(key = key, value = unname(value), unit = unit,
                               definition = definition, source = source)
}
cts <- c(Astro = "astro", Micro = "micro", Oligo = "oligo")

# --- 1. Cohort -----------------------------------------------------------------
src <- "Paper_Source_Data/cohort_prefilter.csv"
pre <- rd(src)
add("n_donors_prefilter", pre$n_donors[pre$cell_type == "All"], "count",
    "PsychAD donors with any astrocytes, microglia or oligodendrocytes in the unfiltered (full) shards, before the age and cell-count filters", src)
for (ct in names(cts)) {
  add(paste0("n_donors_prefilter_", cts[[ct]]), pre$n_donors[pre$cell_type == ct], "count",
      sprintf("PsychAD donors in the unfiltered %s shards", ct), src)
}

src <- "Paper_Source_Data/dataset_statistics.csv"
ds <- rd(src) %>% filter(Dataset == "PsychAD")
for (ct in names(cts)) {
  r <- ds[ds$Cell_Type == ct, ]
  add(paste0("mean_cells_per_donor_", cts[[ct]]), r$Mean_Cells_Per_Donor, "count",
      sprintf("Mean %s cells per donor, filtered cohort", ct), src)
  add(paste0("mean_umis_per_cell_", cts[[ct]]), r$Mean_UMIs_Per_Cell, "count",
      sprintf("Mean UMIs per %s cell, filtered cohort", ct), src)
  add(paste0("mean_umis_per_donor_millions_", cts[[ct]]),
      r$Mean_Cells_Per_Donor * r$Mean_UMIs_Per_Cell / 1e6, "millions",
      sprintf("Mean cells per donor x mean UMIs per cell, %s, in millions", ct), src)
}
add("umis_per_donor_ratio_oligo_micro",
    values$mean_umis_per_donor_millions_oligo$value / values$mean_umis_per_donor_millions_micro$value,
    "ratio", "Mean UMIs per donor, oligodendrocytes divided by microglia", src)

# --- 2. Sex-label validation ---------------------------------------------------
src <- "Sex_Validation/sex_exclusion_combined.csv"
sx <- rd(src)
n_flag <- length(unique(sx$donor_id))
add("n_sex_discordant_donors", n_flag, "count",
    "Distinct donors flagged as sex-discordant in any cell type", src)
add("pct_sex_discordant_donors", 100 * n_flag / values$n_donors_prefilter$value, "percent",
    "Sex-discordant donors as a percentage of n_donors_prefilter", src)

# --- 3. Sex-biased truth genes: autosomal vs sex-chromosome, tiny bin ---------
src <- "_scratch_autosomal_vs_sexchr_summary.csv"
av <- rd(src) %>% filter(grepl("^Tiny", bin))
for (ct in names(cts)) {
  for (g in c("Autosomal", "Sex_chr")) {
    r <- av[av$cell_type == ct & av$chr_group == g, ]
    tag <- if (g == "Autosomal") "auto" else "sexchr"
    add(sprintf("tiny_median_dispersion_%s_%s", tag, cts[[ct]]), r$dispersion_median, "dispersion",
        sprintf("Median DESeq2 dispersion of %s full-cohort truth genes with |LFC| 0.2-1, %s", g, ct), src)
    add(sprintf("tiny_median_basemean_%s_%s", tag, cts[[ct]]), r$baseMean_median, "baseMean",
        sprintf("Median baseMean of %s full-cohort truth genes with |LFC| 0.2-1, %s", g, ct), src)
  }
}

# --- 4. Recall / precision grid (method-pooled, cells mode, all genes) --------
src <- "Paper_Source_Data/recall_precision_grid.csv"
grid <- rd(src) %>%
  filter(method == "all_methods", mode == "cells", gene_category == "Total") %>%
  mutate(cells = suppressWarnings(as.numeric(y_value)),          # NA = all available cells ("Max")
         bin = sub(" .*", "", lfc_bin))
g_oligo <- function(b) filter(grid, cell_type == "Oligo", bin == b)
cell_at <- function(b, d, c, col = "mean_recall") {
  x <- g_oligo(b) %>% filter(n_donors == d, cells == c) %>% pull(.data[[col]])
  stopifnot(length(x) == 1)
  x
}
min_donors <- function(df, col, thr) {
  d <- df %>% group_by(n_donors) %>% summarise(m = max(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    filter(m >= thr)
  if (nrow(d) == 0) NA_real_ else min(d$n_donors)
}
GRID_DEF <- "oligodendrocytes, method-pooled (all_methods), cells-per-donor mode, all truth genes"

add("recall_large_min_donors_80_c1500_oligo",
    min_donors(g_oligo("Large") %>% filter(cells >= 1500), "mean_recall", 0.8), "donors",
    paste("Fewest donors at which large-effect (|LFC| > 3) recall reaches 0.80 at any cells-per-donor level >= 1,500;", GRID_DEF), src)
add("recall_large_d50_c500_oligo", 100 * cell_at("Large", 50, 500), "percent",
    paste("Large-effect recall at 50 donors x 500 cells per donor;", GRID_DEF), src)
add("recall_large_d500_c10_oligo", 100 * cell_at("Large", 500, 10), "percent",
    paste("Large-effect recall at 500 donors x 10 cells per donor;", GRID_DEF), src)
add("recall_medium_min_donors_80_c750_oligo",
    min_donors(g_oligo("Medium") %>% filter(cells >= 750), "mean_recall", 0.8), "donors",
    paste("Fewest donors at which medium-effect (|LFC| 2-3) recall reaches 0.80 at any cells-per-donor level >= 750;", GRID_DEF), src)
add("medium_mean_truth_genes_oligo", mean(g_oligo("Medium")$ref_gt_count), "genes",
    paste("Mean number of medium-effect truth genes per grid cell;", GRID_DEF), src)
add("recall_max_lfc_lt2", 100 * max(grid$mean_recall[grid$bin %in% c("Tiny", "Small")]), "percent",
    "Highest recall of any tiny or small-effect (|LFC| < 2) grid cell, any cell type; method-pooled, cells mode, all truth genes", src)
add("recall_tiny_min_donors_10_oligo", min_donors(g_oligo("Tiny"), "mean_recall", 0.10), "donors",
    paste("Fewest donors at which tiny-effect (|LFC| 0.2-1) recall reaches 0.10 at any cells-per-donor level (including all cells);", GRID_DEF), src)
add("recall_tiny_max_d200_oligo", 100 * max(g_oligo("Tiny")$mean_recall[g_oligo("Tiny")$n_donors == 200]), "percent",
    paste("Highest tiny-effect recall at 200 donors over all cells-per-donor levels;", GRID_DEF), src)
add("recall_tiny_max_oligo", 100 * max(g_oligo("Tiny")$mean_recall), "percent",
    paste("Highest tiny-effect recall anywhere in the grid;", GRID_DEF), src)

lg <- filter(grid, bin == "Large", n_donors > 5)
add("precision_large_min_d_gt5", 100 * min(lg$mean_precision, na.rm = TRUE), "percent",
    "Lowest large-effect precision at more than 5 donors, any cell type; method-pooled, cells mode, all truth genes", src)
add("precision_large_pct_cells_gt90_d_gt5", 100 * mean(lg$mean_precision > 0.9, na.rm = TRUE), "percent",
    "Percentage of large-effect grid cells (more than 5 donors, all three cell types) with precision above 0.90", src)
add("precision_medium_min_donors_80_oligo", min_donors(g_oligo("Medium"), "mean_precision", 0.8), "donors",
    paste("Fewest donors at which medium-effect precision reaches 0.80 at any cells-per-donor level;", GRID_DEF), src)
add("precision_small_min_donors_80_oligo", min_donors(g_oligo("Small"), "mean_precision", 0.8), "donors",
    paste("Fewest donors at which small-effect (|LFC| 1-2) precision reaches 0.80 at any cells-per-donor level;", GRID_DEF), src)
add("precision_tiny_max_oligo", 100 * max(g_oligo("Tiny")$mean_precision, na.rm = TRUE), "percent",
    paste("Highest tiny-effect precision anywhere in the grid;", GRID_DEF), src)
add("precision_tiny_n_cells_80_oligo", sum(g_oligo("Tiny")$mean_precision >= 0.8, na.rm = TRUE), "grid cells",
    paste("Number of tiny-effect grid cells with precision >= 0.80;", GRID_DEF), src)

# --- 5. Depth / cell-count equalisation (astrocytes, |LFC| > 0.2) -------------
src <- "Depth_Eq_DEG_Loss/deg_loss_table.csv"
dt <- rd(src)
if ("Cell type" %in% names(dt)) dt <- dt[dt[["Cell type"]] == "Astrocytes", ]
pct_of_native <- function(config, column) {
  s <- dt[dt$Configuration == config, column]
  stopifnot(length(s) == 1)
  as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", s))
}
for (m in c("DESeq2", "edgeR")) {
  add(paste0("deg_loss_pct_full_eq_", tolower(m)),
      100 - pct_of_native("Depth + cell equalised", sprintf("%s (|LFC| > 0.2)", m)), "percent",
      sprintf("Astrocyte truth DEGs (|LFC| > 0.2) lost under depth + cell-count equalisation, %s, %% of native", m), src)
}
add("deg_limma_pct_native_cell_eq", pct_of_native("Cell-count-equalised", "limma (|LFC| > 0.2)"), "percent",
    "limma-voom astrocyte truth DEGs (|LFC| > 0.2) under cell-count equalisation, % of native", src)
add("deg_limma_pct_native_full_eq", pct_of_native("Depth + cell equalised", "limma (|LFC| > 0.2)"), "percent",
    "limma-voom astrocyte truth DEGs (|LFC| > 0.2) under depth + cell-count equalisation, % of native", src)

src <- "R/depth_equalization/config.R"
cfg <- new.env()
old <- Sys.getenv("DEPTH_EQ_CELL_TYPE", unset = NA)
Sys.setenv(DEPTH_EQ_CELL_TYPE = "Astro")
sys.source(project_path(src), envir = cfg)
if (is.na(old)) Sys.unsetenv("DEPTH_EQ_CELL_TYPE") else Sys.setenv(DEPTH_EQ_CELL_TYPE = old)
add("depth_eq_factor_astro", cfg$DOWNSAMPLE_FRACTION, "proportion",
    "Multinomial downsampling fraction applied to astrocyte UMIs (microglia mean / astrocyte mean UMIs per cell)", src)

# --- 6. XGBoost power model -----------------------------------------------------
xgb <- "Power_Curve_Fitting/Models/xgboost"
src <- file.path(xgb, "model_selection.csv")
ms <- rd(src) %>% arrange(bin_balanced_MAE)
add("xgb_n_configs", nrow(ms), "count", "Hyperparameter configurations searched", src)
for (p in c("alpha", "max_depth", "eta", "min_child_weight")) {
  add(paste0("xgb_winner_", p), ms[[p]][1], "hyperparameter",
      sprintf("%s of the selected configuration (lowest out-of-fold bin-balanced MAE)", p), src)
}
add("xgb_mae_worst", max(ms$bin_balanced_MAE), "MAE",
    "Highest out-of-fold bin-balanced MAE of any configuration", src)

src <- file.path(xgb, "calibration_residuals.csv")
cr <- rd(src)
lo <- as.numeric(sub("^[\\[(]([0-9.]+),.*", "\\1", cr$pred_bin))
hi <- as.numeric(sub(".*,([0-9.]+)\\]$", "\\1", cr$pred_bin))
mid <- !is.na(lo) & lo >= 0.2 - 1e-9 & hi <= 0.8 + 1e-9
add("xgb_within010_pct_bin_0_01", 100 * cr$pct_within_010[cr$pred_bin == "[0,0.1]"], "percent",
    "Selected model: % of out-of-fold predictions in [0, 0.1] within 0.10 of the observed detection rate", src)
add("xgb_within010_pct_bin_09_1", 100 * cr$pct_within_010[cr$pred_bin == "(0.9,1]"], "percent",
    "Selected model: % of out-of-fold predictions in (0.9, 1] within 0.10 of the observed detection rate", src)
add("xgb_within010_pct_mid_min", 100 * min(cr$pct_within_010[mid]), "percent",
    "Selected model: lowest % within 0.10 over the prediction bins between 0.2 and 0.8", src)
add("xgb_within010_pct_mid_max", 100 * max(cr$pct_within_010[mid]), "percent",
    "Selected model: highest % within 0.10 over the prediction bins between 0.2 and 0.8", src)
add("xgb_mae_mid_min", min(cr$mae[mid]), "MAE",
    "Selected model: lowest MAE over the prediction bins between 0.2 and 0.8", src)
add("xgb_mae_mid_max", max(cr$mae[mid]), "MAE",
    "Selected model: highest MAE over the prediction bins between 0.2 and 0.8", src)

src <- file.path(xgb, "comparison", "unified_calibration_table.csv")
uc <- rd(src)
ulo <- as.numeric(sub("^[\\[(]([0-9.]+),.*", "\\1", uc$pred_bin))
uhi <- as.numeric(sub(".*,([0-9.]+)\\]$", "\\1", uc$pred_bin))
umid <- !is.na(ulo) & ulo >= 0.2 - 1e-9 & uhi <= 0.8 + 1e-9 & !is.na(uc$alpha) & uc$alpha > 0
add("xgb_reweighted_max_underprediction_pp", 100 * max(uc$mean_obs[umid] - uc$mean_pred[umid]), "percentage points",
    "Largest gap between observed and predicted detection rate (observed minus predicted) in any prediction bin between 0.2 and 0.8, over all configurations with reweighting exponent alpha > 0", src)

src <- "Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data.qs"
fd <- qs::qread(dl(src))
r <- fd$detection_rate
eps <- 1e-9
for (b in list(c(0.1, 0.9), c(0.2, 0.8), c(0.3, 0.7))) {
  add(sprintf("pct_rows_detection_%02d_%02d", round(10 * b[1]), round(10 * b[2])),
      100 * mean(r > b[1] + eps & r < b[2] - eps), "percent",
      sprintf("%% of power-model training rows (n_cells >= 100) whose observed detection rate lies strictly between %.1f and %.1f", b[1], b[2]),
      src)
}
rm(fd)

# --- 7. Gene groups (Figure 6), oligodendrocytes, 50 donors x 1,000 cells ------
src <- "Gene_Classification/oligo_gene_group_power_summary.csv"
gs <- rd(src)
sets <- list(
  easy = c("Chromatin regulator", "Ubiquitin ligase", "LoF-intolerant"),
  hard = c("Intermediate filament", "Forkhead", "Cytokine")
)
for (s in names(sets)) {
  for (l in c(2, 1)) {
    v <- gs$pct_detectable[gs$group %in% sets[[s]] & abs(gs$abs_lfc - l) < 1e-9]
    stopifnot(length(v) == length(sets[[s]]))
    lbl <- paste(sets[[s]], collapse = ", ")
    add(sprintf("genegroup_%s_lfc%d_pct_detectable_min", s, l), min(v), "percent",
        sprintf("Lowest %% Detectable (power > 0.8) among expressed genes over the groups {%s} at |LFC| = %d", lbl, l), src)
    add(sprintf("genegroup_%s_lfc%d_pct_detectable_max", s, l), max(v), "percent",
        sprintf("Highest %% Detectable (power > 0.8) among expressed genes over the groups {%s} at |LFC| = %d", lbl, l), src)
  }
}
add("genegroup_any_lfc02_pct_detectable_max", max(gs$pct_detectable[abs(gs$abs_lfc - 0.2) < 1e-9], na.rm = TRUE),
    "percent", "Highest % Detectable over all Figure 6 groups at |LFC| = 0.2", src)

# --- write ----------------------------------------------------------------------
out <- bind_rows(values)
out_path <- file.path(OUTPUT_BASE, "paper_values.csv")
write.csv(out, out_path, row.names = FALSE)
message(sprintf("Wrote %d paper values to %s", nrow(out), out_path))
