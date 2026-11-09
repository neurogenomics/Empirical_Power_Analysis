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
#   AD_vs_Control/Full_Cohort_DE/sex_de_results.csv                    (ad_de)
#   Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_all.qs    (pcf_aggregate)
#   Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv          (crossscore)
#   AD_vs_Control/Full_Cohort_DE/ad_vs_control_de_results.csv          (ad_de)
#   Cartesian_Downsampling/<CT>/Roussos_ADvsCtrl/cells/all_methods/Plots/Aggregated_Summary_Data.qs (ad_cartesian)
# Also writes Paper_Source_Data/de_package_agreement.csv (Supplementary Table).
#   R/depth_equalization/config.R (the downsampling factor actually used)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(tidyr)
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
# The 80% frontier for medium effects runs from many donors with few cells to
# fewer donors with many cells (pre-rerun: 300 x 500, 200 x 750, 150 x >= 1,000).
# The text quotes the middle point, so the design is fixed at 750 cells per donor.
med_c750 <- g_oligo("Medium") %>% filter(cells == 750, mean_recall >= 0.8)
add("recall_medium_min_donors_80_c750_oligo",
    if (nrow(med_c750) == 0) NA_real_ else min(med_c750$n_donors), "donors",
    paste("Fewest donors at which medium-effect (|LFC| 2-3) recall reaches 0.80 at exactly 750 cells per donor (middle point of the 80% donors x cells frontier);", GRID_DEF), src)
add("recall_medium_at_min_donors_c750_oligo",
    if (nrow(med_c750) == 0) NA_real_ else 100 * med_c750$mean_recall[which.min(med_c750$n_donors)], "percent",
    paste("Medium-effect recall at that design (fewest donors reaching 0.80 at 750 cells per donor);", GRID_DEF), src)
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
donor_floor <- function(df, col, thr) {
  # Smallest donor count D such that every grid cell with n_donors >= D has col > thr.
  ok <- df %>% group_by(n_donors) %>% summarise(ok = all(.data[[col]] > thr, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n_donors))
  bad <- ok$n_donors[!ok$ok]
  above <- ok$n_donors[ok$ok & (length(bad) == 0 | ok$n_donors > max(c(bad, -Inf)))]
  if (length(above) == 0) NA_real_ else min(above)
}
for (ct in names(cts)) {
  add(paste0("precision_large_floor_90_", cts[[ct]]),
      donor_floor(filter(grid, bin == "Large", cell_type == ct), "mean_precision", 0.9), "donors",
      sprintf("Smallest donor count D such that large-effect (|LFC| > 3) precision is above 0.90 at every design with >= D donors; %s, method-pooled, cells mode, all truth genes", ct), src)
}
add("precision_large_floor_90_all",
    donor_floor(filter(grid, bin == "Large"), "mean_precision", 0.9), "donors",
    "Smallest donor count D such that large-effect (|LFC| > 3) precision is above 0.90 at every design with >= D donors in all three cell types; method-pooled, cells mode, all truth genes", src)
add("precision_medium_min_donors_80_oligo", min_donors(g_oligo("Medium"), "mean_precision", 0.8), "donors",
    paste("Fewest donors at which medium-effect precision reaches 0.80 at any cells-per-donor level;", GRID_DEF), src)
add("precision_small_min_donors_80_oligo", min_donors(g_oligo("Small"), "mean_precision", 0.8), "donors",
    paste("Fewest donors at which small-effect (|LFC| 1-2) precision reaches 0.80 at any cells-per-donor level;", GRID_DEF), src)
add("precision_tiny_max_oligo", 100 * max(g_oligo("Tiny")$mean_precision, na.rm = TRUE), "percent",
    paste("Highest tiny-effect precision anywhere in the grid;", GRID_DEF), src)
tny <- filter(grid, bin == "Tiny")
add("precision_tiny_max_all", 100 * max(tny$mean_precision, na.rm = TRUE), "percent",
    "Highest tiny-effect (|LFC| 0.2-1) precision of any grid cell in any cell type; method-pooled, cells mode, all truth genes", src)
add("precision_tiny_reaches_80_all", as.numeric(any(tny$mean_precision >= 0.8, na.rm = TRUE)), "flag",
    "1 if any tiny-effect grid cell in any cell type has precision >= 0.80, else 0 (the text claims 0); method-pooled, cells mode, all truth genes", src)
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
add("xgb_reweighted_max_overprediction_pp", 100 * max(uc$mean_pred[umid] - uc$mean_obs[umid]), "percentage points",
    "Largest over-prediction (predicted minus observed detection rate) in any prediction bin between 0.2 and 0.8, over all configurations with reweighting exponent alpha > 0", src)

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

# --- 8. DE package differences (sex contrast; DESeq2, edgeR, limma-voom) ------
# Truth-set rule throughout: padj < 0.05 and |LFC| > 0.2 (as GT_LFC_THRESHOLD).
PKG <- c(deseq2 = "DESeq2", edger = "edgeR", limma = "limma-voom")

# 8a. Full cohort: every post-filter donor, all three packages on the same genes.
src <- "AD_vs_Control/Full_Cohort_DE/sex_de_results.csv"
fc <- rd(src) %>% filter(contrast == "sex")
fc_sig <- fc %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.2)
agreement <- list()
flips <- list()
for (ct in names(cts)) {
  sets <- lapply(names(PKG), function(m) unique(fc_sig$gene_id[fc_sig$cell_type == ct & fc_sig$de_method == m]))
  names(sets) <- names(PKG)
  called <- Reduce(union, sets)
  shared3 <- Reduce(intersect, sets)
  only_limma <- setdiff(sets$limma, union(sets$deseq2, sets$edger))
  for (m in names(PKG)) {
    add(sprintf("depkg_fc_n_%s_%s", m, cts[[ct]]), length(sets[[m]]), "genes",
        sprintf("Full-cohort sex DEGs (padj < 0.05, |LFC| > 0.2) called by %s, %s, all post-filter donors", PKG[[m]], ct), src)
  }
  add(paste0("depkg_fc_pct_shared3_", cts[[ct]]), 100 * length(shared3) / length(called), "percent",
      sprintf("Genes called by all three packages as a %% of genes called by at least one (full-cohort sex DE, %s)", ct), src)
  add(paste0("depkg_fc_n_limma_only_", cts[[ct]]), length(only_limma), "genes",
      sprintf("Full-cohort sex DEGs called by limma-voom and by neither DESeq2 nor edgeR, %s", ct), src)
  add(paste0("depkg_fc_median_count_limma_only_", cts[[ct]]),
      median(fc$base_mean[fc$cell_type == ct & fc$de_method == "limma" & fc$gene_id %in% only_limma]), "count",
      sprintf("Median gene count (median raw pseudobulk count across donors) of the limma-voom-only full-cohort DEGs, %s", ct), src)

  # Direction: LFC sign of every package, over genes called by at least one.
  lfc <- fc %>% filter(cell_type == ct, gene_id %in% called) %>%
    select(gene_id, de_method, log2FoldChange) %>%
    tidyr::pivot_wider(names_from = de_method, values_from = log2FoldChange) %>%
    filter(!is.na(deseq2), !is.na(edger), !is.na(limma)) %>%
    mutate(agree = sign(deseq2) == sign(edger) & sign(edger) == sign(limma),
           all_gt02 = pmin(abs(deseq2), abs(edger), abs(limma)) > 0.2,
           shared3 = gene_id %in% shared3)
  add(paste0("depkg_fc_pct_sign_agree_", cts[[ct]]), 100 * mean(lfc$agree), "percent",
      sprintf("Genes called by at least one package whose LFC has the same sign in all three packages, %% (full-cohort sex DE, %s)", ct), src)
  agreement[[ct]] <- data.frame(
    `Cell type` = ct, `DESeq2 DEGs` = length(sets$deseq2), `edgeR DEGs` = length(sets$edger),
    `limma-voom DEGs` = length(sets$limma), `Called by any` = length(called),
    `Called by all three` = length(shared3), `All three / any (%)` = round(100 * length(shared3) / length(called), 1),
    `Same direction in all three (%)` = round(100 * mean(lfc$agree), 1), check.names = FALSE)
  d <- lfc %>% filter(!agree) %>%
    mutate(dissent = ifelse(sign(deseq2) == sign(edger), "limma", ifelse(sign(deseq2) == sign(limma), "edger", "deseq2")),
           min_abs = pmin(abs(deseq2), abs(edger), abs(limma)))
  flips[[ct]] <- list(shared3 = lfc[lfc$shared3, ], all_gt02 = lfc[lfc$all_gt02, ], d = d)
}
fl_d <- bind_rows(lapply(flips, `[[`, "d"))
add("depkg_fc_pct_sign_agree_shared3_min", 100 * min(sapply(flips, function(x) mean(x$shared3$agree))), "percent",
    "Lowest, over the three cell types, % of genes called by all three packages whose LFC has the same sign in all three (full-cohort sex DE; the text claims 100)", src)
add("depkg_fc_pct_sign_agree_all_lfc_gt02_min", 100 * min(sapply(flips, function(x) mean(x$all_gt02$agree))), "percent",
    "Lowest, over the three cell types, % same LFC sign among called genes where all three packages estimate |LFC| > 0.2 (full-cohort sex DE; the text claims 100)", src)
add("depkg_fc_n_sign_flips", nrow(fl_d), "genes",
    "Called genes (any package, any cell type) whose LFC sign differs between packages, full-cohort sex DE", src)
add("depkg_fc_n_sign_flips_limma", sum(fl_d$dissent == "limma"), "genes",
    "Of those sign flips, the number where limma-voom is the package that disagrees with the other two", src)
add("depkg_fc_max_min_abs_lfc_flips", max(fl_d$min_abs), "|LFC|",
    "Over all sign flips, the largest value of the smallest |LFC| any package estimates for that gene (full-cohort sex DE; < 0.2 means every flip involves a package estimating |LFC| below 0.2)", src)

# 8b. Per run: truth sets of the ten ground-truth halves, per package.
src <- "Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_all.qs"
fa <- qs::qread(dl(src))
gt_n <- as.data.frame(fa) %>% filter(mode == "cells") %>%
  distinct(cell_type, de_method, run_id, gene_id) %>%
  count(cell_type, de_method, run_id, name = "n")
rm(fa)
for (ct in names(cts)) {
  for (m in names(PKG)) {
    v <- gt_n$n[gt_n$cell_type == ct & gt_n$de_method == m]
    stopifnot(length(v) == 10)
    add(sprintf("depkg_run_n_min_%s_%s", m, cts[[ct]]), min(v), "genes",
        sprintf("Smallest %s truth set (padj < 0.05, |LFC| > 0.2) over the ten ground-truth halves, %s", PKG[[m]], ct), src)
    add(sprintf("depkg_run_n_max_%s_%s", m, cts[[ct]]), max(v), "genes",
        sprintf("Largest %s truth set (padj < 0.05, |LFC| > 0.2) over the ten ground-truth halves, %s", PKG[[m]], ct), src)
  }
}

# 8c. Cross-scoring: each package's subsample calls against each package's truth
# set, cells mode, one reference design, pooled over effect-size bins
# (sum of numerators / sum of denominators).
src <- "Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv"
XS_DONORS <- 200; XS_CELLS <- 250
XS_DEF <- sprintf("cells mode, %d donors x %d cells per donor, 10 runs", XS_DONORS, XS_CELLS)
xs <- rd(src) %>% filter(mode == "cells", n_donors == XS_DONORS, as.character(y_value) == as.character(XS_CELLS))
stopifnot(nrow(xs) > 0)
xp <- xs %>% group_by(cell_type, metric, subsample_method, truth_method) %>%
  summarise(v = sum(numerator) / sum(denom_raw), .groups = "drop")
rc <- filter(xp, metric == "recall"); pr <- filter(xp, metric == "precision")
stopifnot(nrow(rc) == 27, nrow(pr) == 27)
add("depkg_xs_recall_spread_within_truth_max",
    max((rc %>% group_by(cell_type, truth_method) %>% summarise(s = max(v) - min(v), .groups = "drop"))$s), "recall",
    paste("Largest difference in recall between the three calling packages scored against the same truth set, any cell type, pooled over effect sizes;", XS_DEF), src)
add("depkg_xs_recall_spread_between_truth_max",
    max((rc %>% group_by(cell_type, subsample_method) %>% summarise(s = max(v) - min(v), .groups = "drop"))$s), "recall",
    paste("Largest difference in one calling package's recall between the three truth sets, any cell type, pooled over effect sizes;", XS_DEF), src)
am <- filter(rc, cell_type %in% c("Astro", "Micro"))
add("depkg_xs_recall_vs_limma_min_astro_micro", min(am$v[am$truth_method == "limma"]), "recall",
    paste("Lowest recall of any package against the limma-voom truth set, astrocytes and microglia, pooled over effect sizes;", XS_DEF), src)
add("depkg_xs_recall_vs_limma_max_astro_micro", max(am$v[am$truth_method == "limma"]), "recall",
    paste("Highest recall of any package against the limma-voom truth set, astrocytes and microglia, pooled over effect sizes;", XS_DEF), src)
add("depkg_xs_recall_vs_other_min_astro_micro", min(am$v[am$truth_method != "limma"]), "recall",
    paste("Lowest recall of any package against the DESeq2 or edgeR truth sets, astrocytes and microglia, pooled over effect sizes;", XS_DEF), src)
add("depkg_xs_recall_vs_other_max_astro_micro", max(am$v[am$truth_method != "limma"]), "recall",
    paste("Highest recall of any package against the DESeq2 or edgeR truth sets, astrocytes and microglia, pooled over effect sizes;", XS_DEF), src)
ol <- filter(rc, cell_type == "Oligo")
add("depkg_xs_recall_vs_limma_max_oligo", max(ol$v[ol$truth_method == "limma"]), "recall",
    paste("Highest recall of any package against the limma-voom truth set, oligodendrocytes, pooled over effect sizes;", XS_DEF), src)
# Claim flags (1 = holds): limma-voom truth easiest for every package in Astro
# and Micro, hardest for every package in Oligo; edgeR calls least precise
# against every truth set in every cell type.
rank_flag <- function(d, fn) all(sapply(split(d, d$subsample_method), function(s) s$truth_method[fn(s$v)] == "limma"))
add("depkg_xs_limma_truth_easiest_astro_micro",
    as.numeric(all(sapply(split(am, am$cell_type), rank_flag, fn = which.max))), "flag",
    paste("1 if, in both astrocytes and microglia, every package has its highest recall against the limma-voom truth set;", XS_DEF), src)
add("depkg_xs_limma_truth_hardest_oligo", as.numeric(rank_flag(ol, which.min)), "flag",
    paste("1 if, in oligodendrocytes, every package has its lowest recall against the limma-voom truth set;", XS_DEF), src)
add("depkg_xs_edger_least_precise_all",
    as.numeric(all(sapply(split(pr, list(pr$cell_type, pr$truth_method)),
                          function(s) s$subsample_method[which.min(s$v)] == "edger"))), "flag",
    paste("1 if edgeR's calls have the lowest precision against every truth set in every cell type, pooled over effect sizes;", XS_DEF), src)
lg <- xs %>% filter(lfc_bin == "Large")
stopifnot(nrow(lg) == 54)
add("depkg_xs_recall_large_min", min(lg$value[lg$metric == "recall"], na.rm = TRUE), "recall",
    paste("Lowest large-effect (|LFC| > 3) recall over all nine calling x truth package pairs and three cell types;", XS_DEF), src)
add("depkg_xs_precision_large_min", min(lg$value[lg$metric == "precision"], na.rm = TRUE), "precision",
    paste("Lowest large-effect (|LFC| > 3) precision over all nine calling x truth package pairs and three cell types;", XS_DEF), src)

# Supplementary table: full-cohort agreement between packages.
write.csv(bind_rows(agreement), file.path(OUTPUT_BASE, "de_package_agreement.csv"), row.names = FALSE)


# --- 9. AD versus control (validation that the sex-contrast trends replicate) --
# Full-cohort AD DE: MSSM donors whose only diagnosis is AD versus donors with no
# diagnosis (Roussos_ADvsCtrl shards). "AD DEGs" here = padj < 0.05 with no LFC
# filter, matching the published PsychAD AD lists (Lee et al. 2025) they are
# compared with in the text; the tiny-bin share uses the truth-set rule.
src <- "AD_vs_Control/Full_Cohort_DE/ad_vs_control_de_results.csv"
ad <- rd(src) %>% filter(contrast == "ad_vs_control")
sx_donors <- rd("AD_vs_Control/Full_Cohort_DE/sex_de_results.csv") %>%
  filter(contrast == "sex") %>% distinct(cell_type, n_donors)
ratios <- c()
for (ct in names(cts)) {
  n <- unique(ad$n_donors[ad$cell_type == ct])
  n_sex <- unique(sx_donors$n_donors[sx_donors$cell_type == ct])
  stopifnot(length(n) == 1, length(n_sex) == 1)
  add(sprintf("ad_n_donors_%s", cts[[ct]]), n, "donors",
      sprintf("Donors in the full-cohort AD versus control analysis, %s: Mount Sinai (MSSM) donors whose only diagnosis is AD, plus donors with no diagnosis, after the cohort filters", ct), src)
  ratios[ct] <- n / n_sex
}
add("ad_n_donors_ratio_to_sex_min", min(ratios), "ratio",
    "Smallest ratio, over cell types, of AD versus control donors to sex-contrast donors (full cohort). The ground-truth half is 50% of each, so this is also the ratio of ground-truth sizes", src)
add("ad_n_donors_ratio_to_sex_max", max(ratios), "ratio",
    "Largest ratio, over cell types, of AD versus control donors to sex-contrast donors (full cohort); see ad_n_donors_ratio_to_sex_min", src)

ad_sig <- ad %>% filter(!is.na(padj), padj < 0.05) %>% mutate(a = abs(log2FoldChange))
eff <- ad_sig %>% group_by(cell_type, de_method) %>%
  summarise(med = median(a), lt1 = 100 * mean(a < 1), .groups = "drop")
stopifnot(nrow(eff) == 9)
add("ad_fc_median_abs_lfc_min", min(eff$med), "log2FC",
    "Lowest median |LFC| of full-cohort AD DEGs (padj < 0.05, no LFC filter) over the three DE packages and three cell types", src)
add("ad_fc_median_abs_lfc_max", max(eff$med), "log2FC",
    "Highest median |LFC| of full-cohort AD DEGs (padj < 0.05, no LFC filter) over the three DE packages and three cell types", src)
add("ad_fc_pct_abs_lfc_lt1_min", min(eff$lt1), "percent",
    "Lowest percentage of full-cohort AD DEGs (padj < 0.05, no LFC filter) with |LFC| < 1, over the three DE packages and three cell types", src)
tiny <- ad_sig %>% filter(a > 0.2) %>% group_by(cell_type, de_method) %>%
  summarise(p = 100 * mean(a <= 1), .groups = "drop")
add("ad_fc_pct_tiny_min", min(tiny$p), "percent",
    "Lowest percentage, over the three DE packages and three cell types, of full-cohort AD DEGs passing the truth-set rule (padj < 0.05, |LFC| > 0.2) that fall in the tiny bin (|LFC| 0.2-1)", src)

# AD downsampling grid (method-pooled, cells-per-donor mode): truth-set sizes per
# effect-size bin and the donor trend in tiny-bin recall with all cells.
med_large <- c()
donor_trend <- c()
for (ct in names(cts)) {
  src <- sprintf("Cartesian_Downsampling/%s/Roussos_ADvsCtrl/cells/all_methods/Plots/Aggregated_Summary_Data.qs", ct)
  s <- qs::qread(dl(src))$df_summary %>% filter(Gene_Category == "Total")
  med_large[ct] <- max(s$Ref_GT_Count[s$LFC_Bin %in% c("Medium (2.0-3.0)", "Large (3.0+)")], na.rm = TRUE)
  t <- s %>% filter(LFC_Bin == "Tiny (0.2-1.0)", Cells == "Max", !is.na(Mean_Recall)) %>% arrange(Donors)
  donor_trend[ct] <- all(t$Mean_Recall[nrow(t)] > t$Mean_Recall[-nrow(t)])
}
src <- "Cartesian_Downsampling/<CT>/Roussos_ADvsCtrl/cells/all_methods/Plots/Aggregated_Summary_Data.qs"
add("ad_gt_medlarge_mean_genes_max", max(med_large), "genes",
    "Largest mean number of AD ground-truth genes (averaged over the ten split-half runs, method-pooled) in the medium (|LFC| 2-3) or large (|LFC| > 3) bin, over cell types and design points", src)
add("ad_recall_tiny_max_donors_highest_all", as.numeric(all(donor_trend)), "flag",
    "1 if, in every cell type, AD tiny-bin recall (method-pooled, all available cells per donor) is higher at the largest donor count than at every smaller donor count", src)

# --- write ----------------------------------------------------------------------
out <- bind_rows(values)
out_path <- file.path(OUTPUT_BASE, "paper_values.csv")
write.csv(out, out_path, row.names = FALSE)
message(sprintf("Wrote %d paper values to %s", nrow(out), out_path))
