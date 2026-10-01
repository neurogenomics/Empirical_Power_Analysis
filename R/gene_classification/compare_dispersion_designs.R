# ==============================================================================
# SCRIPT: compare_dispersion_designs.R
# PURPOSE: Empirically settle whether DESeq2 dispersion is meaningfully
#          design-dependent for the Oligo transcriptome -- i.e. does fitting the
#          paper's SEX design (~ covariates + age_scaled + sex) give a materially
#          different per-gene dispersion from an intercept-only (~1) fit?
#
#          This decides whether the power calculator can legitimately ingest a
#          design-FREE dispersion read off any atlas (~1 fit) and feed it to the
#          model that was TRAINED on sex-design dispersion, with negligible error
#          (=> no retrain needed), or whether the two diverge enough to require a
#          design-free retrain of the XGBoost model.
#
# METHOD:
#   1. Reuse the cached sex-design features from classify_oligo_genes.R:
#        Datalink/Gene_Classification/oligo_gene_features_checkpoint.qs
#      (raw basemean + DESeq2 dispersion under ~ covariates + age_scaled + sex).
#   2. Re-aggregate the SAME Oligo/Roussos pseudobulk, apply the SAME
#      MIN_GENE_COUNT filter, and fit an INTERCEPT-ONLY (~1) DESeq2 dispersion
#      (estimateSizeFactors -> estimateDispersions -> dispersions()).
#   3. Join by gene_id; compare on log10 dispersion:
#        Pearson r, Spearman rho, median |log2 ratio|, and per-gene ratio.
#   4. Propagate to POWER: predict at the fixed design (100 donors, 250 cells,
#      5742 reads/cell, |LFC|=1.5, Oligo, deseq2) under BOTH dispersions, classify
#      with the skewed trinary scheme, count class changes.
#
# OUTPUT: Datalink/Gene_Classification/oligo_dispersion_design_comparison.csv
#         (chromosome / TF stratification is done downstream in Python, where the
#          HGNC + Lambert maps already live.)
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(BiocParallel)
  library(xgboost)
})

if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
source(file.path(Sys.getenv("POWER_ANALYSIS_ROOT"),
                 "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

CELL_TYPE    <- "Oligo"
DATASET      <- "Roussos"
FIXED_LFC    <- 1.5
FIXED_DONORS <- 100
FIXED_CELLS  <- 250
FIXED_READS  <- 5742
N_CORES_USE  <- as.integer(Sys.getenv("CLASSIFY_CORES", "16"))

out_dir    <- project_path("Datalink", "Gene_Classification")
ckpt_file  <- file.path(out_dir, "oligo_gene_features_checkpoint.qs")
disp1_ckpt <- file.path(out_dir, "oligo_intercept_dispersion_checkpoint.qs")

log_section("Dispersion design comparison: ~1 vs ~sex (Oligo transcriptome)")

if (!file.exists(ckpt_file)) {
  stop(sprintf("Sex-design feature checkpoint missing: %s\n  Run classify_oligo_genes.R first.", ckpt_file))
}
feat <- qs::qread(ckpt_file)
gene_universe <- feat$gene_universe
basemean_raw  <- feat$basemean_raw          # raw pseudobulk mean (method-neutral)
disp_sex      <- feat$dispersion            # DESeq2 dispersion under SEX design
log_info("Loaded sex-design features: %d genes (%d with finite dispersion)",
         length(gene_universe), sum(is.finite(disp_sex)))

cfg      <- DATASET_CONFIG[[DATASET]]
de_param <- setup_parallel(N_CORES_USE)

# --- Intercept-only (~1) dispersion on the SAME pseudobulk -------------------
if (file.exists(disp1_ckpt) && !identical(tolower(Sys.getenv("CLASSIFY_FORCE_REFIT")), "true")) {
  log_info("[CHECKPOINT] Loading cached ~1 dispersion: %s", disp1_ckpt)
  disp_int <- qs::qread(disp1_ckpt)$dispersion_intercept
} else {
  shard_dir <- get_shard_dir(CELL_TYPE, DATASET)
  log_info("Aggregating raw pseudobulk from: %s", shard_dir)
  vars_requested <- unique(c(cfg$cols$donor, cfg$cols$sex, cfg$cols$age,
                             cfg$covariates, cfg$numeric_covariates))
  pb_sce <- aggregate_shards_parallel(
    shard_dir   = shard_dir,
    design_vars = vars_requested,
    n_cores     = min(N_CORES_USE, 8)
  )
  full_counts <- SummarizedExperiment::assay(pb_sce, "counts")

  # SAME low-count filter as run_deseq_for_dataset (rowSums >= MIN_GENE_COUNT).
  keep <- Matrix::rowSums(full_counts) >= MIN_GENE_COUNT
  log_info("Intercept fit gene filter: %d -> %d (>= %d counts)",
           nrow(full_counts), sum(keep), MIN_GENE_COUNT)
  counts_keep <- as.matrix(full_counts[keep, , drop = FALSE])

  # Intercept-only DESeq2 dispersion: geneEst -> trend fit -> MAP shrinkage.
  # (Same estimator chain as the sex fit; ONLY the design differs -> isolates the
  #  design effect on dispersion.)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_keep,
    colData   = S4Vectors::DataFrame(dummy = factor(rep("a", ncol(counts_keep)))),
    design    = ~ 1
  )
  log_info("Estimating ~1 size factors + dispersions (the slow step)...")
  dds1 <- DESeq2::estimateSizeFactors(dds1)
  dds1 <- DESeq2::estimateDispersions(dds1, fitType = "parametric")
  disp_int_fit <- DESeq2::dispersions(dds1)
  names(disp_int_fit) <- rownames(dds1)

  disp_int <- setNames(rep(NA_real_, length(gene_universe)), gene_universe)
  disp_int[names(disp_int_fit)] <- disp_int_fit
  qs::qsave(list(dispersion_intercept = disp_int, design = "intercept ~1",
                 min_gene_count = MIN_GENE_COUNT), disp1_ckpt)
  log_info("[CHECKPOINT] Saved ~1 dispersion: %s", disp1_ckpt)
}

# --- Assemble comparison frame ----------------------------------------------
cmp <- data.frame(
  gene_id    = gene_universe,
  basemean   = as.numeric(basemean_raw),
  disp_sex   = as.numeric(disp_sex),
  disp_int   = as.numeric(disp_int),
  stringsAsFactors = FALSE
)
both <- is.finite(cmp$disp_sex) & cmp$disp_sex > 0 &
        is.finite(cmp$disp_int) & cmp$disp_int > 0 &
        is.finite(cmp$basemean) & cmp$basemean > 0
log_info("Genes with finite dispersion under BOTH designs: %d", sum(both))

# --- Headline agreement stats (on log10 dispersion) -------------------------
ls_sex <- log10(cmp$disp_sex[both]); ls_int <- log10(cmp$disp_int[both])
pear <- cor(ls_sex, ls_int, method = "pearson")
spear <- cor(cmp$disp_sex[both], cmp$disp_int[both], method = "spearman")
log2ratio <- log2(cmp$disp_int[both] / cmp$disp_sex[both])   # int relative to sex
cmp$log2_ratio_int_over_sex <- NA_real_
cmp$log2_ratio_int_over_sex[both] <- log2ratio

log_section("Dispersion agreement (~1 vs ~sex)")
log_info("Pearson r (log10 disp)        : %.4f", pear)
log_info("Spearman rho (raw disp)       : %.4f", spear)
log_info("median |log2 ratio| (int/sex) : %.4f", median(abs(log2ratio)))
log_info("mean   log2 ratio  (int/sex)  : %.4f (>0 => ~1 disp higher)", mean(log2ratio))
qs_probs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
qv <- quantile(log2ratio, qs_probs)
for (i in seq_along(qs_probs)) log_info("  log2ratio p%02.0f = %+.4f", 100*qs_probs[i], qv[i])
log_info("Genes with |log2 ratio| > 1 (>2x diff): %d (%.2f%%)",
         sum(abs(log2ratio) > 1), 100*mean(abs(log2ratio) > 1))
log_info("Genes with |log2 ratio| > log2(1.5) (>1.5x): %d (%.2f%%)",
         sum(abs(log2ratio) > log2(1.5)), 100*mean(abs(log2ratio) > log2(1.5)))

# --- Propagate to power + trinary class -------------------------------------
bundle_dir  <- project_path("Datalink", "Power_Curve_Fitting", "Exports", "xgboost", "realistic")
source(file.path(bundle_dir, "predict_power.R"))
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")

classify_trinary <- function(p) {
  cls <- rep(NA_character_, length(p))
  cls[!is.na(p) & p <= 0.4]           <- "Undetectable"
  cls[!is.na(p) & p > 0.4 & p <= 0.8] <- "Unsure"
  cls[!is.na(p) & p > 0.8]            <- "Detectable"
  cls
}
predict_at_fixed <- function(disp_vec, idx) {
  predict_power(
    n_donors = FIXED_DONORS, abs_lfc = FIXED_LFC,
    basemean = cmp$basemean[idx], dispersion = disp_vec[idx],
    n_cells = FIXED_CELLS, reads_per_cell = FIXED_READS,
    cell_type = CELL_TYPE, de_method = "deseq2",
    bundle_path = bundle_path, warn_extrapolation = FALSE
  )
}
idx <- which(both)
cmp$pred_power_sex <- NA_real_; cmp$pred_power_int <- NA_real_
cmp$pred_power_sex[idx] <- predict_at_fixed(cmp$disp_sex, idx)
cmp$pred_power_int[idx] <- predict_at_fixed(cmp$disp_int, idx)
cmp$class_sex <- classify_trinary(cmp$pred_power_sex)
cmp$class_int <- classify_trinary(cmp$pred_power_int)

log_section("Power impact of swapping sex-design -> intercept dispersion")
dp <- cmp$pred_power_int[idx] - cmp$pred_power_sex[idx]
log_info("Predicted-power delta (int - sex): mean %+.4f | median %+.4f | p95 |Δ| %.4f",
         mean(dp), median(dp), quantile(abs(dp), 0.95))
changed <- cmp$class_sex[idx] != cmp$class_int[idx]
log_info("Trinary class CHANGES: %d / %d (%.2f%%)",
         sum(changed, na.rm = TRUE), length(idx),
         100*mean(changed, na.rm = TRUE))
ct <- table(sex = cmp$class_sex[idx], intercept = cmp$class_int[idx])
log_info("Confusion (rows=sex design, cols=intercept):")
print(ct)

out_csv <- file.path(out_dir, "oligo_dispersion_design_comparison.csv")
write.csv(cmp, out_csv, row.names = FALSE)
log_info("Wrote %s", out_csv)
log_info("Agreement summary: r=%.4f rho=%.4f med|log2ratio|=%.4f class_change=%.2f%%",
         pear, spear, median(abs(log2ratio)), 100*mean(changed, na.rm = TRUE))
