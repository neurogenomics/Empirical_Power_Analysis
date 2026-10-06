# ==============================================================================
# SCRIPT: classify_oligo_genes.R
# PURPOSE: Classify EVERY oligodendrocyte gene (full transcriptome) into a
#          skewed-trinary detectability call at a FIXED experimental design and a
#          FIXED effect size, using the XGBoost power model.
#
# DESIGN DECISIONS (locked by Chris, 2026-08-21/22):
#   Gene universe   : full Oligo/Roussos transcriptome (all donor shards).
#   Fixed |LFC|     : 1.5 for every gene.
#   Fixed design    : 100 donors, 250 cells/donor, 5742 reads/cell, Oligo.
#   Trinary scheme  : trinary_04_08 -> breaks (0.4, 0.8)
#                       power <= 0.4        -> "Undetectable"  (Definitely Not)
#                       0.4 < power <= 0.8  -> "Unsure"        (Maybe)
#                       power  > 0.8        -> "Detectable"    (Yes)
#
# PER-GENE MODEL INPUTS  ── CORRECTED 2026-08-22 ──────────────────────────────
#   The XGBoost model trained on TWO gene characteristics, defined EXACTLY as:
#     gt_basemean   = baseMean_raw = Matrix::rowMeans(DESeq2::counts(dds))
#                     (raw pseudobulk mean count; identical across de_methods;
#                      cartesian_downsampling/utils.R:499).
#     gt_dispersion = DESeq2::dispersions(dds)  -- the final MAP dispersion from
#                     a full-cohort DESeq2 fit under the paper's SEX design
#                     (~ covariates + age_scaled + sex); 01_Build_Fitting_Dataset.R:214,226.
#   de_method (deseq2/edger/limma) is a SEPARATE one-hot axis: it flags which
#   method's detection outcome is predicted, and does NOT change basemean/
#   dispersion. Every training row (incl. edger/limma) carried the SAME DESeq2
#   dispersion + raw basemean. So the correct per-method classification is:
#   ONE full-cohort DESeq2 fit -> raw basemean + DESeq2 dispersion, then predict
#   three times toggling the de_method one-hot.
#
#   (This replaces the earlier naive method-of-moments dispersion on raw
#    pseudobulk, which over-stated dispersion by treating all covariate
#    structure as noise -- an intercept-only ~1 fit -- and was off-distribution.)
#
#   Dispersion IS design-dependent: DESeq2 estimates it as the residual
#   overdispersion AFTER the design matrix absorbs what it can. It is ~identical
#   across sex/disease/age designs for most genes, differing only for genes
#   linked to the chosen variable. The model trained under the SEX design, so we
#   condition on the SAME sex design to stay on-distribution.
#
#   Genes below the MIN_GENE_COUNT filter get no DESeq2 dispersion; at |LFC|=1.5
#   they are undetectable anyway, so they are labelled "Undetectable".
#
# OUTPUT: Datalink/Gene_Classification/oligo_gene_power_classification.csv
#         (+ oligo_gene_features_checkpoint.qs -- basemean/dispersion cache)
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(Matrix)
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(BiocParallel)
  library(xgboost)
})

# --- 0. Config + shared DE utils (reuse the pipeline's own code) --------------
if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
source(file.path(Sys.getenv("POWER_ANALYSIS_ROOT"),
                 "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))   # aggregate_shards_parallel(), run_deseq_for_dataset()

CELL_TYPE    <- "Oligo"
DATASET      <- "Roussos"
FIXED_LFC    <- 1.5
FIXED_DONORS <- 100
FIXED_CELLS  <- 250
FIXED_READS  <- 5742
DE_METHODS   <- c("deseq2", "edger", "limma")
N_CORES_USE  <- as.integer(Sys.getenv("CLASSIFY_CORES", "16"))

bundle_dir <- project_path("Datalink", "Power_Curve_Fitting",
                           "Exports", "xgboost")
out_dir    <- project_path("Datalink", "Gene_Classification")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ckpt_file  <- file.path(out_dir, "oligo_gene_features_checkpoint.qs")

log_section("Classify Oligo genes: trinary detectability at fixed design")
log_info("Per-gene features: raw basemean + DESeq2 dispersion (full-cohort sex fit)")

cfg      <- DATASET_CONFIG[[DATASET]]
de_param <- setup_parallel(N_CORES_USE)

# --- 1. Full-cohort features (raw basemean + DESeq2 dispersion) ---------------
# Cache the expensive DESeq2 fit so re-prediction / re-annotation is cheap.
if (file.exists(ckpt_file) && !identical(tolower(Sys.getenv("CLASSIFY_FORCE_REFIT")), "true")) {
  log_info("[CHECKPOINT] Loading cached gene features: %s", ckpt_file)
  feat <- qs::qread(ckpt_file)
  gene_universe <- feat$gene_universe
  basemean_raw  <- feat$basemean_raw
  dispersion    <- feat$dispersion
  n_expr        <- feat$n_donors_expressing
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
  full_counts   <- as.matrix(SummarizedExperiment::assay(pb_sce, "counts"))
  gene_universe <- rownames(full_counts)
  n_expr_all    <- rowSums(full_counts > 0)
  log_info("Full transcriptome: %d genes x %d donors", nrow(full_counts), ncol(full_counts))

  # Full-cohort DESeq2 fit under the paper's sex design (SAME as training GT fit).
  log_info("Running full-cohort DESeq2 (sex design) -- this is the slow step...")
  de <- run_deseq_for_dataset(pb_sce = pb_sce, cfg = cfg, de_param = de_param,
                              min_gene_count = MIN_GENE_COUNT)
  dds <- de$dds

  # baseMean_raw = raw pseudobulk mean over the fitted genes (matches training).
  bm_fit   <- as.numeric(Matrix::rowMeans(DESeq2::counts(dds)))
  disp_fit <- DESeq2::dispersions(dds)
  fit_genes <- rownames(dds)
  log_info("Fitted (post-filter) genes: %d ; NA dispersions: %d",
           length(fit_genes), sum(is.na(disp_fit)))

  # Map onto the full universe (filtered-out / NA-dispersion genes stay NA).
  basemean_raw <- setNames(rep(NA_real_, length(gene_universe)), gene_universe)
  dispersion   <- setNames(rep(NA_real_, length(gene_universe)), gene_universe)
  basemean_raw[fit_genes] <- bm_fit
  dispersion[fit_genes]   <- disp_fit
  n_expr <- setNames(as.integer(n_expr_all), gene_universe)

  qs::qsave(list(gene_universe = gene_universe, basemean_raw = basemean_raw,
                 dispersion = dispersion, n_donors_expressing = n_expr,
                 design = "full-cohort sex fit", min_gene_count = MIN_GENE_COUNT),
            ckpt_file)
  log_info("[CHECKPOINT] Saved gene features: %s", ckpt_file)
}

gene_df <- data.frame(
  gene_id             = gene_universe,
  basemean            = as.numeric(basemean_raw),
  dispersion          = as.numeric(dispersion),
  n_donors_expressing = as.integer(n_expr),
  stringsAsFactors    = FALSE
)

# Scorable = has a finite raw basemean>0 AND a finite DESeq2 dispersion>0.
scorable <- is.finite(gene_df$basemean) & gene_df$basemean > 0 &
            is.finite(gene_df$dispersion) & gene_df$dispersion > 0
log_info("Scorable genes (passed DESeq2 filter w/ finite disp): %d / %d",
         sum(scorable), nrow(gene_df))

# --- 2. Predict power per DE method at the fixed design -----------------------
source(file.path(bundle_dir, "predict_power.R"))
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")
bundle <- readRDS(bundle_path)
bm_rng <- bundle$training_ranges$gt_basemean
dp_rng <- bundle$training_ranges$gt_dispersion

gene_df$basemean_in_range   <- gene_df$basemean   >= bm_rng["min"] & gene_df$basemean   <= bm_rng["max"]
gene_df$dispersion_in_range <- gene_df$dispersion >= dp_rng["min"] & gene_df$dispersion <= dp_rng["max"]

classify_trinary <- function(p) {
  cls <- rep(NA_character_, length(p))
  cls[!is.na(p) & p <= 0.4]           <- "Undetectable"
  cls[!is.na(p) & p > 0.4 & p <= 0.8] <- "Unsure"
  cls[!is.na(p) & p > 0.8]            <- "Detectable"
  cls
}

for (m in DE_METHODS) {
  pred <- rep(NA_real_, nrow(gene_df))
  pred[scorable] <- predict_power(
    n_donors       = FIXED_DONORS,
    abs_lfc        = FIXED_LFC,
    basemean       = gene_df$basemean[scorable],
    dispersion     = gene_df$dispersion[scorable],
    n_cells        = FIXED_CELLS,
    reads_per_cell = FIXED_READS,
    cell_type      = CELL_TYPE,
    de_method      = m,
    bundle_path    = bundle_path,
    warn_extrapolation = FALSE
  )
  cls <- classify_trinary(pred)
  # Genes that failed the DESeq2 count filter are undetectable at |LFC|=1.5.
  cls[!scorable] <- "Undetectable"
  gene_df[[paste0("pred_power_", m)]]   <- pred
  gene_df[[paste0("trinary_class_", m)]] <- cls
}

# Primary columns = DESeq2 (keeps the downstream annotate.py join unchanged).
gene_df$pred_power    <- gene_df$pred_power_deseq2
gene_df$trinary_class <- gene_df$trinary_class_deseq2

# --- 3. Write out ------------------------------------------------------------
gene_df <- gene_df[order(-gene_df$pred_power), ]
col_order <- c("gene_id", "basemean", "dispersion", "n_donors_expressing",
               "pred_power", "trinary_class",
               "basemean_in_range", "dispersion_in_range",
               "pred_power_deseq2", "trinary_class_deseq2",
               "pred_power_edger",  "trinary_class_edger",
               "pred_power_limma",  "trinary_class_limma")
out_csv <- file.path(out_dir, "oligo_gene_power_classification.csv")
write.csv(gene_df[, col_order], out_csv, row.names = FALSE)

log_section("Classification complete")
log_info("Fixed design: %d donors, %d cells, %d reads/cell, |LFC|=%.2f, %s",
         FIXED_DONORS, FIXED_CELLS, FIXED_READS, FIXED_LFC, CELL_TYPE)
for (m in DE_METHODS) {
  tab <- table(factor(gene_df[[paste0("trinary_class_", m)]],
                      levels = c("Detectable","Unsure","Undetectable")))
  log_info("  [%s]  Detectable %d | Unsure %d | Undetectable %d",
           m, tab[["Detectable"]], tab[["Unsure"]], tab[["Undetectable"]])
}
log_info("Extrapolated (deseq2 scorable): basemean %d | dispersion %d",
         sum(scorable & !gene_df$basemean_in_range, na.rm = TRUE),
         sum(scorable & !gene_df$dispersion_in_range, na.rm = TRUE))
log_info("Wrote %s", out_csv)
