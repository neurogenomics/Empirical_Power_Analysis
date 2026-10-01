# ==============================================================================
# SCRIPT: classify_celltype_genes.R
# PURPOSE: Same as classify_oligo_genes.R but for Astro + Micro. Classify EVERY
#          gene (full transcriptome) into a skewed-trinary detectability call at
#          a FIXED design + FIXED |LFC|=1.5, using the realistic XGBoost model.
#          Per-gene inputs match training EXACTLY:
#            gt_basemean   = raw pseudobulk mean = rowMeans(counts(dds))
#            gt_dispersion = DESeq2::dispersions(dds)  (full-cohort SEX design fit)
#          de_method one-hot toggled for deseq2/edger/limma.
#
# Fixed design: 100 donors, 250 cells/donor, 5742 reads/cell, |LFC|=1.5.
#   NOTE reads_per_cell defaults to the model export's per-cell-type mean; we keep
#   5742 (the Oligo/export default) to hold design identical across cell types so
#   the ONLY thing that differs between Astro/Micro/Oligo results is the gene's
#   own basemean/dispersion + the cell_type one-hot. (Matches how the Oligo run
#   was done.)
#
# OUTPUT (per cell type, lower-cased):
#   Datalink/Gene_Classification/<ct>_gene_power_classification.csv
#   Datalink/Gene_Classification/<ct>_gene_features_checkpoint.qs
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

if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
source(file.path(Sys.getenv("POWER_ANALYSIS_ROOT"),
                 "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

DATASET      <- "Roussos"
FIXED_LFC    <- 1.5
FIXED_DONORS <- 100
FIXED_CELLS  <- 250
FIXED_READS  <- 5742
DE_METHODS   <- c("deseq2", "edger", "limma")
N_CORES_USE  <- as.integer(Sys.getenv("CLASSIFY_CORES", "16"))
CELL_TYPES   <- c("Astro", "Micro")

bundle_dir <- project_path("Datalink", "Power_Curve_Fitting",
                           "Exports", "xgboost", "realistic")
out_dir    <- project_path("Datalink", "Gene_Classification")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

source(file.path(bundle_dir, "predict_power.R"))
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")
bundle <- readRDS(bundle_path)
bm_rng <- bundle$training_ranges$gt_basemean
dp_rng <- bundle$training_ranges$gt_dispersion

classify_trinary <- function(p) {
  cls <- rep(NA_character_, length(p))
  cls[!is.na(p) & p <= 0.4]           <- "Undetectable"
  cls[!is.na(p) & p > 0.4 & p <= 0.8] <- "Unsure"
  cls[!is.na(p) & p > 0.8]            <- "Detectable"
  cls
}

cfg      <- DATASET_CONFIG[[DATASET]]
de_param <- setup_parallel(N_CORES_USE)

for (CELL_TYPE in CELL_TYPES) {
  ct_lc     <- tolower(CELL_TYPE)
  ckpt_file <- file.path(out_dir, sprintf("%s_gene_features_checkpoint.qs", ct_lc))
  out_csv   <- file.path(out_dir, sprintf("%s_gene_power_classification.csv", ct_lc))

  log_section(sprintf("Classify %s genes: trinary detectability at fixed design", CELL_TYPE))

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

    log_info("Running full-cohort DESeq2 (sex design) -- slow step...")
    de <- run_deseq_for_dataset(pb_sce = pb_sce, cfg = cfg, de_param = de_param,
                                min_gene_count = MIN_GENE_COUNT)
    dds <- de$dds

    bm_fit   <- as.numeric(Matrix::rowMeans(DESeq2::counts(dds)))
    disp_fit <- DESeq2::dispersions(dds)
    fit_genes <- rownames(dds)
    log_info("Fitted (post-filter) genes: %d ; NA dispersions: %d",
             length(fit_genes), sum(is.na(disp_fit)))

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

  scorable <- is.finite(gene_df$basemean) & gene_df$basemean > 0 &
              is.finite(gene_df$dispersion) & gene_df$dispersion > 0
  log_info("Scorable genes: %d / %d", sum(scorable), nrow(gene_df))

  gene_df$basemean_in_range   <- gene_df$basemean   >= bm_rng["min"] & gene_df$basemean   <= bm_rng["max"]
  gene_df$dispersion_in_range <- gene_df$dispersion >= dp_rng["min"] & gene_df$dispersion <= dp_rng["max"]

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
    cls[!scorable] <- "Undetectable"
    gene_df[[paste0("pred_power_", m)]]    <- pred
    gene_df[[paste0("trinary_class_", m)]] <- cls
  }

  gene_df$pred_power    <- gene_df$pred_power_deseq2
  gene_df$trinary_class <- gene_df$trinary_class_deseq2

  gene_df <- gene_df[order(-gene_df$pred_power), ]
  col_order <- c("gene_id", "basemean", "dispersion", "n_donors_expressing",
                 "pred_power", "trinary_class",
                 "basemean_in_range", "dispersion_in_range",
                 "pred_power_deseq2", "trinary_class_deseq2",
                 "pred_power_edger",  "trinary_class_edger",
                 "pred_power_limma",  "trinary_class_limma")
  write.csv(gene_df[, col_order], out_csv, row.names = FALSE)

  for (m in DE_METHODS) {
    tab <- table(factor(gene_df[[paste0("trinary_class_", m)]],
                        levels = c("Detectable","Unsure","Undetectable")))
    log_info("  [%s][%s]  Detectable %d | Unsure %d | Undetectable %d",
             CELL_TYPE, m, tab[["Detectable"]], tab[["Unsure"]], tab[["Undetectable"]])
  }
  log_info("Wrote %s", out_csv)
}

log_section("Astro + Micro classification complete")
