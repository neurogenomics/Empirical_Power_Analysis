# ==============================================================================
# SCRIPT: deg_dispersion_table.R
# PURPOSE: Emit a per-gene table for the Oligo transcriptome carrying the REAL
#          sex-contrast DE result (log2FoldChange, padj, baseMean) alongside the
#          DESeq2 dispersion, so we can ask empirically:
#            do DEGs have systematically different dispersion from non-DEGs,
#            and does any difference survive controlling for expression?
#
#          Reuses the pipeline's own full-cohort sex-design fit
#          (run_deseq_for_dataset), which is exactly the ground-truth fit the
#          power model trained on (~ covariates + age_scaled + sex). res_df
#          gives baseMean/log2FoldChange/padj (apeglm-shrunk); dispersions() on
#          the dds gives the MAP dispersion. Raw pseudobulk mean (baseMean_raw)
#          is added to match the model's expression feature.
#
# OUTPUT: Datalink/Gene_Classification/oligo_deg_dispersion_table.csv
#         cols: gene_id, baseMean_raw, baseMean_deseq2, dispersion,
#               log2FoldChange, padj
#         (chromosome / DEG-definition stratification is done downstream in
#          Python where the HGNC location map already lives.)
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
})

if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
source(file.path(Sys.getenv("POWER_ANALYSIS_ROOT"),
                 "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

CELL_TYPE   <- "Oligo"
DATASET     <- "Roussos"
N_CORES_USE <- as.integer(Sys.getenv("CLASSIFY_CORES", "16"))

out_dir <- project_path("Datalink", "Gene_Classification")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

log_section("DEG vs non-DEG dispersion table (Oligo, sex contrast)")

cfg      <- DATASET_CONFIG[[DATASET]]
de_param <- setup_parallel(N_CORES_USE)

shard_dir <- get_shard_dir(CELL_TYPE, DATASET)
log_info("Aggregating raw pseudobulk from: %s", shard_dir)
vars_requested <- unique(c(cfg$cols$donor, cfg$cols$sex, cfg$cols$age,
                           cfg$covariates, cfg$numeric_covariates))
pb_sce <- aggregate_shards_parallel(
  shard_dir   = shard_dir,
  design_vars = vars_requested,
  n_cores     = min(N_CORES_USE, 8)
)

log_info("Running full-cohort DESeq2 (sex design) -- the slow step...")
de     <- run_deseq_for_dataset(pb_sce = pb_sce, cfg = cfg, de_param = de_param,
                                min_gene_count = MIN_GENE_COUNT)
dds    <- de$dds
res_df <- de$res_df

bm_raw   <- as.numeric(Matrix::rowMeans(DESeq2::counts(dds)))
disp_fit <- DESeq2::dispersions(dds)
fit_genes <- rownames(dds)

# res_df is aligned to dds rows (apeglm shrink of the sex coef); guard by name.
res_df$gene_id <- rownames(res_df)
tab <- data.frame(
  gene_id         = fit_genes,
  baseMean_raw    = bm_raw,
  dispersion      = as.numeric(disp_fit),
  stringsAsFactors = FALSE
)
tab$baseMean_deseq2  <- res_df$baseMean[match(tab$gene_id, res_df$gene_id)]
tab$log2FoldChange   <- res_df$log2FoldChange[match(tab$gene_id, res_df$gene_id)]
tab$padj             <- res_df$padj[match(tab$gene_id, res_df$gene_id)]

log_info("Fitted genes: %d ; finite disp: %d ; finite padj: %d",
         nrow(tab), sum(is.finite(tab$dispersion)), sum(is.finite(tab$padj)))

# --- quick in-R headline (full detail / expression-matching done in Python) --
sig     <- is.finite(tab$padj) & tab$padj < FDR_THRESHOLD
strict  <- sig & abs(tab$log2FoldChange) > 1
relaxed <- sig & abs(tab$log2FoldChange) > 0.2
ok      <- is.finite(tab$dispersion) & tab$dispersion > 0
med <- function(x) median(x, na.rm = TRUE)

log_section("Headline dispersion medians (in-R; Python does the controlled test)")
log_info("Tested genes with finite disp: %d", sum(ok))
log_info("padj<%.2f (any LFC): DEG n=%d med disp=%.4f | nonDEG n=%d med disp=%.4f",
         FDR_THRESHOLD, sum(sig & ok), med(tab$dispersion[sig & ok]),
         sum(!sig & ok), med(tab$dispersion[!sig & ok]))
log_info("strict (padj<%.2f,|LFC|>1): DEG n=%d med disp=%.4f | nonDEG n=%d med disp=%.4f",
         FDR_THRESHOLD, sum(strict & ok), med(tab$dispersion[strict & ok]),
         sum(!strict & ok), med(tab$dispersion[!strict & ok]))
log_info("relaxed(padj<%.2f,|LFC|>0.2): DEG n=%d med disp=%.4f | nonDEG n=%d med disp=%.4f",
         FDR_THRESHOLD, sum(relaxed & ok), med(tab$dispersion[relaxed & ok]),
         sum(!relaxed & ok), med(tab$dispersion[!relaxed & ok]))
log_info("median baseMean_raw: DEG(sig)=%.2f | nonDEG=%.2f (expression confound check)",
         med(tab$baseMean_raw[sig & ok]), med(tab$baseMean_raw[!sig & ok]))

out_csv <- file.path(out_dir, "oligo_deg_dispersion_table.csv")
write.csv(tab, out_csv, row.names = FALSE)
log_info("Wrote %s (%d genes)", out_csv, nrow(tab))
