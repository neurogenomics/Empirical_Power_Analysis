#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: classify_oligo_multi_lfc.R
# PURPOSE: Predicted detection power for every oligodendrocyte gene at three
#          effect sizes, under the fixed study design of Figure 6, plus
#          loss-of-function intolerance (gnomAD pLI). Input table of
#          plot_oligo_expr_disp_panels.py.
#
# DESIGN (Figure 6 caption): 50 donors, 1,000 cells per donor, 5,742 counts per
#   cell (PsychAD oligodendrocyte mean), cell type Oligo, |LFC| in {0.2, 1.0, 2.0},
#   each for DESeq2, edgeR and limma-voom. Power comes from the shipped
#   XGBoost bundle (Exports/xgboost/predict_power.R); the
#   unsuffixed power_/class_ columns repeat the DESeq2 values.
#   Classes: power <= 0.4 Undetectable, 0.4-0.8 Unsure, > 0.8 Detectable; genes
#   without a usable baseMean/dispersion are Undetectable with NA power.
#
# GENE SET: every gene of the Oligo universe (oligo_gene_classification_annotated.csv,
#   in that order), then the Lambert TFs that are not expressed in oligodendrocytes
#   (expressed_in_oligo = False; drawn as grey points in Figure 6).
#
# INPUT : Datalink/Gene_Classification/oligo_gene_classification_annotated.csv
#         reference/DatabaseExtract_v_1.01.csv
#         reference/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
# OUTPUT: Datalink/Gene_Classification/oligo_gene_power_multiLFC.csv
#         Datalink/Gene_Classification/oligo_tf_power_multiLFC.csv  (Lambert TFs only)
#
# gnomAD v2.1.1 constraint is per transcript; per gene the most constrained
# transcript (maximum pLI) is kept; lof_intol = pLI >= 0.9 (missing -> False).
#
# Rewritten 2026-10-01: the previous version of this script used a stale design
# (100 donors, 250 cells, |LFC| 0.75/1.5/2.25); the figure was produced by the
# lost /tmp/oligo_annot/predict_multi_lfc.py, which this reproduces.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("POWER_ANALYSIS_ROOT", unset = "")
if (!nzchar(root)) root <- here::here()
ref_dir <- Sys.getenv("EPA_REFERENCE_DIR", unset = file.path(root, "reference"))
gc_dir  <- file.path(root, "Datalink", "Gene_Classification")
bundle_dir  <- file.path(root, "Datalink", "Power_Curve_Fitting", "Exports", "xgboost")

need <- function(p) {
  if (!file.exists(p)) stop("missing input: ", p, call. = FALSE)
  p
}

CELL_TYPE    <- "Oligo"
FIXED_DONORS <- 50
FIXED_CELLS  <- 1000
FIXED_READS  <- 5742
DE_METHODS   <- c("deseq2", "edger", "limma")
LFCS         <- c("020" = 0.2, "100" = 1.0, "200" = 2.0)

source(need(file.path(bundle_dir, "predict_power.R")))
bundle_path <- need(file.path(bundle_dir, "power_model_bundle.rds"))

ann <- fread(need(file.path(gc_dir, "oligo_gene_classification_annotated.csv")))
gene_df <- ann[, .(gene_id, symbol, basemean, dispersion,
                   expressed_in_oligo = TRUE, is_tf_lambert, lambert_dbd,
                   panther_protein_class)]

# Lambert TFs absent from the oligodendrocyte universe
lam <- fread(need(file.path(ref_dir, "DatabaseExtract_v_1.01.csv")),
             select = c("Ensembl ID", "HGNC symbol", "DBD", "Is TF?"))
setnames(lam, c("gene_id", "symbol", "lambert_dbd", "is_tf"))
lam <- lam[is_tf == "Yes" & !duplicated(gene_id) & !gene_id %in% gene_df$gene_id]
gene_df <- rbind(gene_df,
                 lam[, .(gene_id, symbol, basemean = NA_real_, dispersion = NA_real_,
                         expressed_in_oligo = FALSE, is_tf_lambert = TRUE, lambert_dbd,
                         panther_protein_class = NA_character_)])

scorable <- is.finite(gene_df$basemean) & gene_df$basemean > 0 &
            is.finite(gene_df$dispersion) & gene_df$dispersion > 0
message(sprintf("[multiLFC] %d genes (%d not expressed), %d scorable",
                nrow(gene_df), sum(!gene_df$expressed_in_oligo), sum(scorable)))

classify_trinary <- function(p) {
  cls <- rep("Undetectable", length(p))
  cls[!is.na(p) & p > 0.4 & p <= 0.8] <- "Unsure"
  cls[!is.na(p) & p > 0.8]            <- "Detectable"
  cls
}

per_method <- list()
for (tag in names(LFCS)) for (m in DE_METHODS) {
  pred <- rep(NA_real_, nrow(gene_df))
  pred[scorable] <- predict_power(
    n_donors = FIXED_DONORS, abs_lfc = LFCS[[tag]],
    basemean = gene_df$basemean[scorable], dispersion = gene_df$dispersion[scorable],
    n_cells = FIXED_CELLS, reads_per_cell = FIXED_READS,
    cell_type = CELL_TYPE, de_method = m,
    bundle_path = bundle_path, warn_extrapolation = FALSE)
  per_method[[sprintf("power_%s_lfc%s", m, tag)]] <- pred
  per_method[[sprintf("class_%s_lfc%s", m, tag)]] <- classify_trinary(pred)
}
for (tag in names(LFCS)) gene_df[[sprintf("power_lfc%s", tag)]] <- per_method[[sprintf("power_deseq2_lfc%s", tag)]]
for (tag in names(LFCS)) gene_df[[sprintf("class_lfc%s", tag)]] <- per_method[[sprintf("class_deseq2_lfc%s", tag)]]
for (nm in names(per_method)) gene_df[[nm]] <- per_method[[nm]]

# gnomAD v2.1.1 loss-of-function constraint
gn <- fread(cmd = paste("gzip -dc", shQuote(need(file.path(ref_dir, "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")))),
            select = c("gene_id", "pLI", "exac_pLI"), na.strings = c("NA", ""))
setorderv(gn, c("gene_id", "pLI"), c(1L, -1L), na.last = TRUE)
gn <- gn[!duplicated(gene_id)]
gene_df <- merge(gene_df, gn, by = "gene_id", all.x = TRUE, sort = FALSE)
gene_df <- gene_df[match(c(ann$gene_id, lam$gene_id), gene_id)]
gene_df[, lof_intol := !is.na(pLI) & pLI >= 0.9]
stopifnot(!anyDuplicated(gene_df$gene_id))

# Booleans as True/False, the format the published table used.
for (cn in c("expressed_in_oligo", "is_tf_lambert", "lof_intol")) {
  set(gene_df, j = cn, value = ifelse(gene_df[[cn]], "True", "False"))
}
fwrite(gene_df, file.path(gc_dir, "oligo_gene_power_multiLFC.csv"))
fwrite(gene_df[is_tf_lambert == "True"], file.path(gc_dir, "oligo_tf_power_multiLFC.csv"))

for (tag in names(LFCS)) {
  tab <- table(factor(gene_df[[sprintf("class_lfc%s", tag)]], levels = c("Detectable", "Unsure", "Undetectable")))
  message(sprintf("[multiLFC] |LFC| = %.1f: Detectable %d, Unsure %d, Undetectable %d",
                  LFCS[[tag]], tab[["Detectable"]], tab[["Unsure"]], tab[["Undetectable"]]))
}
message(sprintf("[multiLFC] LoF-intolerant (pLI >= 0.9): %d", sum(gene_df$lof_intol == "True")))
