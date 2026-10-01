# ==============================================================================
# SCRIPT: classify_oligo_multi_lfc.R
# PURPOSE: Re-predict oligodendrocyte gene detectability at THREE fixed effect
#          sizes (hard / medium / easy) from the CACHED full-cohort features,
#          so we can show how the trinary calls shift with |LFC|.
#
#   Hard   : |LFC| = 0.75
#   Medium : |LFC| = 1.50   (the existing setting)
#   Easy   : |LFC| = 2.25
#
# Everything else is held fixed exactly as in classify_oligo_genes.R:
#   100 donors, 250 cells/donor, 5742 reads/cell, Oligo.
#   Per-gene inputs = raw pseudobulk basemean + DESeq2 dispersion (full-cohort
#   sex fit), read straight from oligo_gene_features_checkpoint.qs (NO refit).
#   Predict per DE method by toggling the de_method one-hot (deseq2/edger/limma).
#   Trinary: power<=0.4 Undetectable | 0.4-0.8 Unsure | >0.8 Detectable.
#   Non-scorable genes (no finite basemean/dispersion) -> Undetectable at all LFC.
#
# OUTPUT: Datalink/Gene_Classification/oligo_gene_power_multiLFC.csv
# ==============================================================================

suppressPackageStartupMessages({ library(xgboost) })

if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT")

CELL_TYPE    <- "Oligo"
FIXED_DONORS <- 100
FIXED_CELLS  <- 250
FIXED_READS  <- 5742
DE_METHODS   <- c("deseq2", "edger", "limma")
LFC_SETTINGS <- c(hard = 0.75, medium = 1.50, easy = 2.25)

out_dir    <- file.path(ROOT, "Datalink", "Gene_Classification")
ckpt_file  <- file.path(out_dir, "oligo_gene_features_checkpoint.qs")
bundle_dir <- file.path(ROOT, "Datalink", "Power_Curve_Fitting",
                        "Exports", "xgboost", "realistic")
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")

message("[multiLFC] loading cached features: ", ckpt_file)
feat <- qs::qread(ckpt_file)
gene_df <- data.frame(
  gene_id             = feat$gene_universe,
  basemean            = as.numeric(feat$basemean_raw),
  dispersion          = as.numeric(feat$dispersion),
  n_donors_expressing = as.integer(feat$n_donors_expressing),
  stringsAsFactors    = FALSE
)
scorable <- is.finite(gene_df$basemean) & gene_df$basemean > 0 &
            is.finite(gene_df$dispersion) & gene_df$dispersion > 0
message(sprintf("[multiLFC] scorable %d / %d genes", sum(scorable), nrow(gene_df)))

source(file.path(bundle_dir, "predict_power.R"))

classify_trinary <- function(p) {
  cls <- rep(NA_character_, length(p))
  cls[!is.na(p) & p <= 0.4]           <- "Undetectable"
  cls[!is.na(p) & p > 0.4 & p <= 0.8] <- "Unsure"
  cls[!is.na(p) & p > 0.8]            <- "Detectable"
  cls
}

for (nm in names(LFC_SETTINGS)) {
  lfc <- LFC_SETTINGS[[nm]]
  tag <- sub("\\.", "", sprintf("%.2f", lfc))  # 0.75 -> 075, 1.50 -> 150
  for (m in DE_METHODS) {
    pred <- rep(NA_real_, nrow(gene_df))
    pred[scorable] <- predict_power(
      n_donors       = FIXED_DONORS,
      abs_lfc        = lfc,
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
    gene_df[[sprintf("power_%s_lfc%s", m, tag)]] <- pred
    gene_df[[sprintf("class_%s_lfc%s", m, tag)]] <- cls
  }
  # convenience: deseq2 primary columns per LFC
  gene_df[[sprintf("power_lfc%s", tag)]] <- gene_df[[sprintf("power_deseq2_lfc%s", tag)]]
  gene_df[[sprintf("class_lfc%s", tag)]] <- gene_df[[sprintf("class_deseq2_lfc%s", tag)]]
}

out_csv <- file.path(out_dir, "oligo_gene_power_multiLFC.csv")
write.csv(gene_df, out_csv, row.names = FALSE)
message("[multiLFC] wrote ", out_csv)

for (nm in names(LFC_SETTINGS)) {
  tag <- sub("\\.", "", sprintf("%.2f", LFC_SETTINGS[[nm]]))
  tab <- table(factor(gene_df[[sprintf("class_lfc%s", tag)]],
                      levels = c("Detectable","Unsure","Undetectable")))
  message(sprintf("[multiLFC] %-6s |LFC|=%.2f  Detectable %d | Unsure %d | Undetectable %d",
                  nm, LFC_SETTINGS[[nm]], tab[["Detectable"]], tab[["Unsure"]], tab[["Undetectable"]]))
}
