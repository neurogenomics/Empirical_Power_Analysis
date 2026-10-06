#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  01_Full_Cohort_DE.R
# PURPOSE: Full-cohort pseudobulk DE for one contrast with the SAME DE workers,
#          covariates and ashr shrinkage the grid uses
#          (R/cartesian_downsampling/utils.R run_all_de_methods), so the
#          sex-biased and AD vs control effect distributions are directly
#          comparable in 02_Compare_Effect_Distributions.R.
#
#   EPA_CONTRAST=sex           all donors of Shards/<CT>/Roussos,
#                              ~ sex + age + pmi + Source + ethnicity + disease
#   EPA_CONTRAST=ad_vs_control all donors of Shards/<CT>/Roussos_ADvsCtrl,
#                              ~ ad_group + age + pmi + sex
#   (SHARD_DATASET must match; R/shared/contrast.R enforces it.)
#
# OUTPUT:  Datalink/AD_vs_Control/Full_Cohort_DE/<contrast>/<CT>/<method>.qs  (DE tables)
#          Datalink/AD_vs_Control/Full_Cohort_DE/<contrast>_de_results.csv     (tidy, all genes)
# USAGE:   EPA_CONTRAST=ad_vs_control SHARD_DATASET=Roussos_ADvsCtrl \
#            Rscript R/ad_contrast/01_Full_Cohort_DE.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
})

ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT", ".")
Sys.setenv(DOWNSAMPLE_MODE = "cells")
source(file.path(ROOT, "R", "cartesian_downsampling", "config.R"))
source(file.path(ROOT, "R", "cartesian_downsampling", "utils.R"))

OUT_BASE <- file.path(ROOT, "Datalink", "AD_vs_Control", "Full_Cohort_DE")
CELL_TYPES <- c("Astro", "Micro", "Oligo")
log_section(sprintf("Full-cohort DE: contrast=%s, dataset=%s, methods=%s",
                    EPA_CONTRAST, SHARD_DATASET, paste(DE_METHODS, collapse = ",")))

rows <- list()
for (ct in CELL_TYPES) {
  shard_dir <- file.path(SHARD_BASE, ct, SHARD_DATASET)
  mf  <- qs::qread(file.path(shard_dir, "00_MASTER_MANIFEST.qs"))$manifest
  ids <- as.character(mf$donor_id)
  log_info("%s: aggregating %d donors", ct, length(ids))
  agg <- aggregate_from_shards(ids, shard_dir, mode = "cells")
  cd  <- SummarizedExperiment::colData(agg$sce)
  log_info("%s: %s", ct, paste(names(table(cd[[SEX_COL]])), table(cd[[SEX_COL]]),
                               sep = "=", collapse = ", "))

  od <- file.path(OUT_BASE, EPA_CONTRAST, ct)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  status <- run_all_de_methods(
    agg, methods = DE_METHODS,
    save_path_results = setNames(file.path(od, paste0(DE_METHODS, ".qs")), DE_METHODS),
    save_path_log     = setNames(file.path(od, paste0(DE_METHODS, "_log.qs")), DE_METHODS))
  log_info("%s DE status: %s", ct, paste(names(status), status, sep = "=", collapse = ", "))
  if (any(startsWith(status, "ERROR"))) stop("DE failed for ", ct, ": ", paste(status, collapse = "; "))

  for (m in DE_METHODS) {
    r <- as.data.frame(qs::qread(file.path(od, paste0(m, ".qs"))))
    bm <- if ("baseMean_raw" %in% names(r)) r$baseMean_raw else r$baseMean
    rows[[paste(ct, m)]] <- data.frame(
      contrast = EPA_CONTRAST, cell_type = ct, de_method = m,
      n_donors = length(ids), gene_id = rownames(r),
      base_mean = as.numeric(bm), log2FoldChange = as.numeric(r$log2FoldChange),
      padj = as.numeric(r$padj), stringsAsFactors = FALSE)
  }
  file.remove(file.path(od, paste0(DE_METHODS, "_log.qs")))
}

res <- do.call(rbind, rows)
chr <- get_chr_map(unique(res$gene_id))
res$chromosome <- if (is.null(chr)) NA_character_ else unname(chr[res$gene_id])
res$chr_class  <- ifelse(is.na(res$chromosome), "Unmapped",
                         ifelse(res$chromosome %in% c("X", "Y"), "SexChr", "Autosomal"))
out_csv <- file.path(OUT_BASE, paste0(EPA_CONTRAST, "_de_results.csv"))
write.csv(res, out_csv, row.names = FALSE)
log_info("Wrote %s (%d rows)", out_csv, nrow(res))
