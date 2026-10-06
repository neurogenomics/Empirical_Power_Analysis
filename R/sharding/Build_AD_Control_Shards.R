#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  Build_AD_Control_Shards.R
# PURPOSE: Build the donor shard set for the AD vs control contrast
#          (Shards/<cell type>/Roussos_ADvsCtrl/) from the filtered PsychAD
#          shards (Shards/<cell type>/Roussos/). The grid, ground truth and
#          power-model pipelines then run on it with EPA_CONTRAST=ad_vs_control
#          (R/shared/contrast.R).
#
# DONOR GROUPS (per-donor diagnosis flags in the shard colData):
#   AD       AD_status == "Yes" and no other diagnosis recorded as "Yes", with
#            every other diagnosis evaluated (not "Not Evaluated")
#   Control  AD_status == "No" and no diagnosis recorded as "Yes", with every
#            diagnosis evaluated
#   Diagnoses: AD, DLBD, FTD, Tauopathy, Vascular, ASCVD, Schizophrenia,
#   Bipolar_Disorder, Parkinson_disease, Tardive_dyskinesia. Tauopathy_status is
#   only ever recorded as "Yes" or "Not Evaluated" (never "No"), so for that one
#   flag "not Yes" counts as clear. The dataset's own `disease == "normal"`
#   label is NOT used: it includes donors with a recorded diagnosis.
#   Every other donor is excluded.
#
# Each selected shard is copied unchanged except for one added colData column,
# ad_group ("AD" / "Control"). The manifest is the source manifest restricted
# to the selected donors, with a record of the definition.
#
# INPUT:  Datalink/Shards/{Astro,Micro,Oligo}/Roussos/{donor_*.qs,00_MASTER_MANIFEST.qs}
# OUTPUT: Datalink/Shards/<CT>/Roussos_ADvsCtrl/{donor_*.qs,00_MASTER_MANIFEST.qs,
#                                               ad_control_donors.csv}
#         Datalink/AD_vs_Control/ad_control_donor_counts.csv
# USAGE:  Rscript R/sharding/Build_AD_Control_Shards.R [n_cores]
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
N_CORES <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 8L
ROOT    <- Sys.getenv("POWER_ANALYSIS_ROOT", ".")
SHARDS  <- file.path(ROOT, "Datalink", "Shards")
OUT_SUMMARY_DIR <- file.path(ROOT, "Datalink", "AD_vs_Control")
dir.create(OUT_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES  <- c("Astro", "Micro", "Oligo")
SOURCE_DS   <- "Roussos"
TARGET_DS   <- "Roussos_ADvsCtrl"
DIAGNOSES   <- c("AD_status", "DLBD_status", "FTD_status", "Tauopathy_status",
                 "Vascular_status", "ASCVD_status", "Schizophrenia",
                 "Bipolar_Disorder", "Parkinson_disease", "Tardive_dyskinesia")
NEVER_NO    <- "Tauopathy_status"   # never recorded as "No"
OTHER_DX    <- setdiff(DIAGNOSES, "AD_status")
MUST_EVAL   <- setdiff(DIAGNOSES, NEVER_NO)

assign_group <- function(d) {
  others_clear <- rowSums(d[OTHER_DX] == "Yes") == 0 &
                  rowSums(d[MUST_EVAL] == "Not Evaluated") == 0
  ifelse(d$AD_status == "Yes" & others_clear, "AD",
         ifelse(d$AD_status == "No" & others_clear, "Control", NA_character_))
}

counts <- list()
for (ct in CELL_TYPES) {
  src_dir <- file.path(SHARDS, ct, SOURCE_DS)
  dst_dir <- file.path(SHARDS, ct, TARGET_DS)
  src_mf  <- qs::qread(file.path(src_dir, "00_MASTER_MANIFEST.qs"))
  mf      <- src_mf$manifest
  files   <- file.path(src_dir, basename(mf$file_path))
  stopifnot(all(file.exists(files)))

  meta <- do.call(rbind, mclapply(files, function(f) {
    cd <- SummarizedExperiment::colData(qs::qread(f))
    missing_cols <- setdiff(c("donor_id", "sex", DIAGNOSES), colnames(cd))
    if (length(missing_cols)) stop(f, " lacks colData columns: ", paste(missing_cols, collapse = ", "))
    as.data.frame(cd[1, c("donor_id", "sex", "Source", DIAGNOSES)])
  }, mc.cores = N_CORES))
  meta$ad_group <- assign_group(meta)
  keep <- !is.na(meta$ad_group)
  message(sprintf("%s: %d AD, %d Control of %d donors", ct,
                  sum(meta$ad_group == "AD", na.rm = TRUE),
                  sum(meta$ad_group == "Control", na.rm = TRUE), nrow(meta)))

  unlink(dst_dir, recursive = TRUE)
  dir.create(dst_dir, recursive = TRUE)
  sel <- which(keep)
  ok <- unlist(mclapply(sel, function(i) {
    sce <- qs::qread(files[i])
    sce$ad_group <- meta$ad_group[i]
    qs::qsave(sce, file.path(dst_dir, basename(files[i])))
    TRUE
  }, mc.cores = N_CORES))
  stopifnot(length(ok) == length(sel), all(ok))

  sel_ids <- meta$donor_id[sel]
  new_mf  <- src_mf
  new_mf$manifest <- mf[mf$donor_id %in% sel_ids, , drop = FALSE]
  new_mf$manifest$file_path <- file.path(dst_dir, basename(new_mf$manifest$file_path))
  if (!is.null(new_mf$feasibility)) {
    new_mf$feasibility <- lapply(new_mf$feasibility, function(v) v[v %in% sel_ids])
  }
  new_mf$config$n_donors    <- length(sel_ids)
  new_mf$config$total_cells <- sum(new_mf$manifest$cell_count)
  new_mf$config$source_dir  <- src_dir
  new_mf$config$contrast    <- "ad_vs_control"
  new_mf$config$ad_control_definition <- paste(
    "AD: AD_status Yes, no other diagnosis Yes, all diagnoses evaluated;",
    "Control: no diagnosis Yes, all diagnoses evaluated (Tauopathy_status is never 'No')")
  qs::qsave(new_mf, file.path(dst_dir, "00_MASTER_MANIFEST.qs"))

  donors <- meta[sel, c("donor_id", "ad_group", "sex", "Source")]
  donors <- donors[order(donors$donor_id), ]
  write.csv(donors, file.path(dst_dir, "ad_control_donors.csv"), row.names = FALSE)

  tab <- as.data.frame(table(ad_group = meta$ad_group[sel], sex = meta$sex[sel]))
  tab$cell_type <- ct
  counts[[ct]] <- tab[, c("cell_type", "ad_group", "sex", "Freq")]
}

counts_df <- do.call(rbind, counts)
names(counts_df)[4] <- "n_donors"
write.csv(counts_df, file.path(OUT_SUMMARY_DIR, "ad_control_donor_counts.csv"), row.names = FALSE)
print(counts_df, row.names = FALSE)
message("Done.")
