#!/usr/bin/env Rscript
# ==============================================================================
# build_supplementary_data.R -- package Supplementary Data 1-4.
# ------------------------------------------------------------------------------
# Writes Datalink/Supplementary_Data/:
#   SD1  Supplementary_Data_1_dataset_summary.xlsx
#          sheet Dataset_Characteristics   <- Paper_Source_Data/dataset_statistics.csv
#          sheet Gene_Universe_Composition <- Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv
#   SD2  Supplementary_Data_2_cartesian_downsampling/ (+ .zip)
#          aggregated grid summaries, <cell_type>/<mode>/<method>/Aggregated_Summary_Data.qs,
#          with the per-run *_raw components removed; all_methods for cells mode only.
#          README.txt is copied verbatim from R/supplementary_data/SD2_README.txt.
#   SD3  Supplementary_Data_3_truth_gene_characteristics.csv
#          <- DE_Diagnostics/Shared_Gene_Diagnostics/Full_Truth_Set_Summary.csv
#   SD4  Supplementary_Data_4_xgboost_hyperparameters.csv
#          <- Power_Curve_Fitting/Models/xgboost/model_selection.csv
#
# Replaces Datalink/Supplementary_Data/_scripts/{package_supplementary_data,slim_sd2}.R
# (hardcoded paths; README and zip were made by hand). The zip is built with
# fixed file times and no extra attributes so it is byte-reproducible.
#
# Usage: Rscript R/supplementary_data/build_supplementary_data.R
#        (project root from POWER_ANALYSIS_ROOT, else here::here())
# ==============================================================================
suppressPackageStartupMessages({
  library(qs)
  library(writexl)
})

root <- Sys.getenv("POWER_ANALYSIS_ROOT", unset = "")
if (!nzchar(root)) root <- here::here()
DL  <- file.path(root, "Datalink")
OUT <- file.path(DL, "Supplementary_Data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

need <- function(p) {
  if (!file.exists(p)) stop("missing input: ", p, call. = FALSE)
  p
}

## ---- SD1: two-sheet xlsx ---------------------------------------------------
sd1a <- read.csv(need(file.path(DL, "Paper_Source_Data", "dataset_statistics.csv")), check.names = FALSE)
sd1b <- read.csv(need(file.path(DL, "Full_Ground_Truth", "Gene_Audit_Full_Universe_Combined.csv")),
                 check.names = FALSE)
sd1_path <- file.path(OUT, "Supplementary_Data_1_dataset_summary.xlsx")
write_xlsx(list(Dataset_Characteristics = sd1a, Gene_Universe_Composition = sd1b), path = sd1_path)
cat(sprintf("SD1: %s (%d + %d rows)\n", basename(sd1_path), nrow(sd1a), nrow(sd1b)))

## ---- SD2: slim aggregated grid summaries ------------------------------------
sd2_name <- "Supplementary_Data_2_cartesian_downsampling"
sd2_root <- file.path(OUT, sd2_name)
unlink(sd2_root, recursive = TRUE)
methods_by_mode <- list(cells = c("deseq2", "edger", "limma", "all_methods"),
                        reads = c("deseq2", "edger", "limma"))
n <- 0L
for (ct in c("Astro", "Micro", "Oligo")) for (mode in names(methods_by_mode)) {
  for (meth in methods_by_mode[[mode]]) {
    x <- qread(need(file.path(DL, "Cartesian_Downsampling", ct, "Roussos", mode, meth,
                              "Plots", "Aggregated_Summary_Data.qs")))
    keep <- names(x)[!grepl("_raw$", names(x))]
    od <- file.path(sd2_root, ct, mode, meth)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    qsave(x[keep], file.path(od, "Aggregated_Summary_Data.qs"))
    n <- n + 1L
    cat(sprintf("SD2: %-6s %-6s %-11s kept %d/%d components\n", ct, mode, meth,
                length(keep), length(x)))
  }
}
stopifnot(n == 21L)
file.copy(need(file.path(root, "R", "supplementary_data", "SD2_README.txt")),
          file.path(sd2_root, "README.txt"), overwrite = TRUE)

# Deterministic zip: fixed timestamps on every entry, no extra file attributes.
zip_path <- file.path(OUT, paste0(sd2_name, ".zip"))
unlink(zip_path)
fixed_time <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
entries <- c(sd2_root, list.files(sd2_root, recursive = TRUE, full.names = TRUE, include.dirs = TRUE))
invisible(Sys.setFileTime(entries, fixed_time))
old_wd <- setwd(OUT)
rel <- sort(sub(paste0("^", OUT, "/"), "", entries))
status <- system2("zip", c("-X", "-q", zip_path, rel))
setwd(old_wd)
if (status != 0) stop("zip failed with status ", status)
cat(sprintf("SD2: %s (%d files)\n", basename(zip_path), n))

## ---- SD3 / SD4: plain copies -------------------------------------------------
file.copy(need(file.path(DL, "DE_Diagnostics", "Shared_Gene_Diagnostics", "Full_Truth_Set_Summary.csv")),
          file.path(OUT, "Supplementary_Data_3_truth_gene_characteristics.csv"), overwrite = TRUE)
file.copy(need(file.path(DL, "Power_Curve_Fitting", "Models", "xgboost", "model_selection.csv")),
          file.path(OUT, "Supplementary_Data_4_xgboost_hyperparameters.csv"), overwrite = TRUE)
cat("SD3, SD4 written\n")
