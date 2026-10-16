#!/usr/bin/env Rscript
# ==============================================================================
# build_supplementary_data.R -- package Supplementary Data 1 and 4.
# ------------------------------------------------------------------------------
# Writes Datalink/Supplementary_Data/:
#   SD1  Supplementary_Data_1_dataset_summary.xlsx
#          sheet Dataset_Characteristics   <- Paper_Source_Data/dataset_statistics.csv
#          sheet Gene_Universe_Composition <- Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv
#   SD4  Supplementary_Data_4_xgboost_hyperparameters.csv
#          <- Power_Curve_Fitting/Models/xgboost/model_selection.csv
#
# SD2 and SD3 are not cited in the manuscript and are no longer built. The
# file names keep their original numbers so SD4 matches the manuscript.
#
# Usage: Rscript R/supplementary_data/build_supplementary_data.R
#        (project root from POWER_ANALYSIS_ROOT, else here::here())
# ==============================================================================
suppressPackageStartupMessages({
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

## ---- SD4: plain copy -------------------------------------------------------
invisible(file.copy(need(file.path(DL, "Power_Curve_Fitting", "Models", "xgboost", "model_selection.csv")),
          file.path(OUT, "Supplementary_Data_4_xgboost_hyperparameters.csv"), overwrite = TRUE))
cat("SD4 written\n")
