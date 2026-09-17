#!/usr/bin/env Rscript
# Thin wrapper to run step 4 standalone (sources config first)
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) > 0) {
  SCRIPT_DIR_LOCAL <- dirname(sub("^--file=", "", file_arg[1]))
} else {
  SCRIPT_DIR_LOCAL <- "."
}
source(file.path(SCRIPT_DIR_LOCAL, "config.R"))
source(file.path(SCRIPT_DIR_LOCAL, "04_DEG_Chromosome_Grid.R"))
