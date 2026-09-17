# ==============================================================================
# Capture Reproducibility Environment Snapshot
# ==============================================================================
# Writes a full snapshot of the R/system environment the pipelines were run
# in, so that future replicators can match package versions exactly.
#
# Run with:   Rscript R/shared/capture_environment.R
# Outputs land under:  Reproducibility/
#   - session_info.txt        (sessionInfo() with all project packages loaded)
#   - installed_packages.csv  (every installed package + version)
#   - project_packages.csv    (only project-used packages + version, hand-picked)
#   - system_info.txt         (OS, CPU, memory, R build platform)
# ==============================================================================

suppressPackageStartupMessages({
  here_ok <- requireNamespace("here", quietly = TRUE)
})

# Figure out the project root robustly
if (here_ok) {
  repo_root <- tryCatch(here::here(), error = function(e) getwd())
} else {
  # script may be invoked from anywhere; walk upwards looking for .here or README.md
  repo_root <- getwd()
  for (i in 1:10) {
    if (file.exists(file.path(repo_root, ".here")) ||
        file.exists(file.path(repo_root, "Empirical_Power_Analysis.Rproj"))) break
    parent <- dirname(repo_root)
    if (parent == repo_root) break
    repo_root <- parent
  }
}

out_dir <- file.path(repo_root, "Reproducibility")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[capture_environment] writing to: %s", out_dir))

# --- Project packages (hand-picked from library()/require() calls in R/) -----
project_packages <- c(
  "akima", "AnnotationDbi", "apeglm", "ashr", "BiocParallel", "biomaRt",
  "cowplot", "data.table", "DESeq2", "dplyr", "future",
  "future.apply", "ggnewscale", "ggplot2", "ggrepel", "glue",
  "grid", "httr", "lmtest", "Matrix", "org.Hs.eg.db",
  "parallel", "pheatmap", "qs", "RColorBrewer", "reticulate",
  "rhdf5", "S4Vectors", "scales", "shadowtext", "SingleCellExperiment",
  "stringr", "SummarizedExperiment", "testthat", "tidyr", "tidyverse",
  "tools", "UpSetR", "viridis", "zellkonverter", "here"
)

# --- Load as many as possible (record success/failure) ----------------------
load_results <- data.frame(
  package = character(0),
  loaded  = logical(0),
  version = character(0),
  stringsAsFactors = FALSE
)

for (pkg in project_packages) {
  ok <- suppressWarnings(suppressMessages(
    requireNamespace(pkg, quietly = TRUE)
  ))
  ver <- if (ok) as.character(utils::packageVersion(pkg)) else NA_character_
  if (ok) suppressPackageStartupMessages(try(library(pkg, character.only = TRUE), silent = TRUE))
  load_results <- rbind(load_results, data.frame(
    package = pkg, loaded = ok, version = ver, stringsAsFactors = FALSE
  ))
}

write.csv(load_results,
          file.path(out_dir, "project_packages.csv"),
          row.names = FALSE)

# --- Full installed.packages() snapshot -------------------------------------
ip <- utils::installed.packages()
ip_df <- as.data.frame(ip[, c("Package", "Version", "Priority",
                              "LibPath", "Built", "NeedsCompilation")],
                       stringsAsFactors = FALSE)
write.csv(ip_df,
          file.path(out_dir, "installed_packages.csv"),
          row.names = FALSE)

# --- sessionInfo() ----------------------------------------------------------
si <- sessionInfo()
si_path <- file.path(out_dir, "session_info.txt")
writeLines(capture.output({
  cat("### sessionInfo()\n\n")
  print(si)
  cat("\n\n### Sys.getlocale()\n")
  cat(Sys.getlocale(), "\n")
  cat("\n### .libPaths()\n")
  print(.libPaths())
}), si_path)

# --- System info ------------------------------------------------------------
sys_path <- file.path(out_dir, "system_info.txt")
writeLines(capture.output({
  cat("### R.version\n")
  print(R.version)
  cat("\n### Sys.info()\n")
  print(Sys.info())
  cat("\n### capabilities()\n")
  print(capabilities())
  cat("\n### extSoftVersion()\n")
  print(extSoftVersion())
  cat("\n### La_library(), La_version()\n")
  print(try(La_library(), silent = TRUE))
  print(try(La_version(), silent = TRUE))
  cat("\n### BLAS/LAPACK info\n")
  cat("BLAS:   "); cat(try(extSoftVersion()[["BLAS"]], silent = TRUE), "\n")
  cat("LAPACK: "); cat(try(extSoftVersion()[["LAPACK"]], silent = TRUE), "\n")
  cat("\n### timestamp\n")
  cat(format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
}), sys_path)

message("[capture_environment] Done.")
message(sprintf("  - %s", file.path(out_dir, "project_packages.csv")))
message(sprintf("  - %s", file.path(out_dir, "installed_packages.csv")))
message(sprintf("  - %s", file.path(out_dir, "session_info.txt")))
message(sprintf("  - %s", file.path(out_dir, "system_info.txt")))
