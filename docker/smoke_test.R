# =============================================================================
# docker/smoke_test.R -- verify the container's software environment.
#
# Checks, without any data:
#   1. every package in renv.lock is installed at exactly the locked version;
#   2. the packages the pipeline relies on load with library() (runs .onLoad);
#   3. every Python package pinned in python-requirements.txt is installed at
#      exactly the pinned version, and reticulate can drive the venv;
#   4. the fonts the figure composers need are present.
# Prints "SMOKE TEST PASSED" and exits 0, or lists every failure and exits 1.
#
# Inside the image the lockfile and requirements are copied to
# /usr/local/share/epa/. Run it with:
#   docker run --rm empirical-power-analysis:v1.0.0            (default command)
#   apptainer exec --cleanenv epa.sif Rscript /usr/local/share/epa/smoke_test.R
# docker/build-sif.sh runs it on every .sif it builds.
# =============================================================================

SHARE <- Sys.getenv("EPA_SHARE_DIR", "/usr/local/share/epa")
failures <- character()
fail <- function(...) failures <<- c(failures, sprintf(...))

# 1. renv.lock versions -------------------------------------------------------
lock <- jsonlite::fromJSON(file.path(SHARE, "renv.lock"), simplifyVector = FALSE)
locked <- vapply(lock$Packages, function(p) p$Version, character(1))
installed <- installed.packages()[, "Version"]
for (pkg in names(locked)) {
  if (!pkg %in% names(installed)) {
    fail("R package %s %s is in renv.lock but not installed", pkg, locked[[pkg]])
  } else if (installed[[pkg]] != locked[[pkg]]) {
    fail("R package %s: installed %s, renv.lock %s", pkg, installed[[pkg]], locked[[pkg]])
  }
}
cat(sprintf("renv.lock: %d packages checked\n", length(locked)))

# 2. packages the pipeline loads ----------------------------------------------
critical <- c("DESeq2", "edgeR", "limma", "ashr", "glmGamPoi", "SingleCellExperiment",
              "zellkonverter", "reticulate", "qs", "data.table", "future.apply",
              "xgboost", "ggplot2", "dplyr", "tidyr", "patchwork", "Matrix",
              "rmarkdown", "bookdown", "knitr", "writexl", "png", "here")
for (pkg in critical) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) {
    fail("library(%s) failed: %s", pkg, conditionMessage(e))
    FALSE
  })
  if (ok) cat(sprintf("OK   %-22s %s\n", pkg, packageVersion(pkg)))
}

# 3. Python pins ----------------------------------------------------------------
req <- readLines(file.path(SHARE, "python-requirements.txt"))
pins <- regmatches(req, regexpr("^[A-Za-z0-9_.-]+==[^ \\\\]+", req))
py <- Sys.getenv("RETICULATE_PYTHON", "/opt/python-env/bin/python")
for (pin in pins) {
  name <- sub("==.*", "", pin)
  want <- sub(".*==", "", pin)
  got <- suppressWarnings(system2(py, c("-c", shQuote(sprintf(
    "import importlib.metadata as m; print(m.version('%s'))", name))), stdout = TRUE, stderr = FALSE))
  if (length(got) != 1 || got != want) {
    fail("Python package %s: installed %s, pinned %s", name,
         if (length(got)) got[1] else "none", want)
  }
}
cat(sprintf("python-requirements.txt: %d pins checked\n", length(pins)))
tryCatch({
  reticulate::use_python(py, required = TRUE)
  reticulate::py_run_string("import anndata, PIL, matplotlib")
}, error = function(e) fail("reticulate / anndata import failed: %s", conditionMessage(e)))

# 4. fonts for the figure composers -------------------------------------------
if (!file.exists("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")) {
  fail("DejaVu fonts missing")
}

if (length(failures)) {
  cat("\nSMOKE TEST FAILED:\n", paste0("  ", failures, "\n"), sep = "")
  quit(status = 1)
}
cat("\nSMOKE TEST PASSED\n")
