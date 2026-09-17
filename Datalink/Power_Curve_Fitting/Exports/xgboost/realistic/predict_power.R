# ==============================================================================
# predict_power.R -- Standalone XGBoost Power Predictor
# ==============================================================================
# Exported: 2026-08-11   Pipeline: v4.2.0   xgboost: 3.2.1.1
# Model: XGBoost 'realistic' (gradient-boosted trees; no closed-form equation).
#
# Requirements: R with the `xgboost` package installed. Nothing else.
#
# Usage:
#   source("predict_power.R")   # looks for power_model_bundle.rds alongside it
#   predict_power(n_donors = 100, abs_lfc = 1.5, basemean = 50, dispersion = 0.5)
#
#   # vectorised over any argument:
#   predict_power(n_donors = c(25, 50, 100, 200), abs_lfc = 2)
#
#   # power-vs-donors curve at fixed gene characteristics:
#   power_curve(donors = c(10, 25, 50, 100, 200, 500), abs_lfc = 1.5,
#               basemean = 50, dispersion = 0.5)
#
# Arguments (scalars or vectors; recycled to common length):
#   n_donors        Number of donors.
#   abs_lfc         Absolute log2 fold change of the target gene.
#   basemean        Mean expression (DESeq2 baseMean).
#   dispersion      Gene-level dispersion (DESeq2 dispersion).
#   n_cells         Cells per donor.            (default from bundle)
#   reads_per_cell  Mean reads per cell.        (default from bundle)
#   cell_type       "Astro" | "Micro" | "Oligo" (default from bundle)
#   de_method       "deseq2" | "edger" | "limma" (default from bundle)
#
# Returns: numeric vector of predicted power (detection probability) in [0, 1].
# ==============================================================================

suppressPackageStartupMessages(library(xgboost))

.power_env <- new.env(parent = emptyenv())

.load_power_bundle <- function(bundle_path = NULL) {
  if (!is.null(.power_env$bundle) && is.null(bundle_path)) return(.power_env$bundle)
  if (is.null(bundle_path)) {
    # Resolve power_model_bundle.rds sitting next to this script.
    this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    base_dir <- if (!is.null(this) && nzchar(this)) dirname(normalizePath(this)) else getwd()
    bundle_path <- file.path(base_dir, "power_model_bundle.rds")
    if (!file.exists(bundle_path)) bundle_path <- "power_model_bundle.rds"
  }
  if (!file.exists(bundle_path))
    stop("power_model_bundle.rds not found. Pass bundle_path= explicitly.")
  b <- readRDS(bundle_path)
  b$booster <- xgb.load.raw(b$booster_raw)
  .power_env$bundle <- b
  b
}

.build_design <- function(df, b) {
  cont_mat <- as.matrix(df[, b$cont_features, drop = FALSE])
  oh_blocks <- lapply(b$cat_features, function(f) {
    lv <- b$cat_levels[[f]]
    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
    colnames(m) <- paste0(f, "=", lv)
    m
  })
  X <- do.call(cbind, c(list(cont_mat), oh_blocks))
  X[, b$feature_names, drop = FALSE]   # exact fit-time column order
}

predict_power <- function(n_donors, abs_lfc, basemean, dispersion,
                          n_cells = NULL, reads_per_cell = NULL,
                          cell_type = NULL, de_method = NULL,
                          bundle_path = NULL, warn_extrapolation = TRUE) {
  b <- .load_power_bundle(bundle_path)
  d <- b$defaults
  if (is.null(n_cells))        n_cells        <- d$n_cells
  if (is.null(reads_per_cell)) reads_per_cell <- d$reads_per_cell
  if (is.null(cell_type))      cell_type      <- d$cell_type
  if (is.null(de_method))      de_method      <- d$de_method

  df <- data.frame(
    n_donors = n_donors, n_cells = n_cells, reads_per_cell = reads_per_cell,
    abs_gt_lfc = abs_lfc, gt_basemean = basemean, gt_dispersion = dispersion,
    cell_type = as.character(cell_type), de_method = as.character(de_method),
    stringsAsFactors = FALSE)   # recycles scalars/vectors to common length

  ok_ct <- df$cell_type %in% b$cat_levels$cell_type
  ok_dm <- df$de_method %in% b$cat_levels$de_method
  if (!all(ok_ct)) stop(sprintf("cell_type must be one of {%s}",
                                paste(b$cat_levels$cell_type, collapse = ", ")))
  if (!all(ok_dm)) stop(sprintf("de_method must be one of {%s}",
                                paste(b$cat_levels$de_method, collapse = ", ")))

  if (warn_extrapolation) {
    for (f in b$cont_features) {
      rng <- b$training_ranges[[f]]
      if (any(df[[f]] < rng["min"] | df[[f]] > rng["max"], na.rm = TRUE))
        warning(sprintf(paste0("%s outside training range [%g, %g]; tree models ",
                               "extrapolate as a flat constant beyond the grid."),
                        f, rng["min"], rng["max"]), call. = FALSE)
    }
  }

  X <- .build_design(df, b)
  as.numeric(predict(b$booster, xgb.DMatrix(X)))
}

# Convenience: power vs a vector of donor counts at fixed gene characteristics.
power_curve <- function(donors, abs_lfc, basemean, dispersion, ...) {
  data.frame(n_donors = donors,
             power = predict_power(n_donors = donors, abs_lfc = abs_lfc,
                                   basemean = basemean, dispersion = dispersion, ...))
}

