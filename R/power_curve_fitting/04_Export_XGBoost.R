# ==============================================================================
# SCRIPT: Export XGBoost Power Model (M2 -- portable, shareable predictor)
# PURPOSE: Take the promoted XGBoost booster (the paper's power model) and emit a
#          self-contained bundle plus a dependency-light R wrapper that anyone
#          can source() to compute power for their own designs. This is the
#          XGBoost replacement for the retired logistic 05_Export_Equations.R.
#
#          There is NO closed-form equation for a gradient-boosted tree, so the
#          "exportable power function" IS the serialised booster + a predict
#          wrapper that reproduces the exact design-matrix recipe used at fit
#          time (02_Fit_XGBoost.R / 03_Select_XGBoost.R).
#
# MODEL: the single power model (trained on the n_cells >= 100 design rows,
#   selected over the {tree hyperparameters x alpha} search), promoted by
#   03_Select_XGBoost.R to Models/xgboost/. XGB_MODEL_SUBDIR exports a side
#   experiment from Models/xgboost/<subdir>/ instead.
#
# DETERMINISM: the export date written into the bundle, wrapper and README
#   comes from SOURCE_DATE_EPOCH when set (the reproduction pipeline sets it
#   to the code commit time), so re-running the export gives identical files.
#
# FEATURE HANDLING (design decision): cell_type, de_method and reads_per_cell
#   are KEPT in the model (dropping them would silently answer cell-type- and
#   method-specific questions with a pooled average) but the wrapper gives them
#   sensible PsychAD defaults so a naive caller supplies only the four things
#   they actually vary (donors, |LFC|, baseMean, dispersion). A power user can
#   override any argument.
#
# INPUT:  Models/xgboost/<subdir>/{booster.ubj,feature_spec.qs,
#             feature_importance.csv,shap_summary.csv,cv_metrics.csv,
#             calibration_residuals.csv}
#         Fitting_Data/<fitting_data>.qs  (for training-range extrapolation guards)
# OUTPUT: Exports/xgboost/<subdir>/{power_model_bundle.rds,predict_power.R,
#             example_power_curve.csv,README.md}
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(xgboost)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
    "config.R",
    "../config.R"
  )
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

XGB_MODEL_SUBDIR <- Sys.getenv("XGB_MODEL_SUBDIR", "")
XGB_GENE_SET     <- Sys.getenv("XGB_GENE_SET", "full")
MODEL_LABEL      <- if (nzchar(XGB_MODEL_SUBDIR)) sprintf("XGBoost '%s'", XGB_MODEL_SUBDIR) else "XGBoost power model"
.sde <- Sys.getenv("SOURCE_DATE_EPOCH", "")
EXPORT_DATE <- if (nzchar(.sde)) {
  as.Date(as.POSIXct(as.numeric(.sde), origin = "1970-01-01", tz = "UTC"))
} else Sys.Date()

log_section(sprintf("STEP 4 (XGBoost): Export portable power model [%s]", MODEL_LABEL))

model_base <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Models", "xgboost")
}
if (!dir.exists(model_base)) {
  stop(sprintf("Model dir not found: %s. Run 02_Fit_XGBoost.R + 03_Select_XGBoost.R first.", model_base))
}

export_dir <- if (nzchar(XGB_MODEL_SUBDIR)) {
  file.path(OUTPUT_BASE, "Exports", "xgboost", XGB_MODEL_SUBDIR)
} else {
  file.path(OUTPUT_BASE, "Exports", "xgboost")
}
safe_mkdir(export_dir)

# ==============================================================================
# 2. LOAD PROMOTED BOOSTER + SPEC + DIAGNOSTICS
# ==============================================================================

# Load the booster from the version-stable .ubj (NOT booster.qs): an xgboost
# booster is an external C++ pointer and the qs round-trip is unreliable across
# the xgboost 2.x -> 3.x rewrite (can rehydrate as a null handle). This mirrors
# the load discipline in 03_Select_XGBoost.R.
booster_ubj <- file.path(model_base, "booster.ubj")
if (!file.exists(booster_ubj)) stop(sprintf("booster.ubj not found: %s", booster_ubj))
booster <- xgb.load(booster_ubj)

feature_spec <- qs::qread(file.path(model_base, "feature_spec.qs"))
CONT_FEATURES <- feature_spec$cont_features
CAT_FEATURES  <- feature_spec$cat_features
cat_levels    <- feature_spec$cat_levels
feature_names <- feature_spec$feature_names

log_info("Loaded booster (%d rounds, objective=%s) + feature_spec (%d features)",
         feature_spec$best_nrounds, feature_spec$objective, length(feature_names))

read_csv_if <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
feature_importance <- read_csv_if(file.path(model_base, "feature_importance.csv"))
shap_summary       <- read_csv_if(file.path(model_base, "shap_summary.csv"))
cv_metrics         <- read_csv_if(file.path(model_base, "cv_metrics.csv"))
calibration_resid  <- read_csv_if(file.path(model_base, "calibration_residuals.csv"))

# ==============================================================================
# 3. TRAINING RANGES (for extrapolation warnings in the wrapper)
# ==============================================================================
# The wrapper warns when a caller supplies continuous inputs outside the grid
# the model was actually trained on, since tree models extrapolate as a flat
# constant beyond the training range (monotone but uninformative).

fitting_data_name <- if (!is.null(feature_spec$fitting_data)) {
  feature_spec$fitting_data
} else {
  "aggregated_fitting_data.qs"
}
agg_path <- file.path(OUTPUT_BASE, "Fitting_Data", fitting_data_name)
if (!file.exists(agg_path)) {
  stop(sprintf("Fitting data not found: %s (needed for training-range guards).", agg_path))
}
dat <- qs::qread(agg_path)
if ("gene_set" %in% colnames(dat)) dat <- dat %>% filter(gene_set == XGB_GENE_SET)
dat <- dat %>%
  filter(is.finite(detection_rate), detection_rate >= 0, detection_rate <= 1,
         n_trials > 0,
         is.finite(n_donors), is.finite(n_cells), is.finite(reads_per_cell),
         is.finite(abs_gt_lfc), is.finite(gt_basemean), is.finite(gt_dispersion))
log_info("Training-range reference: %s (%d rows)", fitting_data_name, nrow(dat))

training_ranges <- lapply(CONT_FEATURES, function(f) {
  v <- dat[[f]]
  c(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))
})
names(training_ranges) <- CONT_FEATURES

# ==============================================================================
# 4. DEFAULTS (PsychAD-representative)
# ==============================================================================
# Conservative reference cell type = Oligo (least depth-censored truth set).
# reads_per_cell default matches that cell type's PsychAD mean; n_cells default
# is a typical realistic per-donor yield. All overridable by the caller.
DEFAULT_CELL_TYPE <- if ("Oligo" %in% cat_levels$cell_type) "Oligo" else cat_levels$cell_type[1]
DEFAULT_DE_METHOD <- if ("deseq2" %in% cat_levels$de_method) "deseq2" else cat_levels$de_method[1]
DEFAULT_READS_PER_CELL <- if (!is.null(MEAN_READS_PER_CELL[[DEFAULT_CELL_TYPE]])) {
  as.numeric(MEAN_READS_PER_CELL[[DEFAULT_CELL_TYPE]])
} else {
  round(median(dat$reads_per_cell, na.rm = TRUE))
}
DEFAULT_N_CELLS <- 250

model_defaults <- list(
  n_cells        = DEFAULT_N_CELLS,
  reads_per_cell = DEFAULT_READS_PER_CELL,
  cell_type      = DEFAULT_CELL_TYPE,
  de_method      = DEFAULT_DE_METHOD
)
log_info("Defaults: cell_type=%s de_method=%s n_cells=%g reads_per_cell=%g",
         DEFAULT_CELL_TYPE, DEFAULT_DE_METHOD, DEFAULT_N_CELLS, DEFAULT_READS_PER_CELL)

# ==============================================================================
# 5. WRITE SELF-CONTAINED BUNDLE
# ==============================================================================
# Booster serialised as raw UBJ bytes (xgb.save.raw): survives the 2.x->3.x
# rewrite and travels inside the .rds so there is no external .ubj to lose.

bundle <- list(
  booster_raw        = xgb.save.raw(booster, raw_format = "ubj"),
  feature_spec       = feature_spec,
  cont_features      = CONT_FEATURES,
  cat_features       = CAT_FEATURES,
  cat_levels         = cat_levels,
  feature_names      = feature_names,
  defaults           = model_defaults,
  training_ranges    = training_ranges,
  feature_importance = feature_importance,
  shap_summary       = shap_summary,
  cv_metrics         = cv_metrics,
  calibration_resid  = calibration_resid,
  metadata = list(
    model_subdir     = XGB_MODEL_SUBDIR,
    fitting_data     = fitting_data_name,
    gene_set         = XGB_GENE_SET,
    objective        = feature_spec$objective,
    best_nrounds     = feature_spec$best_nrounds,
    best_params      = feature_spec$best_params,
    alpha            = feature_spec$alpha,
    pipeline_version = PIPELINE_VERSION,
    exported_on      = as.character(EXPORT_DATE),
    xgboost_version  = as.character(utils::packageVersion("xgboost"))
  )
)

bundle_path <- file.path(export_dir, "power_model_bundle.rds")
saveRDS(bundle, bundle_path)
log_info("Wrote model bundle -> %s (%.1f KB)",
         bundle_path, file.info(bundle_path)$size / 1024)

# ==============================================================================
# 6. WRITE STANDALONE predict_power.R WRAPPER
# ==============================================================================
# Dependency-light: base R + xgboost only. Reads the bundle at load time so it
# never drifts from the serialised booster. Vectorised, validates categoricals,
# warns on extrapolation beyond the training grid.

wrapper_lines <- c(
  '# ==============================================================================',
  '# predict_power.R -- Standalone XGBoost Power Predictor',
  '# ==============================================================================',
  sprintf('# Exported: %s   Pipeline: v%s   xgboost: %s',
          EXPORT_DATE, PIPELINE_VERSION, as.character(utils::packageVersion("xgboost"))),
  sprintf('# Model: %s (gradient-boosted trees; no closed-form equation).', MODEL_LABEL),
  '#',
  '# Requirements: R with the `xgboost` package installed. Nothing else.',
  '#',
  '# Usage:',
  '#   source("predict_power.R")   # looks for power_model_bundle.rds alongside it',
  '#   predict_power(n_donors = 100, abs_lfc = 1.5, basemean = 50, dispersion = 0.5)',
  '#',
  '#   # vectorised over any argument:',
  '#   predict_power(n_donors = c(25, 50, 100, 200), abs_lfc = 2)',
  '#',
  '#   # power-vs-donors curve at fixed gene characteristics:',
  '#   power_curve(donors = c(10, 25, 50, 100, 200, 500), abs_lfc = 1.5,',
  '#               basemean = 50, dispersion = 0.5)',
  '#',
  '# Arguments (scalars or vectors; recycled to common length):',
  '#   n_donors        Number of donors.',
  '#   abs_lfc         Absolute log2 fold change of the target gene.',
  '#   basemean        Mean expression (DESeq2 baseMean).',
  '#   dispersion      Gene-level dispersion (DESeq2 dispersion).',
  '#   n_cells         Cells per donor.            (default from bundle)',
  '#   reads_per_cell  Mean reads per cell.        (default from bundle)',
  '#   cell_type       "Astro" | "Micro" | "Oligo" (default from bundle)',
  '#   de_method       "deseq2" | "edger" | "limma" (default from bundle)',
  '#',
  '# Returns: numeric vector of predicted power (detection probability) in [0, 1].',
  '# ==============================================================================',
  '',
  'suppressPackageStartupMessages(library(xgboost))',
  '',
  '.power_env <- new.env(parent = emptyenv())',
  '',
  '.load_power_bundle <- function(bundle_path = NULL) {',
  '  if (!is.null(.power_env$bundle) && is.null(bundle_path)) return(.power_env$bundle)',
  '  if (is.null(bundle_path)) {',
  '    # Resolve power_model_bundle.rds sitting next to this script.',
  '    this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)',
  '    base_dir <- if (!is.null(this) && nzchar(this)) dirname(normalizePath(this)) else getwd()',
  '    bundle_path <- file.path(base_dir, "power_model_bundle.rds")',
  '    if (!file.exists(bundle_path)) bundle_path <- "power_model_bundle.rds"',
  '  }',
  '  if (!file.exists(bundle_path))',
  '    stop("power_model_bundle.rds not found. Pass bundle_path= explicitly.")',
  '  b <- readRDS(bundle_path)',
  '  b$booster <- xgb.load.raw(b$booster_raw)',
  '  .power_env$bundle <- b',
  '  b',
  '}',
  '',
  '.build_design <- function(df, b) {',
  '  cont_mat <- as.matrix(df[, b$cont_features, drop = FALSE])',
  '  oh_blocks <- lapply(b$cat_features, function(f) {',
  '    lv <- b$cat_levels[[f]]',
  '    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))',
  '    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))',
  '    colnames(m) <- paste0(f, "=", lv)',
  '    m',
  '  })',
  '  X <- do.call(cbind, c(list(cont_mat), oh_blocks))',
  '  X[, b$feature_names, drop = FALSE]   # exact fit-time column order',
  '}',
  '',
  'predict_power <- function(n_donors, abs_lfc, basemean, dispersion,',
  '                          n_cells = NULL, reads_per_cell = NULL,',
  '                          cell_type = NULL, de_method = NULL,',
  '                          bundle_path = NULL, warn_extrapolation = TRUE) {',
  '  b <- .load_power_bundle(bundle_path)',
  '  d <- b$defaults',
  '  if (is.null(n_cells))        n_cells        <- d$n_cells',
  '  if (is.null(reads_per_cell)) reads_per_cell <- d$reads_per_cell',
  '  if (is.null(cell_type))      cell_type      <- d$cell_type',
  '  if (is.null(de_method))      de_method      <- d$de_method',
  '',
  '  df <- data.frame(',
  '    n_donors = n_donors, n_cells = n_cells, reads_per_cell = reads_per_cell,',
  '    abs_gt_lfc = abs_lfc, gt_basemean = basemean, gt_dispersion = dispersion,',
  '    cell_type = as.character(cell_type), de_method = as.character(de_method),',
  '    stringsAsFactors = FALSE)   # recycles scalars/vectors to common length',
  '',
  '  ok_ct <- df$cell_type %in% b$cat_levels$cell_type',
  '  ok_dm <- df$de_method %in% b$cat_levels$de_method',
  '  if (!all(ok_ct)) stop(sprintf("cell_type must be one of {%s}",',
  '                                paste(b$cat_levels$cell_type, collapse = ", ")))',
  '  if (!all(ok_dm)) stop(sprintf("de_method must be one of {%s}",',
  '                                paste(b$cat_levels$de_method, collapse = ", ")))',
  '',
  '  if (warn_extrapolation) {',
  '    for (f in b$cont_features) {',
  '      rng <- b$training_ranges[[f]]',
  '      if (any(df[[f]] < rng["min"] | df[[f]] > rng["max"], na.rm = TRUE))',
  '        warning(sprintf(paste0("%s outside training range [%g, %g]; tree models ",',
  '                               "extrapolate as a flat constant beyond the grid."),',
  '                        f, rng["min"], rng["max"]), call. = FALSE)',
  '    }',
  '  }',
  '',
  '  X <- .build_design(df, b)',
  '  as.numeric(predict(b$booster, xgb.DMatrix(X)))',
  '}',
  '',
  '# Convenience: power vs a vector of donor counts at fixed gene characteristics.',
  'power_curve <- function(donors, abs_lfc, basemean, dispersion, ...) {',
  '  data.frame(n_donors = donors,',
  '             power = predict_power(n_donors = donors, abs_lfc = abs_lfc,',
  '                                   basemean = basemean, dispersion = dispersion, ...))',
  '}',
  ''
)

wrapper_path <- file.path(export_dir, "predict_power.R")
writeLines(wrapper_lines, wrapper_path)
log_info("Wrote standalone wrapper -> %s", wrapper_path)

# ==============================================================================
# 7. WORKED EXAMPLE + README
# ==============================================================================
# Reconstruct the booster the same way the wrapper does and emit a power-vs-
# donors curve at a moderate-effect gene, so the recipient has a runnable check.

booster_rt <- xgb.load.raw(bundle$booster_raw)
build_design_local <- function(df) {
  cont_mat <- as.matrix(df[, CONT_FEATURES, drop = FALSE])
  oh_blocks <- lapply(CAT_FEATURES, function(f) {
    lv <- cat_levels[[f]]
    m <- vapply(lv, function(l) as.numeric(as.character(df[[f]]) == l), numeric(nrow(df)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
    colnames(m) <- paste0(f, "=", lv)
    m
  })
  X <- do.call(cbind, c(list(cont_mat), oh_blocks))
  X[, feature_names, drop = FALSE]
}
example_donors <- c(10, 25, 50, 100, 200, 300, 500, 600)
example_df <- data.frame(
  n_donors = example_donors, n_cells = DEFAULT_N_CELLS,
  reads_per_cell = DEFAULT_READS_PER_CELL,
  abs_gt_lfc = 1.5, gt_basemean = 50, gt_dispersion = 0.5,
  cell_type = DEFAULT_CELL_TYPE, de_method = DEFAULT_DE_METHOD,
  stringsAsFactors = FALSE)
example_df$power <- as.numeric(predict(booster_rt, xgb.DMatrix(build_design_local(example_df))))
write.csv(example_df, file.path(export_dir, "example_power_curve.csv"), row.names = FALSE)
log_info("Wrote worked example -> %s", file.path(export_dir, "example_power_curve.csv"))

readme_lines <- c(
  "# Portable XGBoost Power Model",
  "",
  sprintf("Exported %s from the %s (pipeline v%s, xgboost %s).",
          EXPORT_DATE, MODEL_LABEL,
          PIPELINE_VERSION, as.character(utils::packageVersion("xgboost"))),
  "",
  "## Files",
  "- `power_model_bundle.rds` -- self-contained: serialised booster + feature spec + defaults + training ranges + diagnostics.",
  "- `predict_power.R` -- source() this; needs only base R + the `xgboost` package.",
  "- `example_power_curve.csv` -- a worked power-vs-donors curve to check against.",
  "",
  "## Quick start",
  "```r",
  'source("predict_power.R")',
  "predict_power(n_donors = 100, abs_lfc = 1.5, basemean = 50, dispersion = 0.5)",
  "power_curve(donors = c(10, 25, 50, 100, 200, 500),",
  "            abs_lfc = 1.5, basemean = 50, dispersion = 0.5)",
  "```",
  "",
  "## Inputs",
  "Required: `n_donors`, `abs_lfc`, `basemean`, `dispersion`.",
  sprintf("Optional (defaults in brackets): `n_cells` [%g], `reads_per_cell` [%g], `cell_type` [%s], `de_method` [%s].",
          DEFAULT_N_CELLS, DEFAULT_READS_PER_CELL, DEFAULT_CELL_TYPE, DEFAULT_DE_METHOD),
  "All arguments accept vectors and are recycled to a common length.",
  "",
  "## Notes",
  "- Output is predicted power (detection probability) in [0, 1]; there is no closed-form equation (gradient-boosted trees).",
  "- `cell_type`, `de_method` and `reads_per_cell` are kept in the model but default to PsychAD/Oligo settings; override to vary them.",
  "- Predictions outside the training grid warn: tree models extrapolate as a flat constant beyond the observed range.",
  sprintf("- Reference cell type is Oligo (least depth-censored truth set). Winner: %s.",
          if (!is.null(cv_metrics)) sprintf("config %s, alpha=%s, CV logloss=%.4f",
                                            cv_metrics$config_index[1], cv_metrics$alpha[1],
                                            cv_metrics$cv_logloss[1]) else "see cv_metrics.csv")
)
writeLines(readme_lines, file.path(export_dir, "README.md"))
log_info("Wrote README -> %s", file.path(export_dir, "README.md"))

log_section("STEP 4 (XGBoost) EXPORT COMPLETE")
log_info("Portable power model exported to: %s", export_dir)
log_info("Example power (|LFC|=1.5, baseMean=50, disp=0.5, %s): %s",
         DEFAULT_CELL_TYPE,
         paste(sprintf("N%d=%.3f", example_df$n_donors, example_df$power), collapse = " "))
