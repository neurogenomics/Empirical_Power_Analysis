# ==============================================================================
# SCRIPT: Aggregate Fitting Dataset -> per-grid-point detection rate (M2 / XGBoost)
# PURPOSE: Collapse the per-iteration detection table produced by Step 1 into a
#          per-grid-point table whose target is a detection RATE in [0,1], ready
#          for the XGBoost power model (Step 2 replacement, 02_Fit_XGBoost.R).
#
# Aggregation unit (Open decision A = iterations-only; runs kept separate):
#   (gene_id, cell_type, de_method, mode, n_donors, y_value, run_id)
#
# For each unit, collapse the ~N_ITERATIONS iterations into ONE row:
#   n_detected      = sum(detected)
#   n_trials        = number of iterations present
#   detection_rate  = n_detected / n_trials          <- XGBoost training target
#   n_cells         = mean(n_cells)                  (actual measured depth)
#   reads_per_cell  = mean(reads_per_cell)
#   gene characteristics (abs_gt_lfc, gt_basemean, gt_dispersion) carried through
#
# de_method axis:
#   The per-iteration Step-1 output is written per DE method (selected via the
#   CARTESIAN_DE_METHOD env var when Step 1 runs). This script reads whichever
#   method-tagged per-iteration files it can find under Fitting_Data/ and stacks
#   them, tagging each with de_method. For the FIRST TRIAL we aggregate whatever
#   is present; with only the default (deseq2) per-iteration files on disk this
#   yields a single-method aggregated table, and re-running Step 1 with
#   CARTESIAN_DE_METHOD=edger / =limma populates the other two slices.
#
#   Method-tagged inputs are discovered as:
#     Fitting_Data/{model_name}/unified_fitting_data_per_iteration{_shared}.qs           (deseq2, legacy default)
#     Fitting_Data/{model_name}/unified_fitting_data_per_iteration__{method}{_shared}.qs (method-namespaced)
#   To produce method-namespaced files, run Step 1 with CARTESIAN_DE_METHOD set
#   and PCF_METHOD_NAMESPACE=true (see 01_Build_Fitting_Dataset note); otherwise
#   the legacy filename is treated as de_method = CARTESIAN_DE_METHOD.
#
# INPUT:  Fitting_Data/{model_name}/unified_fitting_data_per_iteration*.qs (from Step 1)
# OUTPUT: Power_Curve_Fitting/Fitting_Data/{model_name}/aggregated_fitting_data.qs
#         Power_Curve_Fitting/Fitting_Data/aggregated_fitting_data_all.qs   (stacked)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
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

log_section("STEP 1b: Aggregate Fitting Dataset -> per-grid-point detection rate")

# ==============================================================================
# 2. HELPERS
# ==============================================================================

#' Parse a per-iteration filename into (de_method, gene_set)
#'
#' Recognises:
#'   unified_fitting_data_per_iteration.qs                 -> (default method, full)
#'   unified_fitting_data_per_iteration_shared.qs          -> (default method, shared)
#'   unified_fitting_data_per_iteration__edger.qs          -> (edger, full)
#'   unified_fitting_data_per_iteration__limma_shared.qs   -> (limma, shared)
parse_per_iter_filename <- function(fname, default_method) {
  stem <- sub("\\.qs$", "", fname)
  stem <- sub("^unified_fitting_data_per_iteration", "", stem)
  gene_set <- if (grepl("_shared$", stem)) "shared" else "full"
  stem <- sub("_shared$", "", stem)
  method <- if (grepl("^__", stem)) sub("^__", "", stem) else default_method
  if (!nzchar(method)) method <- default_method
  list(de_method = method, gene_set = gene_set)
}

#' Aggregate a per-iteration data frame to per-grid-point detection rate.
#'
#' Open decision A (default): keep each run_id as its own row (aggregate
#' iterations only). Set POOL_RUNS=true to also pool across runs.
aggregate_per_grid_point <- function(df, pool_runs = FALSE) {
  group_cols <- c("gene_id", "cell_type", "de_method", "mode",
                  "n_donors", "y_value")
  if (!pool_runs) group_cols <- c(group_cols, "run_id")

  # Gene characteristics are constant per (gene_id, run_id); take the mean so
  # the row is well-defined whether or not runs are pooled.
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_detected     = sum(detected, na.rm = TRUE),
      n_trials       = dplyr::n(),
      detection_rate = sum(detected, na.rm = TRUE) / dplyr::n(),
      n_cells        = mean(n_cells, na.rm = TRUE),
      reads_per_cell = mean(reads_per_cell, na.rm = TRUE),
      abs_gt_lfc     = mean(abs_gt_lfc, na.rm = TRUE),
      gt_basemean    = mean(gt_basemean, na.rm = TRUE),
      gt_dispersion  = mean(gt_dispersion, na.rm = TRUE),
      gt_lfc         = mean(gt_lfc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log_donors         = log(n_donors),
      log_n_cells        = log(n_cells),
      log_reads_per_cell = log(reads_per_cell),
      log_abs_gt_lfc     = log(abs_gt_lfc),
      log_gt_basemean    = log(gt_basemean),
      log_gt_dispersion  = log(gt_dispersion)
    )
}

POOL_RUNS <- tolower(Sys.getenv("PCF_POOL_RUNS", "false")) %in% c("true", "1", "yes")
log_info("Run pooling (Open decision A): %s", if (POOL_RUNS) "POOL runs" else "keep runs separate [default]")

# ==============================================================================
# 3. DISCOVER + AGGREGATE PER MODEL NAME
# ==============================================================================

all_aggregated <- list()

for (ds_config in DATASET_CONFIGS) {
  for (cell_type in ds_config$cell_types) {

    model_name <- get_dataset_model_name(cell_type, ds_config)
    paths <- get_paths_for_dataset(cell_type, ds_config)
    fitting_dir <- paths$fitting_dir

    if (!dir.exists(fitting_dir)) {
      log_warning("No fitting dir for %s: %s", model_name, fitting_dir)
      next
    }

    per_iter_files <- list.files(
      fitting_dir,
      pattern = "^unified_fitting_data_per_iteration.*\\.qs$",
      full.names = FALSE
    )
    if (length(per_iter_files) == 0) {
      log_warning("No per-iteration files for %s -- run Step 1 first", model_name)
      next
    }

    model_aggs <- list()

    for (fname in per_iter_files) {
      meta <- parse_per_iter_filename(fname, default_method = CARTESIAN_DE_METHOD)

      # We only model the FULL gene set (shared dropped by decision, 2026-07-15).
      if (meta$gene_set != "full") {
        log_info("  %-45s set=%-6s : skipped (only 'full' is modelled)", fname, meta$gene_set)
        next
      }

      df <- qs::qread(file.path(fitting_dir, fname))

      # Step-1 output uses model_name in the cell_type column already.
      if (!"cell_type" %in% colnames(df)) df$cell_type <- model_name
      df$de_method <- meta$de_method
      df$gene_set  <- meta$gene_set

      if (!"y_value" %in% colnames(df)) {
        log_error(paste0(
          "%s is missing the y_value column. Re-run Step 1 (01_Build_Fitting_Dataset.R) ",
          "after the y_value patch so grid points can be identified."), fname)
        next
      }

      agg <- aggregate_per_grid_point(df, pool_runs = POOL_RUNS)
      agg$gene_set <- meta$gene_set

      log_info("  %-45s method=%-7s set=%-6s : %d per-iter rows -> %d grid-point rows",
               fname, meta$de_method, meta$gene_set, nrow(df), nrow(agg))

      model_aggs[[fname]] <- agg
    }

    if (length(model_aggs) == 0) next

    model_agg <- bind_rows(model_aggs)

    out_path <- file.path(fitting_dir, "aggregated_fitting_data.qs")
    qs::qsave(model_agg, out_path)
    log_info("Saved %s: %d rows (%d genes, methods: %s)",
             out_path, nrow(model_agg), n_distinct(model_agg$gene_id),
             paste(sort(unique(model_agg$de_method)), collapse = ", "))

    all_aggregated[[model_name]] <- model_agg
  }
}

# ==============================================================================
# 4. STACK ALL MODELS + SANITY REPORT
# ==============================================================================

if (length(all_aggregated) == 0) {
  stop("No aggregated data produced. Run Step 1 first (per-iteration files missing).")
}

combined <- bind_rows(all_aggregated)

combined_path <- file.path(OUTPUT_BASE, "Fitting_Data", "aggregated_fitting_data_all.qs")
safe_mkdir(dirname(combined_path))
qs::qsave(combined, combined_path)

log_section("Aggregation summary")
log_info("Combined aggregated table: %d rows", nrow(combined))
log_info("  cell types : %s", paste(sort(unique(combined$cell_type)), collapse = ", "))
log_info("  de methods : %s", paste(sort(unique(combined$de_method)), collapse = ", "))
log_info("  gene sets  : %s", paste(sort(unique(combined$gene_set)), collapse = ", "))
log_info("  genes      : %d", n_distinct(combined$gene_id))
log_info("  n_trials   : median=%.0f, min=%d, max=%d",
         median(combined$n_trials), min(combined$n_trials), max(combined$n_trials))
log_info("  detection_rate: mean=%.3f (0s: %d, 1s: %d)",
         mean(combined$detection_rate),
         sum(combined$detection_rate == 0), sum(combined$detection_rate == 1))
log_info("Saved combined table: %s", combined_path)

log_section("STEP 1b COMPLETE")
