# ==============================================================================
# SCRIPT: Cross-Method Aggregation (pool DE methods, keep cell types separate)
# PURPOSE: Produce the paper-facing "all-methods" summaries by trial-weighted
#          pooling of the per-(cell_type, method) Step 3 checkpoints. Cell types
#          stay separate; the three DE methods (deseq2/edger/limma) are pooled.
# INPUT:   Per-method Aggregated_Summary_Data.qs checkpoints (Step 3 output),
#          specifically their raw per-iteration frames df_raw (LFC) and
#          df_finecross_raw (fine LFC x BaseMean cross).
# OUTPUT:  An Aggregated_Summary_Data.qs checkpoint in the CROSS_METHOD_LABEL
#          ("all_methods") plots directory per cell type, shaped so 04/07/09 can
#          render it exactly like a per-method checkpoint.
#
# POOLING SEMANTICS (decision, handover §8.1): TRIAL-WEIGHTED. For each grid
# cell x bin x gene-category, every (method, run, iteration) observation is
# pooled by its trial count, so the pooled rate equals total detections over
# total trials:
#   recall            = sum(tp_rec)      / sum(GT_Count)       = wmean(Recall,         GT_Count)
#   recall_harm       = sum(tp_rec_harm) / sum(GT_Harm_Count)  = wmean(Recall_Harm,    GT_Harm_Count)
#   precision         = sum(tp_prec)     / sum(Sub_Count)      = wmean(Precision,      Sub_Count)
#   precision_harm    = sum(tp_prec_harm)/ sum(Sub_Harm_Count) = wmean(Precision_Harm, Sub_Harm_Count)
# This holds because Recall*GT_Count = tp_rec, Precision*Sub_Count = tp_prec, etc.
# (see 03_Aggregate_Results.R), so the weighted mean reconstructs the pooled
# count ratio exactly. This matches how 01b already pools iterations.
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",
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
  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 8: Cross-Method Aggregation (all_methods)")

y_col <- get_y_col_name()
log_info("Mode: %s | Y-column: %s | Pooling: trial-weighted across %s",
         DOWNSAMPLE_MODE, y_col, paste(DE_METHODS, collapse = "/"))

# ==============================================================================
# 2. POOLING HELPERS
# ==============================================================================

# Safe trial-weighted mean: sum(value * weight) / sum(weight). Rows with NA
# value or NA/zero weight contribute nothing. Returns na_empty when the total
# weight is zero (no trials), matching the per-method convention where recall
# defaults to 0 and precision to NA when there is nothing to score.
weighted_rate <- function(value, weight, na_empty = NA_real_) {
  ok <- !is.na(value) & !is.na(weight) & weight > 0
  if (!any(ok)) return(na_empty)
  sum(value[ok] * weight[ok]) / sum(weight[ok])
}

# Pool one raw per-iteration frame (stacked across methods) into a summary with
# the exact column set 04/07/09 expect. `group_cols` are the facet/grid keys.
pool_raw <- function(df_raw, group_cols) {
  df_raw %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      Mean_Precision      = weighted_rate(Precision,      Sub_Count),
      Mean_Recall         = weighted_rate(Recall,         GT_Count,      na_empty = 0),
      Mean_Precision_Harm = weighted_rate(Precision_Harm, Sub_Harm_Count),
      Mean_Recall_Harm    = weighted_rate(Recall_Harm,    GT_Harm_Count, na_empty = 0),
      Ref_GT_Count        = mean(GT_Count, na.rm = TRUE),
      Ref_Harm_Count      = mean(GT_Harm_Count, na.rm = TRUE),
      Mean_Sub_Count      = mean(Sub_Count, na.rm = TRUE),
      Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
      GT_Count_String     = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
      Harm_Count_String   = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
      CV_GT = {
        vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
        if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
      },
      CV_Harm = {
        vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
        if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
      },
      .groups = "drop"
    )
}

# ==============================================================================
# 3. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {

  log_section(sprintf("Cell Type: %s | pooling methods", cell_type))

  # Load every available per-method checkpoint for this cell type.
  raw_lfc_parts       <- list()
  raw_finecross_parts <- list()
  xist_parts          <- list()
  reads_summary_ref   <- NULL
  methods_used        <- character(0)

  for (method in DE_METHODS) {
    paths <- get_paths(cell_type, method = method)
    checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
    if (!file.exists(checkpoint_file)) {
      log_warning("No Step 3 checkpoint for %s/%s — skipping this method", cell_type, method)
      next
    }
    cp <- qs::qread(checkpoint_file)

    if (is.null(cp$df_raw)) {
      log_warning("Checkpoint for %s/%s has no df_raw — skipping", cell_type, method)
      next
    }
    if (!"Sub_Harm_Count" %in% names(cp$df_raw)) {
      stop(sprintf(paste0("Checkpoint for %s/%s predates the Sub_Harm_Count column. ",
                          "Re-run Step 3 (03_Aggregate_Results.R) before cross-method pooling."),
                   cell_type, method))
    }

    methods_used <- c(methods_used, method)
    raw_lfc_parts[[method]] <- cp$df_raw
    if (!is.null(cp$df_finecross_raw)) raw_finecross_parts[[method]] <- cp$df_finecross_raw
    if (!is.null(cp$df_xist_raw))      xist_parts[[method]] <- cp$df_xist_raw
    if (is.null(reads_summary_ref) && !is.null(cp$df_reads_summary)) {
      # Total reads per matrix are method-independent (same downsampling), so a
      # single method's reads summary is representative for the pool.
      reads_summary_ref <- cp$df_reads_summary
    }
  }

  if (length(raw_lfc_parts) == 0) {
    log_warning("No usable per-method checkpoints for %s — skipping cross-method pool", cell_type)
    next
  }
  log_info("Pooling %d method(s): %s", length(methods_used), paste(methods_used, collapse = ", "))

  # ----- Pool the standard LFC frame -----
  df_lfc_all <- do.call(rbind, raw_lfc_parts)
  df_summary <- pool_raw(df_lfc_all, c("Donors", y_col, "LFC_Bin", "Gene_Category"))
  log_info("Cross-method LFC summary: %d rows", nrow(df_summary))

  # ----- Pool the fine LFC x BaseMean cross frame (Figure T / S source) -----
  df_finecross_summary <- NULL
  if (length(raw_finecross_parts) > 0) {
    df_fc_all <- do.call(rbind, raw_finecross_parts)
    df_finecross_summary <- pool_raw(
      df_fc_all, c("Donors", y_col, "LFC_Bin", "BaseMean_Bin", "Gene_Category")
    )
    log_info("Cross-method fine-cross summary: %d rows", nrow(df_finecross_summary))
  } else {
    log_warning("No df_finecross_raw in any method checkpoint for %s — Figure T/S all_methods will be skipped", cell_type)
  }

  # ----- Pool XIST detection (unweighted mean across methods x runs x iters) --
  df_xist_summary <- NULL
  if (length(xist_parts) > 0) {
    df_xist_all <- do.call(rbind, xist_parts)
    df_xist_summary <- df_xist_all %>%
      group_by(Donors, !!sym(y_col), LFC_Bin) %>%
      summarise(Mean_Detection = mean(Detected, na.rm = TRUE), .groups = "drop")
  }

  # ==============================================================================
  # 4. SAVE all_methods CHECKPOINT (shaped like a per-method Step 3 checkpoint)
  # ==============================================================================
  # bm / cart / exprfilt summaries are intentionally NULL: the cross-method
  # deliverable is the standard LFC heatmaps + Figure T/S only (handover §8.5).
  # 04/07 guard NULL for those, so they render only the LFC set + XIST + reads.

  out_paths <- get_paths(cell_type, method = CROSS_METHOD_LABEL)
  if (!dir.exists(out_paths$plots_dir)) {
    dir.create(out_paths$plots_dir, recursive = TRUE)
    log_info("Created all_methods plots directory: %s", out_paths$plots_dir)
  }
  checkpoint_file <- file.path(out_paths$plots_dir, "Aggregated_Summary_Data.qs")

  qs::qsave(list(
    df_summary           = df_summary,
    df_xist_summary      = df_xist_summary,
    df_reads_summary     = reads_summary_ref,
    df_bm_summary        = NULL,
    df_cart_summary      = NULL,
    df_exprfilt_summary  = NULL,
    df_finecross_summary = df_finecross_summary,
    pooled_methods       = methods_used,
    pooling              = "trial-weighted"
  ), checkpoint_file)

  log_info("Saved all_methods checkpoint -> %s", checkpoint_file)
}

log_section("STEP 8 COMPLETE")
log_info("Cross-method aggregation complete for all cell types")
