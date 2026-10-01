# ==============================================================================
# STEP 08: Tidy Long-Format Recall / Precision Grid
# ==============================================================================
# Flattens the per-(cell type, mode, method) Aggregated_Summary_Data.qs
# checkpoints into a single tidy long-format CSV, one row per
#   (cell_type, mode, method, n_donors, y_value, lfc_bin, gene_category)
# carrying mean_recall, mean_precision and the mean ground-truth gene count.
#
# This is the queryable source the paper's numeric-reference binding system
# (_paper_values.R / val()) needs for medium-tier claims that quote a specific
# grid cell (e.g. "in the tiny-effect bin at 300 donors x 100 cells: microglia
# 20% recall, astrocytes 15%, oligodendrocytes 14%"). The grids otherwise live
# only as .qs checkpoints and rendered heatmap PNGs, which are not queryable.
#
# Standalone, lightweight (qs + dplyr only), no HPC. Filesystem-driven: it
# exports exactly the checkpoints that exist, so absent methods/modes are simply
# skipped rather than erroring.
#
# Inputs:
#   - Datalink/Cartesian_Downsampling/{CellType}/{Dataset}/{mode}/{method}/Plots/
#       Aggregated_Summary_Data.qs   (produced by cartesian_downsampling/03)
#     method in {deseq2, edger, limma, all_methods}; mode in {cells, reads}
#
# Outputs:
#   - Datalink/Paper_Source_Data/recall_precision_grid.csv
#
# NA in y_value denotes the "Max" step (all available cells / reads per donor).
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 08: TIDY RECALL / PRECISION GRID ###\n")

cartesian_base <- project_path("Datalink", "Cartesian_Downsampling")
out_path       <- project_path("Datalink", "Paper_Source_Data", "recall_precision_grid.csv")

# ------------------------------------------------------------------------------
# Discover checkpoints
#   path = .../Cartesian_Downsampling/<Cell>/<Dataset>/<mode>/<method>/Plots/
#          Aggregated_Summary_Data.qs
# ------------------------------------------------------------------------------
ckpts <- Sys.glob(file.path(cartesian_base, "*", "*", "*", "*", "Plots",
                            "Aggregated_Summary_Data.qs"))

if (length(ckpts) == 0) {
  stop(sprintf("No Aggregated_Summary_Data.qs checkpoints found under %s.\n  Run cartesian_downsampling/03 first.",
               cartesian_base))
}

message(sprintf("Found %d checkpoint(s).", length(ckpts)))

parse_meta <- function(p) {
  parts <- strsplit(p, .Platform$file.sep, fixed = TRUE)[[1]]
  n <- length(parts)
  # ... <Cell> <Dataset> <mode> <method> Plots Aggregated_Summary_Data.qs
  list(
    cell_type = parts[n - 5],
    dataset   = parts[n - 4],
    mode      = parts[n - 3],
    method    = parts[n - 2]
  )
}

parts_list <- list()

for (p in ckpts) {
  meta <- parse_meta(p)
  y_col <- if (meta$mode == "reads") "Target_Reads" else "Cells"

  cp <- tryCatch(qs::qread(p), error = function(e) {
    message(sprintf("  [skip] unreadable checkpoint: %s (%s)", p, conditionMessage(e)))
    NULL
  })
  if (is.null(cp)) next

  df <- cp$df_summary
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("  [skip] empty df_summary: %s/%s/%s", meta$cell_type, meta$mode, meta$method))
    next
  }

  required <- c("Donors", y_col, "LFC_Bin", "Gene_Category", "Mean_Recall", "Mean_Precision")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    message(sprintf("  [skip] %s/%s/%s missing columns: %s",
                    meta$cell_type, meta$mode, meta$method, paste(missing, collapse = ", ")))
    next
  }

  tidy <- df %>%
    transmute(
      cell_type      = meta$cell_type,
      dataset        = meta$dataset,
      mode           = meta$mode,
      method         = meta$method,
      n_donors       = Donors,
      y_value        = .data[[y_col]],
      y_type         = y_col,
      lfc_bin        = as.character(LFC_Bin),
      gene_category  = as.character(Gene_Category),
      mean_recall    = Mean_Recall,
      mean_precision = Mean_Precision,
      ref_gt_count   = if ("Ref_GT_Count" %in% names(df)) Ref_GT_Count else NA_real_
    )

  parts_list[[p]] <- tidy
  message(sprintf("  [ok] %-6s %-6s %-11s -> %d rows",
                  meta$cell_type, meta$mode, meta$method, nrow(tidy)))
}

if (length(parts_list) == 0) {
  stop("No usable checkpoints produced any rows.")
}

out <- bind_rows(parts_list) %>%
  arrange(cell_type, mode, method, n_donors, y_value, lfc_bin, gene_category)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)
message(sprintf("\n[OK] Wrote %s  (%d rows, %d checkpoints)",
                out_path, nrow(out), length(parts_list)))

# ------------------------------------------------------------------------------
# Verification: reproduce the paper's tiny-effect-bin example so the numbers can
# be checked against the draft's "[TODO: Check actual numbers]" placeholder.
# Draft claim: tiny bin, 300 donors x 100 cells -> micro 20%, astro 15%, oligo 14%.
# ------------------------------------------------------------------------------
message("\n--- Verification: tiny-effect bin, 300 donors x 100 cells (recall %) ---")
check <- out %>%
  filter(mode == "cells",
         n_donors == 300,
         y_value == 100,
         grepl("^Tiny", lfc_bin)) %>%
  arrange(gene_category, method, cell_type)

if (nrow(check) == 0) {
  message("  (no matching grid cell found)")
} else {
  for (gc in unique(check$gene_category)) {
    message(sprintf("  gene_category = %s", gc))
    sub <- check[check$gene_category == gc, ]
    for (i in seq_len(nrow(sub))) {
      message(sprintf("    %-11s %-6s recall=%5.1f%%  precision=%5.1f%%  (gt_genes=%.1f)",
                      sub$method[i], sub$cell_type[i],
                      100 * sub$mean_recall[i], 100 * sub$mean_precision[i],
                      sub$ref_gt_count[i]))
    }
  }
}

message("\n[COMPLETE] Step 08")
