# ==============================================================================
# STEP 2b: Regenerate FDP + IoU Heatmaps from Cached CSVs
# ==============================================================================
# Fast-path viz-only script. Reads CSVs persisted by Step 2 and redraws the
# three heatmaps (raw IoU per cell type, universe-controlled IoU per cell type,
# combined FDP across cell types) using the current formatting conventions.
#
# Run this after editing pheatmap parameters (fonts, labels, column orders)
# without rerunning the full DE/IoU/FDP computation loop in Step 2.
#
# Inputs (per cell type, from Step 2):
#   {OUTPUT_BASE}/{CellType}/Diagnostic_Results/Expanded_IoU_Matrix.csv
#   {OUTPUT_BASE}/{CellType}/Diagnostic_Results/Expanded_IoU_Matrix_Universe_Controlled.csv
#   {OUTPUT_BASE}/{CellType}/Diagnostic_Results/FDP_per_source.csv
#
# Outputs:
#   {OUTPUT_BASE}/{CellType}/Diagnostic_Results/Expanded_IoU_Matrix.png
#   {OUTPUT_BASE}/{CellType}/Diagnostic_Results/Expanded_IoU_Matrix_Universe_Controlled.png
#   {OUTPUT_BASE}/FDP_Heatmap.png
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- PROJECT_ROOT
}

source(file.path(PROJECT_ROOT, "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

suppressPackageStartupMessages({
  library(pheatmap)
  library(glue)
  library(grid)
})

log_section("STEP 2b: Regenerate Heatmaps from Cached CSVs")
log_info("Output base: %s", OUTPUT_BASE)
log_info("Cell types: %s", paste(CELL_TYPES, collapse = ", "))

# ------------------------------------------------------------------------------
# Helper: classify IoU matrix columns into groups by name pattern
# ------------------------------------------------------------------------------
classify_iou_columns <- function(col_names) {
  grp <- character(length(col_names))
  for (i in seq_along(col_names)) {
    nm <- col_names[i]
    if (grepl("^Roussos_Half_", nm)) {
      grp[i] <- "Split-Half"
    } else if (grepl("^Null_", nm)) {
      grp[i] <- "Null"
    } else {
      grp[i] <- "Cross-Dataset"
    }
  }
  grp
}

# Shared annotation colours (match Step 2)
ann_colors <- list(
  Group = c(
    "Split-Half"    = "#2166AC",
    "Cross-Dataset" = "#4DAF4A",
    "Null"          = "#E41A1C"
  )
)

# ------------------------------------------------------------------------------
# Helper: render a single IoU heatmap
# ------------------------------------------------------------------------------
render_iou_heatmap <- function(iou_csv, out_png, main_title,
                               number_fontsize = 14) {
  if (!file.exists(iou_csv)) {
    log_warn("  Missing: %s", iou_csv)
    return(invisible(FALSE))
  }

  iou_df <- read.csv(iou_csv, row.names = 1, check.names = FALSE)
  iou_mat <- as.matrix(iou_df)

  col_groups <- classify_iou_columns(colnames(iou_mat))

  # Reorder columns/rows: Split-Half, Cross-Dataset, Null (preserve within-group order)
  group_priority <- setNames(seq_along(c("Split-Half", "Cross-Dataset", "Null")),
                             c("Split-Half", "Cross-Dataset", "Null"))
  col_idx <- order(group_priority[col_groups], seq_along(col_groups))
  iou_mat <- iou_mat[col_idx, col_idx, drop = FALSE]
  col_groups <- col_groups[col_idx]

  n_upper <- sum(col_groups == "Split-Half")
  n_datasets <- sum(col_groups == "Cross-Dataset")

  col_annotation <- data.frame(
    Group = factor(col_groups,
                   levels = c("Split-Half", "Cross-Dataset", "Null")),
    row.names = colnames(iou_mat)
  )

  # Prettify labels for display
  iou_mat_display <- iou_mat
  rownames(iou_mat_display) <- prettify_label(rownames(iou_mat_display))
  colnames(iou_mat_display) <- prettify_label(colnames(iou_mat_display))

  col_annotation_display <- col_annotation
  rownames(col_annotation_display) <- prettify_label(rownames(col_annotation_display))

  png(out_png, width = 2200, height = 1800, res = 120)
  pheatmap(iou_mat_display,
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = number_fontsize,
           number_color = "black",
           main = main_title,
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           breaks = seq(0, 1, length.out = 101),
           annotation_col = col_annotation_display,
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_col = c(n_upper, n_upper + n_datasets),
           gaps_row = c(n_upper, n_upper + n_datasets),
           fontsize = 15,
           fontsize_row = 14, fontsize_col = 14,
           show_rownames = FALSE,
           show_colnames = FALSE,
           annotation_names_col = FALSE,
           annotation_legend = FALSE,
           legend = TRUE)
  dev.off()

  log_info("  Saved %s", basename(out_png))
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# 1. Per-cell-type IoU heatmaps (raw + universe-controlled)
# ------------------------------------------------------------------------------
fdp_collector <- list()

for (cell_type in CELL_TYPES) {
  log_section(sprintf("Regenerating heatmaps: %s", cell_type))

  paths <- get_paths(cell_type)
  diag_dir <- paths$diagnostic_dir

  if (!dir.exists(diag_dir)) {
    log_warn("%s: Diagnostic directory missing, skipping.", cell_type)
    next
  }

  # Raw IoU
  render_iou_heatmap(
    iou_csv    = file.path(diag_dir, "Expanded_IoU_Matrix.csv"),
    out_png    = file.path(diag_dir, "Expanded_IoU_Matrix.png"),
    main_title = glue("IoU Without Universe Control ({REFERENCE_DATASET} Significant Genes): {prettify_label(cell_type)}"),
    number_fontsize = 14
  )

  # Universe-controlled IoU
  render_iou_heatmap(
    iou_csv    = file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.csv"),
    out_png    = file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.png"),
    main_title = glue("Universe-Controlled IoU ({REFERENCE_DATASET} Significant Genes): {prettify_label(cell_type)}"),
    number_fontsize = 13
  )

  # Collect FDP CSV for post-loop combined heatmap
  fdp_csv <- file.path(diag_dir, "FDP_per_source.csv")
  if (file.exists(fdp_csv)) {
    fdp_collector[[cell_type]] <- read.csv(fdp_csv, stringsAsFactors = FALSE)
  } else {
    log_warn("  Missing: %s", fdp_csv)
  }
}

# ------------------------------------------------------------------------------
# 2. Combined FDP heatmap across cell types
# ------------------------------------------------------------------------------
log_section("FDP Heatmap: Combining across cell types")

tryCatch({
  if (length(fdp_collector) > 0) {

    combined_fdp <- do.call(rbind, fdp_collector)

    # Save combined CSV
    combined_csv <- file.path(OUTPUT_BASE, "FDP_Summary.csv")
    write.csv(combined_fdp, combined_csv, row.names = FALSE, quote = TRUE)
    log_info("Saved combined FDP summary: %s", combined_csv)

    # Pivot to wide format: rows = cell types, columns = sources
    fdp_wide <- reshape(
      combined_fdp[, c("CellType", "Source", "FDP")],
      idvar     = "CellType",
      timevar   = "Source",
      direction = "wide"
    )
    rownames(fdp_wide) <- fdp_wide$CellType
    fdp_wide$CellType <- NULL
    colnames(fdp_wide) <- gsub("^FDP\\.", "", colnames(fdp_wide))

    # Group lookup (from per-cell-type CSV which already has Group column)
    first_ct <- fdp_collector[[1]]
    source_order <- first_ct$Source
    group_order  <- first_ct$Group

    available_cols <- intersect(source_order, colnames(fdp_wide))
    fdp_mat <- as.matrix(fdp_wide[, available_cols, drop = FALSE])

    row_order <- intersect(CELL_TYPES, rownames(fdp_mat))
    fdp_mat <- fdp_mat[row_order, , drop = FALSE]

    col_groups <- group_order[match(available_cols, source_order)]

    # Reorder columns to Split-Half, Cross-Dataset, Null
    group_priority <- setNames(seq_along(c("Split-Half", "Cross-Dataset", "Null")),
                               c("Split-Half", "Cross-Dataset", "Null"))
    col_order_idx <- order(group_priority[col_groups], seq_along(col_groups))
    fdp_mat <- fdp_mat[, col_order_idx, drop = FALSE]
    available_cols <- available_cols[col_order_idx]
    col_groups <- col_groups[col_order_idx]

    # Abbreviate row labels to A/M/O
    pretty_rows <- prettify_label(rownames(fdp_mat))
    rownames(fdp_mat) <- substr(pretty_rows, 1, 1)
    colnames(fdp_mat) <- prettify_label(colnames(fdp_mat))

    col_annotation <- data.frame(
      Group = factor(col_groups,
                     levels = c("Split-Half", "Cross-Dataset", "Null")),
      row.names = colnames(fdp_mat)
    )

    group_runs <- rle(as.character(col_groups))
    gap_positions <- cumsum(group_runs$lengths[-length(group_runs$lengths)])

    png(file.path(OUTPUT_BASE, "FDP_Heatmap.png"),
        width = 1800, height = 700, res = 150)
    pheatmap(fdp_mat,
             display_numbers = TRUE,
             number_format   = "%.2f",
             fontsize_number = 13,
             main = glue("False Discovery Proportion vs {REFERENCE_DATASET} Full Ground Truth"),
             color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#FC8D59", "#D73027"))(100),
             breaks = seq(0, 1, length.out = 101),
             annotation_col  = col_annotation,
             annotation_colors = ann_colors,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             gaps_col = gap_positions,
             angle_col = 315,
             fontsize = 15,
             fontsize_row = 15, fontsize_col = 14,
             annotation_legend = TRUE,
             annotation_names_col = FALSE)
    dev.off()
    log_info("Saved FDP_Heatmap.png")

  } else {
    log_warn("No FDP data collected. Skipping FDP heatmap.")
  }
}, error = function(e) {
  log_error("FDP heatmap generation failed: %s", e$message)
  try(dev.off(), silent = TRUE)
})

log_section("Step 2b Complete")
log_info("Regenerated: 3 IoU raw + 3 IoU universe-controlled + 1 FDP heatmap")
