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
  PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
    Sys.getenv("POWER_ANALYSIS_ROOT")
  } else {
    here::here()
  }
}

source(file.path(PROJECT_ROOT, "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
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
group_colors <- c("Split-Half" = "#2166AC", "Cross-Dataset" = "#4DAF4A", "Null" = "#E41A1C")

# ------------------------------------------------------------------------------
# Helper: render a single IoU heatmap
# ------------------------------------------------------------------------------
render_iou_heatmap <- function(iou_csv, out_png, main_title,
                               number_fontsize = 0.5,
                               label_fontsize = 3.5,
                               show_labels = TRUE) {
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

  n <- nrow(iou_mat)

  # Prettify labels
  pretty_names <- prettify_label(colnames(iou_mat))

  # Build long-format data for lower triangle only (row > col, no diagonal)
  tiles <- do.call(rbind, lapply(2:n, function(i) {
    data.frame(col = 1:(i - 1), row = i, value = iou_mat[i, 1:(i - 1)],
               col_group = col_groups[1:(i - 1)],
               row_group = col_groups[i])
  }))
  # Null-region numbers rendered white so they pop on the blue/grey background
  tiles$text_colour <- ifelse(tiles$col_group == "Null" | tiles$row_group == "Null",
                              "white", "grey25")

  # Group annotation bar — thin strip at row position just above first data row
  group_colors <- c("Split-Half" = "#2166AC", "Cross-Dataset" = "#4DAF4A", "Null" = "#E41A1C")
  ann_df <- data.frame(
    col = seq_len(n - 1), row = 1.15, group = col_groups[1:(n - 1)],
    stringsAsFactors = FALSE
  )

  # Group labels centered on each group's span in the annotation bar
  grp_runs <- rle(col_groups[1:(n - 1)])
  grp_ends <- cumsum(grp_runs$lengths)
  grp_starts <- c(1, grp_ends[-length(grp_ends)] + 1)
  grp_label_df <- data.frame(
    x = (grp_starts + grp_ends) / 2,
    y = 1.15,
    label = grp_runs$values,
    stringsAsFactors = FALSE
  )

  # Row labels: positioned just right of the diagonal cell for each row
  row_labels <- data.frame(
    x = 2:n - 1 + 0.6,
    y = 2:n,
    label = pretty_names[2:n],
    stringsAsFactors = FALSE
  )

  # Colour scale matching the original blue-white-red
  fill_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

  # Scale figure dimensions to matrix size
  cell_size <- 0.55
  fig_w <- n * cell_size + 0.8  # extra for labels (legend inside plot)
  fig_h <- n * cell_size + 2    # extra for title + x labels

  p <- ggplot() +
    # Heatmap tiles
    geom_tile(data = tiles, aes(x = col, y = row, fill = value),
              colour = "white", linewidth = 0.6) +
    # Value labels: non-null region (dark text)
    geom_text(data = tiles[tiles$text_colour != "white", ],
              aes(x = col, y = row,
                  label = ifelse(is.na(value), "",
                                 sprintf("%d", round(value * 100)))),
              size = number_fontsize, colour = "grey25", fontface = "bold") +
    # Value labels: null region (white text)
    geom_text(data = tiles[tiles$text_colour == "white", ],
              aes(x = col, y = row,
                  label = ifelse(is.na(value), "",
                                 sprintf("%d", round(value * 100)))),
              size = number_fontsize, colour = "white", fontface = "bold") +
    # Group annotation bar (no internal borders — continuous coloured strip)
    geom_tile(data = ann_df, aes(x = col, y = row),
              fill = group_colors[ann_df$group],
              height = 0.3, colour = NA) +
    # Group labels on the annotation bar
    geom_text(data = grp_label_df, aes(x = x, y = y, label = label),
              size = 2.8, colour = "white", fontface = "bold") +
    # Row labels right of diagonal
    geom_text(data = row_labels, aes(x = x, y = y, label = label),
              hjust = 0, size = label_fontsize, fontface = "bold") +
    # Fill scale
    scale_fill_gradientn(
      colours = fill_colors,
      limits = c(0, 1),
      name = "IoU"
    ) +
    # X-axis: column labels (only 1 to n-1)
    scale_x_continuous(
      breaks = seq_len(n - 1),
      labels = pretty_names[1:(n - 1)],
      expand = expansion(add = c(0.3, 0.3)),
      position = "bottom"
    ) +
    scale_y_reverse(
      breaks = NULL,
      expand = expansion(add = c(0.3, 0.3))
    ) +
    coord_cartesian(
      clip = "off",
      xlim = c(0.4, n + max(nchar(pretty_names)) * 0.09),
      ylim = c(n + 0.5, 0.7)
    ) +
    labs(title = main_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 10, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.margin = margin(5, 10, 5, 10),
      legend.position = "inside",
      legend.position.inside = c(0.92, 0.08),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
      legend.margin = margin(2, 4, 2, 4),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.35, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  ggsave(out_png, plot = p, width = fig_w, height = fig_h, dpi = 300, bg = "white")

  log_info("  Saved %s (%gx%g in)", basename(out_png), fig_w, fig_h)
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

  # Raw IoU — show row/col labels (rotated down-right) for panel B of FigRepro
  render_iou_heatmap(
    iou_csv    = file.path(diag_dir, "Expanded_IoU_Matrix.csv"),
    out_png    = file.path(diag_dir, "Expanded_IoU_Matrix.png"),
    main_title = glue("IoU Without Universe Control ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
    number_fontsize = 4.4,
    show_labels = TRUE
  )

  # Universe-controlled IoU — also show row/col labels so it lines up with the
  # raw IoU panel when stacked in the supplementary.
  render_iou_heatmap(
    iou_csv    = file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.csv"),
    out_png    = file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.png"),
    main_title = glue("Universe-Controlled IoU ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
    number_fontsize = 4.4,
    show_labels = TRUE
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
  # Fallback: read from combined FDP_Summary.csv if per-cell-type files were missing
  if (length(fdp_collector) == 0) {
    combined_csv <- file.path(OUTPUT_BASE, "FDP_Summary.csv")
    if (file.exists(combined_csv)) {
      log_info("Per-cell-type FDP CSVs missing; reading from %s", combined_csv)
      combined_fdp <- read.csv(combined_csv, stringsAsFactors = FALSE)
    }
  } else {
    combined_fdp <- do.call(rbind, fdp_collector)
    combined_csv <- file.path(OUTPUT_BASE, "FDP_Summary.csv")
    write.csv(combined_fdp, combined_csv, row.names = FALSE, quote = TRUE)
    log_info("Saved combined FDP summary: %s", combined_csv)
  }

  if (exists("combined_fdp") && nrow(combined_fdp) > 0) {

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

    # Group lookup (from first cell type's data — either collected or from combined CSV)
    if (length(fdp_collector) > 0) {
      first_ct <- fdp_collector[[1]]
    } else {
      first_ct <- combined_fdp[combined_fdp$CellType == combined_fdp$CellType[1], ]
    }
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

    # Build FDP heatmap with ggplot2 for full control over group labels
    fdp_long <- do.call(rbind, lapply(seq_len(nrow(fdp_mat)), function(i) {
      data.frame(
        row = rownames(fdp_mat)[i],
        col = colnames(fdp_mat),
        col_idx = seq_len(ncol(fdp_mat)),
        row_idx = i,
        value = fdp_mat[i, ],
        stringsAsFactors = FALSE
      )
    }))

    # Group annotation bar data
    n_fdp_cols <- ncol(fdp_mat)
    fdp_ann_df <- data.frame(
      col_idx = seq_len(n_fdp_cols),
      y = 0.4,
      group = col_groups,
      stringsAsFactors = FALSE
    )

    # Group labels centered on each group's span
    grp_ends_fdp <- cumsum(group_runs$lengths)
    grp_starts_fdp <- c(1, grp_ends_fdp[-length(grp_ends_fdp)] + 1)
    fdp_grp_labels <- data.frame(
      x = (grp_starts_fdp + grp_ends_fdp) / 2,
      y = 0.4,
      label = group_runs$values,
      stringsAsFactors = FALSE
    )

    fdp_fill <- colorRampPalette(c("#F7F7F7", "#FEE08B", "#FC8D59", "#D73027"))(100)

    # Row order: A, M, O (from top to bottom = row_idx 1, 2, 3)
    row_labels_fdp <- rownames(fdp_mat)

    p_fdp <- ggplot() +
      # Heatmap tiles
      geom_tile(data = fdp_long, aes(x = col_idx, y = row_idx, fill = value),
                colour = "white", linewidth = 0.8) +
      # Value labels
      geom_text(data = fdp_long, aes(x = col_idx, y = row_idx,
                                      label = sprintf("%.2f", value)),
                size = 4, colour = "grey15") +
      # Group annotation bar above the heatmap (no internal borders)
      geom_tile(data = fdp_ann_df, aes(x = col_idx, y = y),
                fill = group_colors[fdp_ann_df$group],
                height = 0.35, colour = NA) +
      # Group labels on the bar
      geom_text(data = fdp_grp_labels, aes(x = x, y = y, label = label),
                size = 3, colour = "white", fontface = "bold") +
      # Fill scale
      scale_fill_gradientn(
        colours = fdp_fill,
        limits = c(0, 1),
        name = "FDP"
      ) +
      # X-axis: dataset names
      scale_x_continuous(
        breaks = seq_len(n_fdp_cols),
        labels = colnames(fdp_mat),
        expand = expansion(add = 0.3)
      ) +
      # Y-axis: cell type abbreviations
      scale_y_reverse(
        breaks = seq_along(row_labels_fdp),
        labels = row_labels_fdp,
        expand = expansion(add = c(0.3, 0.3))
      ) +
      coord_cartesian(clip = "off", ylim = c(nrow(fdp_mat) + 0.5, 0.1)) +
      labs(
        title = glue("False Discovery Proportion vs {prettify_label(REFERENCE_DATASET)} Full Ground Truth"),
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 10, 5, 10),
        legend.position = "right",
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )

    fdp_w <- n_fdp_cols * 0.65 + 2.5
    fdp_h <- nrow(fdp_mat) * 0.8 + 3.5
    ggsave(file.path(OUTPUT_BASE, "FDP_Heatmap.png"), plot = p_fdp,
           width = fdp_w, height = fdp_h, dpi = 300, bg = "white")
    log_info("Saved FDP_Heatmap.png (%gx%g in)", fdp_w, fdp_h)

  } else {
    log_warn("No FDP data collected. Skipping FDP heatmap.")
  }
}, error = function(e) {
  log_error("FDP heatmap generation failed: %s", e$message)
  try(dev.off(), silent = TRUE)
})

log_section("Step 2b Complete")
log_info("Regenerated: 3 IoU raw + 3 IoU universe-controlled + 1 FDP heatmap")
