# ==============================================================================
# SCRIPT: Figure T / Figure S — small-effect "figure of 4"
# PURPOSE: For each of the two target small-effect LFC bins, show how the metric
#          (recall / precision) varies across the donors x cells|reads grid,
#          faceted by the four BaseMean expression categories (a 2x2 "figure of
#          four"). Reveals whether power dips as expression decreases within a
#          fixed small-effect LFC band.
# INPUT:   df_finecross_summary inside Aggregated_Summary_Data.qs (Step 3 for
#          per-method; Step 8 for the pooled all_methods checkpoint).
# OUTPUT:  PNGs under <plots_dir>/Figure_TS/ for the current DOWNSAMPLE_MODE.
#
# Layout (decision, handover §8.2): separate precision + recall figures (4 per
# view = 2 LFC bins x 2 metrics); both donors x cells and donors x reads views
# are produced by running this once per DOWNSAMPLE_MODE. Rendered for each
# per-(cell_type, method) checkpoint AND for the cross-method all_methods pool.
# Plotting conventions mirror 04_Visualise_Power.R (magma fill 0-1, 0.8
# threshold overlay, adaptive text colour, 14x10).
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)

  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required. Install: install.packages('ggnewscale')")
  }
  library(ggnewscale)
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
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 9: Figure T / Figure S (small-effect figure of 4)")

y_col <- get_y_col_name()
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Counts per Donor" else "Cells per Donor"
clean_y_factor <- if (DOWNSAMPLE_MODE == "reads") clean_reads_factor else clean_cells_factor
log_info("Mode: %s | Y-column: %s", DOWNSAMPLE_MODE, y_col)

# Figure label per target LFC bin: "Tiny (0.2-1.0)" -> Figure T, etc.
figure_tag_for_bin <- function(bin_name) {
  first_word <- toupper(substr(trimws(bin_name), 1, 1))
  sprintf("Figure_%s", first_word)
}

# ==============================================================================
# 2. PLOTTING FUNCTION (self-contained; mirrors 04_Visualise_Power.R)
# ==============================================================================

#' 2x2 BaseMean-faceted metric heatmap for a fixed LFC bin.
#'
#' @param data   df with Donors, y_col, BaseMean_Bin, and val_col (factors on
#'               the two axes already applied, grid completed).
#' @param val_col Metric column ("Mean_Recall" or "Mean_Precision").
plot_ts_facet <- function(data, val_col, title_str, y_col_name = y_col,
                          y_axis_label = y_label) {

  # Facet label is the bare BaseMean category. Gene counts are intentionally
  # omitted: paper convention is NO "(n = x genes)" annotations on these grids.
  data <- data %>%
    mutate(Facet_Display_Label = as.character(BaseMean_Bin))

  # Order facets by BASEMEAN_INTERVALS definition order.
  ordered_labels <- unique(
    data$Facet_Display_Label[order(match(data$BaseMean_Bin, names(BASEMEAN_INTERVALS)))]
  )
  data$Facet_Display_Label <- factor(data$Facet_Display_Label, levels = ordered_labels)

  # Grey out entire facets that carry essentially no ground-truth genes (mean
  # gene count below threshold): metrics for a near-empty bin (e.g. High
  # expression within a small-LFC bin) are misleading, so blank the whole facet
  # to a uniform grey rather than showing spurious values.
  if ("Ref_GT_Count" %in% names(data)) {
    facet_gene_means <- tapply(data$Ref_GT_Count, data$Facet_Display_Label,
                               function(v) mean(v, na.rm = TRUE))
    empty_facets <- names(facet_gene_means)[is.nan(facet_gene_means) |
                                              facet_gene_means < EMPTY_FACET_GENE_THRESHOLD]
    if (length(empty_facets) > 0) {
      data[[val_col]][as.character(data$Facet_Display_Label) %in% empty_facets] <- NA
    }
  }

  metric_label <- if (grepl("Recall", val_col)) "Recall" else "Precision"

  p <- ggplot(data, aes(x = Donors, y = !!sym(y_col_name))) +
    geom_tile(aes(fill = !!sym(val_col)), color = "white") +
    facet_wrap(~Facet_Display_Label, ncol = 2) +
    labs(title = title_str, x = "Donors", y = y_axis_label) +
    theme_publication(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 13),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "right",
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(1.3, "cm"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    ) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50",
                         name = sub("^Mean_", "Mean ", val_col),
                         guide = guide_colorbar(order = 1))

  # 0.8 threshold overlay (cyan border on cells that reach 80%).
  data <- data %>%
    mutate(Threshold_Met = factor(
      ifelse(!is.na(!!sym(val_col)) & !!sym(val_col) >= 0.8,
             paste0("\u2265 80% ", metric_label),
             paste0("< 80% ", metric_label)),
      levels = c(paste0("\u2265 80% ", metric_label), paste0("< 80% ", metric_label))
    ))
  any_threshold_met <- any(!is.na(data[[val_col]]) & data[[val_col]] >= 0.8)

  p <- p +
    new_scale_color() +
    geom_tile(data = data, aes(color = Threshold_Met), fill = NA, linewidth = 0.6,
              width = 0.85, height = 0.85) +
    scale_color_manual(
      name = "Threshold",
      values = setNames(c("#00FFFF", "grey70"),
                        c(paste0("\u2265 80% ", metric_label),
                          paste0("< 80% ", metric_label))),
      limits = c(paste0("\u2265 80% ", metric_label), paste0("< 80% ", metric_label)),
      drop = FALSE,
      guide = if (any_threshold_met) {
        guide_legend(order = 2, override.aes = list(fill = c(NA, "grey85"),
                     color = c("#00FFFF", "grey70"), linewidth = 1.2))
      } else {
        "none"
      }
    )

  # Adaptive text colour (dark on light tiles, white on dark).
  data$.text_col <- ifelse(is.na(data[[val_col]]) | data[[val_col]] > 0.45, "black", "white")
  for (tc in unique(data$.text_col)) {
    sub <- data[data$.text_col == tc, ]
    if (nrow(sub) == 0) next
    p <- p + geom_text(data = sub,
                       aes(label = ifelse(is.na(!!sym(val_col)), "",
                                          round(!!sym(val_col) * 100))),
                       color = tc, size = 3.0, show.legend = FALSE)
  }

  p
}

# ==============================================================================
# 3. PROCESS EACH CELL TYPE x METHOD (per-method + all_methods)
# ==============================================================================

cell_types <- get_cell_types()
viz_methods <- get_viz_methods()
log_info("Cell types: %s | methods: %s",
         paste(cell_types, collapse = ", "), paste(viz_methods, collapse = ", "))

metrics_to_plot <- list(
  list(col = "Mean_Recall",    desc = "Recall"),
  list(col = "Mean_Precision", desc = "Precision")
)

for (cell_type in cell_types) {
  for (method in viz_methods) {

    paths <- get_paths(cell_type, method = method)
    checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
    if (!file.exists(checkpoint_file)) {
      log_warning("No checkpoint for %s/%s — skipping Figure T/S", cell_type, method)
      next
    }

    cp <- qs::qread(checkpoint_file)
    df_fc <- cp$df_finecross_summary
    if (is.null(df_fc) || nrow(df_fc) == 0) {
      log_warning("No df_finecross_summary for %s/%s — skipping (re-run Step 3/8)", cell_type, method)
      next
    }

    # Paper figure uses the Total gene category.
    df_fc <- df_fc %>% filter(Gene_Category == "Total")
    if (nrow(df_fc) == 0) next

    # Axis factors + complete the grid so missing cells render grey.
    df_fc$Donors <- clean_donors_factor(df_fc$Donors)
    df_fc[[y_col]] <- clean_y_factor(df_fc[[y_col]])
    df_fc <- df_fc %>%
      complete(Donors, !!sym(y_col), LFC_Bin, BaseMean_Bin)

    fig_dir <- file.path(paths$plots_dir, "Figure_TS")
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

    for (fc_lfc in FINECROSS_LFC_BINS) {
      df_lfc <- df_fc %>% filter(LFC_Bin == fc_lfc)
      if (nrow(df_lfc) == 0) {
        log_warning("No fine-cross rows for LFC bin '%s' (%s/%s)", fc_lfc, cell_type, method)
        next
      }
      fig_tag <- figure_tag_for_bin(fc_lfc)

      for (m in metrics_to_plot) {
        title_str <- sprintf("%s across expression (%s): %s [%s]",
                             m$desc, fc_lfc, prettify_label(cell_type), method)
        p <- plot_ts_facet(df_lfc, m$col, title_str)

        fname <- sprintf("%s_%s_%s.png", fig_tag, m$desc, DOWNSAMPLE_MODE)
        ggsave(file.path(fig_dir, fname), p, width = 14, height = 10, dpi = 150)
        log_info("Saved: %s/%s", basename(fig_dir), fname)
      }
    }
  }
}

log_section("STEP 9 COMPLETE")
log_info("Figure T / Figure S generation complete (%s view)", DOWNSAMPLE_MODE)
