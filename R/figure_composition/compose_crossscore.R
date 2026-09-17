# =============================================================================
# Cross-method 3x3 truth-scoring visualiser.
#
# Renders the 3x3 cross-scoring matrix produced by
#   R/cartesian_downsampling/10_CrossScore_Methods.R
# Rows = subsample method (i), cols = truth method (j). Diagonal (i == j)
# reproduces self-scored numbers; off-diagonal asks "do method i's subsamples
# recover method j's truth?".
#
# Style matches R/cartesian_downsampling/04_Visualise_Power.R (discrete) and
# 07_Visualise_Continuous_Grid.R (continuous) so cross-score panels compose
# cleanly with the rest of the paper's figures.
#
# INPUT : Datalink/Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv
# OUTPUT: R/figure_composition/output_crossscore/<name>.png (+ .txt caption)
# PUBLISH: Datalink/Composed_Figures/Cross_Score/<name>.png (+ .txt)
#
# FIGURES  (recall and precision share the SAME layouts — precision is now
#           per-LFC-bin, call-anchored, exactly like every other precision
#           figure in the paper; see 10_CrossScore_Methods.R header):
#   mega   : outer 3x3 (subsample x truth), inner 2x2 of LFC bins.
#            No per-cell labels; cyan tile outline on cells >= 80% (discrete)
#            or 80% cyan contour line (continuous). Explicit row/col strips.
#   perbin : plain 3x3 of method pairs for one LFC bin. Adaptive per-cell
#            value labels (percent) on discrete; 80% cyan contour on
#            continuous. Cyan tile outline on cells >= 80% (discrete).
#
# All figures are TITLE-FREE. Each <name>.png ships with a sibling <name>.txt
# containing a paper-ready caption.
#
# Run:
#   Rscript R/figure_composition/compose_crossscore.R                         # everything, cells view
#   Rscript R/figure_composition/compose_crossscore.R --metric recall --figset mega
#   Rscript R/figure_composition/compose_crossscore.R --mode reads --style continuous
#   Rscript R/figure_composition/compose_crossscore.R --celltypes Oligo
#   Rscript R/figure_composition/compose_crossscore.R --datalink /path/to/Datalink
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  if (requireNamespace("ggnewscale", quietly = TRUE)) library(ggnewscale)
  if (requireNamespace("shadowtext", quietly = TRUE)) library(shadowtext)
})

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------

HERE <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE))
if (is.na(HERE) || HERE == "") {
  HERE <- normalizePath(".", mustWork = FALSE)
}
PROJECT_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)
DEFAULT_DATALINK <- file.path(PROJECT_ROOT, "Datalink")
OUT_ROOT <- file.path(HERE, "output_crossscore")

source(file.path(HERE, "interp_helpers.R"))

# Reuse the paper's theme_publication if available; otherwise fall back.
.prettify_path <- file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R")
if (file.exists(.prettify_path)) source(.prettify_path)
if (!exists("theme_publication")) {
  theme_publication <- function(base_size = 12) theme_bw(base_size = base_size)
}

METHODS      <- c("deseq2", "edger", "limma")
METHOD_LABEL <- c(deseq2 = "DESeq2", edger = "edgeR", limma = "limma-voom")
CELLTYPES    <- c("Astro", "Micro", "Oligo")
CELL_LABEL   <- c(Astro = "Astrocytes", Micro = "Microglia",
                  Oligo = "Oligodendrocytes")

# Bin order top-to-bottom on facets: biggest effect first (matches 04/07).
LFC_BINS  <- c("Large", "Medium", "Small", "Tiny")
LFC_LABEL <- c(Tiny   = "|LFC| 0.2\u20131",
               Small  = "|LFC| 1\u20132",
               Medium = "|LFC| 2\u20133",
               Large  = "|LFC| \u2265 3")

DIAG_COLOR      <- "#00c2ff"    # panel outline for self-scored (i == j)
THRESHOLD_COLOR <- "#00FFFF"    # cell outline / contour for >= 80%
CMAP_OPT        <- "magma"      # viridis option name

REQUIRED_COLS <- c("subsample_method", "truth_method", "cell_type", "mode",
                   "n_donors", "y_value", "lfc_bin", "metric",
                   "value")

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------

parse_args <- function(argv) {
  defaults <- list(
    datalink   = DEFAULT_DATALINK,
    metric     = "both",       # recall | precision | both
    mode       = "cells",      # cells | reads
    style      = "both",       # discrete | continuous | both
    value      = "raw",        # retained for CLI compatibility; raw only
    figset     = "both",       # mega | perbin | both
    celltypes  = paste(CELLTYPES, collapse = ",")
  )
  it <- 1
  while (it <= length(argv)) {
    a <- argv[it]
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    key <- kv[1]
    if (length(kv) == 2) {
      val <- kv[2]
    } else if (it < length(argv)) {
      val <- argv[it + 1]
      it  <- it + 1
    } else {
      stop("missing value for ", key)
    }
    nm <- sub("^--", "", key)
    if (!nm %in% names(defaults)) stop("unknown argument: ", key)
    defaults[[nm]] <- val
    it <- it + 1
  }
  validate <- function(name, allowed) {
    if (!defaults[[name]] %in% allowed) {
      stop(sprintf("--%s must be one of %s", name,
                   paste(allowed, collapse = "/")))
    }
  }
  validate("metric", c("recall", "precision", "both"))
  validate("mode",   c("cells", "reads"))
  validate("style",  c("discrete", "continuous", "both"))
  validate("value",  c("raw"))
  validate("figset", c("mega", "perbin", "both"))
  defaults$datalink  <- normalizePath(defaults$datalink, mustWork = FALSE)
  defaults$celltypes <- trimws(strsplit(defaults$celltypes, ",",
                                        fixed = TRUE)[[1]])
  bad <- setdiff(defaults$celltypes, CELLTYPES)
  if (length(bad) > 0) stop("unknown cell types: ", paste(bad, collapse = ","))
  defaults
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

load_matrix <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop("cross_score_matrix.csv not found at ", csv_path,
         " - run 10_CrossScore_Methods.R first (cartesian cross-score phase).")
  }
  cat(sprintf("[load] reading %s\n", csv_path))
  dt <- fread(csv_path)
  missing <- setdiff(REQUIRED_COLS, names(dt))
  if (length(missing) > 0) {
    stop("cross_score_matrix.csv missing columns: ",
         paste(missing, collapse = ","))
  }
  dt <- dt[metric %in% c("recall", "precision")]
  dt[, n_donors := as.numeric(n_donors)]
  dt[, y_value  := as.numeric(y_value)]
  dt
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

y_axis_label <- function(mode) {
  if (mode == "reads") "Target Counts per Donor" else "Cells per Donor"
}

y_labels_fn <- function(mode) {
  if (mode == "reads") scientific_power_labels else waiver()
}

pick_fill_col <- function(sub, value_col) {
  "value"
}

# Write a paper-ready caption sibling for a figure.
write_caption <- function(png_path, caption_text) {
  txt_path <- sub("\\.png$", ".txt", png_path)
  writeLines(caption_text, txt_path)
  invisible(txt_path)
}

# Save a plot to output_crossscore/ and mirror to Datalink/Composed_Figures/Cross_Score/.
finalise <- function(plt, out_name, width, height, datalink,
                     caption_text, dpi = 150) {
  dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(OUT_ROOT, out_name)
  ggsave(out_path, plt, width = width, height = height, dpi = dpi,
         bg = "white", limitsize = FALSE)
  write_caption(out_path, caption_text)
  pub_dir <- file.path(datalink, "Composed_Figures", "Cross_Score")
  dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(out_path, file.path(pub_dir, basename(out_path)), overwrite = TRUE)
  file.copy(sub("\\.png$", ".txt", out_path),
            file.path(pub_dir, sub("\\.png$", ".txt", basename(out_path))),
            overwrite = TRUE)
  out_path
}

# ---------------------------------------------------------------------------
# Discrete panel: uniform tile grid matching 04_ style.
# `sub`         : data.table filtered to one panel's rows (or a facetable set)
# `value_col`   : retained for CLI compatibility; recall and precision both use value
# `facet_by`    : NULL | "lfc_bin"; adds facet_wrap when set
# `show_labels` : TRUE for perbin/precision panels; FALSE for mega inner cells
# `mode`        : "cells" or "reads"; drives Y label formatter
# =============================================================================
build_discrete_panel <- function(sub, value_col, facet_by, show_labels, mode,
                                 base_size = 12, tile_lw = 0.4,
                                 outline_lw = 0.6, text_size = 3.0) {
  fill_col <- pick_fill_col(sub, value_col)
  sub[, .fill := as.numeric(get(fill_col))]
  sub[, .donor_f := factor(n_donors, levels = sort(unique(n_donors)))]
  sub[, .y_f     := factor(y_value,  levels = sort(unique(y_value)))]

  p <- ggplot(sub, aes(x = .donor_f, y = .y_f)) +
    geom_tile(aes(fill = .fill), colour = "white", linewidth = tile_lw) +
    scale_fill_viridis_c(option = CMAP_OPT, limits = c(0, 1),
                         na.value = "grey50", name = NULL,
                         guide = guide_colorbar(order = 1))

  if (!is.null(facet_by) && facet_by %in% names(sub)) {
    sub[, .facet := factor(get(facet_by), levels = LFC_BINS,
                           labels = LFC_LABEL[LFC_BINS])]
    p <- p + facet_wrap(~ .facet, ncol = 2, drop = FALSE)
  }

  # Cyan outline for tiles at or above the 80% threshold.
  hi <- sub[!is.na(.fill) & .fill >= 0.8]
  if (nrow(hi) > 0) {
    p <- p +
      (if (exists("new_scale_color")) new_scale_color() else NULL) +
      geom_tile(data = hi, aes(x = .donor_f, y = .y_f),
                fill = NA, colour = THRESHOLD_COLOR,
                width = 0.90, height = 0.90, linewidth = outline_lw,
                inherit.aes = FALSE)
  }

  if (show_labels) {
    sub[, .text_col := ifelse(is.na(.fill) | .fill > 0.45, "black", "white")]
    for (tc in unique(sub$.text_col)) {
      s <- sub[.text_col == tc]
      if (nrow(s) == 0) next
      p <- p + geom_text(data = s,
        aes(x = .donor_f, y = .y_f,
            label = ifelse(is.na(.fill), "", round(.fill * 100))),
        colour = tc, size = text_size, show.legend = FALSE, inherit.aes = FALSE)
    }
  }

  p +
    theme_publication(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.15, "lines"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right"
    )
}

# ---------------------------------------------------------------------------
# Continuous panel: bilinear-interpolated surface with concave frontier mask
# and optional 80% contour. Matches 07_.
# =============================================================================
build_continuous_panel <- function(sub, value_col, facet_by, mode,
                                   base_size = 12, contour_lw = 0.9) {
  fill_col <- pick_fill_col(sub, value_col)
  sub_local <- copy(sub)
  sub_local[, .fill := as.numeric(get(fill_col))]

  # Interpolate per-facet if facet_by set, else once.
  interp_list <- if (!is.null(facet_by) && facet_by %in% names(sub_local)) {
    split(sub_local, sub_local[[facet_by]], drop = TRUE)
  } else {
    list(all = sub_local)
  }

  all_interp <- rbindlist(lapply(names(interp_list), function(nm) {
    part <- interp_list[[nm]]
    interp <- interpolate_crossscore_panel(part, val_col = ".fill",
                                           n_interp = 180)
    if (is.null(interp)) return(NULL)
    if (!is.null(facet_by)) interp[[facet_by]] <- nm
    interp
  }), fill = TRUE)

  if (nrow(all_interp) == 0) return(NULL)

  if (!is.null(facet_by)) {
    all_interp$.facet <- factor(all_interp[[facet_by]],
                                levels = LFC_BINS,
                                labels = LFC_LABEL[LFC_BINS])
  }

  p <- ggplot(all_interp,
              aes(x = Donors, y = Y_Value, fill = .value)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_viridis_c(option = CMAP_OPT, limits = c(0, 1),
                         na.value = "grey50", name = NULL,
                         guide = guide_colorbar(order = 1)) +
    coord_cartesian(expand = FALSE) +
    scale_y_continuous(labels = y_labels_fn(mode))

  # 80% contour where the surface actually crosses the threshold.
  crosses_threshold <- function(df) {
    v <- df$.value[is.finite(df$.value)]
    length(v) > 0 && min(v) < 0.8 && max(v) > 0.8
  }
  if (!is.null(facet_by)) {
    contour_data <- all_interp[, if (crosses_threshold(.SD)) .SD else NULL,
                               by = .facet]
  } else {
    contour_data <- if (crosses_threshold(all_interp)) all_interp else all_interp[0]
  }
  if (nrow(contour_data) > 0) {
    p <- p + geom_contour(data = contour_data,
                          aes(x = Donors, y = Y_Value, z = .value),
                          breaks = c(0.8), colour = THRESHOLD_COLOR,
                          linewidth = contour_lw, na.rm = TRUE,
                          inherit.aes = FALSE)
  }

  if (!is.null(facet_by)) {
    p <- p + facet_wrap(~ .facet, ncol = 2, drop = FALSE)
  }

  p +
    theme_publication(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.6, "lines"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right"
    )
}

# ---------------------------------------------------------------------------
# Row / column label strip helpers (used by mega layout).
# ---------------------------------------------------------------------------

label_panel <- function(text, angle = 0, size = 5.5) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = text,
             angle = angle, size = size, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------

# MEGA: outer 3x3 (subsample x truth) of method pairs, each inner 2x2 (LFC bins).
render_mega <- function(dt, cell, mode, style, value_col, datalink,
                        metric_sel = "recall") {
  mode_sel <- mode
  sub <- dt[metric == metric_sel & cell_type == cell & mode == mode_sel]
  if (nrow(sub) == 0) {
    return(list(path = NA, info = sprintf("no %s grids for %s/%s",
                                          metric_sel, cell, mode)))
  }
  is_prec <- metric_sel == "precision"

  build_inner <- function(cell_dt) {
    if (nrow(cell_dt) == 0) {
      return(ggplot() +
               annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = -Inf, ymax = Inf, fill = "grey85") +
               theme_void())
    }
    if (style == "continuous") {
      p <- build_continuous_panel(cell_dt, value_col,
                                  facet_by = "lfc_bin", mode = mode,
                                  base_size = 8, contour_lw = 0.6)
      if (is.null(p)) {
        return(ggplot() +
                 annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = -Inf, ymax = Inf, fill = "grey85") +
                 theme_void())
      }
    } else {
      p <- build_discrete_panel(cell_dt, value_col,
                                facet_by = "lfc_bin", show_labels = FALSE,
                                mode = mode, base_size = 8,
                                tile_lw = 0.15, outline_lw = 0.45)
    }
    p +
      theme(
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 6.5,
                                  margin = margin(0.6, 0, 0.6, 0)),
        legend.position = "none",
        plot.margin = margin(1, 1, 1, 1)
      )
  }

  # Build the 9 inner panels, marking the diagonal.
  inners <- vector("list", 9)
  k <- 1
  for (i_sub in METHODS) {
    for (j_truth in METHODS) {
      cell_dt <- sub[subsample_method == i_sub & truth_method == j_truth]
      p <- build_inner(cell_dt)
      if (i_sub == j_truth) {
        p <- p + theme(plot.background = element_rect(colour = DIAG_COLOR,
                                                      fill = NA, linewidth = 1.4))
      }
      inners[[k]] <- p
      k <- k + 1
    }
  }

  # Label strips.
  corner <- label_panel("", angle = 0)
  col_lbls <- lapply(METHODS, function(m)
    label_panel(paste0("truth: ", METHOD_LABEL[[m]]), angle = 0, size = 5.5))
  row_lbls <- lapply(METHODS, function(m)
    label_panel(paste0("subsample: ", METHOD_LABEL[[m]]), angle = 90,
                size = 5.5))

  # 4x4 layout: row 0 = corner + col labels; rows 1..3 = row label + 3 panels.
  # Pass elements in reading order (row-major) so patchwork's ncol/byrow places
  # them correctly, avoiding any ambiguity with `design` letter mapping.
  els <- list(
    corner,       col_lbls[[1]], col_lbls[[2]], col_lbls[[3]],
    row_lbls[[1]], inners[[1]],  inners[[2]],  inners[[3]],
    row_lbls[[2]], inners[[4]],  inners[[5]],  inners[[6]],
    row_lbls[[3]], inners[[7]],  inners[[8]],  inners[[9]]
  )
  combined <- wrap_plots(els, ncol = 4, nrow = 4, byrow = TRUE) +
    plot_layout(widths  = c(0.10, 1, 1, 1),
                heights = c(0.10, 1, 1, 1))

  metric_word <- if (is_prec) "precision" else "recall"
  out_name <- sprintf("crossscore_%s_mega_%s_%s_%s_%s.png",
                      metric_word, cell, mode, style, value_col)
  bin_note <- if (is_prec) {
    "Precision is per |LFC| bin, call-anchored: called genes are binned by their OWN estimated |LFC| and scored against the full truth set (all bins), identical to the standard per-bin precision used elsewhere in the paper. "
  } else {
    ""
  }
  cap <- sprintf(
    "Cross-method %s in %s (%s panels; %s). Each outer cell shows one (subsample method, truth method) pair: rows are the DE method whose downsampled calls are being scored, columns are the DE method providing the full-data truth set. Each inner 2x2 shows the four |LFC| bins (top row: |LFC| >= 3 and |LFC| 2-3; bottom row: |LFC| 1-2 and |LFC| 0.2-1). %sColour scale is %s on [0, 1] (magma). Cyan panel outlines mark self-scored diagonals. Cyan %s marks the 80%% %s threshold. Y axis: %s; X axis: number of donors. Grey cells: missing or non-computable grid points.",
    metric_word, CELL_LABEL[[cell]], style,
    if (mode == "reads") "target reads per donor" else "cells per donor",
    bin_note, metric_word,
    if (style == "continuous") "contour" else "cell outlines",
    metric_word,
    y_axis_label(mode))
  path <- finalise(combined, out_name, width = 16, height = 16, datalink, cap)
  list(path = path, info = out_name)
}

# Plain 3x3 grid via facet_grid (used for perbin recall and precision).
render_grid_3x3 <- function(sub, cell, mode, style, value_col, out_name,
                            caption, datalink, is_precision = FALSE) {
  if (nrow(sub) == 0) {
    return(list(path = NA, info = sprintf("no data for %s", out_name)))
  }
  sub[, subsample_method := factor(subsample_method, levels = METHODS,
                                    labels = METHOD_LABEL[METHODS])]
  sub[, truth_method     := factor(truth_method, levels = METHODS,
                                    labels = METHOD_LABEL[METHODS])]

  fill_col <- pick_fill_col(sub, value_col)
  sub[, .fill := as.numeric(get(fill_col))]

  if (style == "continuous") {
    parts <- split(sub, list(sub$subsample_method, sub$truth_method),
                   drop = TRUE, sep = "||")
    all_interp <- rbindlist(lapply(names(parts), function(nm) {
      part <- parts[[nm]]
      interp <- interpolate_crossscore_panel(part, val_col = ".fill",
                                             n_interp = 180)
      if (is.null(interp)) return(NULL)
      pieces <- strsplit(nm, "||", fixed = TRUE)[[1]]
      interp$subsample_method <- factor(pieces[1],
                                        levels = METHOD_LABEL[METHODS])
      interp$truth_method     <- factor(pieces[2],
                                        levels = METHOD_LABEL[METHODS])
      interp
    }), fill = TRUE)

    if (nrow(all_interp) == 0) {
      return(list(path = NA, info = sprintf("no interpolable data for %s",
                                            out_name)))
    }

    p <- ggplot(all_interp, aes(x = Donors, y = Y_Value, fill = .value)) +
      geom_raster(interpolate = FALSE) +
      scale_fill_viridis_c(option = CMAP_OPT, limits = c(0, 1),
                           na.value = "grey50",
                           name = if (is_precision) "Precision" else "Recall") +
      coord_cartesian(expand = FALSE) +
      scale_y_continuous(labels = y_labels_fn(mode)) +
      facet_grid(subsample_method ~ truth_method, switch = "y")

    # 80% contour per facet where surface crosses threshold.
    fac <- all_interp[, list(cross = min(.value, na.rm = TRUE) < 0.8 &
                                       max(.value, na.rm = TRUE) > 0.8),
                      by = list(subsample_method, truth_method)]
    contour_data <- merge(all_interp,
                          fac[cross == TRUE,
                              list(subsample_method, truth_method)],
                          by = c("subsample_method", "truth_method"))
    if (nrow(contour_data) > 0) {
      p <- p + geom_contour(data = contour_data,
                            aes(x = Donors, y = Y_Value, z = .value),
                            breaks = c(0.8), colour = THRESHOLD_COLOR,
                            linewidth = 0.9, na.rm = TRUE, inherit.aes = FALSE)
    }
  } else {
    sub[, .donor_f := factor(n_donors, levels = sort(unique(n_donors)))]
    sub[, .y_f     := factor(y_value,  levels = sort(unique(y_value)))]

    p <- ggplot(sub, aes(x = .donor_f, y = .y_f)) +
      geom_tile(aes(fill = .fill), colour = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = CMAP_OPT, limits = c(0, 1),
                           na.value = "grey50",
                           name = if (is_precision) "Precision" else "Recall") +
      facet_grid(subsample_method ~ truth_method, switch = "y")

    hi <- sub[!is.na(.fill) & .fill >= 0.8]
    if (nrow(hi) > 0) {
      p <- p +
        (if (exists("new_scale_color")) new_scale_color() else NULL) +
        geom_tile(data = hi,
                  aes(x = .donor_f, y = .y_f),
                  fill = NA, colour = THRESHOLD_COLOR,
                  width = 0.90, height = 0.90, linewidth = 0.5,
                  inherit.aes = FALSE)
    }

    sub[, .text_col := ifelse(is.na(.fill) | .fill > 0.45, "black", "white")]
    for (tc in unique(sub$.text_col)) {
      s <- sub[.text_col == tc]
      if (nrow(s) == 0) next
      p <- p + geom_text(data = s,
        aes(x = .donor_f, y = .y_f,
            label = ifelse(is.na(.fill), "", round(.fill * 100))),
        colour = tc, size = 2.6, show.legend = FALSE, inherit.aes = FALSE)
    }
  }

  # Cyan panel outline on the diagonal (data columns route via facet_grid).
  diag_df <- data.frame(
    subsample_method = factor(METHOD_LABEL[METHODS],
                              levels = METHOD_LABEL[METHODS]),
    truth_method     = factor(METHOD_LABEL[METHODS],
                              levels = METHOD_LABEL[METHODS])
  )
  p <- p + geom_rect(data = diag_df,
                     xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                     fill = NA, colour = DIAG_COLOR, linewidth = 1.2,
                     inherit.aes = FALSE)

  p <- p +
    labs(x = "Number of Donors", y = y_axis_label(mode)) +
    theme_publication(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      strip.placement = "outside",
      panel.spacing = unit(0.15, "lines"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right"
    )

  path <- finalise(p, out_name, width = 14, height = 14, datalink, caption)
  list(path = path, info = out_name)
}

render_perbin <- function(dt, cell, mode, lfc, style, value_col, datalink,
                          metric_sel = "recall") {
  mode_sel <- mode
  sub <- dt[metric == metric_sel & cell_type == cell & mode == mode_sel
            & lfc_bin == lfc]
  is_prec <- metric_sel == "precision"
  metric_word <- if (is_prec) "precision" else "recall"
  out_name <- sprintf("crossscore_%s_%s_%s_%s_%s_%s.png",
                      metric_word, lfc, cell, mode, style, value_col)
  anchor_note <- if (is_prec) {
    sprintf("Cell values (discrete panels) are precision percentages: of the genes each method calls with its OWN estimated |LFC| in the %s band, the fraction that lie in the column method's full truth set (all bins) — the standard per-bin precision used elsewhere in the paper.",
            LFC_LABEL[[lfc]])
  } else {
    "Cell values (discrete panels) are recall percentages."
  }
  self_note <- if (is_prec) {
    "Diagonal panels (cyan outline) are self-scored."
  } else {
    "Diagonal panels (cyan outline) are self-scored and reproduce the standard within-method recall grid."
  }
  cap <- sprintf(
    "Cross-method %s for %s genes (%s, %s panels). Rows are the DE method whose downsampled calls are being scored, columns are the DE method providing the full-data truth set. %s %s marks the 80%% threshold. %s Y axis: %s.",
    metric_word, LFC_LABEL[[lfc]], CELL_LABEL[[cell]], style,
    anchor_note,
    if (style == "continuous") "The cyan contour" else "Cyan tile outlines",
    self_note,
    y_axis_label(mode))
  render_grid_3x3(sub, cell, mode, style, value_col, out_name, cap,
                  datalink, is_precision = is_prec)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  argv <- argv[argv != "--rebuild-cache"]  # CLI parity with previous versions
  opts <- parse_args(argv)

  csv_path <- file.path(opts$datalink, "Cartesian_Downsampling", "Cross_Score",
                        "cross_score_matrix.csv")
  dt <- load_matrix(csv_path)

  cells  <- opts$celltypes
  styles <- if (opts$style == "both") c("discrete", "continuous") else opts$style
  do_rec <- opts$metric %in% c("recall", "both")
  do_prc <- opts$metric %in% c("precision", "both")

  built <- 0L
  skipped <- 0L
  log_lines <- character(0)
  record <- function(res) {
    if (is.na(res$path)) {
      skipped <<- skipped + 1L
      msg <- paste0("[SKIP] ", res$info)
      message(msg)
      log_lines <<- c(log_lines, msg)
    } else {
      built <<- built + 1L
      kb <- file.info(res$path)$size %/% 1024
      msg <- sprintf("[OK] %s (%d KB)", basename(res$path), kb)
      cat(msg, "\n", sep = "")
      log_lines <<- c(log_lines, msg)
    }
  }

  for (cell in cells) {
    for (style in styles) {
      if (do_rec) {
        if (opts$figset %in% c("mega", "both")) {
          record(render_mega(dt, cell, opts$mode, style, opts$value,
                             opts$datalink, metric_sel = "recall"))
        }
        if (opts$figset %in% c("perbin", "both")) {
          for (b in c("Large", "Medium", "Small", "Tiny")) {
            record(render_perbin(dt, cell, opts$mode, b, style, opts$value,
                                 opts$datalink, metric_sel = "recall"))
          }
        }
      }
      if (do_prc) {
        if (opts$figset %in% c("mega", "both")) {
          record(render_mega(dt, cell, opts$mode, style, opts$value,
                             opts$datalink, metric_sel = "precision"))
        }
        if (opts$figset %in% c("perbin", "both")) {
          for (b in c("Large", "Medium", "Small", "Tiny")) {
            record(render_perbin(dt, cell, opts$mode, b, style, opts$value,
                                 opts$datalink, metric_sel = "precision"))
          }
        }
      }
    }
  }

  dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
  writeLines(log_lines, file.path(OUT_ROOT, "build_log.txt"))
  cat(sprintf("\nDone: %d built, %d skipped (mode=%s, metric=%s, style=%s, value=%s).\n",
              built, skipped, opts$mode, opts$metric, opts$style, opts$value))
}

if (sys.nframe() == 0) main()
