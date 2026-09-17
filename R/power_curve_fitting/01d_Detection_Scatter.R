# ==============================================================================
# SCRIPT: Tiny-bin Detection Scatter (dispersion x baseMean x detection)  — 01d
# PURPOSE: For the Tiny effect-size bin (|LFC| 0.2-1.0), plot every truth gene as
#          a point in (baseMean_raw, dispersion) space, coloured by its detection
#          probability (viridis). The hypothesis: high-dispersion / low-baseMean
#          genes have near-zero power, and power increases as you move along both
#          axes. Rendered at three fixed experiment designs:
#              A. 600 donors x 100 cells
#              B. 100 donors x 3000 cells
#              C. 600 donors x MAX cells (all available cells per donor)
#
# RULERS (shared across methods, so the axes are the same physical quantity for
#         deseq2 / edger / limma):
#   x = gt_basemean = baseMean_raw  (raw pseudobulk mean counts; identical
#       definition across all three DE workers — see utils.R baseMean_raw)
#   y = gt_dispersion = DESeq2 full-data GT dispersion (the only dispersion the
#       pipeline estimates; one GT_DDS per cell type per run, method-independent)
#   colour = detection probability = fraction of subsample iterations (at the
#       fixed grid point, pooled over robustness runs) in which the gene was
#       called significant (padj < DETECT_PADJ_THRESHOLD)
#
# DETECTION SEMANTICS mirror 01_Build_Fitting_Dataset.R exactly: truth genes are
# defined per (method, run) from that method's GT; a gene's detection is measured
# over the iterations at the grid point within the runs where it is a truth gene;
# probabilities pool trials across those runs. Gene characteristics are averaged
# over the contributing runs.
#
# GT_DDS DEPENDENCY: dispersion comes from Power_Curve_Fitting/GT_DDS_Checkpoints
# (produced by Step 00, method-independent). Inside Stage A these exist on local
# scratch only when the aggregation block ran (NEED_AGG=true). If they are
# absent (e.g. a resume that skipped aggregation), this script logs a warning and
# writes nothing — it is BEST-EFFORT and must never fail the job.
#
# READS ONLY the cartesian tree + GT_DDS checkpoints. Writes into
#   Datalink/Power_Curve_Fitting/Cross_Score/Detection_Scatter/
# so it rides the existing cross-score stage-out tar (no new PBS tar/gate).
#
# SCOPING (env vars — for cheap interactive probes):
#   SCATTER_METHODS    default "deseq2,edger,limma"
#   SCATTER_CELLTYPES  default CELL_TYPES
#   SCATTER_RUNS       default "1..N_ROBUSTNESS_RUNS"
#   SCATTER_MAX_ITERS  default "" (all iterations per grid point)
#   SCATTER_RENDER     default "true" (set false to write data only, no PNGs)
#
# OUTPUT: Cross_Score/Detection_Scatter/
#           detection_scatter.qs / .csv   (tidy long: one row per gene x method
#                                           x cell_type x design point)
#           Plots/scatter_<design>.png    (viridis scatter, facet cell_type x method)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------------------------

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

log_section("STEP 1d: Tiny-bin Detection Scatter (dispersion x baseMean x detection)")

# ------------------------------------------------------------------------------
# 2. RESOLVE SCOPE
# ------------------------------------------------------------------------------

.env_list <- function(name, default_vec) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default_vec)
  out <- trimws(strsplit(v, "[,[:space:]]+")[[1]])
  out[nzchar(out)]
}

SC_METHODS   <- .env_list("SCATTER_METHODS",   c("deseq2", "edger", "limma"))
SC_CELLTYPES <- .env_list("SCATTER_CELLTYPES", CELL_TYPES)
SC_RUNS      <- as.integer(.env_list("SCATTER_RUNS", as.character(seq_len(N_ROBUSTNESS_RUNS))))
SC_RUNS      <- SC_RUNS[!is.na(SC_RUNS)]
.max_iters_raw <- Sys.getenv("SCATTER_MAX_ITERS")
SC_MAX_ITERS <- if (nzchar(.max_iters_raw)) as.integer(.max_iters_raw) else NA_integer_
SC_RENDER    <- tolower(Sys.getenv("SCATTER_RENDER", "true")) %in% c("true", "1", "yes")

stopifnot(length(SC_METHODS) >= 1, length(SC_CELLTYPES) >= 1, length(SC_RUNS) >= 1)

# Tiny LFC band (matches the paper's "Tiny (0.2-1.0)" bin).
TINY_LFC_LO <- GT_LFC_THRESHOLD   # 0.2
TINY_LFC_HI <- 1.0

# all_methods pooled panel: trial-weighted pool across the DE methods (sum det +
# trials per gene), added as an extra facet column alongside the per-method ones.
# On by default when >=2 methods; disable with SCATTER_ALLMETHODS=false.
SC_ALLMETHODS_LABEL <- if (exists("CROSS_METHOD_LABEL")) CROSS_METHOD_LABEL else "all_methods"
SC_ADD_ALLMETHODS   <- (tolower(Sys.getenv("SCATTER_ALLMETHODS", "true")) %in%
                          c("true", "1", "yes")) && length(SC_METHODS) >= 2L

# The three fixed experiment designs (cells mode). cells = NA_real_ means the
# MAX-cells step (all available cells per donor), stored as the CNA iteration
# files (y_value parsed to NA by parse_iteration_filename).
DESIGN_POINTS <- list(
  list(label = "600don_100cell",  donors = 600L, cells = 100),
  list(label = "100don_3000cell", donors = 100L, cells = 3000),
  list(label = "600don_maxcell",  donors = 600L, cells = NA_real_)
)

log_info("methods    : %s", paste(SC_METHODS, collapse = ", "))
log_info("cell types : %s", paste(SC_CELLTYPES, collapse = ", "))
log_info("runs       : %s", paste(SC_RUNS, collapse = ", "))
log_info("max iters  : %s", if (is.na(SC_MAX_ITERS)) "all" else as.character(SC_MAX_ITERS))
log_info("designs    : %s", paste(vapply(DESIGN_POINTS, function(d) d$label, character(1)),
                                   collapse = ", "))

# ------------------------------------------------------------------------------
# 3. ACCUMULATE PER-GENE DETECTION AT THE THREE DESIGN POINTS
# ------------------------------------------------------------------------------
# det_parts : per (method, cell_type, design, gene) integer sums n_det, n_trials
# char_parts: per (method, cell_type, gene) gene characteristics per contributing
#             run (averaged at the end)
det_parts  <- list()
char_parts <- list()

# Match iteration-file rows to a design point (cells mode).
match_design <- function(iter_files, donors, cells) {
  if (is.na(cells)) {
    iter_files[iter_files$n_donors == donors & is.na(iter_files$y_value), , drop = FALSE]
  } else {
    iter_files[iter_files$n_donors == donors &
                 !is.na(iter_files$y_value) & iter_files$y_value == cells, , drop = FALSE]
  }
}

for (ct in SC_CELLTYPES) {
  for (run in SC_RUNS) {

    # --- DESeq2 full-data GT dispersion (shared y-axis ruler) ----------------
    dds_path <- file.path(OUTPUT_BASE, "GT_DDS_Checkpoints", ct,
                          sprintf("GT_DDS_Run_%d.qs", run))
    if (!file.exists(dds_path)) {
      log_warning("GT_DDS missing for %s Run %d (%s) — skipping run (no dispersion ruler)",
                  ct, run, dds_path)
      next
    }
    gt_dds <- tryCatch(qs::qread(dds_path), error = function(e) NULL)
    if (is.null(gt_dds)) { log_warning("GT_DDS unreadable: %s", dds_path); next }
    gt_disp <- DESeq2::dispersions(gt_dds)
    names(gt_disp) <- rownames(gt_dds)

    for (m in SC_METHODS) {

      # --- method m's Tiny truth genes (this run) with characteristics -------
      gt_path <- file.path(CARTESIAN_BASE, ct, SHARD_DATASET, "cells", m,
                           paste0("Run_", run), "Ground_Truth.qs")
      if (!file.exists(gt_path)) {
        log_warning("truth missing: %s/%s Run %d — skipping", ct, m, run)
        next
      }
      gt <- tryCatch(as.data.frame(qs::qread(gt_path)), error = function(e) NULL)
      if (is.null(gt) || nrow(gt) == 0) next
      gt$gene_id <- rownames(gt)
      # Shared expression ruler: baseMean_raw (fall back to baseMean pre-rerun).
      bm_raw <- if (!is.null(gt$baseMean_raw)) gt$baseMean_raw else gt$baseMean

      keep <- !is.na(gt$padj) & gt$padj < GT_PADJ_THRESHOLD &
        !is.na(gt$log2FoldChange) &
        abs(gt$log2FoldChange) >= TINY_LFC_LO & abs(gt$log2FoldChange) < TINY_LFC_HI
      tiny <- data.table(
        gene_id     = gt$gene_id[keep],
        gt_basemean = as.numeric(bm_raw[keep]),
        abs_gt_lfc  = abs(gt$log2FoldChange[keep])
      )
      if (nrow(tiny) == 0) next
      tiny[, gt_dispersion := gt_disp[gene_id]]
      tiny <- tiny[!is.na(gt_dispersion) & gt_dispersion > 0 & gt_basemean > 0]
      if (nrow(tiny) == 0) next

      # gene characteristics for this (method, run) — averaged across runs later
      char_parts[[length(char_parts) + 1L]] <- data.table(
        method = m, cell_type = ct,
        gene_id = tiny$gene_id,
        gt_basemean = tiny$gt_basemean,
        gt_dispersion = tiny$gt_dispersion,
        abs_gt_lfc = tiny$abs_gt_lfc
      )

      truth_vec <- tiny$gene_id

      run_dir <- file.path(CARTESIAN_BASE, ct, SHARD_DATASET, "cells", m,
                           paste0("Run_", run))
      if (!dir.exists(run_dir)) next
      iter_files <- list_iteration_files(run_dir, "cells")
      if (is.null(iter_files) || nrow(iter_files) == 0) next

      for (dp in DESIGN_POINTS) {
        gp_files <- match_design(iter_files, dp$donors, dp$cells)
        if (nrow(gp_files) == 0) next
        if (!is.na(SC_MAX_ITERS) && nrow(gp_files) > SC_MAX_ITERS) {
          gp_files <- gp_files[seq_len(SC_MAX_ITERS), , drop = FALSE]
        }

        n_det <- integer(length(truth_vec))   # per truth gene, this design+run
        n_trials <- 0L
        for (r in seq_len(nrow(gp_files))) {
          res <- tryCatch(as.data.frame(qs::qread(gp_files$file_path[r])),
                          error = function(e) NULL)
          if (is.null(res) || nrow(res) == 0) next
          det_genes <- rownames(res)[!is.na(res$padj) &
                                       res$padj < DETECT_PADJ_THRESHOLD]
          n_det <- n_det + as.integer(truth_vec %in% det_genes)
          n_trials <- n_trials + 1L
        }
        if (n_trials == 0L) next

        det_parts[[length(det_parts) + 1L]] <- data.table(
          method = m, cell_type = ct, design = dp$label,
          gene_id = truth_vec, n_det = n_det, n_trials = n_trials
        )
      }

      log_info("scanned  ct=%s run=%d method=%s  (%d tiny truth genes)",
               ct, run, m, length(truth_vec))
    }
  }
}

if (length(det_parts) == 0) {
  log_warning("No detection rows produced (no inputs / GT_DDS in scope). Nothing written.")
  log_section("STEP 1d COMPLETE (empty)")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# 4. AGGREGATE ACROSS RUNS -> detection probability + averaged characteristics
# ------------------------------------------------------------------------------

det_all <- rbindlist(det_parts)
det_all <- det_all[, .(n_det = sum(n_det), n_trials = sum(n_trials)),
                   by = .(method, cell_type, design, gene_id)]

char_all <- rbindlist(char_parts)
char_all <- char_all[, .(
  gt_basemean   = mean(gt_basemean,   na.rm = TRUE),
  gt_dispersion = mean(gt_dispersion, na.rm = TRUE),
  abs_gt_lfc    = mean(abs_gt_lfc,    na.rm = TRUE)
), by = .(method, cell_type, gene_id)]

# all_methods pooled panel: sum detections + trials across methods (trial-
# weighted), average the (near-identical) shared characteristics across methods.
if (SC_ADD_ALLMETHODS) {
  det_am <- det_all[, .(n_det = sum(n_det), n_trials = sum(n_trials)),
                    by = .(cell_type, design, gene_id)]
  det_am[, method := SC_ALLMETHODS_LABEL]
  det_all <- rbind(det_all, det_am, use.names = TRUE)

  char_am <- char_all[, .(
    gt_basemean   = mean(gt_basemean,   na.rm = TRUE),
    gt_dispersion = mean(gt_dispersion, na.rm = TRUE),
    abs_gt_lfc    = mean(abs_gt_lfc,    na.rm = TRUE)
  ), by = .(cell_type, gene_id)]
  char_am[, method := SC_ALLMETHODS_LABEL]
  char_all <- rbind(char_all, char_am, use.names = TRUE)
}

det_all[, detection_prob := fifelse(n_trials > 0, n_det / n_trials, NA_real_)]

scatter <- merge(det_all, char_all, by = c("method", "cell_type", "gene_id"),
                 all.x = TRUE)
scatter <- scatter[!is.na(gt_basemean) & !is.na(gt_dispersion)]

# ------------------------------------------------------------------------------
# 5. WRITE TIDY DATA
# ------------------------------------------------------------------------------

out_dir <- file.path(OUTPUT_BASE, "Cross_Score", "Detection_Scatter")
safe_mkdir(out_dir)
scatter_df <- as.data.frame(scatter)
qs::qsave(scatter_df, file.path(out_dir, "detection_scatter.qs"))
data.table::fwrite(scatter, file.path(out_dir, "detection_scatter.csv"))
log_info("Wrote %d scatter rows -> %s", nrow(scatter_df),
         file.path(out_dir, "detection_scatter.csv"))

# ------------------------------------------------------------------------------
# 6. RENDER (best-effort; a plotting failure never fails the job) —------------
#    viridis scatter, one PNG per design, facet cell_type x method.
# ------------------------------------------------------------------------------

if (SC_RENDER) {
  render_ok <- tryCatch({
    suppressPackageStartupMessages({
      library(ggplot2)
      library(viridis)
    })
    plots_dir <- file.path(out_dir, "Plots")
    safe_mkdir(plots_dir)

    design_title <- c(
      "600don_100cell"  = "600 donors x 100 cells",
      "100don_3000cell" = "100 donors x 3000 cells",
      "600don_maxcell"  = "600 donors x max cells"
    )

    for (dp in DESIGN_POINTS) {
      d <- scatter[design == dp$label]
      if (nrow(d) == 0) {
        log_warning("no rows for design %s — skipping plot", dp$label)
        next
      }
      d <- as.data.frame(d)
      plot_method_levels <- if (SC_ADD_ALLMETHODS) c(SC_METHODS, SC_ALLMETHODS_LABEL) else SC_METHODS
      d$method    <- factor(d$method, levels = plot_method_levels)
      d$cell_type <- factor(d$cell_type, levels = SC_CELLTYPES)

      ttl <- if (!is.na(design_title[dp$label])) design_title[dp$label] else dp$label
      p <- ggplot(d, aes(x = gt_basemean, y = gt_dispersion, colour = detection_prob)) +
        geom_point(size = 1.4, alpha = 0.85) +
        facet_grid(cell_type ~ method) +
        scale_x_log10() +
        scale_y_log10() +
        scale_colour_viridis_c(option = "viridis", limits = c(0, 1),
                               name = "Detection\nprobability") +
        labs(
          title = sprintf("Tiny-bin (|LFC| 0.2-1.0) detection: %s", ttl),
          subtitle = "Each point a truth gene; colour = detection probability",
          x = "Mean expression (baseMean_raw, log10)",
          y = "Dispersion (DESeq2 GT, log10)"
        ) +
        theme_bw(base_size = 13) +
        theme(
          strip.text = element_text(face = "bold"),
          legend.key.height = unit(1.2, "cm")
        )

      fname <- file.path(plots_dir, sprintf("scatter_%s.png", dp$label))
      ggsave(fname, p, width = 12, height = 10, dpi = 150)
      log_info("Saved: %s", fname)
    }
    TRUE
  }, error = function(e) {
    log_warning("rendering failed (data already written): %s", conditionMessage(e))
    FALSE
  })
  if (!render_ok) log_warning("Detection-scatter PNGs not produced; CSV/qs are intact.")
}

log_section("STEP 1d COMPLETE")
