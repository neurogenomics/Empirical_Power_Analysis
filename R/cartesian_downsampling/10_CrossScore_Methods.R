# ==============================================================================
# SCRIPT: Cross-Method Truth Scoring (3x3) - 10
# PURPOSE: Score each DE method's downsampled calls against EVERY method's
#          full-data truth set, to separate cross-method detection
#          from self-consistency inflation.
#
#   Self-scoring (the existing pipeline) compares method i's subsamples against
#   method i's OWN full-data truth. That conflates (a) true detection ability
#   with (b) how self-reproducible a method's own gene ranking is under
#   subsampling. A method with strong variance shrinkage (e.g. limma-voom's
#   empirical-Bayes moderation) can look better on BOTH precision and recall
#   simply because its own truth is a more stable target, an artefact, not a
#   real advantage.
#
#   This script fills the 3x3 matrix:  rows = subsample method (i),
#   cols = truth method (j). The DIAGONAL (i == j) reproduces the existing
#   self-scored numbers (a built-in sanity check); the OFF-DIAGONAL asks "do
#   method i's subsamples recover method j's truth genes?". If method i only
#   shines on its own diagonal, its advantage was self-consistency; if it stays
#   high across columns, it is detecting shared signal across packages.
#
# METRICS  (per subsample_method i x truth_method j x cell_type x mode x
#           n_donors x y_value, aggregated trial-weighted over iterations+runs):
#   recall      = sum|truth_j(bin) ∩ detected_i| / sum|truth_j(bin)|
#                 (per LFC bin: Tiny 0.2-1, Small 1-2, Medium 2-3, Large 3+;
#                  truth-anchored — bin defined by truth gene's |LFC|)
#   precision   = sum|detected_i(bin) ∩ truth_j(all bins)| / sum|detected_i(bin)|
#                 (per LFC bin, call-anchored — detected_i(bin) = called genes
#                  whose OWN estimated |LFC| falls in that bin, scored against
#                  method j's full truth set across all bins. This is the exact
#                  per-bin precision definition used everywhere else in the
#                  paper — see cartesian_downsampling/03_Aggregate_Results.R
#                  (sub_in_bin ∩ gt_all_sig / |sub_in_bin|) — so the 3x3 matrix
#                  is consistent with every other precision figure.)
#
#   detected_i         : genes with padj < DETECT_PADJ_THRESHOLD in subsample i
#   truth_j            : genes with padj < GT_PADJ_THRESHOLD and
#                        |log2FoldChange| > GT_LFC_THRESHOLD in method j's GT
#   (padj / log2FoldChange are standard column names written by all three DE
#    workers, so the same code path is valid for DESeq2 / edgeR / limma.)
#
# READS ONLY from the per-run cartesian grid tree. Does NOT modify any existing
# per-run pipeline output. Intended to run where the raw runs 1-10 tree is
# present, inside Stage A on node-local scratch (the raw Run_* dirs only exist
# there), or standalone in an interactive session against the same staged tree.
#
# SCOPING (env vars - for cheap interactive probes against one slice):
#   XSCORE_METHODS    default "deseq2,edger,limma"
#   XSCORE_CELLTYPES  default "Astro,Micro,Oligo"
#   XSCORE_MODES      default "cells,reads"
#   XSCORE_RUNS       default "1..N_ROBUSTNESS_RUNS"  (comma list, e.g. "1,2")
#   XSCORE_MAX_ITERS  default "" (all iterations per grid point)
#
# OUTPUT: Datalink/Cartesian_Downsampling/Cross_Score/
#           cross_score_matrix.qs   (tidy long data.frame)
#           cross_score_matrix.csv  (same, for quick inspection)
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
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 10: Cross-Method Truth Scoring (3x3)")

# ------------------------------------------------------------------------------
# 2. RESOLVE SCOPE
# ------------------------------------------------------------------------------

.env_list <- function(name, default_vec) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default_vec)
  out <- trimws(strsplit(v, "[,[:space:]]+")[[1]])
  out[nzchar(out)]
}

CARTESIAN_BASE <- OUTPUT_BASE
GT_PADJ_THRESHOLD <- FDR_THRESHOLD
GT_LFC_THRESHOLD <- LFC_INTERVALS[["Tiny (0.2-1.0)"]][1]
DETECT_PADJ_THRESHOLD <- FDR_THRESHOLD
XS_DEFAULT_CELLTYPES <- c("Astro", "Micro", "Oligo")
XS_DEFAULT_MODES <- c("cells", "reads")

XS_METHODS   <- .env_list("XSCORE_METHODS",   DE_METHODS)
XS_CELLTYPES <- .env_list("XSCORE_CELLTYPES", XS_DEFAULT_CELLTYPES)
XS_MODES     <- .env_list("XSCORE_MODES",     XS_DEFAULT_MODES)
XS_RUNS      <- as.integer(.env_list("XSCORE_RUNS", as.character(seq_len(N_ROBUSTNESS_RUNS))))
XS_RUNS      <- XS_RUNS[!is.na(XS_RUNS)]
.max_iters_raw <- Sys.getenv("XSCORE_MAX_ITERS")
XS_MAX_ITERS <- if (nzchar(.max_iters_raw)) as.integer(.max_iters_raw) else NA_integer_

stopifnot(length(XS_METHODS) >= 1, length(XS_CELLTYPES) >= 1,
          length(XS_MODES) >= 1, length(XS_RUNS) >= 1)

log_info("methods    : %s", paste(XS_METHODS, collapse = ", "))
log_info("cell types : %s", paste(XS_CELLTYPES, collapse = ", "))
log_info("modes      : %s", paste(XS_MODES, collapse = ", "))
log_info("runs       : %s", paste(XS_RUNS, collapse = ", "))
log_info("max iters  : %s", if (is.na(XS_MAX_ITERS)) "all" else as.character(XS_MAX_ITERS))

# LFC bin edges (match the paper truth-set bands; truth set is |LFC| > 0.2).
LFC_BREAKS <- c(GT_LFC_THRESHOLD, 1.0, 2.0, 3.0, Inf)
LFC_LABELS <- c("Tiny", "Small", "Medium", "Large")
lfc_bin_of <- function(abs_lfc) {
  as.character(cut(abs_lfc, breaks = LFC_BREAKS, labels = LFC_LABELS,
                   right = FALSE, include.lowest = TRUE))
}

# BaseMean expression bins for the filtering-diagnostic metrics below. These
# mirror BASEMEAN_INTERVALS in cartesian_downsampling/config.R. Keep in sync if
# the cartesian bins ever change.
BASEMEAN_BREAKS <- c(0, 10, 50, 200, Inf)
BASEMEAN_LABELS <- c("Very Low (<10)", "Low (10-50)", "Medium (50-200)", "High (200+)")
# Tiny LFC band for Metric 2 (matches LFC_INTERVALS "Tiny (0.2-1.0)").
TINY_LFC_LO <- GT_LFC_THRESHOLD   # 0.2
TINY_LFC_HI <- 1.0

# ------------------------------------------------------------------------------
# 3. CROSS-SCORE
# ------------------------------------------------------------------------------
# Accumulate INTEGER count sums (never rates) so runs can be summed and the
# rate computed once at the end. Two accumulators:
#   recall_acc    keyed by (i, j, cell_type, mode, n_donors, y_value, lfc_bin)
#                 sums: n_truth, n_recovered
#   prec_acc      keyed by (i, j, cell_type, mode, n_donors, y_value)
#                 sums: n_tp, n_detected
# Per (cell_type, run, mode, i) chunk we collapse iterations immediately, so the
# master tables stay small (one row per grid point x bin, not per iteration).
recall_parts <- list()
prec_parts   <- list()
# Filtering-diagnostic accumulators (properties of subsample method i alone, not
# cross i x j). filt = mean testable-gene count; tiny = mean Tiny-band DEG count
# split by expression bin. Both keyed on subsample method + grid point.
filt_parts   <- list()
tiny_parts   <- list()

for (ct in XS_CELLTYPES) {
  for (run in XS_RUNS) {

    # --- load all methods' truth sets for this (cell_type, run) once ---------
    truth <- list()
    gt_median <- list()
    for (m in XS_METHODS) {
      gt_path <- file.path(CARTESIAN_BASE, ct, SHARD_DATASET, "cells", m,
                           paste0("Run_", run), "Ground_Truth.qs")
      if (!file.exists(gt_path)) {
        log_warning("truth missing: %s/%s Run %d (%s) — skipping this truth",
                    ct, m, run, gt_path)
        next
      }
      gt <- tryCatch(as.data.frame(qs::qread(gt_path)), error = function(e) NULL)
      if (is.null(gt) || nrow(gt) == 0) next
      gt$gene_id <- rownames(gt)
      # GT-half median raw count per gene, used to bin method m's subsample calls.
      gt_median[[m]] <- setNames(if (!is.null(gt$median_raw)) gt$median_raw else gt$baseMean,
                                 gt$gene_id)
      tj <- gt[!is.na(gt$padj) &
                 gt$padj < GT_PADJ_THRESHOLD &
                 !is.na(gt$log2FoldChange) &
                 abs(gt$log2FoldChange) > GT_LFC_THRESHOLD, , drop = FALSE]
      if (nrow(tj) == 0) next
      truth[[m]] <- data.table(
        gene_id = tj$gene_id,
        lfc_bin = lfc_bin_of(abs(tj$log2FoldChange))
      )
    }
    if (length(truth) == 0) {
      log_warning("no truth sets for %s Run %d — skipping run", ct, run)
      next
    }

    # per-(j) truth-set-size by bin, and full truth gene vector (for precision)
    truth_bin_n  <- lapply(truth, function(t) t[, .N, by = lfc_bin])
    truth_gene_v <- lapply(truth, function(t) t$gene_id)
    truth_bin_v  <- lapply(truth, function(t) split(t$gene_id, t$lfc_bin))

    for (mode in XS_MODES) {
      for (i_method in XS_METHODS) {

        run_dir <- file.path(CARTESIAN_BASE, ct, SHARD_DATASET, mode, i_method,
                             paste0("Run_", run))
        if (!dir.exists(run_dir)) next

        iter_files <- list_iteration_files(run_dir, mode)
        if (is.null(iter_files) || nrow(iter_files) == 0) next
        iter_files <- iter_files[!is.na(iter_files$y_value), , drop = FALSE]
        if (nrow(iter_files) == 0) next

        # per-chunk iteration-level accumulators (collapsed before append)
        rec_rows  <- vector("list", nrow(iter_files) * length(truth) * length(LFC_LABELS))
        prc_rows  <- vector("list", nrow(iter_files) * length(truth) * length(LFC_LABELS))
        flt_rows  <- vector("list", nrow(iter_files))
        tny_rows  <- vector("list", nrow(iter_files) * length(BASEMEAN_LABELS))
        rk <- 0L; pk <- 0L; fk <- 0L; tk <- 0L

        # optional per-grid iteration cap (for cheap probes)
        if (!is.na(XS_MAX_ITERS)) {
          iter_files <- as.data.frame(
            data.table(iter_files)[, head(.SD, XS_MAX_ITERS),
                                   by = .(n_donors, y_value)]
          )
        }

        for (r in seq_len(nrow(iter_files))) {
          res <- tryCatch(as.data.frame(qs::qread(iter_files$file_path[r])),
                          error = function(e) NULL)
          if (is.null(res) || nrow(res) == 0) next
          gene_id <- rownames(res)
          # Expression ruler: each called gene's median raw count in method
          # i_method's ground-truth half (0 if the gene failed the GT min-count
          # filter), so a gene's expression bin does not change with the grid point.
          res_bm_raw <- unname(gt_median[[i_method]][gene_id])
          res_bm_raw[is.na(res_bm_raw)] <- 0
          has_padj <- !is.na(res$padj)
          universe_i <- gene_id[has_padj]
          detected_i <- gene_id[has_padj & res$padj < DETECT_PADJ_THRESHOLD]
          n_det <- length(detected_i)

          # Split the called genes by their OWN estimated |LFC| bin, for the
          # per-bin (call-anchored) precision below. Same lfc_bin_of() breaks as
          # the truth bins, so precision bins align with recall bins. Genes with
          # |LFC| < 0.2 or NA fall outside all bins and are dropped (identical to
          # 03_Aggregate_Results.R, where sub_in_bin requires |LFC| >= bound).
          det_by_bin <- if (n_det > 0) {
            split(detected_i,
                  factor(lfc_bin_of(abs(res$log2FoldChange[has_padj &
                           res$padj < DETECT_PADJ_THRESHOLD])),
                         levels = LFC_LABELS))
          } else {
            setNames(vector("list", length(LFC_LABELS)), LFC_LABELS)
          }

          # --- METRIC 1: testable-gene count for subsample i (independent-
          #     filtering signal). Unified "testable = non-NA padj": for DESeq2
          #     this drops independent-filtered / Cook's-outlier / all-zero genes
          #     (present as NA-padj rows); for edgeR/limma the min-count-filtered
          #     genes are absent rows, so this equals the present-and-tested set.
          #     Same rule catches "untestable" correctly for all three methods.
          fk <- fk + 1L
          flt_rows[[fk]] <- list(
            i = i_method,
            n_donors = iter_files$n_donors[r], y_value = iter_files$y_value[r],
            n_testable = length(universe_i), n_obs = 1L
          )

          # --- METRIC 2: Tiny-band DEG count (padj < 0.05 AND |LFC| in Tiny),
          #     split by GT-half median expression bin (shared cross-method ruler).
          #     Counts ALL significant Tiny calls regardless of truth membership
          #     (true+false positives), using the subsample's own LFC.
          #     Zeros recorded for
          #     every bin so the per-iteration mean stays honest.
          tiny_counts <- setNames(integer(length(BASEMEAN_LABELS)), BASEMEAN_LABELS)
          if (n_det > 0) {
            det_idx  <- has_padj & res$padj < DETECT_PADJ_THRESHOLD
            det_lfc  <- abs(res$log2FoldChange[det_idx])
            det_bm   <- res_bm_raw[det_idx]
            tiny_sel <- !is.na(det_lfc) & det_lfc >= TINY_LFC_LO & det_lfc < TINY_LFC_HI
            if (any(tiny_sel)) {
              bm_bin <- cut(det_bm[tiny_sel], breaks = BASEMEAN_BREAKS,
                            labels = BASEMEAN_LABELS, right = FALSE,
                            include.lowest = TRUE)
              tc <- table(bm_bin)
              tiny_counts[names(tc)] <- as.integer(tc)
            }
          }
          for (bn in BASEMEAN_LABELS) {
            tk <- tk + 1L
            tny_rows[[tk]] <- list(
              i = i_method,
              n_donors = iter_files$n_donors[r], y_value = iter_files$y_value[r],
              basemean_bin = bn, n_deg = tiny_counts[[bn]], n_obs = 1L
            )
          }

          for (j_method in names(truth)) {
            tg_all <- truth_gene_v[[j_method]]
            # precision per LFC bin (call-anchored): within each bin of called
            # genes (binned on their own |LFC|), count how many are in method j's
            # full truth set. Identical to the standard per-bin precision in
            # 03_Aggregate_Results.R (sub_in_bin ∩ gt_all_sig / |sub_in_bin|).
            for (b in LFC_LABELS) {
              det_b   <- det_by_bin[[b]]
              n_det_b <- length(det_b)
              n_tp_b  <- if (n_det_b > 0) length(intersect(det_b, tg_all)) else 0L
              pk <- pk + 1L
              prc_rows[[pk]] <- list(
                i = i_method, j = j_method,
                n_donors = iter_files$n_donors[r], y_value = iter_files$y_value[r],
                lfc_bin = b,
                n_tp = n_tp_b, n_detected = n_det_b
              )
            }

            # recall per LFC bin (truth-anchored)
            tbin <- truth_bin_v[[j_method]]
            for (b in names(tbin)) {
              tb <- tbin[[b]]
              n_truth_raw  <- length(tb)
              n_recovered  <- length(intersect(tb, detected_i))
              rk <- rk + 1L
              rec_rows[[rk]] <- list(
                i = i_method, j = j_method,
                n_donors = iter_files$n_donors[r], y_value = iter_files$y_value[r],
                lfc_bin = b,
                n_truth_raw = n_truth_raw,
                n_recovered = n_recovered
              )
            }
          }
        }

        if (rk == 0L) next

        # collapse this chunk's iterations -> per grid-point x bin sums
        rec_dt <- rbindlist(rec_rows[seq_len(rk)])
        rec_dt <- rec_dt[, .(
          n_truth_raw  = sum(n_truth_raw),
          n_recovered  = sum(n_recovered)
        ), by = .(i, j, n_donors, y_value, lfc_bin)]
        rec_dt[, `:=`(cell_type = ct, mode = mode)]
        recall_parts[[length(recall_parts) + 1L]] <- rec_dt

        prc_dt <- rbindlist(prc_rows[seq_len(pk)])
        prc_dt <- prc_dt[, .(
          n_tp       = sum(n_tp),
          n_detected = sum(n_detected)
        ), by = .(i, j, n_donors, y_value, lfc_bin)]
        prc_dt[, `:=`(cell_type = ct, mode = mode)]
        prec_parts[[length(prec_parts) + 1L]] <- prc_dt

        # collapse filtering-diagnostic chunks -> per grid-point sums
        if (fk > 0L) {
          flt_dt <- rbindlist(flt_rows[seq_len(fk)])
          flt_dt <- flt_dt[, .(n_testable = sum(n_testable), n_obs = sum(n_obs)),
                           by = .(i, n_donors, y_value)]
          flt_dt[, `:=`(cell_type = ct, mode = mode)]
          filt_parts[[length(filt_parts) + 1L]] <- flt_dt
        }
        if (tk > 0L) {
          tny_dt <- rbindlist(tny_rows[seq_len(tk)])
          tny_dt <- tny_dt[, .(n_deg = sum(n_deg), n_obs = sum(n_obs)),
                           by = .(i, n_donors, y_value, basemean_bin)]
          tny_dt[, `:=`(cell_type = ct, mode = mode)]
          tiny_parts[[length(tiny_parts) + 1L]] <- tny_dt
        }

        log_info("scored  ct=%s run=%d mode=%s i=%s  (%d grid-point rows)",
                 ct, run, mode, i_method, nrow(rec_dt))
      }
    }
  }
}

if (length(recall_parts) == 0) {
  log_warning("No cross-score rows produced (no inputs found in scope). Nothing written.")
  log_section("STEP 10 COMPLETE (empty)")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# 4. AGGREGATE ACROSS RUNS + COMPUTE RATES
# ------------------------------------------------------------------------------

recall_all <- rbindlist(recall_parts)
recall_all <- recall_all[, .(
  n_truth_raw  = sum(n_truth_raw),
  n_recovered  = sum(n_recovered)
), by = .(i, j, cell_type, mode, n_donors, y_value, lfc_bin)]
recall_all[, `:=`(
  recall = fifelse(n_truth_raw > 0, n_recovered / n_truth_raw, NA_real_),
  metric = "recall"
)]

prec_all <- rbindlist(prec_parts)
prec_all <- prec_all[, .(
  n_tp       = sum(n_tp),
  n_detected = sum(n_detected)
), by = .(i, j, cell_type, mode, n_donors, y_value, lfc_bin)]
prec_all[, `:=`(
  precision = fifelse(n_detected > 0, n_tp / n_detected, NA_real_),
  metric    = "precision"
)]

# tidy long form: one column set covering both metrics
recall_long <- recall_all[, .(
  subsample_method = i, truth_method = j, cell_type, mode,
  n_donors, y_value, lfc_bin, metric,
  value     = recall,
  numerator = n_recovered,
  denom_raw = n_truth_raw
)]
prec_long <- prec_all[, .(
  subsample_method = i, truth_method = j, cell_type, mode,
  n_donors, y_value, lfc_bin, metric,
  value     = precision,
  numerator = n_tp,
  denom_raw = n_detected
)]

cross_score <- rbindlist(list(recall_long, prec_long), use.names = TRUE)
cross_score[, self := (subsample_method == truth_method)]

# ------------------------------------------------------------------------------
# 5. WRITE
# ------------------------------------------------------------------------------

out_dir <- file.path(OUTPUT_BASE, "Cross_Score")
safe_mkdir(out_dir)
out_qs  <- file.path(out_dir, "cross_score_matrix.qs")
out_csv <- file.path(out_dir, "cross_score_matrix.csv")

cross_df <- as.data.frame(cross_score)
qs::qsave(cross_df, out_qs)
data.table::fwrite(cross_score, out_csv)

log_info("Wrote %d cross-score rows", nrow(cross_df))
log_info("  qs : %s", out_qs)
log_info("  csv: %s", out_csv)

# quick collapsed 3x3 preview (Tiny-bin recall, cells) for the log
preview <- cross_score[metric == "recall" & lfc_bin == "Tiny" & mode == "cells",
  .(recall = sum(numerator) / sum(denom_raw)),
  by = .(cell_type, subsample_method, truth_method)]
if (nrow(preview) > 0) {
  log_section("Preview: Tiny-bin recall (cells) - rows=subsample, cols=truth")
  for (ctp in unique(preview$cell_type)) {
    log_info("[%s]", ctp)
    sub <- preview[cell_type == ctp]
    wide <- dcast(sub, subsample_method ~ truth_method, value.var = "recall")
    tmp <- capture.output(print(wide))
    for (ln in tmp) log_info("  %s", ln)
  }
}

# ------------------------------------------------------------------------------
# 6. FILTERING-DIAGNOSTIC METRICS (per subsample method; independent-filtering
#    exploration). Written alongside the cross-score matrix so they ride the
#    same Cross_Score stage-out tar. Means are over iterations + robustness runs.
# ------------------------------------------------------------------------------

if (length(filt_parts) > 0) {
  filt_all <- rbindlist(filt_parts)
  filt_all <- filt_all[, .(n_testable = sum(n_testable), n_obs = sum(n_obs)),
                       by = .(i, cell_type, mode, n_donors, y_value)]
  filt_all[, mean_testable := fifelse(n_obs > 0, n_testable / n_obs, NA_real_)]
  filt_out <- filt_all[, .(
    method = i, cell_type, mode, n_donors, y_value,
    mean_testable, n_obs
  )]
  filt_csv <- file.path(out_dir, "filtering_testable_genes.csv")
  qs::qsave(as.data.frame(filt_out), file.path(out_dir, "filtering_testable_genes.qs"))
  data.table::fwrite(filt_out, filt_csv)
  log_info("Wrote %d testable-gene rows -> %s", nrow(filt_out), filt_csv)
} else {
  log_warning("No testable-gene rows produced.")
}

if (length(tiny_parts) > 0) {
  tiny_all <- rbindlist(tiny_parts)
  tiny_all <- tiny_all[, .(n_deg = sum(n_deg), n_obs = sum(n_obs)),
                       by = .(i, cell_type, mode, n_donors, y_value, basemean_bin)]
  tiny_all[, mean_deg := fifelse(n_obs > 0, n_deg / n_obs, NA_real_)]
  tiny_out <- tiny_all[, .(
    method = i, cell_type, mode, n_donors, y_value, basemean_bin,
    mean_deg, n_obs
  )]
  tiny_csv <- file.path(out_dir, "filtering_tiny_deg_by_expr.csv")
  qs::qsave(as.data.frame(tiny_out), file.path(out_dir, "filtering_tiny_deg_by_expr.qs"))
  data.table::fwrite(tiny_out, tiny_csv)
  log_info("Wrote %d Tiny-DEG-by-expression rows -> %s", nrow(tiny_out), tiny_csv)
} else {
  log_warning("No Tiny-DEG-by-expression rows produced.")
}

log_section("STEP 10 COMPLETE")
