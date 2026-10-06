# ==============================================================================
# SCRIPT: rank_oligo_power.R
# PURPOSE: Rank EVERY oligodendrocyte gene by predicted detectability, for the
#          threshold-free GO enrichment (prerank GSEA) in go_prerank_oligo.py.
#
#   This is the R half of the two-step GO pipeline (Chris, "rank in R, prerank in
#   Python"). Ranking is done in R so it reuses the SAME renv-pinned xgboost and
#   the SAME shipped `predict_power.R` bundle wrapper that every other power
#   number in the paper uses. Python then only consumes the ranked lists, so no
#   model ever loads twice / in two runtimes.
#
#   REPLACES the ranking half of the old over-representation analysis (ORA),
#   which took the top-200 easiest + bottom-200 hardest genes as discrete hit
#   lists. Nathan flagged that cutoff as arbitrary. Here we emit the FULL ranked
#   list; go_prerank_oligo.py runs pre-ranked GSEA over all of it (no cutoff,
#   no background).
#
# RANKING METRIC (Chris, confirmed):
#   Per gene, predict power at THREE effect sizes |LFC| in {0.5, 1.5, 2.5} at the
#   fixed PsychAD/Oligo design (100 donors, 250 cells/donor, 5742 reads/cell,
#   de_method = deseq2), then rank by the MEAN of the three. Predicted power
#   saturates near 1.0 for the easy majority; the |LFC| = 0.5 setting is far from
#   saturation and spreads those ties so the ranking has real resolution.
#   Non-scorable genes (no finite basemean/dispersion) -> score 0 (bottom).
#
# TWO GENE SETS (both requested):
#   allgene : the full oligo transcriptome (every gene with a symbol).
#   pc      : protein-coding only  (locus_group == "protein-coding gene").
#
# INPUT:
#   Datalink/Gene_Classification/oligo_gene_classification_annotated.csv
#       -> gene_id, symbol, locus_group, basemean, dispersion. basemean and
#          dispersion are the exact checkpoint values used by classify_oligo_*
#          (verified identical to oligo_gene_features_checkpoint.qs), so no refit.
#   Datalink/Power_Curve_Fitting/Exports/xgboost/{predict_power.R,
#          power_model_bundle.rds}  -> the promoted power model.
#
# OUTPUT (Datalink/Gene_Classification/GO_Enrichment/):
#   oligo_<set>_prerank_input.rnk   two-column (symbol<TAB>score), no header,
#                                   descending, duplicate symbols collapsed by
#                                   MAX. This is the GSEA .rnk fed to Python.
#   oligo_<set>_ranked_genes.csv    same ranking with the per-LFC powers kept,
#                                   for inspection/provenance.
#       where <set> in {allgene, pc}
#
# RUN (on iclneuro or inside the Docker image, POWER_ANALYSIS_ROOT set):
#   Rscript R/gene_classification/rank_oligo_power.R
#   then:  python R/gene_classification/go_prerank_oligo.py
# ==============================================================================

suppressPackageStartupMessages({ library(xgboost) })

if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}
ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT")
if (!nzchar(ROOT)) stop("POWER_ANALYSIS_ROOT is not set and here::here() unavailable")

# ------------------------------------------------------------------ PARAMETERS
CELL_TYPE    <- "Oligo"
DE_METHOD    <- "deseq2"
FIXED_DONORS <- 100
FIXED_CELLS  <- 250
FIXED_READS  <- 5742
LFC_SETTINGS <- c(0.5, 1.5, 2.5)          # rank score = mean power over these
PC_LOCUS_GROUP <- "protein-coding gene"
SELFCHECK_GENES <- c("GSK3B", "APC", "ITCH")

gc_dir     <- file.path(ROOT, "Datalink", "Gene_Classification")
out_dir    <- file.path(gc_dir, "GO_Enrichment")
anno_csv   <- file.path(gc_dir, "oligo_gene_classification_annotated.csv")
bundle_dir <- file.path(ROOT, "Datalink", "Power_Curve_Fitting",
                        "Exports", "xgboost")
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------- LOAD
message("[rank] loading annotated genes: ", anno_csv)
df <- read.csv(anno_csv, stringsAsFactors = FALSE, check.names = FALSE)
for (req in c("gene_id", "symbol", "locus_group", "basemean", "dispersion")) {
  if (!req %in% names(df)) stop(sprintf("missing required column '%s' in %s", req, anno_csv))
}
df$basemean   <- suppressWarnings(as.numeric(df$basemean))
df$dispersion <- suppressWarnings(as.numeric(df$dispersion))
scorable <- is.finite(df$basemean) & df$basemean > 0 &
            is.finite(df$dispersion) & df$dispersion > 0
message(sprintf("[rank] %d genes; protein-coding = %d; scorable = %d",
                nrow(df), sum(df$locus_group == PC_LOCUS_GROUP, na.rm = TRUE),
                sum(scorable)))

source(file.path(bundle_dir, "predict_power.R"))

# --------------------------------------------------------------------- PREDICT
# Predict power at each |LFC| for scorable genes; non-scorable -> 0 (bottom).
lfc_tags <- sprintf("power_lfc%03d", as.integer(round(LFC_SETTINGS * 100)))
for (j in seq_along(LFC_SETTINGS)) {
  lfc <- LFC_SETTINGS[[j]]
  p <- rep(0, nrow(df))
  p[scorable] <- predict_power(
    n_donors       = FIXED_DONORS,
    abs_lfc        = lfc,
    basemean       = df$basemean[scorable],
    dispersion     = df$dispersion[scorable],
    n_cells        = FIXED_CELLS,
    reads_per_cell = FIXED_READS,
    cell_type      = CELL_TYPE,
    de_method      = DE_METHOD,
    bundle_path    = bundle_path,
    warn_extrapolation = FALSE
  )
  df[[lfc_tags[j]]] <- p
  message(sprintf("[rank] |LFC|=%.1f: mean power (scorable) = %.4f",
                  lfc, mean(p[scorable])))
}
df$score <- rowMeans(df[, lfc_tags])       # rank metric = mean over the 3 LFC

# Self-check: print per-gene power so it can be diffed against classify_oligo_*.
sc <- df[df$symbol %in% SELFCHECK_GENES, , drop = FALSE]
if (nrow(sc) > 0) {
  message("[rank] self-check (compare to classify_oligo_multi_lfc.R):")
  for (i in seq_len(nrow(sc))) {
    r <- sc[i, ]
    message(sprintf("  %-8s bm=%.3f disp=%.4f p05=%.4f p15=%.4f p25=%.4f mean=%.4f",
                    r$symbol, r$basemean, r$dispersion,
                    r[[lfc_tags[1]]], r[[lfc_tags[2]]], r[[lfc_tags[3]]], r$score))
  }
}

# --------------------------------------------------------------------- EMIT
# For each set, collapse duplicate symbols to the MAX score (a symbol is as
# detectable as its easiest gene id), drop blank symbols, sort descending, and
# write the .rnk (Python's prerank input) plus a provenance CSV.
emit_set <- function(sub, tag) {
  sub <- sub[!is.na(sub$symbol) & nzchar(sub$symbol), , drop = FALSE]
  ord <- order(sub$score, decreasing = TRUE)
  sub <- sub[ord, , drop = FALSE]
  sub <- sub[!duplicated(sub$symbol), , drop = FALSE]     # keep max (first after sort)

  rnk <- data.frame(symbol = sub$symbol, score = sub$score)
  rnk_path <- file.path(out_dir, sprintf("oligo_%s_prerank_input.rnk", tag))
  write.table(rnk, rnk_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  keep <- c("gene_id", "symbol", "locus_group", "basemean", "dispersion",
            lfc_tags, "score")
  csv_path <- file.path(out_dir, sprintf("oligo_%s_ranked_genes.csv", tag))
  write.csv(sub[, keep], csv_path, row.names = FALSE)

  message(sprintf("[rank] [%s] %d unique symbols -> %s", tag, nrow(rnk),
                  basename(rnk_path)))
}

emit_set(df, "allgene")
emit_set(df[df$locus_group == PC_LOCUS_GROUP, , drop = FALSE], "pc")
message("[rank] done.")
