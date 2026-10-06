# ==============================================================================
# SCRIPT: group_analysis_oligo.R
# PURPOSE: Compare functional gene GROUPS (not individual genes) on a nominal,
#          junk-free footing. For each group, take the MEDIAN (+ IQR) of the two
#          model characteristics (raw basemean, DESeq2 dispersion) over the
#          expressed genes in that group, then push that single REPRESENTATIVE
#          entry through the XGBoost power model at the FIXED design.
#          This answers "are TFs actually lower-expressed than other groups?"
#          without pseudogene/lncRNA junk diluting any comparison, because it
#          compares group central tendencies rather than counting genes.
#
# WHY THIS DESIGN (Chris, 2026-08-23):
#   The headline "TFs are highly detectable" was an artefact of comparing a clean
#   curated TF list against a background half-full of undetectable junk. The only
#   fair comparison is group-vs-group on matched summary stats. So we:
#     1. keep only expressed genes (finite raw basemean>0 AND DESeq2 dispersion>0),
#     2. group by PANTHER protein class / Lambert DBD family / TF-vs-nonTF,
#     3. take median + IQR of basemean & dispersion per group,
#     4. predict a single representative power per group (+ a power band from the
#        basemean IQR, holding dispersion at its median),
#     5. classify the representative with the same skewed trinary (0.4 / 0.8).
#
# NOTE (marginal medians): median(basemean) and median(dispersion) are taken
# independently, so the representative row decouples the per-gene basemean-
# dispersion pairing. That is the intended "nominal representative" per Chris.
#
# FIXED design: 100 donors, 250 cells/donor, 5742 reads/cell, |LFC|=1.5, Oligo.
# Trinary: power<=0.4 Undetectable | 0.4-0.8 Unsure | >0.8 Detectable.
#
# INPUT : Datalink/Gene_Classification/oligo_gene_classification_annotated.csv
# OUTPUT: Datalink/Gene_Classification/oligo_group_comparison_{panther,dbd,tf}.csv
# ==============================================================================

suppressPackageStartupMessages({ library(xgboost) })

ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT")
if (!nzchar(ROOT) && requireNamespace("here", quietly = TRUE)) ROOT <- here::here()

gc_dir     <- file.path(ROOT, "Datalink", "Gene_Classification")
bundle_dir <- file.path(ROOT, "Datalink", "Power_Curve_Fitting",
                        "Exports", "xgboost")
source(file.path(bundle_dir, "predict_power.R"))
bundle_path <- file.path(bundle_dir, "power_model_bundle.rds")

FIXED_LFC <- 1.5; FIXED_DONORS <- 100; FIXED_CELLS <- 250
FIXED_READS <- 5742; CELL_TYPE <- "Oligo"; DE_METHOD <- "deseq2"
MIN_N <- 15   # minimum expressed genes for a group to be reported

classify_trinary <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p <= 0.4, "Undetectable",
                ifelse(p <= 0.8, "Unsure", "Detectable")))
}

df <- read.csv(file.path(gc_dir, "oligo_gene_classification_annotated.csv"),
               stringsAsFactors = FALSE)
tf_flag <- df$is_tf_lambert %in% c(TRUE, "True", "TRUE", "true")

# Useful / expressed universe: finite basemean>0 AND finite dispersion>0.
keep <- is.finite(df$basemean) & df$basemean > 0 &
        is.finite(df$dispersion) & df$dispersion > 0
df <- df[keep, ]; tf_flag <- tf_flag[keep]
cat(sprintf("Expressed genes for group analysis: %d\n", nrow(df)))

summarise_groups <- function(data, group_col) {
  gv <- data[[group_col]]
  parts <- split(seq_len(nrow(data)), gv)
  rows <- lapply(names(parts), function(g) {
    idx <- parts[[g]]
    if (length(idx) < MIN_N || is.na(g) || g == "") return(NULL)
    bm <- data$basemean[idx]; dp <- data$dispersion[idx]
    data.frame(
      group          = g,
      n              = length(idx),
      basemean_med   = median(bm),
      basemean_q1    = as.numeric(quantile(bm, 0.25)),
      basemean_q3    = as.numeric(quantile(bm, 0.75)),
      disp_med       = median(dp),
      disp_q1        = as.numeric(quantile(dp, 0.25)),
      disp_q3        = as.numeric(quantile(dp, 0.75)),
      pct_detectable = 100 * mean(data$trinary_class_deseq2[idx] == "Detectable"),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

predict_rep <- function(tbl) {
  if (is.null(tbl)) return(NULL)
  pw <- function(bm, dp) predict_power(
    n_donors = FIXED_DONORS, abs_lfc = FIXED_LFC, basemean = bm, dispersion = dp,
    n_cells = FIXED_CELLS, reads_per_cell = FIXED_READS, cell_type = CELL_TYPE,
    de_method = DE_METHOD, bundle_path = bundle_path, warn_extrapolation = FALSE)
  tbl$power_rep <- pw(tbl$basemean_med, tbl$disp_med)         # median|median
  tbl$power_lo  <- pw(tbl$basemean_q1,  tbl$disp_med)         # Q1 basemean band
  tbl$power_hi  <- pw(tbl$basemean_q3,  tbl$disp_med)         # Q3 basemean band
  tbl$class_rep <- classify_trinary(tbl$power_rep)
  tbl[order(-tbl$basemean_med), ]
}

# --- 1. PANTHER protein class (main comparison) ------------------------------
panther <- predict_rep(summarise_groups(df, "panther_protein_class"))
write.csv(panther, file.path(gc_dir, "oligo_group_comparison_panther.csv"),
          row.names = FALSE)

# --- 2. Lambert DBD family (TF subtypes) -------------------------------------
dbd <- predict_rep(summarise_groups(df[tf_flag, ], "lambert_dbd"))
write.csv(dbd, file.path(gc_dir, "oligo_group_comparison_dbd.csv"),
          row.names = FALSE)

# --- 3. Coarse TF vs protein-coding non-TF vs other --------------------------
df$tf_group <- ifelse(tf_flag, "Lambert TF",
                ifelse(df$locus_group == "protein-coding gene",
                       "Protein-coding non-TF", "Other / non-coding"))
tf <- predict_rep(summarise_groups(df, "tf_group"))
write.csv(tf, file.path(gc_dir, "oligo_group_comparison_tf.csv"),
          row.names = FALSE)

cat("\n=== PANTHER protein class (top/bottom by median basemean) ===\n")
print(utils::head(panther[, c("group","n","basemean_med","disp_med","power_rep","class_rep","pct_detectable")], 12))
print(utils::tail(panther[, c("group","n","basemean_med","disp_med","power_rep","class_rep","pct_detectable")], 8))
cat("\n=== Lambert DBD family ===\n")
print(dbd[, c("group","n","basemean_med","disp_med","power_rep","power_lo","power_hi","class_rep")])
cat("\n=== TF vs non-TF ===\n")
print(tf[, c("group","n","basemean_med","disp_med","power_rep","class_rep","pct_detectable")])
cat("\nDone.\n")
