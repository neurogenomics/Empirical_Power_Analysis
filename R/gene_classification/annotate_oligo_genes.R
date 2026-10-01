#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: annotate_oligo_genes.R
# PURPOSE: Annotate every oligodendrocyte gene with its symbol (HGNC), PANTHER
#          protein class/family and Lambert et al. (2018) transcription-factor
#          status, for the gene-class panels of Figure 6.
#
# INPUT : Datalink/Gene_Classification/oligo_gene_power_classification.csv
#           (classify_oligo_genes.R)
#         reference/ (offline bundle, see pipeline/reference/README.md):
#           hgnc_complete_set_2026-09-04.txt       HGNC complete set, monthly archive
#           PTHR19.0_human                         PANTHER 19.0 sequence classification
#           DatabaseExtract_v_1.01.csv             Lambert 2018 human TFs v1.01
# OUTPUT: Datalink/Gene_Classification/oligo_gene_classification_annotated.csv
#         Datalink/Gene_Classification/oligo_lambert_tfs.csv   (TF rows only)
#
# JOINS (all left joins onto the classification table, one row per gene):
#   HGNC    ensembl_gene_id == gene_id     -> symbol, hgnc_id, locus_group
#           (genes without an HGNC record get locus_group "unannotated (no HGNC)")
#   PANTHER HGNC=<id> in column 1 == hgnc_id (falling back to the symbol in
#           column 3) -> panther_family (family name, column 5) and
#           panther_protein_class (first class in column 10, "#PCxxxxx" removed)
#   Lambert "Ensembl ID" == gene_id       -> lambert_dbd ("DBD", every Lambert
#           row) and is_tf_lambert ("Is TF?" == "Yes")
#
# Reconstructed 2026-10-01 from the record of the lost /tmp/oligo_annot/annotate.py
# (diary 2026-08-25..09-03) and checked against its output.
#
# Usage: Rscript R/gene_classification/annotate_oligo_genes.R
#        env EPA_REFERENCE_DIR (default <root>/reference)
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("POWER_ANALYSIS_ROOT", unset = "")
if (!nzchar(root)) root <- here::here()
ref_dir <- Sys.getenv("EPA_REFERENCE_DIR", unset = file.path(root, "reference"))
gc_dir  <- file.path(root, "Datalink", "Gene_Classification")

need <- function(p) {
  if (!file.exists(p)) stop("missing input: ", p, call. = FALSE)
  p
}

cls <- fread(need(file.path(gc_dir, "oligo_gene_power_classification.csv")))

# ---- HGNC ---------------------------------------------------------------------
hgnc <- fread(need(file.path(ref_dir, Sys.getenv("EPA_HGNC_FILE", "hgnc_complete_set_2026-09-04.txt"))),
              select = c("hgnc_id", "symbol", "locus_group", "ensembl_gene_id"),
              quote = "", na.strings = "")
hgnc <- hgnc[!is.na(ensembl_gene_id) & nzchar(ensembl_gene_id)]
hgnc <- hgnc[!duplicated(ensembl_gene_id)]

# ---- PANTHER 19.0 -------------------------------------------------------------
pth <- fread(need(file.path(ref_dir, "PTHR19.0_human")), sep = "\t", header = FALSE,
             quote = "", na.strings = "", fill = TRUE)
pth[, hgnc_id := ifelse(grepl("\\|HGNC=[0-9]+\\|", V1),
                        sub(".*\\|HGNC=([0-9]+)\\|.*", "HGNC:\\1", V1), NA_character_)]
pth[, symbol := V3]
pth[, panther_family := V5]
# Column 10 lists the protein class first, then its parent classes
# ("child#PC1;parent#PC2"); keep the most specific (first) class.
pth[, panther_protein_class := trimws(sub(";.*$", "", gsub("#PC[0-9]+", "", V10)))]
pth[panther_protein_class == "", panther_protein_class := NA_character_]
# Several UniProt entries can map to one gene: prefer an annotated entry.
pth[, n_ann := (!is.na(panther_family)) + (!is.na(panther_protein_class))]
setorderv(pth, "n_ann", -1L)
pth_by_id  <- pth[!is.na(hgnc_id) & !duplicated(hgnc_id), .(hgnc_id, panther_family, panther_protein_class)]
pth_by_sym <- pth[!is.na(symbol)  & !duplicated(symbol),  .(symbol,  panther_family, panther_protein_class)]

# ---- Lambert 2018 TFs -----------------------------------------------------------
lam <- fread(need(file.path(ref_dir, "DatabaseExtract_v_1.01.csv")),
             select = c("Ensembl ID", "DBD", "Is TF?"))
setnames(lam, c("gene_id", "lambert_dbd", "is_tf"))
lam <- lam[!duplicated(gene_id)]
lam[, is_tf_lambert := is_tf == "Yes"]

# ---- join -------------------------------------------------------------------------
ann <- merge(cls, hgnc, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
ann[is.na(locus_group), locus_group := "unannotated (no HGNC)"]
ann <- merge(ann, pth_by_id, by = "hgnc_id", all.x = TRUE, sort = FALSE)
# Fallback for genes whose PANTHER entry is keyed to another HGNC id (column 1
# sometimes carries a read-through locus, e.g. BCL2L2): match on the symbol.
miss <- is.na(ann$panther_family) & is.na(ann$panther_protein_class) & !is.na(ann$symbol)
fb <- pth_by_sym[match(ann$symbol[miss], pth_by_sym$symbol)]
ann[miss, `:=`(panther_family = fb$panther_family, panther_protein_class = fb$panther_protein_class)]
ann <- merge(ann, lam[, .(gene_id, lambert_dbd, is_tf_lambert)], by = "gene_id",
             all.x = TRUE, sort = FALSE)
ann[is.na(is_tf_lambert), is_tf_lambert := FALSE]
ann <- ann[match(cls$gene_id, gene_id)]          # keep the classification row order

per_method <- grep("^(pred_power|trinary_class)_", names(cls), value = TRUE)
core <- setdiff(names(cls), c("gene_id", per_method))
setcolorder(ann, c("gene_id", "symbol", "hgnc_id", "locus_group", core,
                   "panther_protein_class", "panther_family", "is_tf_lambert",
                   "lambert_dbd", per_method))
stopifnot(nrow(ann) == nrow(cls), !anyDuplicated(ann$gene_id))

fwrite(ann, file.path(gc_dir, "oligo_gene_classification_annotated.csv"))
fwrite(ann[is_tf_lambert == TRUE], file.path(gc_dir, "oligo_lambert_tfs.csv"))
cat(sprintf("annotated %d genes: %d with HGNC symbol, %d PANTHER class, %d Lambert TFs\n",
            nrow(ann), sum(!is.na(ann$symbol)), sum(!is.na(ann$panther_protein_class)),
            sum(ann$is_tf_lambert)))
