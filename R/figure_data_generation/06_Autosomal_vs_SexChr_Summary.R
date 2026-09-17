# ==============================================================================
# STEP 06: Autosomal vs Sex-Chromosome Truth Gene Properties
# ==============================================================================
# Diagnostic: are autosomal truth genes lower-expressed (baseMean) or more
# variable (dispersion) than sex chromosome truth genes, within the same
# |LFC| bin?
#
# Purpose: justify that sex chromosome genes are an optimistic bound for
# power analysis, i.e. they succeed in the Tiny and Small bins because they
# have systematically different distributional properties, not because the
# power equation generalises there.
#
# Inputs:
#   - Datalink/Power_Curve_Fitting/GT_DDS_Checkpoints/{CellType}/GT_DDS_Run_*.qs
#     (produced by power_curve_fitting/00_Generate_GT_DDS_Checkpoints.R)
#   - Datalink/Full_Ground_Truth/{CellType}/Roussos_DE_results.csv
#     (produced by full_ground_truth/01_Run_DE_All_Datasets.R)
#   - Datalink/Full_Ground_Truth/gene_audit_biomart_cache.csv
#     (produced by full_ground_truth/04_Gene_Audit_Matrix.R)
#
# Outputs:
#   - Datalink/_scratch_autosomal_vs_sexchr_summary.csv       (main table)
#   - Datalink/_scratch_autosomal_vs_sexchr_summary_pivot.csv (wide ratio view)
#
# The main CSV is consumed directly by paper_final_draft.rmd via
#   here::here("Datalink", "_scratch_autosomal_vs_sexchr_summary.csv")
# so the filename is preserved verbatim for paper compatibility.
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 06: AUTOSOMAL vs SEX-CHROMOSOME TRUTH GENE PROPERTIES ###\n")

# ==============================================================================
# Paths (resolved via PROJECT_ROOT, set in config.R)
# ==============================================================================

gt_dds_base    <- project_path("Datalink", "Power_Curve_Fitting", "GT_DDS_Checkpoints")
de_results_dir <- project_path("Datalink", "Full_Ground_Truth")
biomart_cache  <- project_path("Datalink", "Full_Ground_Truth", "gene_audit_biomart_cache.csv")
out_path       <- project_path("Datalink", "_scratch_autosomal_vs_sexchr_summary.csv")

cell_types <- CELL_TYPES
runs       <- 1:6

# ==============================================================================
# Helpers
# ==============================================================================

classify_chr <- function(chr) {
  dplyr::case_when(
    as.character(chr) == "X"                    ~ "Sex_chr",
    as.character(chr) == "Y"                    ~ "Sex_chr",
    as.character(chr) %in% c("MT", "M")         ~ "Mito",
    as.character(chr) %in% as.character(1:22)   ~ "Autosomal",
    TRUE                                        ~ "Other"
  )
}

bin_lfc <- function(lfc) {
  a <- abs(lfc)
  dplyr::case_when(
    a >= 0.2 & a < 1.0 ~ "Tiny (0.2-1)",
    a >= 1.0 & a < 2.0 ~ "Small (1-2)",
    a >= 2.0 & a < 3.0 ~ "Medium (2-3)",
    a >= 3.0           ~ "Large (3+)",
    TRUE               ~ "Below"
  )
}

# Safe summary: handles empty groups
summary_row <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) {
    return(list(n = 0L, mean = NA, median = NA,
                q25 = NA, q75 = NA, min = NA, max = NA))
  }
  list(
    n = length(x), mean = mean(x), median = median(x),
    q25 = as.numeric(quantile(x, 0.25)),
    q75 = as.numeric(quantile(x, 0.75)),
    min = min(x), max = max(x)
  )
}

# ==============================================================================
# Load annotations
# ==============================================================================

if (!file.exists(biomart_cache)) {
  stop(sprintf("BioMart cache not found: %s\n  Run full_ground_truth/04_Gene_Audit_Matrix.R first.",
               biomart_cache))
}

message(sprintf("[%s] Loading BioMart cache...", format(Sys.time(), "%H:%M:%S")))
cache <- read.csv(biomart_cache, stringsAsFactors = FALSE)
cache$chr_group <- classify_chr(cache$chromosome_name)
chr_lookup <- setNames(cache$chr_group, cache$ensembl_gene_id)

# ==============================================================================
# Assemble gene-level table per cell type
# ==============================================================================

all_rows <- list()

for (ct in cell_types) {
  message(sprintf("\n=== %s ===", ct))

  de_path <- file.path(de_results_dir, ct, "Roussos_DE_results.csv")
  if (!file.exists(de_path)) {
    message(sprintf("  [skip] missing DE results: %s", de_path))
    next
  }
  de <- read.csv(de_path, stringsAsFactors = FALSE)
  # Column "gene" = ensembl id, baseMean, log2FoldChange, padj
  de$chr_group <- unname(chr_lookup[de$gene])
  de$chr_group[is.na(de$chr_group)] <- "Unannotated"
  de$bin <- bin_lfc(de$log2FoldChange)

  # Pool dispersions across all 6 runs (median per gene for stability).
  disp_list <- list()
  for (r in runs) {
    dds_path <- file.path(gt_dds_base, ct, sprintf("GT_DDS_Run_%d.qs", r))
    if (!file.exists(dds_path)) next
    dds <- qs::qread(dds_path)
    d <- DESeq2::dispersions(dds)
    names(d) <- rownames(dds)
    disp_list[[as.character(r)]] <- data.frame(
      gene = names(d), dispersion = as.numeric(d),
      run = r, stringsAsFactors = FALSE
    )
  }
  if (length(disp_list) == 0) {
    message("  [skip] no GT_DDS checkpoints found")
    next
  }
  disp_df <- dplyr::bind_rows(disp_list) %>%
    dplyr::filter(is.finite(dispersion), dispersion > 0) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dispersion = median(dispersion), .groups = "drop")

  gene_tab <- de %>%
    dplyr::left_join(disp_df, by = "gene") %>%
    dplyr::filter(
      !is.na(padj), padj < 0.05,
      !is.na(log2FoldChange), abs(log2FoldChange) >= 0.2,
      !is.na(dispersion), dispersion > 0,
      !is.na(baseMean), baseMean > 0,
      chr_group %in% c("Autosomal", "Sex_chr")
    )

  message(sprintf("  Truth genes (padj<0.05 & |LFC|>=0.2): %d (auto=%d, sex=%d)",
                  nrow(gene_tab),
                  sum(gene_tab$chr_group == "Autosomal"),
                  sum(gene_tab$chr_group == "Sex_chr")))

  # Per-bin / per-group summary
  for (b in unique(gene_tab$bin)) {
    for (g in c("Autosomal", "Sex_chr")) {
      sub <- gene_tab[gene_tab$bin == b & gene_tab$chr_group == g, ]
      if (nrow(sub) == 0) next
      bm  <- summary_row(sub$baseMean)
      ds  <- summary_row(sub$dispersion)
      lfc <- summary_row(abs(sub$log2FoldChange))
      row <- data.frame(
        cell_type         = ct,
        bin               = b,
        chr_group         = g,
        n                 = bm$n,
        baseMean_mean     = bm$mean,   baseMean_median   = bm$median,
        baseMean_q25      = bm$q25,    baseMean_q75      = bm$q75,
        dispersion_mean   = ds$mean,   dispersion_median = ds$median,
        dispersion_q25    = ds$q25,    dispersion_q75    = ds$q75,
        absLFC_mean       = lfc$mean,  absLFC_median     = lfc$median,
        stringsAsFactors  = FALSE
      )
      all_rows[[length(all_rows) + 1]] <- row
    }
  }
}

if (length(all_rows) == 0) {
  stop("No rows produced. Check that GT_DDS checkpoints and DE results exist for at least one cell type.")
}

agg <- dplyr::bind_rows(all_rows) %>%
  dplyr::mutate(bin = factor(bin, levels = c("Tiny (0.2-1)", "Small (1-2)",
                                             "Medium (2-3)", "Large (3+)"))) %>%
  dplyr::arrange(cell_type, bin, chr_group)

write.csv(agg, out_path, row.names = FALSE)
message(sprintf("\n[OK] Wrote %s", out_path))

# ==============================================================================
# Ratio view: Autosomal vs Sex_chr within each (cell, bin)
# ==============================================================================

pivot <- agg %>%
  dplyr::select(cell_type, bin, chr_group, baseMean_median, dispersion_median, n) %>%
  tidyr::pivot_wider(names_from = chr_group,
                     values_from = c(baseMean_median, dispersion_median, n))

if (nrow(pivot) > 0) {
  if (all(c("baseMean_median_Autosomal", "baseMean_median_Sex_chr") %in% names(pivot))) {
    pivot <- pivot %>%
      dplyr::mutate(
        baseMean_ratio_auto_over_sex   = baseMean_median_Autosomal   / baseMean_median_Sex_chr,
        dispersion_ratio_auto_over_sex = dispersion_median_Autosomal / dispersion_median_Sex_chr
      )
  }
  pivot_path <- sub("\\.csv$", "_pivot.csv", out_path)
  write.csv(pivot, pivot_path, row.names = FALSE)
  message(sprintf("[OK] Wrote %s", pivot_path))
}

message("\n[COMPLETE] Step 06")
