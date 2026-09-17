# ==============================================================================
# STEP 2: Cross-Dataset Diagnostic Analyses
# ==============================================================================
# Performs cross-dataset diagnostics on DE results from Step 1.
# Generates UpSet plots, overlap heatmaps, concordance tables, LFC correlation
# matrices, sign agreement matrices, forest plots, XIST detection summaries,
# and expanded LFC agreement matrices with upper/lower bounds from cartesian
# downsampling split-half comparisons and MI-minimised null label flips.
#
# Input:  {OUTPUT_BASE}/{CellType}/{Dataset}_DE_results.csv  (from Step 1)
#         {CARTESIAN_BASE}/{CellType}/Roussos/reads/Run_{1-3}/  (for upper bound)
# Output: {OUTPUT_BASE}/{CellType}/Diagnostic_Results/*.png, *.csv
# ==============================================================================

# ==============================================================================
# 0. SOURCE PIPELINE INFRASTRUCTURE
# ==============================================================================

# Define PROJECT_ROOT if not already defined
if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- PROJECT_ROOT
}

source(file.path(PROJECT_ROOT, "R", "full_ground_truth", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

# Source shared prettify functions
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(viridis)
  library(UpSetR)
  library(RColorBrewer)
  library(ggrepel)
  library(glue)
  library(grid)
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
  library(Matrix)
})

log_section("STEP 2: Cross-Dataset Diagnostic Analyses")
log_info("Output base: %s", OUTPUT_BASE)
log_info("Cell types: %s", paste(CELL_TYPES, collapse = ", "))
log_info("Reference dataset: %s", REFERENCE_DATASET)
log_info("FDR threshold: %g, LFC threshold: %g", FDR_THRESHOLD, LFC_THRESHOLD)

# ==============================================================================
# 2. MAIN LOOP OVER CELL TYPES
# ==============================================================================

completion_summary <- list()
fdp_collector <- list()

for (cell_type in CELL_TYPES) {

  log_section(sprintf("Processing: %s", cell_type))

  # --------------------------------------------------------------------------
  # 2a. Setup paths
  # --------------------------------------------------------------------------
  paths <- get_paths(cell_type)
  diag_dir <- paths$diagnostic_dir
  if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

  # --------------------------------------------------------------------------
  # 2b. Load all DE results
  # --------------------------------------------------------------------------
  all_results <- load_all_de_results(paths$output_dir)

  if (length(all_results) < 2) {
    log_warn("%s: Only %d dataset(s) found. Need >= 2 for diagnostics. Skipping.",
             cell_type, length(all_results))
    completion_summary[[cell_type]] <- "SKIPPED (insufficient datasets)"
    next
  }

  dataset_names <- names(all_results)
  log_info("Loaded %d datasets: %s", length(dataset_names),
           paste(dataset_names, collapse = ", "))

  # --------------------------------------------------------------------------
  # 2c. Extract significant genes
  # --------------------------------------------------------------------------
  sig_list <- extract_significant_genes(all_results,
                                        padj_cutoff = FDR_THRESHOLD,
                                        lfc_cutoff  = LFC_THRESHOLD)

  for (ds in names(sig_list)) {
    log_info("  %s: %d significant genes", ds, length(sig_list[[ds]]))
  }

  union_sig <- unique(unlist(sig_list))
  n_sets    <- length(sig_list)

  # ==========================================================================
  # DIAGNOSTIC 1: UpSet Plot (overall significant genes)
  # ==========================================================================
  log_info("[1/16] UpSet Plot (overall)...")

  tryCatch({
    if (length(union_sig) > 0 && n_sets > 1) {
      binary_mat <- matrix(0L, nrow = length(union_sig), ncol = n_sets)
      rownames(binary_mat) <- union_sig
      colnames(binary_mat) <- prettify_label(names(sig_list))
      for (nm in names(sig_list)) {
        binary_mat[sig_list[[nm]], prettify_label(nm)] <- 1L
      }

      png(file.path(diag_dir, "Diagnostic_UpSet_Plot.png"),
          width = 1600, height = 900, res = 150)
      print(upset(as.data.frame(binary_mat),
                  sets = colnames(binary_mat), order.by = "freq",
                  mainbar.y.label = "Intersection Size",
                  sets.x.label = "Set Size"))
      dev.off()
      log_info("  Saved Diagnostic_UpSet_Plot.png")
    } else {
      log_warn("  Not enough significant genes for UpSet plot")
    }
  }, error = function(e) {
    log_error("  UpSet plot (overall) failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 2: UpSet Plot (male-biased genes, LFC > lfc_cutoff)
  # ==========================================================================
  log_info("[2/16] UpSet Plot (male-biased)...")

  tryCatch({
    male_list <- lapply(all_results, function(df) {
      sig <- df[!is.na(df$padj) & df$padj < FDR_THRESHOLD &
                  df$log2FoldChange > LFC_THRESHOLD, ]
      rownames(sig)
    })

    union_male <- unique(unlist(male_list))
    if (length(union_male) > 0 && n_sets > 1) {
      bmat_m <- matrix(0L, nrow = length(union_male), ncol = n_sets)
      rownames(bmat_m) <- union_male
      colnames(bmat_m) <- prettify_label(names(male_list))
      for (nm in names(male_list)) bmat_m[male_list[[nm]], prettify_label(nm)] <- 1L

      png(file.path(diag_dir, "UpSet_Male_Biased_Genes.png"),
          width = 1600, height = 900, res = 150)
      print(upset(as.data.frame(bmat_m),
                  sets = colnames(bmat_m), order.by = "freq",
                  mainbar.y.label = "Intersection Size",
                  sets.x.label = "Set Size"))
      dev.off()
      log_info("  Saved UpSet_Male_Biased_Genes.png")
    } else {
      log_warn("  No male-biased genes found for UpSet plot")
    }
  }, error = function(e) {
    log_error("  UpSet plot (male-biased) failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 3: UpSet Plot (female-biased genes, LFC < -lfc_cutoff)
  # ==========================================================================
  log_info("[3/16] UpSet Plot (female-biased)...")

  tryCatch({
    female_list <- lapply(all_results, function(df) {
      sig <- df[!is.na(df$padj) & df$padj < FDR_THRESHOLD &
                  df$log2FoldChange < -LFC_THRESHOLD, ]
      rownames(sig)
    })

    union_female <- unique(unlist(female_list))
    if (length(union_female) > 0 && n_sets > 1) {
      bmat_f <- matrix(0L, nrow = length(union_female), ncol = n_sets)
      rownames(bmat_f) <- union_female
      colnames(bmat_f) <- prettify_label(names(female_list))
      for (nm in names(female_list)) bmat_f[female_list[[nm]], prettify_label(nm)] <- 1L

      png(file.path(diag_dir, "UpSet_Female_Biased_Genes.png"),
          width = 1600, height = 900, res = 150)
      print(upset(as.data.frame(bmat_f),
                  sets = colnames(bmat_f), order.by = "freq",
                  mainbar.y.label = "Intersection Size",
                  sets.x.label = "Set Size"))
      dev.off()
      log_info("  Saved UpSet_Female_Biased_Genes.png")
    } else {
      log_warn("  No female-biased genes found for UpSet plot")
    }
  }, error = function(e) {
    log_error("  UpSet plot (female-biased) failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 4: Overlap Heatmap
  # ==========================================================================
  log_info("[4/16] Overlap Heatmap...")

  tryCatch({
    overlap_mat <- build_overlap_matrix(sig_list)

      rownames(overlap_mat) <- prettify_label(rownames(overlap_mat))
      colnames(overlap_mat) <- prettify_label(colnames(overlap_mat))

      png(file.path(diag_dir, "Diagnostic_Overlap_Heatmap.png"),
          width = 1000, height = 900, res = 120)
      pheatmap(overlap_mat,
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%d",
               main = paste0("Significant Gene Overlap: ", prettify_label(cell_type)),
               color = colorRampPalette(brewer.pal(9, "Blues"))(100),
               fontsize = 14, fontsize_number = 12,
               fontsize_row = 13, fontsize_col = 13)
    dev.off()
    log_info("  Saved Diagnostic_Overlap_Heatmap.png")
  }, error = function(e) {
    log_error("  Overlap heatmap failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 5: Concordance Analysis vs Reference
  # ==========================================================================
  log_info("[5/16] Concordance Analysis vs %s...", REFERENCE_DATASET)

  tryCatch({
    concordance_df <- compute_concordance(sig_list, reference = REFERENCE_DATASET)

    if (nrow(concordance_df) > 0) {
      conc_file <- file.path(diag_dir,
                             paste0("Concordance_Analysis_vs_", REFERENCE_DATASET, ".csv"))
      write.csv(concordance_df, conc_file, row.names = FALSE, quote = TRUE)
      log_info("  Saved Concordance_Analysis_vs_%s.csv", REFERENCE_DATASET)

      for (i in seq_len(nrow(concordance_df))) {
        log_info("    %s: shared=%d, concordance=%.3f",
                 concordance_df$Dataset[i],
                 concordance_df$Shared[i],
                 concordance_df$Concordance_Rate[i])
      }
    } else {
      log_warn("  Reference dataset '%s' not found in results", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  Concordance analysis failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 6: LFC Correlation (Top N genes from reference by baseMean)
  # ==========================================================================
  log_info("[6/16] LFC Correlation (Top %d genes)...", N_TOP_GENES)

  ref_name <- grep(paste0("^", REFERENCE_DATASET),
                   names(all_results), value = TRUE, ignore.case = TRUE)[1]

  tryCatch({
    if (!is.na(ref_name)) {
      ref_df <- all_results[[ref_name]]

      # Select top genes: significant, sorted by baseMean
      ref_sig <- ref_df[!is.na(ref_df$padj) & ref_df$padj < FDR_THRESHOLD &
                          abs(ref_df$log2FoldChange) > LFC_THRESHOLD, ]

      if ("baseMean" %in% colnames(ref_sig)) {
        ref_sig <- ref_sig[order(-ref_sig$baseMean), ]
      } else {
        ref_sig <- ref_sig[order(ref_sig$padj), ]
      }

      top_genes <- head(ref_sig$gene_id, N_TOP_GENES)

      if (length(top_genes) > 0) {
        lfc_mat <- build_lfc_matrix(top_genes, all_results)
        cor_mat <- cor(lfc_mat, use = "pairwise.complete.obs")
        cor_mat[is.na(cor_mat)] <- 0

        # Prettify axis labels
        rownames(cor_mat) <- prettify_label(rownames(cor_mat))
        colnames(cor_mat) <- prettify_label(colnames(cor_mat))

        # Plot
        png(file.path(diag_dir, "LFC_Correlation_Top100.png"),
            width = 1200, height = 900, res = 120)
        pheatmap(cor_mat,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = glue("LFC Correlation (Top {N_TOP_GENES} {prettify_label(ref_name)} Genes): {prettify_label(cell_type)}"),
                 color = colorRampPalette(brewer.pal(11, "RdBu"))(100),
                 fontsize = 14, fontsize_number = 12,
                 fontsize_row = 13, fontsize_col = 13)
        dev.off()
        log_info("  Saved LFC_Correlation_Top100.png")

        # CSV
        write.csv(cor_mat,
                  file.path(diag_dir, "LFC_Correlation_Matrix.csv"),
                  quote = TRUE)
        log_info("  Saved LFC_Correlation_Matrix.csv")
      } else {
        log_warn("  No significant genes in reference for LFC correlation")
      }
    } else {
      log_warn("  Reference dataset '%s' not found", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  LFC correlation failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 7: LFC Agreement Matrix (all reference significant genes)
  # ==========================================================================
  log_info("[7/16] LFC Agreement Matrix (%s sig genes)...", REFERENCE_DATASET)

  tryCatch({
    if (!is.na(ref_name) && ref_name %in% names(sig_list)) {
      ref_sig_genes <- sig_list[[ref_name]]

      if (length(ref_sig_genes) > 0) {
        lfc_mat_all <- build_lfc_matrix(ref_sig_genes, all_results)
        lfc_cor_all <- cor(lfc_mat_all, use = "pairwise.complete.obs")
        lfc_cor_all[is.na(lfc_cor_all)] <- 0

        rownames(lfc_cor_all) <- prettify_label(rownames(lfc_cor_all))
        colnames(lfc_cor_all) <- prettify_label(colnames(lfc_cor_all))

        png(file.path(diag_dir,
                      glue("LFC_Agreement_Matrix_{REFERENCE_DATASET}_Sig_Genes.png")),
            width = 1200, height = 900, res = 120)
        pheatmap(lfc_cor_all,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = glue("LFC Pearson Correlation ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
                 color = viridis(100),
                 breaks = seq(0.5, 1, length.out = 101),
                 fontsize = 14, fontsize_number = 12,
                 fontsize_row = 13, fontsize_col = 13)
        dev.off()
        log_info("  Saved LFC_Agreement_Matrix_%s_Sig_Genes.png", REFERENCE_DATASET)
      } else {
        log_warn("  No significant genes in reference for LFC agreement")
      }
    } else {
      log_warn("  Reference dataset '%s' not found for LFC agreement", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  LFC agreement matrix failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 8: Sign Agreement Matrix
  # ==========================================================================
  log_info("[8/16] Sign Agreement Matrix (%s sig genes)...", REFERENCE_DATASET)

  tryCatch({
    if (!is.na(ref_name) && ref_name %in% names(sig_list)) {
      ref_sig_genes <- sig_list[[ref_name]]

      if (length(ref_sig_genes) > 0) {
        lfc_mat_sign <- build_lfc_matrix(ref_sig_genes, all_results)
        sign_agree   <- build_sign_agreement(lfc_mat_sign)

        # Replace NAs with 0 for pheatmap clustering (NA cells shown as grey)
        sign_agree_plot <- sign_agree
        has_na <- any(is.na(sign_agree_plot))
        sign_agree_plot[is.na(sign_agree_plot)] <- 0

        rownames(sign_agree_plot) <- prettify_label(rownames(sign_agree_plot))
        colnames(sign_agree_plot) <- prettify_label(colnames(sign_agree_plot))

        png(file.path(diag_dir,
                      glue("Sign_Agreement_Matrix_{REFERENCE_DATASET}_Sig_Genes.png")),
            width = 1200, height = 900, res = 120)
        pheatmap(sign_agree_plot,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = glue("Sign Agreement ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
                 color = viridis(100),
                 breaks = seq(0, 1, length.out = 101),
                 na_col = "grey80",
                 fontsize = 14, fontsize_number = 12,
                 fontsize_row = 13, fontsize_col = 13)
        dev.off()
        log_info("  Saved Sign_Agreement_Matrix_%s_Sig_Genes.png", REFERENCE_DATASET)
        if (has_na) log_warn("  Note: some dataset pairs had no shared genes (shown as 0.00)")
      } else {
        log_warn("  No significant genes in reference for sign agreement")
      }
    } else {
      log_warn("  Reference dataset '%s' not found for sign agreement", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  Sign agreement matrix failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 9: Forest Plots (Top 20 genes from reference)
  # ==========================================================================
  log_info("[9/16] Forest Plots (top 20 genes)...")

  tryCatch({
    if (!is.na(ref_name)) {
      ref_df <- all_results[[ref_name]]

      # Top 20 significant genes by baseMean
      ref_sig <- ref_df[!is.na(ref_df$padj) & ref_df$padj < FDR_THRESHOLD &
                          abs(ref_df$log2FoldChange) > LFC_THRESHOLD, ]

      if ("baseMean" %in% colnames(ref_sig)) {
        ref_sig <- ref_sig[order(-ref_sig$baseMean), ]
      } else {
        ref_sig <- ref_sig[order(ref_sig$padj), ]
      }

      top20_genes <- head(ref_sig$gene_id, 20)

      if (length(top20_genes) > 0) {
        # Build long-format data for ggplot
        forest_data <- data.frame()

        for (ds in names(all_results)) {
          ds_df <- all_results[[ds]]
          for (g in top20_genes) {
            if (g %in% rownames(ds_df)) {
              row <- ds_df[g, ]
              lfc <- row$log2FoldChange
              se  <- if ("lfcSE" %in% colnames(row)) row$lfcSE else NA_real_

              forest_data <- rbind(forest_data, data.frame(
                Gene    = g,
                Dataset = prettify_label(ds),
                LFC     = lfc,
                SE      = se,
                LFC_lo  = lfc - 1.96 * ifelse(is.na(se), 0, se),
                LFC_hi  = lfc + 1.96 * ifelse(is.na(se), 0, se),
                stringsAsFactors = FALSE
              ))
            }
          }
        }

        if (nrow(forest_data) > 0) {
          # Maintain gene order from reference ranking
          forest_data$Gene <- factor(forest_data$Gene, levels = rev(top20_genes))

          p <- ggplot(forest_data, aes(x = LFC, y = Dataset,
                                       xmin = LFC_lo, xmax = LFC_hi)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            geom_errorbarh(height = 0.2, color = "grey40") +
            geom_point(aes(color = Dataset), size = 2.5) +
            facet_wrap(~ Gene, ncol = 4, scales = "free_x") +
            labs(title = glue("Forest Plot: Top 20 {prettify_label(REFERENCE_DATASET)} Genes ({prettify_label(cell_type)})"),
                 x = "log2 Fold Change (Male vs Female)",
                 y = NULL) +
            theme_publication(base_size = 16) +
            theme(legend.position = "bottom",
                  strip.text = element_text(size = 10),
                  axis.text.y = element_text(size = 10))

          ggsave(file.path(diag_dir, "Forest_Plot_Top20_Genes.png"),
                 plot = p, width = 16, height = 12, dpi = 150)
          log_info("  Saved Forest_Plot_Top20_Genes.png")
        } else {
          log_warn("  No forest plot data could be assembled")
        }
      } else {
        log_warn("  No significant genes in reference for forest plot")
      }
    } else {
      log_warn("  Reference dataset not found for forest plot")
    }
  }, error = function(e) {
    log_error("  Forest plot failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 10: XIST Detection Summary
  # ==========================================================================
  log_info("[10/16] XIST Detection Summary...")

  tryCatch({
    xist_rows <- data.frame()

    for (ds in names(all_results)) {
      ds_df <- all_results[[ds]]
      idx   <- which(ds_df$gene_id == TARGET_GENE_ID)

      if (length(idx) > 0) {
        row <- ds_df[idx[1], ]
        xist_rows <- rbind(xist_rows, data.frame(
          Dataset        = ds,
          Found          = TRUE,
          log2FoldChange = row$log2FoldChange,
          padj           = row$padj,
          baseMean       = if ("baseMean" %in% colnames(row)) row$baseMean else NA_real_,
          stringsAsFactors = FALSE
        ))
      } else {
        xist_rows <- rbind(xist_rows, data.frame(
          Dataset        = ds,
          Found          = FALSE,
          log2FoldChange = NA_real_,
          padj           = NA_real_,
          baseMean       = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
    }

    xist_file <- file.path(diag_dir, "XIST_Detection_Summary.csv")
    write.csv(xist_rows, xist_file, row.names = FALSE, quote = TRUE)
    log_info("  Saved XIST_Detection_Summary.csv")

    for (i in seq_len(nrow(xist_rows))) {
      r <- xist_rows[i, ]
      if (r$Found) {
        log_info("    %s: %s found (LFC=%.3f, padj=%.2e)",
                 r$Dataset, TARGET_GENE_NAME, r$log2FoldChange, r$padj)
      } else {
        log_warn("    %s: %s NOT found", r$Dataset, TARGET_GENE_NAME)
      }
    }
  }, error = function(e) {
    log_error("  XIST detection failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 11: Significant Gene Counts
  # ==========================================================================
  log_info("[11/16] Significant Gene Counts...")

  tryCatch({
    sig_counts <- data.frame(
      Dataset          = names(sig_list),
      Significant_Genes = sapply(sig_list, length),
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
    sig_counts <- sig_counts[order(-sig_counts$Significant_Genes), ]

    sig_file <- file.path(diag_dir, "Significant_Gene_Counts.csv")
    write.csv(sig_counts, sig_file, row.names = FALSE, quote = TRUE)
    log_info("  Saved Significant_Gene_Counts.csv")

    for (i in seq_len(nrow(sig_counts))) {
      log_info("    %s: %d genes", sig_counts$Dataset[i], sig_counts$Significant_Genes[i])
    }
  }, error = function(e) {
    log_error("  Significant gene counts failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 12: Upper Bound - Split-Half DE from Cartesian Runs 1-3
  # ==========================================================================
  log_info("[12/16] Upper Bound: Split-half DE (Runs 1-%d)...", N_UPPER_BOUND_RUNS)

  upper_bound_results <- list()
  de_param_local <- setup_parallel(min(N_CORES_DE, 6))

  tryCatch({

    for (run_id in seq_len(N_UPPER_BOUND_RUNS)) {
      log_info("  --- Upper Bound Run %d ---", run_id)

      # Get donor split from cartesian run
      split_info <- get_cartesian_split_donors(cell_type, run_id)
      if (is.null(split_info)) {
        log_warn("  Could not load split info for Run %d. Skipping.", run_id)
        next
      }

      # Half A: Load existing cartesian Ground_Truth results (GT half)
      half_a_name <- sprintf("Roussos_Half_A_R%d", run_id)
      half_a_df <- load_cartesian_gt_results(cell_type, run_id)
      if (!is.null(half_a_df)) {
        upper_bound_results[[half_a_name]] <- half_a_df
        log_info("  %s: %d genes loaded from cartesian GT", half_a_name, nrow(half_a_df))
      }

      # Half B: Run fresh DESeq2 on validation donors
      half_b_name <- sprintf("Roussos_Half_B_R%d", run_id)
      half_b_cache <- file.path(diag_dir,
                                 sprintf("Expanded_HalfB_R%d_DE_results.qs", run_id))

      half_b_df <- run_deseq_for_roussos_subset(
        donor_ids    = split_info$val_donors,
        cell_type    = cell_type,
        de_param     = de_param_local,
        cache_path   = half_b_cache,
        force_recalc = FORCE_RECALC
      )

      if (!is.null(half_b_df)) {
        upper_bound_results[[half_b_name]] <- half_b_df
        log_info("  %s: %d genes", half_b_name, nrow(half_b_df))
      }

      gc(verbose = FALSE)
    }

    log_info("  Upper bound: %d result sets loaded", length(upper_bound_results))
  }, error = function(e) {
    log_error("  Upper bound computation failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 13: Lower Bound - MI-Minimised Null Label Flips
  # ==========================================================================
  log_info("[13/16] Lower Bound: MI-minimised null (N=%d runs)...", N_NULL_RUNS)

  null_results <- list()

  tryCatch({
    # Get all Roussos donor IDs and their sex
    all_roussos_ids <- get_all_roussos_donors(cell_type)

    if (!is.null(all_roussos_ids) && length(all_roussos_ids) > 0) {
      # Get sex info (needed for MI-optimal flip calculation)
      roussos_sex_df <- get_roussos_donor_sex(cell_type, all_roussos_ids)
      cfg_roussos <- DATASET_CONFIG[["Roussos"]]
      male_val  <- cfg_roussos$sex_values$male
      female_val <- cfg_roussos$sex_values$female

      n_male   <- sum(roussos_sex_df$sex == male_val)
      n_female <- sum(roussos_sex_df$sex == female_val)
      log_info("  Roussos donors: %d total (%d M, %d F)", nrow(roussos_sex_df), n_male, n_female)

      for (null_id in seq_len(N_NULL_RUNS)) {
        log_info("  --- Null Run %d ---", null_id)

        null_seed <- SEED_BASE + NULL_SEED_OFFSET + null_id
        null_name <- sprintf("Null_%d", null_id)
        null_cache <- file.path(diag_dir,
                                 sprintf("Expanded_Null%d_DE_results.qs", null_id))

        # Create label flip function for this null run
        # Use local() to capture current iteration values (avoid closure pitfall)
        flip_fn <- local({
          .seed      <- null_seed
          .sex_col   <- cfg_roussos$cols$sex
          .male_val  <- male_val
          .female_val <- female_val
          function(meta_df) {
            flip_sex_labels_mi_minimised(
              meta_df    = meta_df,
              sex_col    = .sex_col,
              male_val   = .male_val,
              female_val = .female_val,
              seed       = .seed
            )
          }
        })

        null_df <- run_deseq_for_roussos_subset(
          donor_ids      = all_roussos_ids,
          cell_type      = cell_type,
          label_flip_fn  = flip_fn,
          de_param       = de_param_local,
          cache_path     = null_cache,
          force_recalc   = FORCE_RECALC
        )

        if (!is.null(null_df)) {
          null_results[[null_name]] <- null_df
          log_info("  %s: %d genes", null_name, nrow(null_df))
        }

        gc(verbose = FALSE)
      }
    } else {
      log_warn("  No Roussos donors found for null computation")
    }

    log_info("  Lower bound: %d null result sets", length(null_results))
  }, error = function(e) {
    log_error("  Lower bound computation failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 14 & 15: Expanded LFC Agreement Matrices (Pearson & Spearman)
  # ==========================================================================
  log_info("[14-15/16] Expanded LFC Agreement Matrices...")

  tryCatch({
    if (!is.na(ref_name) && ref_name %in% names(sig_list) &&
        (length(upper_bound_results) > 0 || length(null_results) > 0)) {

      ref_sig_genes <- sig_list[[ref_name]]

      if (length(ref_sig_genes) > 0) {
        # Combine all result sets: upper bound + datasets + null
        combined_results <- c(upper_bound_results, all_results, null_results)

        # Build expanded LFC matrix
        expanded_lfc_mat <- build_lfc_matrix(ref_sig_genes, combined_results)

        # Define column groups for annotation
        n_upper <- length(upper_bound_results)
        n_datasets <- length(all_results)
        n_null <- length(null_results)

        group_labels <- c(
          rep("Upper Bound", n_upper),
          rep("Cross-Dataset", n_datasets),
          rep("Lower Bound", n_null)
        )

        # Create annotation data frame for pheatmap
        col_annotation <- data.frame(
          Group = factor(group_labels,
                         levels = c("Upper Bound",
                                    "Cross-Dataset",
                                    "Lower Bound")),
          row.names = colnames(expanded_lfc_mat)
        )

        # Annotation colours
        ann_colors <- list(
          Group = c(
            "Upper Bound" = "#2166AC",
            "Cross-Dataset" = "#4DAF4A",
            "Lower Bound" = "#E41A1C"
          )
        )

        # --- Pearson correlation ---
        log_info("  Computing Pearson correlation matrix...")
        pearson_cor <- cor(expanded_lfc_mat, use = "pairwise.complete.obs",
                           method = "pearson")
        pearson_cor[is.na(pearson_cor)] <- 0

        # Save CSV
        write.csv(pearson_cor,
                  file.path(diag_dir, "Expanded_LFC_Pearson_Matrix.csv"),
                  quote = TRUE)
        log_info("  Saved Expanded_LFC_Pearson_Matrix.csv")

        # Prettify axis labels for display
        rownames(pearson_cor) <- prettify_label(rownames(pearson_cor))
        colnames(pearson_cor) <- prettify_label(colnames(pearson_cor))
        col_annotation_display <- col_annotation
        rownames(col_annotation_display) <- prettify_label(rownames(col_annotation_display))

        # Plot Pearson
        png(file.path(diag_dir, "Expanded_LFC_Pearson_Matrix.png"),
            width = 1800, height = 1400, res = 120)
        pheatmap(pearson_cor,
                 display_numbers = TRUE, number_format = "%.2f",
                 fontsize_number = 10,
                 main = glue("LFC Pearson Correlation ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
                 color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
                 breaks = seq(-0.5, 1, length.out = 101),
                 annotation_col = col_annotation_display,
                 annotation_colors = ann_colors,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 gaps_col = c(n_upper, n_upper + n_datasets),
                 gaps_row = c(n_upper, n_upper + n_datasets),
                 fontsize = 14,
                 fontsize_row = 13, fontsize_col = 13)
        dev.off()
        log_info("  Saved Expanded_LFC_Pearson_Matrix.png")

        # --- Spearman correlation ---
        log_info("  Computing Spearman correlation matrix...")
        spearman_cor <- cor(expanded_lfc_mat, use = "pairwise.complete.obs",
                             method = "spearman")
        spearman_cor[is.na(spearman_cor)] <- 0

        # Save CSV
        write.csv(spearman_cor,
                  file.path(diag_dir, "Expanded_LFC_Spearman_Matrix.csv"),
                  quote = TRUE)
        log_info("  Saved Expanded_LFC_Spearman_Matrix.csv")

        # Prettify axis labels for display
        rownames(spearman_cor) <- prettify_label(rownames(spearman_cor))
        colnames(spearman_cor) <- prettify_label(colnames(spearman_cor))

        # Plot Spearman
        png(file.path(diag_dir, "Expanded_LFC_Spearman_Matrix.png"),
            width = 1800, height = 1400, res = 120)
        pheatmap(spearman_cor,
                 display_numbers = TRUE, number_format = "%.2f",
                 fontsize_number = 10,
                 main = glue("LFC Spearman Correlation ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
                 color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
                 breaks = seq(-0.5, 1, length.out = 101),
                 annotation_col = col_annotation_display,
                 annotation_colors = ann_colors,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 gaps_col = c(n_upper, n_upper + n_datasets),
                 gaps_row = c(n_upper, n_upper + n_datasets),
                 fontsize = 14,
                 fontsize_row = 13, fontsize_col = 13)
        dev.off()
        log_info("  Saved Expanded_LFC_Spearman_Matrix.png")

        # --- Summary statistics ---
        log_info("  === Expanded Matrix Summary ===")

        # Upper bound: within-run half correlations
        for (run_id in seq_len(N_UPPER_BOUND_RUNS)) {
          ha <- sprintf("Roussos_Half_A_R%d", run_id)
          hb <- sprintf("Roussos_Half_B_R%d", run_id)
          if (ha %in% colnames(pearson_cor) && hb %in% colnames(pearson_cor)) {
            log_info("    Run %d split-half Pearson:  %.3f", run_id, pearson_cor[ha, hb])
            log_info("    Run %d split-half Spearman: %.3f", run_id, spearman_cor[ha, hb])
          }
        }

        # Cross-dataset range
        ds_names <- names(all_results)
        ds_pearson <- pearson_cor[ds_names, ds_names]
        ds_off_diag <- ds_pearson[lower.tri(ds_pearson)]
        log_info("    Cross-dataset Pearson range:  [%.3f, %.3f] (median %.3f)",
                 min(ds_off_diag), max(ds_off_diag), median(ds_off_diag))

        ds_spearman <- spearman_cor[ds_names, ds_names]
        ds_off_diag_s <- ds_spearman[lower.tri(ds_spearman)]
        log_info("    Cross-dataset Spearman range: [%.3f, %.3f] (median %.3f)",
                 min(ds_off_diag_s), max(ds_off_diag_s), median(ds_off_diag_s))

        # Null range
        null_names <- names(null_results)
        if (length(null_names) > 0 && length(ds_names) > 0) {
          null_vs_ds <- pearson_cor[null_names, ds_names, drop = FALSE]
          log_info("    Null vs datasets Pearson range: [%.3f, %.3f] (median %.3f)",
                   min(null_vs_ds), max(null_vs_ds), median(null_vs_ds))
        }

      } else {
        log_warn("  No significant genes in reference for expanded matrix")
      }
    } else {
      log_warn("  Skipping expanded matrix: missing reference or bound data")
    }
  }, error = function(e) {
    log_error("  Expanded LFC agreement matrix failed: %s", e$message)
    try(dev.off(), silent = TRUE)
  })

  # ==========================================================================
  # DIAGNOSTIC 16: False Discovery Proportion (FDP) vs Roussos Full GT
  # ==========================================================================
  # For each query source (cross-dataset, split-half, null), compute:
  #   FDP = |Q_sig ∩ Roussos_insig| / |Q_sig ∩ Roussos_universe|
  # where Roussos_insig = genes Roussos tested but deemed insignificant.
  # ==========================================================================
  log_info("[16/16] False Discovery Proportion (FDP) vs Roussos Full GT...")

  tryCatch({
    if (!is.na(ref_name) && ref_name %in% names(all_results)) {

      ref_df <- all_results[[ref_name]]
      roussos_universe <- rownames(ref_df)
      roussos_sig <- rownames(ref_df[!is.na(ref_df$padj) &
                                       ref_df$padj < FDR_THRESHOLD &
                                       abs(ref_df$log2FoldChange) > LFC_THRESHOLD, ])
      roussos_insig <- setdiff(roussos_universe, roussos_sig)

      log_info("  Roussos universe: %d genes, sig: %d, insig: %d",
               length(roussos_universe), length(roussos_sig), length(roussos_insig))

      # Helper: compute FDP for a single query DE result data.frame
      compute_fdp <- function(query_df) {
        q_sig <- rownames(query_df[!is.na(query_df$padj) &
                                     query_df$padj < FDR_THRESHOLD &
                                     abs(query_df$log2FoldChange) > LFC_THRESHOLD, ])
        q_sig_in_universe <- intersect(q_sig, roussos_universe)
        q_sig_disagreed   <- intersect(q_sig_in_universe, roussos_insig)

        fdp <- if (length(q_sig_in_universe) > 0) {
          length(q_sig_disagreed) / length(q_sig_in_universe)
        } else {
          NA_real_
        }

        list(
          fdp           = fdp,
          n_sig         = length(q_sig),
          n_in_universe = length(q_sig_in_universe),
          n_disagreed   = length(q_sig_disagreed)
        )
      }

      # Assemble all query sources with group labels
      query_sources <- list()
      query_groups  <- character(0)

      # Cross-dataset sources (exclude Roussos itself)
      for (ds in setdiff(names(all_results), ref_name)) {
        query_sources[[ds]] <- all_results[[ds]]
        query_groups <- c(query_groups, "Cross-Dataset")
      }

      # Split-half sources (upper bound)
      for (src_name in names(upper_bound_results)) {
        query_sources[[src_name]] <- upper_bound_results[[src_name]]
        query_groups <- c(query_groups, "Split-Half")
      }

      # Null sources (lower bound)
      for (src_name in names(null_results)) {
        query_sources[[src_name]] <- null_results[[src_name]]
        query_groups <- c(query_groups, "Null")
      }

      # Compute FDP for each source
      fdp_rows <- list()
      for (i in seq_along(query_sources)) {
        src_name <- names(query_sources)[i]
        res <- compute_fdp(query_sources[[src_name]])

        fdp_rows[[i]] <- data.frame(
          CellType      = cell_type,
          Source        = src_name,
          Group         = query_groups[i],
          FDP           = res$fdp,
          N_Sig         = res$n_sig,
          N_In_Universe = res$n_in_universe,
          N_Disagreed   = res$n_disagreed,
          stringsAsFactors = FALSE
        )

        log_info("    %s [%s]: FDP = %.3f (%d/%d disagreed, %d total sig)",
                 src_name, query_groups[i],
                 ifelse(is.na(res$fdp), NA_real_, res$fdp),
                 res$n_disagreed, res$n_in_universe, res$n_sig)
      }

      ct_fdp <- do.call(rbind, fdp_rows)

      # Save per-cell-type CSV
      fdp_csv <- file.path(diag_dir, "FDP_per_source.csv")
      write.csv(ct_fdp, fdp_csv, row.names = FALSE, quote = TRUE)
      log_info("  Saved FDP_per_source.csv (%d sources)", nrow(ct_fdp))

      # Collect for post-loop heatmap
      fdp_collector[[cell_type]] <- ct_fdp

    } else {
      log_warn("  Reference dataset '%s' not found. Skipping FDP.", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  FDP computation failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 17: Intersection over Union (IoU) Analysis
  # ==========================================================================
  log_info("[17/17] Intersection over Union (IoU) Analysis...")

  tryCatch({
    if (!is.na(ref_name) && ref_name %in% names(sig_list)) {
      
      # Build IoU matrix using all available result sets
      # Combine: upper bound + datasets + null
      iou_input_list <- c(upper_bound_results, all_results, null_results)
      
      if (length(iou_input_list) >= 2) {
        log_info("  Computing IoU matrix for %d result sets...", length(iou_input_list))
        
        # Compute IoU
        iou_res <- build_iou_matrix(
          combined_results = iou_input_list,
          ref_name = ref_name,
          fdr_threshold = FDR_THRESHOLD,
          lfc_threshold = LFC_THRESHOLD
        )
        
        if (!is.null(iou_res$iou_matrix)) {
          # Add group annotations
          upper_names <- names(upper_bound_results)
          null_names_vec <- names(null_results)
          
          iou_res$summary_df <- annotate_iou_groups(
            iou_res$summary_df,
            upper_bound_names = upper_names,
            null_names = null_names_vec
          )
          
          # Save IoU matrix CSV
          iou_csv <- file.path(diag_dir, "Expanded_IoU_Matrix.csv")
          write.csv(iou_res$iou_matrix, iou_csv, quote = TRUE)
          log_info("  Saved Expanded_IoU_Matrix.csv")
          
          # Save detailed summary CSV
          summary_csv <- file.path(diag_dir, "IoU_Summary_Statistics.csv")
          write.csv(iou_res$summary_df, summary_csv, row.names = FALSE, quote = TRUE)
          log_info("  Saved IoU_Summary_Statistics.csv")
          
          # Create column annotations for heatmap
          n_upper <- length(upper_bound_results)
          n_datasets <- length(all_results)
          n_null <- length(null_results)
          
          group_labels <- c(
            rep("Split-Half", n_upper),
            rep("Cross-Dataset", n_datasets),
            rep("Null", n_null)
          )

          col_annotation <- data.frame(
            Group = factor(group_labels,
                           levels = c("Split-Half",
                                      "Cross-Dataset",
                                      "Null")),
            row.names = colnames(iou_res$iou_matrix)
          )

          ann_colors <- list(
            Group = c(
              "Split-Half"    = "#2166AC",
              "Cross-Dataset" = "#4DAF4A",
              "Null"          = "#E41A1C"
            )
          )
          
          # Prettify row and column labels for display
          iou_mat_display <- iou_res$iou_matrix
          rownames(iou_mat_display) <- prettify_label(rownames(iou_mat_display))
          colnames(iou_mat_display) <- prettify_label(colnames(iou_mat_display))
          
          # Prettify annotation row names too
          col_annotation_display <- col_annotation
          rownames(col_annotation_display) <- prettify_label(rownames(col_annotation_display))
          
          # Plot IoU heatmap
          png(file.path(diag_dir, "Expanded_IoU_Matrix.png"),
              width = 2200, height = 1800, res = 120)
          pheatmap(iou_mat_display,
                   display_numbers = TRUE, number_format = "%.2f",
                   fontsize_number = 14,
                   number_color = "black",
                   main = glue("IoU Without Universe Control ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
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
          log_info("  Saved Expanded_IoU_Matrix.png")
          
          # Log summary statistics
          null_vs_real <- iou_res$summary_df[iou_res$summary_df$Group == "Lower Bound" &
                                               iou_res$summary_df$Set1 == ref_name, ]
          if (nrow(null_vs_real) > 0) {
            log_info("  IoU (null vs %s): range [%.3f, %.3f], median %.3f",
                     ref_name,
                     min(null_vs_real$IoU, na.rm = TRUE),
                     max(null_vs_real$IoU, na.rm = TRUE),
                     median(null_vs_real$IoU, na.rm = TRUE))
          }
          
          real_vs_real <- iou_res$summary_df[iou_res$summary_df$Group == "Cross-Dataset" &
                                               iou_res$summary_df$Set1 == ref_name, ]
          if (nrow(real_vs_real) > 0) {
            log_info("  IoU (real datasets vs %s): range [%.3f, %.3f], median %.3f",
                     ref_name,
                     min(real_vs_real$IoU, na.rm = TRUE),
                     max(real_vs_real$IoU, na.rm = TRUE),
                     median(real_vs_real$IoU, na.rm = TRUE))
          }
        }
      } else {
        log_warn("  Insufficient result sets for IoU computation (need >= 2)")
      }
    } else {
      log_warn("  Reference dataset '%s' not found. Skipping IoU analysis.", REFERENCE_DATASET)
    }
  }, error = function(e) {
    log_error("  IoU analysis failed: %s", e$message)
  })

  # ==========================================================================
  # DIAGNOSTIC 18: Universe-Controlled IoU Analysis
  # ==========================================================================
  log_info("[18/18] Universe-Controlled IoU Analysis...")

  tryCatch({
    if (length(iou_input_list) >= 2) {
      log_info("  Computing universe-controlled IoU matrix for %d result sets...", length(iou_input_list))
      
      # Compute universe-controlled IoU
      iou_uc_res <- build_iou_matrix_universe_controlled(
        combined_results = iou_input_list,
        fdr_threshold = FDR_THRESHOLD,
        lfc_threshold = LFC_THRESHOLD
      )
      
      if (!is.null(iou_uc_res$iou_matrix)) {
        # Add group annotations (reusing same grouping logic)
        iou_uc_res$summary_df <- annotate_iou_groups(
          iou_uc_res$summary_df,
          upper_bound_names = upper_names,
          null_names = null_names_vec
        )
        
        # Save universe-controlled IoU matrix CSV
        iou_uc_csv <- file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.csv")
        write.csv(iou_uc_res$iou_matrix, iou_uc_csv, quote = TRUE)
        log_info("  Saved Expanded_IoU_Matrix_Universe_Controlled.csv")
        
        # Save detailed summary CSV
        summary_uc_csv <- file.path(diag_dir, "IoU_Summary_Statistics_Universe_Controlled.csv")
        write.csv(iou_uc_res$summary_df, summary_uc_csv, row.names = FALSE, quote = TRUE)
        log_info("  Saved IoU_Summary_Statistics_Universe_Controlled.csv")
        
        # Create column annotations for heatmap (same as standard IoU)
        col_annotation_uc <- col_annotation  # Reuse from standard IoU
        
        # Prettify row and column labels for display
        iou_uc_mat_display <- iou_uc_res$iou_matrix
        rownames(iou_uc_mat_display) <- prettify_label(rownames(iou_uc_mat_display))
        colnames(iou_uc_mat_display) <- prettify_label(colnames(iou_uc_mat_display))
        
        # Prettify annotation row names too
        col_annotation_uc_display <- col_annotation_uc
        rownames(col_annotation_uc_display) <- prettify_label(rownames(col_annotation_uc_display))
        
        # Plot universe-controlled IoU heatmap
        png(file.path(diag_dir, "Expanded_IoU_Matrix_Universe_Controlled.png"),
            width = 2200, height = 1800, res = 120)
        pheatmap(iou_uc_mat_display,
                 display_numbers = TRUE, number_format = "%.2f",
                 fontsize_number = 13,
                 number_color = "black",
                 main = glue("Universe-Controlled IoU ({prettify_label(REFERENCE_DATASET)} Significant Genes): {prettify_label(cell_type)}"),
                 color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
                 breaks = seq(0, 1, length.out = 101),
                 annotation_col = col_annotation_uc_display,
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
        log_info("  Saved Expanded_IoU_Matrix_Universe_Controlled.png")
        
        # Log summary statistics
        null_vs_real_uc <- iou_uc_res$summary_df[iou_uc_res$summary_df$Group == "Lower Bound" &
                                                   iou_uc_res$summary_df$Set1 == ref_name, ]
        if (nrow(null_vs_real_uc) > 0) {
          log_info("  Universe-controlled IoU (null vs %s): range [%.3f, %.3f], median %.3f",
                   ref_name,
                   min(null_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE),
                   max(null_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE),
                   median(null_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE))
        }
        
        real_vs_real_uc <- iou_uc_res$summary_df[iou_uc_res$summary_df$Group == "Cross-Dataset" &
                                                   iou_uc_res$summary_df$Set1 == ref_name, ]
        if (nrow(real_vs_real_uc) > 0) {
          log_info("  Universe-controlled IoU (real datasets vs %s): range [%.3f, %.3f], median %.3f",
                   ref_name,
                   min(real_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE),
                   max(real_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE),
                   median(real_vs_real_uc$IoU_Universe_Controlled, na.rm = TRUE))
        }
      }
    } else {
      log_warn("  Insufficient result sets for universe-controlled IoU computation")
    }
  }, error = function(e) {
    log_error("  Universe-controlled IoU analysis failed: %s", e$message)
  })

  # --------------------------------------------------------------------------
  # Cell type complete
  # --------------------------------------------------------------------------
  n_diag_total <- 18
  completion_summary[[cell_type]] <- sprintf("OK (%d datasets, %d diagnostics, %d upper bound, %d null)",
                                              length(dataset_names), n_diag_total,
                                              length(upper_bound_results),
                                              length(null_results))
  log_info("%s diagnostics complete. Output: %s", cell_type, diag_dir)
}

# ==============================================================================
# 3. POST-LOOP: FDP Heatmap (all cell types combined)
# ==============================================================================

log_section("FDP Heatmap: Combining across cell types")

tryCatch({
  if (length(fdp_collector) > 0) {

    combined_fdp <- do.call(rbind, fdp_collector)

    # Save combined CSV
    combined_csv <- file.path(OUTPUT_BASE, "FDP_Summary.csv")
    write.csv(combined_fdp, combined_csv, row.names = FALSE, quote = TRUE)
    log_info("Saved combined FDP summary: %s", combined_csv)

    # Build wide matrix for heatmap: rows = cell types, columns = sources
    # Pivot to wide format
    fdp_wide <- reshape(
      combined_fdp[, c("CellType", "Source", "FDP")],
      idvar     = "CellType",
      timevar   = "Source",
      direction = "wide"
    )
    rownames(fdp_wide) <- fdp_wide$CellType
    fdp_wide$CellType <- NULL
    colnames(fdp_wide) <- gsub("^FDP\\.", "", colnames(fdp_wide))

    # Ensure consistent column ordering: Cross-Dataset, Split-Half, Null
    # Get group assignments from the first cell type's data
    first_ct <- fdp_collector[[1]]
    source_order <- first_ct$Source
    group_order  <- first_ct$Group

    # Reorder columns to match group ordering
    available_cols <- intersect(source_order, colnames(fdp_wide))
    fdp_mat <- as.matrix(fdp_wide[, available_cols, drop = FALSE])

    # Row ordering: Astro, Micro, Oligo
    row_order <- intersect(CELL_TYPES, rownames(fdp_mat))
    fdp_mat <- fdp_mat[row_order, , drop = FALSE]

    # Compute column group assignments in original source order
    col_groups <- group_order[match(available_cols, source_order)]

    # Reorder columns so groups appear in the same order as the IoU panels:
    # Split-Half, Cross-Dataset, Null. Within a group, preserve source order.
    group_priority <- setNames(seq_along(c("Split-Half", "Cross-Dataset", "Null")),
                               c("Split-Half", "Cross-Dataset", "Null"))
    col_order_idx <- order(group_priority[col_groups], seq_along(col_groups))
    fdp_mat <- fdp_mat[, col_order_idx, drop = FALSE]
    available_cols <- available_cols[col_order_idx]
    col_groups <- col_groups[col_order_idx]

    # Abbreviate row labels to single letters (A, M, O) for compact display
    pretty_rows <- prettify_label(rownames(fdp_mat))
    rownames(fdp_mat) <- substr(pretty_rows, 1, 1)
    colnames(fdp_mat) <- prettify_label(colnames(fdp_mat))

    # Column annotation for pheatmap (factor level order matches reordered columns)
    col_annotation <- data.frame(
      Group = factor(col_groups,
                     levels = c("Split-Half", "Cross-Dataset", "Null")),
      row.names = colnames(fdp_mat)
    )

    ann_colors <- list(
      Group = c(
        "Split-Half"    = "#2166AC",
        "Cross-Dataset" = "#4DAF4A",
        "Null"          = "#E41A1C"
      )
    )

    # Compute gap positions between groups (after reordering)
    group_runs <- rle(as.character(col_groups))
    gap_positions <- cumsum(group_runs$lengths[-length(group_runs$lengths)])

    # Generate heatmap
    png(file.path(OUTPUT_BASE, "FDP_Heatmap.png"),
        width = 1800, height = 700, res = 150)
    pheatmap(fdp_mat,
             display_numbers = TRUE,
             number_format   = "%.2f",
             fontsize_number = 13,
             main = glue("False Discovery Proportion vs {prettify_label(REFERENCE_DATASET)} Full Ground Truth"),
             color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#FC8D59", "#D73027"))(100),
             breaks = seq(0, 1, length.out = 101),
             annotation_col  = col_annotation,
             annotation_colors = ann_colors,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             gaps_col = gap_positions,
             angle_col = 45,
             fontsize = 15,
             fontsize_row = 15, fontsize_col = 14,
             annotation_legend = TRUE,
             annotation_names_col = FALSE)
    dev.off()
    log_info("Saved FDP_Heatmap.png")

    # Summary statistics
    log_info("=== FDP Summary ===")
    for (grp in c("Cross-Dataset", "Split-Half", "Null")) {
      grp_data <- combined_fdp[combined_fdp$Group == grp, ]
      if (nrow(grp_data) > 0) {
        log_info("  %s: mean FDP = %.3f, range = [%.3f, %.3f]",
                 grp,
                 mean(grp_data$FDP, na.rm = TRUE),
                 min(grp_data$FDP, na.rm = TRUE),
                 max(grp_data$FDP, na.rm = TRUE))
      }
    }

  } else {
    log_warn("No FDP data collected. Skipping heatmap.")
  }
}, error = function(e) {
  log_error("FDP heatmap generation failed: %s", e$message)
  try(dev.off(), silent = TRUE)
})

# ==============================================================================
# 4. COMPLETION SUMMARY
# ==============================================================================

log_section("Step 2 Complete: Cross-Dataset Diagnostics")

for (ct in names(completion_summary)) {
  log_info("  %s: %s", ct, completion_summary[[ct]])
}

log_info("All diagnostic results saved under: %s/{CellType}/Diagnostic_Results/",
         OUTPUT_BASE)
