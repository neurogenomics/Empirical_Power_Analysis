# ==============================================================================
# SCRIPT: Aggregate Downsampling Results
# PURPOSE: Combine results across runs and compute precision/recall metrics
# INPUT: DESeq2 results from Step 2
# OUTPUT: Aggregated summary data for visualization
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(parallel)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",  # From project root
    "config.R",                            # From script directory
    "../config.R"                          # From tests/ or toy/
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

log_section("STEP 3: Aggregate Results")

# Mode-aware configuration
file_patterns <- get_file_pattern()
y_col <- get_y_col_name()
log_info("Mode: %s | Y-column: %s | File pattern: %s", DOWNSAMPLE_MODE, y_col, file_patterns$pattern)

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {

  for (method in DE_METHODS) {
  
  log_section(sprintf("Cell Type: %s | Method: %s", cell_type, method))
  
  # Get paths
  paths <- get_paths(cell_type, method = method)
  
  # Create plots directory
  if (!dir.exists(paths$plots_dir)) {
    dir.create(paths$plots_dir, recursive = TRUE)
    log_info("Created plots directory: %s", paths$plots_dir)
  }
  
  # Checkpoint file path
  checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
  
  # ==============================================================================
  # 3. CHECK FOR EXISTING CHECKPOINT
  # ==============================================================================
  
  if (file.exists(checkpoint_file) && !FORCE_RECALC) {
    log_info("Found checkpoint file, loading pre-computed summaries...")
    checkpoint_data <- qs::qread(checkpoint_file)
    log_info("Loaded summaries with %d rows", nrow(checkpoint_data$df_summary))
    next
  }
  
  # ==============================================================================
  # 4. PROCESS ALL RUNS IN PARALLEL
  # ==============================================================================
  
  log_info("No checkpoint found (or recalc forced). Running aggregation...")
  
  # Find run directories
  run_folders <- list.dirs(paths$results_dir, recursive = FALSE, full.names = TRUE)
  run_folders <- run_folders[grepl("Run_", basename(run_folders))]
  
  if (length(run_folders) == 0) {
    log_warning("No run folders found for %s", cell_type)
    next
  }
  
  log_info("Found %d run folders", length(run_folders))
  
  # Process runs in parallel
  parallel_results <- mclapply(run_folders, function(rf) {
    
    run_id <- basename(rf)
    local_results <- list()
    local_xist <- list()
    local_reads <- list()
    local_bm_results <- list()
    local_cart_results <- list()
    local_exprfilt_results <- list()
    local_finecross_results <- list()
    local_pval_q_results <- list()

    # Load ground truth
    gt_path <- file.path(rf, "Ground_Truth.qs")
    if (!file.exists(gt_path)) {
      return(NULL)
    }
    
    gt_res <- as.data.frame(qs::qread(gt_path))
    all_genes_in_run <- rownames(gt_res)

    # Canonical expression ruler: shared raw pseudobulk mean (baseMean_raw), so the
    # expression axis is identical across DESeq2/edgeR/limma. Fall back to the
    # normalised baseMean only for pre-rerun outputs that predate baseMean_raw.
    gt_res$bm_ruler <- if ("baseMean_raw" %in% colnames(gt_res)) gt_res$baseMean_raw else gt_res$baseMean

    # Genes with low ground-truth expression (baseMean < 50) for filtered DEG counts
    genes_lowexpr <- rownames(gt_res)[!is.na(gt_res$bm_ruler) & gt_res$bm_ruler < 50]
    
    # Get chromosome mapping
    chr_map <- get_chr_map(all_genes_in_run)
    if (!is.null(chr_map)) {
      genes_autosomal <- names(chr_map)[chr_map %in% as.character(1:22)]
      genes_sex <- names(chr_map)[chr_map %in% c("X", "Y")]
    } else {
      genes_autosomal <- character(0)
      genes_sex <- character(0)
    }
    
    # Ground truth significant genes
    gt_all_sig <- rownames(gt_res)[!is.na(gt_res$padj) & gt_res$padj < FDR_THRESHOLD]

    # P-value quartile boundaries: computed once per run on GT-sig genes only.
    # Q1 = most significant (smallest padj), Q4 = least significant (still passing).
    # Each run has slightly different GT donors, so break values vary across runs;
    # the bin labels (Q1–Q4) are stable and semantically consistent.
    # Use logical indexing: gt_all_sig is a character vector so gt_res$padj[gt_all_sig]
    # would return NAs (unnamed numeric vector indexed by names). Logical mask is correct.
    gt_sig_padj_vals <- gt_res$padj[!is.na(gt_res$padj) & gt_res$padj < FDR_THRESHOLD &
                                    !is.na(gt_res$log2FoldChange) & abs(gt_res$log2FoldChange) >= 0.2]
    # Inner quartile breaks only (0.25/0.5/0.75); Q1's lower bound is the abstract
    # 0 and Q4's upper bound is FDR_THRESHOLD. Membership is identical to using the
    # true min/max extremes (all sig genes lie within [0, FDR) by definition).
    pval_q_intervals <- if (length(gt_sig_padj_vals) >= 4) {
      pq_brks <- quantile(gt_sig_padj_vals, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      list(
        "Q1 (most sig)"  = c(0,             pq_brks[[1]]),
        "Q2"             = c(pq_brks[[1]],  pq_brks[[2]]),
        "Q3"             = c(pq_brks[[2]],  pq_brks[[3]]),
        "Q4 (least sig)" = c(pq_brks[[3]],  FDR_THRESHOLD)
      )
    } else {
      NULL
    }

    # Find subsample files (mode-aware pattern)
    sub_files <- list.files(rf, pattern = file_patterns$pattern, full.names = TRUE)
    
    for (sf in sub_files) {
      
      fname <- basename(sf)
      meta <- str_match(fname, file_patterns$parse_regex)
      
      if (is.na(meta[1])) next
      
      n_donors <- as.numeric(meta[2])
      y_value <- meta[3]  # Cells or Target_Reads depending on mode
      iter <- as.numeric(meta[4])
      
      # Load subsample results
      sub_res <- tryCatch({
        as.data.frame(qs::qread(sf))
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(sub_res)) next

      # Same shared raw ruler on the subsample table (see gt_res above).
      sub_res$bm_ruler <- if ("baseMean_raw" %in% colnames(sub_res)) sub_res$baseMean_raw else sub_res$baseMean

      sub_all_genes <- rownames(sub_res)
      sub_all_sig <- rownames(sub_res)[!is.na(sub_res$padj) & sub_res$padj < FDR_THRESHOLD]
      
      # Total reads (baseMean sum)
      total_reads_val <- sum(sub_res$bm_ruler, na.rm = TRUE)
      reads_row <- data.frame(
        Donors = n_donors,
        Y_Value = y_value,
        Iteration = iter,
        Total_Reads = total_reads_val,
        stringsAsFactors = FALSE
      )
      names(reads_row)[names(reads_row) == "Y_Value"] <- y_col
      local_reads[[length(local_reads) + 1]] <- reads_row
      
      # Process each LFC bin
      for (bin_name in names(LFC_INTERVALS)) {
        
        bounds <- LFC_INTERVALS[[bin_name]]
        
        # Ground truth genes in this bin
        gt_in_bin <- rownames(gt_res)[
          !is.na(gt_res$padj) & 
          gt_res$padj < FDR_THRESHOLD & 
          abs(gt_res$log2FoldChange) >= bounds[1] & 
          abs(gt_res$log2FoldChange) < bounds[2]
        ]
        
        # Subsample genes in this bin
        sub_in_bin <- rownames(sub_res)[
          !is.na(sub_res$padj) & 
          sub_res$padj < FDR_THRESHOLD & 
          abs(sub_res$log2FoldChange) >= bounds[1] & 
          abs(sub_res$log2FoldChange) < bounds[2]
        ]
        
        # Process each chromosome category
        for (cat in CHR_CATEGORIES) {
          
          if (cat == "Total") {
            gt_target <- gt_in_bin
            sub_target <- sub_in_bin
            gt_sig_ref <- gt_all_sig
            sub_sig_ref <- sub_all_sig
          } else if (cat == "Autosomal") {
            gt_target <- intersect(gt_in_bin, genes_autosomal)
            sub_target <- intersect(sub_in_bin, genes_autosomal)
            gt_sig_ref <- intersect(gt_all_sig, genes_autosomal)
            sub_sig_ref <- intersect(sub_all_sig, genes_autosomal)
          } else if (cat == "SexChr") {
            gt_target <- intersect(gt_in_bin, genes_sex)
            sub_target <- intersect(sub_in_bin, genes_sex)
            gt_sig_ref <- intersect(gt_all_sig, genes_sex)
            sub_sig_ref <- intersect(sub_all_sig, genes_sex)
          }
          
          # Raw and harmonized counts
          n_gt_raw <- length(gt_target)
          gt_target_harm <- intersect(gt_target, sub_all_genes)
          n_gt_harm <- length(gt_target_harm)
          
          # Precision: TP / (TP + FP) - of what we called significant, how many are in GT?
          tp_prec <- length(intersect(sub_target, gt_sig_ref))
          fp_prec <- length(setdiff(sub_target, gt_sig_ref))
          precision <- if ((tp_prec + fp_prec) > 0) tp_prec / (tp_prec + fp_prec) else NA_real_
          
          # Recall: TP / (TP + FN) - of GT significant, how many did we find?
          tp_rec <- length(intersect(gt_target, sub_sig_ref))
          fn_rec <- length(setdiff(gt_target, sub_sig_ref))
          recall <- if ((tp_rec + fn_rec) > 0) tp_rec / (tp_rec + fn_rec) else 0
          
          # Harmonized recall (only genes present in both)
          tp_rec_harm <- length(intersect(gt_target_harm, sub_sig_ref))
          fn_rec_harm <- length(setdiff(gt_target_harm, sub_sig_ref))
          recall_harm <- if ((tp_rec_harm + fn_rec_harm) > 0) tp_rec_harm / (tp_rec_harm + fn_rec_harm) else 0
          
          # Harmonized precision
          sub_target_harm <- intersect(sub_target, all_genes_in_run)
          tp_prec_harm <- length(intersect(sub_target_harm, gt_sig_ref))
          fp_prec_harm <- length(setdiff(sub_target_harm, gt_sig_ref))
          precision_harm <- if ((tp_prec_harm + fp_prec_harm) > 0) tp_prec_harm / (tp_prec_harm + fp_prec_harm) else NA_real_
          
          res_row <- data.frame(
            Run = run_id,
            Donors = n_donors,
            Y_Value = y_value,
            Iteration = iter,
            LFC_Bin = bin_name,
            Gene_Category = cat,
            Precision = precision,
            Recall = recall,
            Precision_Harm = precision_harm,
            Recall_Harm = recall_harm,
            GT_Count = n_gt_raw,
            GT_Harm_Count = n_gt_harm,
            Sub_Count = length(sub_target),
            Sub_Harm_Count = length(sub_target_harm),
            Sub_Count_LowExpr = length(intersect(sub_target, genes_lowexpr)),
            stringsAsFactors = FALSE
          )
          names(res_row)[names(res_row) == "Y_Value"] <- y_col
          local_results[[length(local_results) + 1]] <- res_row
        }
        
        # XIST detection check
        xist_detected <- 0
        if (GENE_OF_INTEREST %in% rownames(sub_res)) {
          x_padj <- sub_res[GENE_OF_INTEREST, "padj"]
          x_lfc <- abs(sub_res[GENE_OF_INTEREST, "log2FoldChange"])
          if (!is.na(x_padj) && x_padj < FDR_THRESHOLD && 
              x_lfc >= bounds[1] && x_lfc < bounds[2]) {
            xist_detected <- 1
          }
        }
        
        xist_row <- data.frame(
          Run = run_id,
          Donors = n_donors,
          Y_Value = y_value,
          Iteration = iter,
          LFC_Bin = bin_name,
          Detected = xist_detected,
          stringsAsFactors = FALSE
        )
        names(xist_row)[names(xist_row) == "Y_Value"] <- y_col
        local_xist[[length(local_xist) + 1]] <- xist_row
      }

      # --- BaseMean bins (expression-level stratification) ---
      for (bm_bin_name in names(BASEMEAN_INTERVALS)) {
        bm_bounds <- BASEMEAN_INTERVALS[[bm_bin_name]]

        gt_in_bm <- rownames(gt_res)[
          !is.na(gt_res$padj) &
          gt_res$padj < FDR_THRESHOLD &
          gt_res$bm_ruler >= bm_bounds[1] &
          gt_res$bm_ruler < bm_bounds[2]
        ]

        sub_in_bm <- rownames(sub_res)[
          !is.na(sub_res$padj) &
          sub_res$padj < FDR_THRESHOLD &
          sub_res$bm_ruler >= bm_bounds[1] &
          sub_res$bm_ruler < bm_bounds[2]
        ]

        for (cat in CHR_CATEGORIES) {
          if (cat == "Total") {
            gt_target <- gt_in_bm; sub_target <- sub_in_bm
            gt_sig_ref <- gt_all_sig; sub_sig_ref <- sub_all_sig
          } else if (cat == "Autosomal") {
            gt_target <- intersect(gt_in_bm, genes_autosomal)
            sub_target <- intersect(sub_in_bm, genes_autosomal)
            gt_sig_ref <- intersect(gt_all_sig, genes_autosomal)
            sub_sig_ref <- intersect(sub_all_sig, genes_autosomal)
          } else if (cat == "SexChr") {
            gt_target <- intersect(gt_in_bm, genes_sex)
            sub_target <- intersect(sub_in_bm, genes_sex)
            gt_sig_ref <- intersect(gt_all_sig, genes_sex)
            sub_sig_ref <- intersect(sub_all_sig, genes_sex)
          }

          n_gt_raw <- length(gt_target)
          gt_target_harm <- intersect(gt_target, sub_all_genes)
          n_gt_harm <- length(gt_target_harm)

          tp_prec <- length(intersect(sub_target, gt_sig_ref))
          fp_prec <- length(setdiff(sub_target, gt_sig_ref))
          precision <- if ((tp_prec + fp_prec) > 0) tp_prec / (tp_prec + fp_prec) else NA_real_

          tp_rec <- length(intersect(gt_target, sub_sig_ref))
          fn_rec <- length(setdiff(gt_target, sub_sig_ref))
          recall <- if ((tp_rec + fn_rec) > 0) tp_rec / (tp_rec + fn_rec) else 0

          tp_rec_harm <- length(intersect(gt_target_harm, sub_sig_ref))
          fn_rec_harm <- length(setdiff(gt_target_harm, sub_sig_ref))
          recall_harm <- if ((tp_rec_harm + fn_rec_harm) > 0) tp_rec_harm / (tp_rec_harm + fn_rec_harm) else 0

          sub_target_harm <- intersect(sub_target, all_genes_in_run)
          tp_prec_harm <- length(intersect(sub_target_harm, gt_sig_ref))
          fp_prec_harm <- length(setdiff(sub_target_harm, gt_sig_ref))
          precision_harm <- if ((tp_prec_harm + fp_prec_harm) > 0) tp_prec_harm / (tp_prec_harm + fp_prec_harm) else NA_real_

          bm_row <- data.frame(
            Run = run_id, Donors = n_donors, Y_Value = y_value,
            Iteration = iter, BaseMean_Bin = bm_bin_name, Gene_Category = cat,
            Precision = precision, Recall = recall,
            Precision_Harm = precision_harm, Recall_Harm = recall_harm,
            GT_Count = n_gt_raw, GT_Harm_Count = n_gt_harm,
            Sub_Count = length(sub_target),
            Sub_Count_LowExpr = length(intersect(sub_target, genes_lowexpr)),
            stringsAsFactors = FALSE
          )
          names(bm_row)[names(bm_row) == "Y_Value"] <- y_col
          local_bm_results[[length(local_bm_results) + 1]] <- bm_row
        }
      }

      # --- Cartesian product: LFC_COARSE x BASEMEAN_COARSE ---
      for (lfc_c_name in names(LFC_COARSE)) {
        lfc_c_bounds <- LFC_COARSE[[lfc_c_name]]
        for (bm_c_name in names(BASEMEAN_COARSE)) {
          bm_c_bounds <- BASEMEAN_COARSE[[bm_c_name]]

          cart_label <- paste0(lfc_c_name, " / ", bm_c_name)

          gt_in_cart <- rownames(gt_res)[
            !is.na(gt_res$padj) &
            gt_res$padj < FDR_THRESHOLD &
            abs(gt_res$log2FoldChange) >= lfc_c_bounds[1] &
            abs(gt_res$log2FoldChange) < lfc_c_bounds[2] &
            gt_res$bm_ruler >= bm_c_bounds[1] &
            gt_res$bm_ruler < bm_c_bounds[2]
          ]

          sub_in_cart <- rownames(sub_res)[
            !is.na(sub_res$padj) &
            sub_res$padj < FDR_THRESHOLD &
            abs(sub_res$log2FoldChange) >= lfc_c_bounds[1] &
            abs(sub_res$log2FoldChange) < lfc_c_bounds[2] &
            sub_res$bm_ruler >= bm_c_bounds[1] &
            sub_res$bm_ruler < bm_c_bounds[2]
          ]

          for (cat in CHR_CATEGORIES) {
            if (cat == "Total") {
              gt_target <- gt_in_cart; sub_target <- sub_in_cart
              gt_sig_ref <- gt_all_sig; sub_sig_ref <- sub_all_sig
            } else if (cat == "Autosomal") {
              gt_target <- intersect(gt_in_cart, genes_autosomal)
              sub_target <- intersect(sub_in_cart, genes_autosomal)
              gt_sig_ref <- intersect(gt_all_sig, genes_autosomal)
              sub_sig_ref <- intersect(sub_all_sig, genes_autosomal)
            } else if (cat == "SexChr") {
              gt_target <- intersect(gt_in_cart, genes_sex)
              sub_target <- intersect(sub_in_cart, genes_sex)
              gt_sig_ref <- intersect(gt_all_sig, genes_sex)
              sub_sig_ref <- intersect(sub_all_sig, genes_sex)
            }

            n_gt_raw <- length(gt_target)
            gt_target_harm <- intersect(gt_target, sub_all_genes)
            n_gt_harm <- length(gt_target_harm)

            tp_prec <- length(intersect(sub_target, gt_sig_ref))
            fp_prec <- length(setdiff(sub_target, gt_sig_ref))
            precision <- if ((tp_prec + fp_prec) > 0) tp_prec / (tp_prec + fp_prec) else NA_real_

            tp_rec <- length(intersect(gt_target, sub_sig_ref))
            fn_rec <- length(setdiff(gt_target, sub_sig_ref))
            recall <- if ((tp_rec + fn_rec) > 0) tp_rec / (tp_rec + fn_rec) else 0

            tp_rec_harm <- length(intersect(gt_target_harm, sub_sig_ref))
            fn_rec_harm <- length(setdiff(gt_target_harm, sub_sig_ref))
            recall_harm <- if ((tp_rec_harm + fn_rec_harm) > 0) tp_rec_harm / (tp_rec_harm + fn_rec_harm) else 0

            sub_target_harm <- intersect(sub_target, all_genes_in_run)
            tp_prec_harm <- length(intersect(sub_target_harm, gt_sig_ref))
            fp_prec_harm <- length(setdiff(sub_target_harm, gt_sig_ref))
            precision_harm <- if ((tp_prec_harm + fp_prec_harm) > 0) tp_prec_harm / (tp_prec_harm + fp_prec_harm) else NA_real_

            cart_row <- data.frame(
              Run = run_id, Donors = n_donors, Y_Value = y_value,
              Iteration = iter, Cartesian_Bin = cart_label, Gene_Category = cat,
              Precision = precision, Recall = recall,
              Precision_Harm = precision_harm, Recall_Harm = recall_harm,
              GT_Count = n_gt_raw, GT_Harm_Count = n_gt_harm,
              Sub_Count = length(sub_target),
              Sub_Count_LowExpr = length(intersect(sub_target, genes_lowexpr)),
              stringsAsFactors = FALSE
            )
            names(cart_row)[names(cart_row) == "Y_Value"] <- y_col
            local_cart_results[[length(local_cart_results) + 1]] <- cart_row
          }
        }
      }

      # --- Expression-filtered LFC bins (baseMean > EXPR_FILTER_THRESHOLD) ---
      # For recall:    only GT genes with GT baseMean > threshold
      # For precision: only predicted genes with subsample baseMean > threshold
      for (bin_name in names(LFC_INTERVALS)) {
        bounds <- LFC_INTERVALS[[bin_name]]

        # GT genes in this LFC bin AND high expression in GT
        gt_in_bin_filt <- rownames(gt_res)[
          !is.na(gt_res$padj) &
          gt_res$padj < FDR_THRESHOLD &
          abs(gt_res$log2FoldChange) >= bounds[1] &
          abs(gt_res$log2FoldChange) < bounds[2] &
          gt_res$bm_ruler > EXPR_FILTER_THRESHOLD
        ]

        # Sub genes in this LFC bin AND high expression in subsample
        sub_in_bin_filt <- rownames(sub_res)[
          !is.na(sub_res$padj) &
          sub_res$padj < FDR_THRESHOLD &
          abs(sub_res$log2FoldChange) >= bounds[1] &
          abs(sub_res$log2FoldChange) < bounds[2] &
          sub_res$bm_ruler > EXPR_FILTER_THRESHOLD
        ]

        # Reference sets also filtered by expression
        gt_all_sig_filt <- rownames(gt_res)[
          !is.na(gt_res$padj) & gt_res$padj < FDR_THRESHOLD &
          gt_res$bm_ruler > EXPR_FILTER_THRESHOLD
        ]
        sub_all_sig_filt <- rownames(sub_res)[
          !is.na(sub_res$padj) & sub_res$padj < FDR_THRESHOLD &
          sub_res$bm_ruler > EXPR_FILTER_THRESHOLD
        ]

        for (cat in CHR_CATEGORIES) {
          if (cat == "Total") {
            gt_target <- gt_in_bin_filt; sub_target <- sub_in_bin_filt
            gt_sig_ref <- gt_all_sig_filt; sub_sig_ref <- sub_all_sig_filt
          } else if (cat == "Autosomal") {
            gt_target <- intersect(gt_in_bin_filt, genes_autosomal)
            sub_target <- intersect(sub_in_bin_filt, genes_autosomal)
            gt_sig_ref <- intersect(gt_all_sig_filt, genes_autosomal)
            sub_sig_ref <- intersect(sub_all_sig_filt, genes_autosomal)
          } else if (cat == "SexChr") {
            gt_target <- intersect(gt_in_bin_filt, genes_sex)
            sub_target <- intersect(sub_in_bin_filt, genes_sex)
            gt_sig_ref <- intersect(gt_all_sig_filt, genes_sex)
            sub_sig_ref <- intersect(sub_all_sig_filt, genes_sex)
          }

          n_gt_raw <- length(gt_target)
          gt_target_harm <- intersect(gt_target, sub_all_genes)
          n_gt_harm <- length(gt_target_harm)

          tp_prec <- length(intersect(sub_target, gt_sig_ref))
          fp_prec <- length(setdiff(sub_target, gt_sig_ref))
          precision <- if ((tp_prec + fp_prec) > 0) tp_prec / (tp_prec + fp_prec) else NA_real_

          tp_rec <- length(intersect(gt_target, sub_sig_ref))
          fn_rec <- length(setdiff(gt_target, sub_sig_ref))
          recall <- if ((tp_rec + fn_rec) > 0) tp_rec / (tp_rec + fn_rec) else 0

          tp_rec_harm <- length(intersect(gt_target_harm, sub_sig_ref))
          fn_rec_harm <- length(setdiff(gt_target_harm, sub_sig_ref))
          recall_harm <- if ((tp_rec_harm + fn_rec_harm) > 0) tp_rec_harm / (tp_rec_harm + fn_rec_harm) else 0

          sub_target_harm <- intersect(sub_target, all_genes_in_run)
          tp_prec_harm <- length(intersect(sub_target_harm, gt_sig_ref))
          fp_prec_harm <- length(setdiff(sub_target_harm, gt_sig_ref))
          precision_harm <- if ((tp_prec_harm + fp_prec_harm) > 0) tp_prec_harm / (tp_prec_harm + fp_prec_harm) else NA_real_

          ef_row <- data.frame(
            Run = run_id, Donors = n_donors, Y_Value = y_value,
            Iteration = iter, LFC_Bin = bin_name, Gene_Category = cat,
            Precision = precision, Recall = recall,
            Precision_Harm = precision_harm, Recall_Harm = recall_harm,
            GT_Count = n_gt_raw, GT_Harm_Count = n_gt_harm,
            Sub_Count = length(sub_target),
            Sub_Count_LowExpr = length(intersect(sub_target, genes_lowexpr)),
            stringsAsFactors = FALSE
          )
          names(ef_row)[names(ef_row) == "Y_Value"] <- y_col
          local_exprfilt_results[[length(local_exprfilt_results) + 1]] <- ef_row
        }
      }

      # --- Fine cross: fixed small-effect LFC bin x fine BaseMean category ---
      # Feeds the "figure of 4" (Figure T / Figure S): for each of the two
      # target small-effect LFC bins, stratify GT / predicted genes into the
      # four BASEMEAN_INTERVALS categories and compute recall/precision per
      # (donors, cells/reads) grid cell. Only the FINECROSS_LFC_BINS (Tiny,
      # Small) are crossed to keep the output compact.
      for (fc_lfc_name in FINECROSS_LFC_BINS) {
        fc_bounds <- LFC_INTERVALS[[fc_lfc_name]]
        for (fc_bm_name in names(BASEMEAN_INTERVALS)) {
          fc_bm_bounds <- BASEMEAN_INTERVALS[[fc_bm_name]]

          gt_in_fc <- rownames(gt_res)[
            !is.na(gt_res$padj) &
            gt_res$padj < FDR_THRESHOLD &
            abs(gt_res$log2FoldChange) >= fc_bounds[1] &
            abs(gt_res$log2FoldChange) < fc_bounds[2] &
            gt_res$bm_ruler >= fc_bm_bounds[1] &
            gt_res$bm_ruler < fc_bm_bounds[2]
          ]

          sub_in_fc <- rownames(sub_res)[
            !is.na(sub_res$padj) &
            sub_res$padj < FDR_THRESHOLD &
            abs(sub_res$log2FoldChange) >= fc_bounds[1] &
            abs(sub_res$log2FoldChange) < fc_bounds[2] &
            sub_res$bm_ruler >= fc_bm_bounds[1] &
            sub_res$bm_ruler < fc_bm_bounds[2]
          ]

          for (cat in CHR_CATEGORIES) {
            if (cat == "Total") {
              gt_target <- gt_in_fc; sub_target <- sub_in_fc
              gt_sig_ref <- gt_all_sig; sub_sig_ref <- sub_all_sig
            } else if (cat == "Autosomal") {
              gt_target <- intersect(gt_in_fc, genes_autosomal)
              sub_target <- intersect(sub_in_fc, genes_autosomal)
              gt_sig_ref <- intersect(gt_all_sig, genes_autosomal)
              sub_sig_ref <- intersect(sub_all_sig, genes_autosomal)
            } else if (cat == "SexChr") {
              gt_target <- intersect(gt_in_fc, genes_sex)
              sub_target <- intersect(sub_in_fc, genes_sex)
              gt_sig_ref <- intersect(gt_all_sig, genes_sex)
              sub_sig_ref <- intersect(sub_all_sig, genes_sex)
            }

            n_gt_raw <- length(gt_target)
            gt_target_harm <- intersect(gt_target, sub_all_genes)
            n_gt_harm <- length(gt_target_harm)

            tp_prec <- length(intersect(sub_target, gt_sig_ref))
            fp_prec <- length(setdiff(sub_target, gt_sig_ref))
            precision <- if ((tp_prec + fp_prec) > 0) tp_prec / (tp_prec + fp_prec) else NA_real_

            tp_rec <- length(intersect(gt_target, sub_sig_ref))
            fn_rec <- length(setdiff(gt_target, sub_sig_ref))
            recall <- if ((tp_rec + fn_rec) > 0) tp_rec / (tp_rec + fn_rec) else 0

            tp_rec_harm <- length(intersect(gt_target_harm, sub_sig_ref))
            fn_rec_harm <- length(setdiff(gt_target_harm, sub_sig_ref))
            recall_harm <- if ((tp_rec_harm + fn_rec_harm) > 0) tp_rec_harm / (tp_rec_harm + fn_rec_harm) else 0

            sub_target_harm <- intersect(sub_target, all_genes_in_run)
            tp_prec_harm <- length(intersect(sub_target_harm, gt_sig_ref))
            fp_prec_harm <- length(setdiff(sub_target_harm, gt_sig_ref))
            precision_harm <- if ((tp_prec_harm + fp_prec_harm) > 0) tp_prec_harm / (tp_prec_harm + fp_prec_harm) else NA_real_

            fc_row <- data.frame(
              Run = run_id, Donors = n_donors, Y_Value = y_value,
              Iteration = iter, LFC_Bin = fc_lfc_name, BaseMean_Bin = fc_bm_name,
              Gene_Category = cat,
              Precision = precision, Recall = recall,
              Precision_Harm = precision_harm, Recall_Harm = recall_harm,
              GT_Count = n_gt_raw, GT_Harm_Count = n_gt_harm,
              Sub_Count = length(sub_target),
              Sub_Harm_Count = length(sub_target_harm),
              Sub_Count_LowExpr = length(intersect(sub_target, genes_lowexpr)),
              stringsAsFactors = FALSE
            )
            names(fc_row)[names(fc_row) == "Y_Value"] <- y_col
            local_finecross_results[[length(local_finecross_results) + 1]] <- fc_row
          }
        }
      }
      # --- P-value quartile bins (sensitivity check vs LFC binning) ---
      # Bins defined by GT padj quartiles computed at the top of this run's block.
      # Both GT and sub genes are assigned by their own padj vs the same cutoffs,
      # mirroring the LFC bin approach (bin boundaries from GT, applied to both sides).
      # P-value quartile bins (Total only — gene universe is fixed within a run,
      # so Autosomal/SexChr breakdown and harmonisation add no value here).
      # Recall = proportion of GT Q_k genes detected by sub as significant.
      # Precision = proportion of all sub sig genes that are in GT Q_k
      #             (shows how much of sub's signal is contributed by each quartile).
      if (!is.null(pval_q_intervals)) {
        for (pq_name in names(pval_q_intervals)) {
          pq_bounds <- pval_q_intervals[[pq_name]]

          gt_in_pq <- rownames(gt_res)[
            !is.na(gt_res$padj) &
            gt_res$padj < FDR_THRESHOLD &
            !is.na(gt_res$log2FoldChange) &
            abs(gt_res$log2FoldChange) >= 0.2 &
            gt_res$padj >= pq_bounds[1] &
            gt_res$padj <  pq_bounds[2]
          ]

          # Recall: how many GT Q_k genes does the sub run detect at all?
          tp_rec <- length(intersect(gt_in_pq, sub_all_sig))
          recall  <- if (length(gt_in_pq) > 0) tp_rec / length(gt_in_pq) else 0

          # Precision: of all sub sig genes, what fraction are in GT Q_k?
          precision <- if (length(sub_all_sig) > 0) tp_rec / length(sub_all_sig) else NA_real_

          pq_row <- data.frame(
            Run = run_id, Donors = n_donors, Y_Value = y_value,
            Iteration = iter, PVal_Quartile = pq_name,
            Precision = precision, Recall = recall,
            GT_Count = length(gt_in_pq),
            Sub_Count = length(sub_all_sig),
            PQ_Lo = pq_bounds[1],
            PQ_Hi = pq_bounds[2],
            stringsAsFactors = FALSE
          )
          names(pq_row)[names(pq_row) == "Y_Value"] <- y_col
          local_pval_q_results[[length(local_pval_q_results) + 1]] <- pq_row
        }
      }
    }

    return(list(
      res = do.call(rbind, local_results),
      xist = do.call(rbind, local_xist),
      reads = do.call(rbind, local_reads),
      bm = if (length(local_bm_results) > 0) do.call(rbind, local_bm_results) else NULL,
      cart = if (length(local_cart_results) > 0) do.call(rbind, local_cart_results) else NULL,
      exprfilt = if (length(local_exprfilt_results) > 0) do.call(rbind, local_exprfilt_results) else NULL,
      finecross = if (length(local_finecross_results) > 0) do.call(rbind, local_finecross_results) else NULL,
      pval_q = if (length(local_pval_q_results) > 0) do.call(rbind, local_pval_q_results) else NULL
    ))

  }, mc.cores = min(N_CORES, length(run_folders)))
  
  # ==============================================================================
  # 5. COMBINE AND SUMMARIZE
  # ==============================================================================
  
  # Remove NULL results
  parallel_results <- parallel_results[!sapply(parallel_results, is.null)]
  
  if (length(parallel_results) == 0) {
    log_warning("No valid results for %s", cell_type)
    next
  }
  
  # Combine all results
  df <- do.call(rbind, lapply(parallel_results, `[[`, "res"))
  df_xist <- do.call(rbind, lapply(parallel_results, `[[`, "xist"))
  df_reads <- do.call(rbind, lapply(parallel_results, `[[`, "reads"))

  # Combine baseMean and cartesian results
  bm_parts <- Filter(Negate(is.null), lapply(parallel_results, `[[`, "bm"))
  df_bm <- if (length(bm_parts) > 0) do.call(rbind, bm_parts) else NULL

  cart_parts <- Filter(Negate(is.null), lapply(parallel_results, `[[`, "cart"))
  df_cart <- if (length(cart_parts) > 0) do.call(rbind, cart_parts) else NULL

  exprfilt_parts <- Filter(Negate(is.null), lapply(parallel_results, `[[`, "exprfilt"))
  df_exprfilt <- if (length(exprfilt_parts) > 0) do.call(rbind, exprfilt_parts) else NULL

  finecross_parts <- Filter(Negate(is.null), lapply(parallel_results, `[[`, "finecross"))
  df_finecross <- if (length(finecross_parts) > 0) do.call(rbind, finecross_parts) else NULL

  pval_q_parts <- Filter(Negate(is.null), lapply(parallel_results, `[[`, "pval_q"))
  df_pval_q <- if (length(pval_q_parts) > 0) do.call(rbind, pval_q_parts) else NULL

  # Clean up Y-axis values (Cells or Target_Reads depending on mode)
  df[[y_col]] <- ifelse(df[[y_col]] == "NA", "Max", df[[y_col]])
  df_xist[[y_col]] <- ifelse(df_xist[[y_col]] == "NA", "Max", df_xist[[y_col]])
  df_reads[[y_col]] <- ifelse(df_reads[[y_col]] == "NA", "Max", df_reads[[y_col]])
  if (!is.null(df_bm)) df_bm[[y_col]] <- ifelse(df_bm[[y_col]] == "NA", "Max", df_bm[[y_col]])
  if (!is.null(df_cart)) df_cart[[y_col]] <- ifelse(df_cart[[y_col]] == "NA", "Max", df_cart[[y_col]])
  if (!is.null(df_exprfilt)) df_exprfilt[[y_col]] <- ifelse(df_exprfilt[[y_col]] == "NA", "Max", df_exprfilt[[y_col]])
  if (!is.null(df_finecross)) df_finecross[[y_col]] <- ifelse(df_finecross[[y_col]] == "NA", "Max", df_finecross[[y_col]])
  if (!is.null(df_pval_q)) df_pval_q[[y_col]] <- ifelse(df_pval_q[[y_col]] == "NA", "Max", df_pval_q[[y_col]])

  log_info("Combined %d result rows, %d XIST rows, %d read rows",
           nrow(df), nrow(df_xist), nrow(df_reads))
  if (!is.null(df_bm)) log_info("BaseMean-binned rows: %d", nrow(df_bm))
  if (!is.null(df_cart)) log_info("Cartesian-binned rows: %d", nrow(df_cart))
  if (!is.null(df_exprfilt)) log_info("Expression-filtered LFC rows: %d", nrow(df_exprfilt))
  if (!is.null(df_pval_q)) log_info("P-value quartile rows: %d", nrow(df_pval_q))

  # ==============================================================================
  # 6. CREATE SUMMARIES
  # ==============================================================================
  
  log_info("Computing summaries...")
  
  df_summary <- df %>%
    group_by(Donors, !!sym(y_col), LFC_Bin, Gene_Category) %>%
    summarise(
      Mean_Precision = mean(Precision, na.rm = TRUE),
      Mean_Recall = mean(Recall, na.rm = TRUE),
      Mean_Precision_Harm = mean(Precision_Harm, na.rm = TRUE),
      Mean_Recall_Harm = mean(Recall_Harm, na.rm = TRUE),
      Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
      Ref_Harm_Count = mean(GT_Harm_Count, na.rm = TRUE),
      Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
      Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
      GT_Count_String = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
      Harm_Count_String = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
      CV_GT = {
        vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
        if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
      },
      CV_Harm = {
        vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
        if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
      },
      .groups = "drop"
    )
  
  df_xist_summary <- df_xist %>%
    group_by(Donors, !!sym(y_col), LFC_Bin) %>%
    summarise(Mean_Detection = mean(Detected, na.rm = TRUE), .groups = "drop")
  
  df_reads_summary <- df_reads %>%
    group_by(Donors, !!sym(y_col)) %>%
    summarise(Mean_Total_Reads = mean(Total_Reads, na.rm = TRUE), .groups = "drop")

  # BaseMean-binned summary
  df_bm_summary <- NULL
  if (!is.null(df_bm) && nrow(df_bm) > 0) {
    df_bm_summary <- df_bm %>%
      group_by(Donors, !!sym(y_col), BaseMean_Bin, Gene_Category) %>%
      summarise(
        Mean_Precision = mean(Precision, na.rm = TRUE),
        Mean_Recall = mean(Recall, na.rm = TRUE),
        Mean_Precision_Harm = mean(Precision_Harm, na.rm = TRUE),
        Mean_Recall_Harm = mean(Recall_Harm, na.rm = TRUE),
        Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
        Ref_Harm_Count = mean(GT_Harm_Count, na.rm = TRUE),
        Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
      Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
        GT_Count_String = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        Harm_Count_String = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        CV_GT = {
          vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        CV_Harm = {
          vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        .groups = "drop"
      )
    log_info("BaseMean summary: %d rows", nrow(df_bm_summary))
  }

  # Cartesian-binned summary
  df_cart_summary <- NULL
  if (!is.null(df_cart) && nrow(df_cart) > 0) {
    df_cart_summary <- df_cart %>%
      group_by(Donors, !!sym(y_col), Cartesian_Bin, Gene_Category) %>%
      summarise(
        Mean_Precision = mean(Precision, na.rm = TRUE),
        Mean_Recall = mean(Recall, na.rm = TRUE),
        Mean_Precision_Harm = mean(Precision_Harm, na.rm = TRUE),
        Mean_Recall_Harm = mean(Recall_Harm, na.rm = TRUE),
        Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
        Ref_Harm_Count = mean(GT_Harm_Count, na.rm = TRUE),
        Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
      Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
        GT_Count_String = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        Harm_Count_String = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        CV_GT = {
          vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        CV_Harm = {
          vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        .groups = "drop"
      )
    log_info("Cartesian summary: %d rows", nrow(df_cart_summary))
  }

  # Expression-filtered LFC summary (baseMean > EXPR_FILTER_THRESHOLD)
  df_exprfilt_summary <- NULL
  if (!is.null(df_exprfilt) && nrow(df_exprfilt) > 0) {
    df_exprfilt_summary <- df_exprfilt %>%
      group_by(Donors, !!sym(y_col), LFC_Bin, Gene_Category) %>%
      summarise(
        Mean_Precision = mean(Precision, na.rm = TRUE),
        Mean_Recall = mean(Recall, na.rm = TRUE),
        Mean_Precision_Harm = mean(Precision_Harm, na.rm = TRUE),
        Mean_Recall_Harm = mean(Recall_Harm, na.rm = TRUE),
        Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
        Ref_Harm_Count = mean(GT_Harm_Count, na.rm = TRUE),
        Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
      Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
        GT_Count_String = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        Harm_Count_String = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        CV_GT = {
          vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        CV_Harm = {
          vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        .groups = "drop"
      )
    log_info("Expression-filtered LFC summary: %d rows", nrow(df_exprfilt_summary))
  }

  # Fine LFC x BaseMean cross summary (Figure T / Figure S source)
  df_finecross_summary <- NULL
  if (!is.null(df_finecross) && nrow(df_finecross) > 0) {
    df_finecross_summary <- df_finecross %>%
      group_by(Donors, !!sym(y_col), LFC_Bin, BaseMean_Bin, Gene_Category) %>%
      summarise(
        Mean_Precision = mean(Precision, na.rm = TRUE),
        Mean_Recall = mean(Recall, na.rm = TRUE),
        Mean_Precision_Harm = mean(Precision_Harm, na.rm = TRUE),
        Mean_Recall_Harm = mean(Recall_Harm, na.rm = TRUE),
        Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
        Ref_Harm_Count = mean(GT_Harm_Count, na.rm = TRUE),
        Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
        Mean_Sub_Count_LowExpr = mean(Sub_Count_LowExpr, na.rm = TRUE),
        GT_Count_String = paste(round(tapply(GT_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        Harm_Count_String = paste(round(tapply(GT_Harm_Count, Run, mean, na.rm = TRUE), 1), collapse = "/"),
        CV_GT = {
          vals <- tapply(GT_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        CV_Harm = {
          vals <- tapply(GT_Harm_Count, Run, mean, na.rm = TRUE)
          if (mean(vals, na.rm = TRUE) > 0) sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE) else 0
        },
        .groups = "drop"
      )
    log_info("Fine-cross (LFC x BaseMean) summary: %d rows", nrow(df_finecross_summary))
  }

  # P-value quartile summary (Total only, no Harm)
  df_pval_q_summary <- NULL
  if (!is.null(df_pval_q) && nrow(df_pval_q) > 0) {
    df_pval_q_summary <- df_pval_q %>%
      group_by(Donors, !!sym(y_col), PVal_Quartile) %>%
      summarise(
        Mean_Precision = mean(Precision, na.rm = TRUE),
        Mean_Recall    = mean(Recall, na.rm = TRUE),
        Ref_GT_Count   = mean(GT_Count, na.rm = TRUE),
        Mean_Sub_Count = mean(Sub_Count, na.rm = TRUE),
        PQ_Lo          = mean(PQ_Lo, na.rm = TRUE),
        PQ_Hi          = mean(PQ_Hi, na.rm = TRUE),
        .groups = "drop"
      )
    log_info("P-value quartile summary: %d rows", nrow(df_pval_q_summary))
  }

  # ==============================================================================
  # 7. SAVE CHECKPOINT
  # ==============================================================================

  log_info("Saving checkpoint...")

  qs::qsave(list(
    df_summary = df_summary,
    df_xist_summary = df_xist_summary,
    df_reads_summary = df_reads_summary,
    df_bm_summary = df_bm_summary,
    df_cart_summary = df_cart_summary,
    df_exprfilt_summary = df_exprfilt_summary,
    df_finecross_summary = df_finecross_summary,
    df_pval_q_summary = df_pval_q_summary,
    df_raw = df,
    df_xist_raw = df_xist,
    df_reads_raw = df_reads,
    df_bm_raw = df_bm,
    df_cart_raw = df_cart,
    df_exprfilt_raw = df_exprfilt,
    df_finecross_raw = df_finecross,
    df_pval_q_raw = df_pval_q
  ), checkpoint_file)

  log_info("Saved checkpoint to %s", checkpoint_file)
  log_info("Summary: %d rows | XIST: %d rows | Reads: %d rows",
           nrow(df_summary), nrow(df_xist_summary), nrow(df_reads_summary))
  if (!is.null(df_bm_summary)) log_info("BaseMean summary: %d rows", nrow(df_bm_summary))
  if (!is.null(df_cart_summary)) log_info("Cartesian summary: %d rows", nrow(df_cart_summary))
  if (!is.null(df_exprfilt_summary)) log_info("Expression-filtered summary: %d rows", nrow(df_exprfilt_summary))
  if (!is.null(df_pval_q_summary)) log_info("P-value quartile summary: %d rows", nrow(df_pval_q_summary))

  } # End method loop

} # End cell type loop

log_section("STEP 3 COMPLETE")
log_info("Aggregation complete for all cell types")
