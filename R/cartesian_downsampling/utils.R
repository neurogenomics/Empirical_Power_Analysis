# ==============================================================================
# UTILS: Cartesian Downsampling Pipeline - Shared Utility Functions
# ==============================================================================
# This file contains all shared utility functions for the Cartesian downsampling
# pipeline. Source this file after config.R.
# ==============================================================================

# ==============================================================================
# 1. AGE PARSING FUNCTIONS
# ==============================================================================

#' Check if a string represents an unknown/missing value
#' 
#' @param x Character vector to check
#' @return Logical vector indicating unknown values
is_unknown_string <- function(x) {
  s <- tolower(trimws(as.character(x)))
  grepl("^(unknown|not available|na|n/a|nan|)$", s)
}

#' Parse age stage strings to numeric years
#' 
#' Handles formats like "65 year-old stage", "6 month-old stage", etc.
#' 
#' @param x Character vector of age strings
#' @return Numeric vector of ages in years
parse_age_stage_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
  
  # Extract numeric value
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))
  
  # Detect time units
  is_year    <- grepl("year[ -]?old", s)
  is_month   <- grepl("month[ -]?old", s)
  is_week    <- grepl("week[ -]?old", s)
  is_day     <- grepl("day[ -]?old", s)
  is_unknown <- is_unknown_string(s)
  
  # Convert to years

y <- rep(NA_real_, length(s))
  y[is_year]    <- n[is_year]
  y[is_month]   <- n[is_month] / 12
  y[is_week]    <- n[is_week] / 52.1429
  y[is_day]     <- n[is_day] / 365.25
  y[is_unknown] <- NA_real_
  
  # Fallback: if it's just a number with no unit, assume years
  just_number <- !is_year & !is_month & !is_week & !is_day & !is_unknown & !is.na(n)
  y[just_number] <- n[just_number]
  
  return(y)
}

#' Standardize age values (mean 0, variance 1)
#' 
#' @param age_vec Numeric vector of ages
#' @param impute_na Whether to impute NA values with mean
#' @return Standardized numeric vector
standardize_age <- function(age_vec, impute_na = TRUE) {
  if (impute_na && any(is.na(age_vec))) {
    mean_age <- mean(age_vec, na.rm = TRUE)
    if (is.na(mean_age)) mean_age <- 50  
    age_vec[is.na(age_vec)] <- mean_age
    log_warning("Imputed %d NA age values with mean (%.1f)", 
                sum(is.na(age_vec)), mean_age)
  }
  
  as.numeric(scale(age_vec, center = TRUE, scale = TRUE))
}

# ==============================================================================
# 2. SHARD AGGREGATION FUNCTIONS
# ==============================================================================

#' Aggregate counts from donor shards into pseudobulk
#' 
#' @param donor_ids Character vector of donor IDs to aggregate
#' @param shard_dir Path to directory containing donor shards
#' @param n_cells_to_sample Number of cells to sample per donor (NA = all) - used in cells mode
#' @param target_reads Target total reads per donor (NA = all) - used in reads mode
#' @param mode Downsampling mode: "cells" or "reads"
#' @param seed Random seed for cell sampling
#' @param verbose Whether to print progress messages
#' @return List with 'sce' (SingleCellExperiment) and 'sampling_log' (cell indices)
aggregate_from_shards <- function(donor_ids, shard_dir, 
                                   n_cells_to_sample = NA,
                                   target_reads = NA,
                                   mode = "cells",
                                   seed = 42, verbose = FALSE) {
  
  results_list <- lapply(donor_ids, function(id) {
    # Construct shard path
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
    
    # Check file exists
    if (!file.exists(shard_path)) {
      log_error("Missing shard: %s", shard_path)
      return(NULL)
    }
    
    # Load shard
    sce_shard <- tryCatch({
      qs::qread(shard_path)
    }, error = function(e) {
      log_error("Failed to load shard %s: %s", shard_path, e$message)
      return(NULL)
    })
    
    if (is.null(sce_shard)) return(NULL)
    
    available <- ncol(sce_shard)
    
    # Mode-specific sampling logic
    if (mode == "reads") {
      # --- READS MODE: Sample cells until target read count is reached ---
      if (!is.na(target_reads)) {
        # Calculate reads per cell
        reads_per_cell <- Matrix::colSums(SummarizedExperiment::assay(sce_shard, "counts"))
        total_available_reads <- sum(reads_per_cell)
        
        # Skip donor if insufficient reads (strict less than, not <=)
        if (total_available_reads < target_reads) {
          if (verbose) log_info("Donor %s has only %.0f reads (need %.0f), skipping",
                                id, total_available_reads, target_reads)
          return(NULL)
        }
        
        # Edge case: if target equals total, just use all cells
        if (total_available_reads == target_reads) {
          chosen_idx <- seq_len(available)
          # Continue to aggregation (no downsampling needed)
        } else {
          # Normal case: target < total, need to sample
          
          # 1. Shuffle cells randomly (using per-job seed, not fixed seed)
          set.seed(seed)
          shuffled_indices <- sample(seq_len(available))
          shuffled_reads <- reads_per_cell[shuffled_indices]
          
          # 2. Cumulative sum to find cutoff
          cum_reads <- cumsum(shuffled_reads)
          cutoff_idx <- which(cum_reads >= target_reads)[1]
          
          # 3. Select cells up to cutoff
          final_indices <- shuffled_indices[1:cutoff_idx]
          sce_shard <- sce_shard[, final_indices]
          
          # 4. Multinomial downsample the overflow cell
          current_total <- cum_reads[cutoff_idx]
          excess_reads <- current_total - target_reads
          
          if (excess_reads > 0) {
            last_cell_pos <- length(final_indices)
            last_cell_counts <- SummarizedExperiment::assay(sce_shard, "counts")[, last_cell_pos]
            original_cell_total <- sum(last_cell_counts)
            new_cell_total <- original_cell_total - excess_reads
            
            if (new_cell_total > 0) {
              # Multinomial downsampling: preserve gene distribution
              probs <- last_cell_counts / original_cell_total
              downsampled_counts <- as.numeric(rmultinom(1, size = new_cell_total, prob = probs))
              
              # Replace counts in SCE (keep as sparse matrix)
              new_counts_matrix <- SummarizedExperiment::assay(sce_shard, "counts")
              new_counts_matrix[, last_cell_pos] <- downsampled_counts
              SummarizedExperiment::assay(sce_shard, "counts") <- new_counts_matrix
            } else {
              # Edge case: overflow cell needs 0 reads, drop it
              sce_shard <- sce_shard[, -last_cell_pos]
              final_indices <- final_indices[-last_cell_pos]
            }
          }
          
          chosen_idx <- final_indices
        }
      } else {
        # NA target_reads: use all cells
        chosen_idx <- seq_len(available)
      }
      
    } else {
      # --- CELLS MODE: Random cell selection (original behavior) ---
      if (!is.na(n_cells_to_sample)) {
        if (available < n_cells_to_sample) {
          if (verbose) log_info("Donor %s has only %d cells (need %d), skipping", 
                                id, available, n_cells_to_sample)
          return(NULL)
        }
        
        set.seed(seed)
        chosen_idx <- sort(sample(seq_len(available), n_cells_to_sample))
        sce_shard <- sce_shard[, chosen_idx]
      } else {
        chosen_idx <- seq_len(available)
      }
    }
    
    # Aggregate to pseudobulk (sum counts across cells)
    counts_vec <- Matrix::rowSums(SummarizedExperiment::assay(sce_shard, "counts"))
    meta_row <- as.data.frame(SummarizedExperiment::colData(sce_shard)[1, , drop = FALSE])
    
    return(list(
      counts = counts_vec,
      meta = meta_row,
      genes = names(counts_vec),
      donor_id = id,
      indices = chosen_idx
    ))
  })
  
  # Remove NULL results (failed donors)
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) == 0) {
    log_warning("No valid donors after aggregation")
    return(NULL)
  }
  
  # Combine counts (ensure gene alignment)
  all_genes <- results_list[[1]]$genes
  
  combined_counts <- do.call(cbind, lapply(results_list, function(x) {
    if (!identical(x$genes, all_genes)) {
      # Reorder genes to match reference
      return(x$counts[match(all_genes, x$genes)])
    }
    return(x$counts)
  }))
  
  # Combine metadata
  combined_meta <- do.call(rbind, lapply(results_list, `[[`, "meta"))
  
  # Set row/column names
  rownames(combined_counts) <- all_genes
  colnames(combined_counts) <- as.character(combined_meta[[DONOR_ID_COL]])
  rownames(combined_meta) <- colnames(combined_counts)
  
  # Create sampling log
  sampling_log <- lapply(results_list, `[[`, "indices")
  names(sampling_log) <- as.character(combined_meta[[DONOR_ID_COL]])
  
  # Create SingleCellExperiment (pseudobulk)
  sce_out <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = combined_counts),
    colData = S4Vectors::DataFrame(combined_meta)
  )
  
  return(list(sce = sce_out, sampling_log = sampling_log))
}

# ==============================================================================
# 3. DIFFERENTIAL EXPRESSION WORKERS
# ==============================================================================

prepare_de_colData <- function(agg_result,
                               sex_col = SEX_COL,
                               age_col = AGE_COL,
                               cat_vars = CATEGORICAL_VARS,
                               num_vars = if (exists("NUMERIC_VARS")) NUMERIC_VARS else character()) {
  pb_sce <- agg_result$sce
  cts <- as.data.frame(SummarizedExperiment::colData(pb_sce))

  if (!sex_col %in% colnames(cts)) {
    stop(sprintf("Sex column '%s' not present in colData", sex_col))
  }
  cts[[sex_col]] <- factor(as.character(cts[[sex_col]]),
                           levels = c(FEMALE_VAL, MALE_VAL))
  if (any(is.na(cts[[sex_col]]))) {
    stop(sprintf("Sex column '%s' contains values outside %s/%s",
                 sex_col, FEMALE_VAL, MALE_VAL))
  }
  if (nlevels(droplevels(cts[[sex_col]])) < 2) {
    stop("DE design needs both male and female donors")
  }

  raw_age <- as.character(cts[[age_col]])
  num_age <- parse_age_stage_to_years(raw_age)
  age_has_variance <- length(unique(num_age[!is.na(num_age)])) > 1

  if (age_has_variance) {
    cts$age_scaled <- standardize_age(num_age, impute_na = TRUE)
  } else {
    log_warning("Age has no variance, excluding from model")
  }

  valid_vars <- character()
  for (v in cat_vars) {
    if (v %in% colnames(cts) && length(unique(cts[[v]])) > 1) {
      cts[[v]] <- factor(cts[[v]])
      sex_var_table <- table(cts[[sex_col]], cts[[v]])
      if (all(rowSums(sex_var_table > 0) > 1) &&
          all(colSums(sex_var_table > 0) > 1)) {
        valid_vars <- c(valid_vars, v)
      } else {
        log_warning("Covariate '%s' appears collinear with sex, excluding", v)
      }
    }
  }

  numeric_terms <- character()
  for (nv in num_vars) {
    if (!nv %in% colnames(cts)) {
      log_warning("Numeric covariate '%s' not present in colData. Skipping.", nv)
      next
    }
    raw_val <- suppressWarnings(as.numeric(as.character(cts[[nv]])))
    if (all(is.na(raw_val))) {
      log_warning("Numeric covariate '%s' is all NA. Skipping.", nv)
      next
    }
    if (any(is.na(raw_val))) {
      raw_val[is.na(raw_val)] <- mean(raw_val, na.rm = TRUE)
    }
    if (sd(raw_val, na.rm = TRUE) == 0) {
      log_warning("Numeric covariate '%s' has zero variance. Skipping.", nv)
      next
    }
    scaled_name <- paste0(nv, "_scaled")
    cts[[scaled_name]] <- as.numeric(scale(raw_val))
    numeric_terms <- c(numeric_terms, scaled_name)
  }

  formula_parts <- c(sex_col)
  if (age_has_variance) formula_parts <- c(formula_parts, "age_scaled")
  formula_parts <- c(formula_parts, numeric_terms, valid_vars)
  formula_str <- paste0("~ ", paste(formula_parts, collapse = " + "))

  design_matrix <- tryCatch({
    stats::model.matrix(stats::as.formula(formula_str), data = cts)
  }, error = function(e) {
    log_warning("Model '%s' failed: %s. Trying sex-only model.", formula_str, e$message)
    NULL
  })

  if (is.null(design_matrix) || qr(design_matrix)$rank < ncol(design_matrix)) {
    if (!is.null(design_matrix)) {
      log_warning("Model '%s' is rank-deficient. Trying sex-only model.", formula_str)
    }
    formula_parts <- c(sex_col)
    formula_str <- paste0("~ ", sex_col)
    design_matrix <- stats::model.matrix(stats::as.formula(formula_str), data = cts)
  }

  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop(sprintf("Sex-only model is rank-deficient for column '%s'", sex_col))
  }

  sex_coef_name <- paste0(sex_col, MALE_VAL)
  if (!sex_coef_name %in% colnames(design_matrix)) {
    candidates <- grep(paste0("^", sex_col, ".*", MALE_VAL, "$"),
                       colnames(design_matrix), value = TRUE)
    if (length(candidates) == 1) {
      sex_coef_name <- candidates
    } else {
      stop(sprintf("Could not identify sex coefficient in design: %s",
                   paste(colnames(design_matrix), collapse = ", ")))
    }
  }

  list(
    pb_sce = pb_sce,
    counts = SummarizedExperiment::assay(pb_sce, "counts"),
    col_data = cts,
    design = stats::as.formula(formula_str),
    design_matrix = design_matrix,
    formula_str = formula_str,
    sex_col = sex_col,
    sex_coef_name = sex_coef_name,
    deseq_coef_name = paste0(sex_col, "_", MALE_VAL, "_vs_", FEMALE_VAL)
  )
}

sex_only_prepped <- function(prepped) {
  formula_str <- paste0("~ ", prepped$sex_col)
  design_matrix <- stats::model.matrix(stats::as.formula(formula_str),
                                       data = prepped$col_data)
  prepped$design <- stats::as.formula(formula_str)
  prepped$design_matrix <- design_matrix
  prepped$formula_str <- formula_str
  prepped$sex_coef_name <- paste0(prepped$sex_col, MALE_VAL)
  prepped
}

save_de_inputs <- function(agg_result, save_path_log) {
  if (!is.null(save_path_log)) {
    dir.create(dirname(save_path_log), recursive = TRUE, showWarnings = FALSE)
    qs::qsave(agg_result$sampling_log, save_path_log)
  }
}

filtered_counts <- function(prepped, min_gene_count) {
  keep <- Matrix::rowSums(prepped$counts) >= min_gene_count
  counts <- prepped$counts[keep, , drop = FALSE]
  if (nrow(counts) == 0) {
    stop(sprintf("No genes remain after min count filter (%s)", min_gene_count))
  }
  counts
}

#' Apply ashr adaptive shrinkage to a vector of log2 fold changes.
#'
#' Provides a uniform LFC-shrinkage step across all three DE methods so effect
#' sizes are directly comparable. ashr (Stephens 2017) takes (estimate, SE)
#' pairs and returns posterior means under an adaptive unimodal prior. The
#' `mixcompdist = "normal"`, `method = "shrink"` settings match the call
#' DESeq2's own `lfcShrink(type = "ashr")` makes, so all three methods use an
#' identical estimator. Genes with a non-finite or non-positive standard error
#' (ashr requires SE > 0) keep their unshrunk log2FC.
#'
#' @param log2fc Numeric vector of log2 fold changes (the point estimates).
#' @param se     Numeric vector of standard errors on the same (log2) scale.
#' @return Numeric vector of ashr posterior-mean log2 fold changes.
ashr_shrink_logfc <- function(log2fc, se) {
  log2fc <- as.numeric(log2fc)
  se <- as.numeric(se)
  ok <- is.finite(log2fc) & is.finite(se) & se > 0
  out <- log2fc
  if (any(ok)) {
    fit <- ashr::ash(log2fc[ok], se[ok], mixcompdist = "normal", method = "shrink")
    out[ok] <- ashr::get_pm(fit)
  }
  out
}

#' Back out a per-gene standard error from edgeR's quasi-likelihood F statistic.
#'
#' For a single tested coefficient the QL F-test has 1 numerator df, so the
#' signed t-ratio is t = sign(logFC) * sqrt(F). Because the estimate/SE ratio is
#' scale-invariant, SE = |logFC| / sqrt(F) on the log2 scale. Non-finite F,
#' F <= 0, or logFC == 0 yield NA (handled downstream by keeping the unshrunk
#' logFC), so ashr only ever receives strictly positive finite SEs.
#'
#' @param log2fc Numeric vector of log2 fold changes from edgeR topTags.
#' @param fstat  Numeric vector of QL F statistics (edgeR `F` column).
#' @return Numeric vector of approximate standard errors on the log2 scale.
se_from_edger_f <- function(log2fc, fstat) {
  log2fc <- as.numeric(log2fc)
  fstat <- as.numeric(fstat)
  se <- abs(log2fc) / sqrt(fstat)
  se[!is.finite(fstat) | fstat <= 0 | log2fc == 0] <- NA_real_
  se
}

canonical_de_table <- function(tab,
                               counts,
                               logfc_col,
                               padj_col,
                               pvalue_col = NULL,
                               method) {
  common <- intersect(rownames(tab), rownames(counts))
  tab <- tab[common, , drop = FALSE]
  counts <- counts[common, , drop = FALSE]

  # `median_raw` is the shared cross-method expression ruler: the median raw
  # pseudobulk count across donors, computed identically for all three DE methods
  # (see run_deseq_worker). Every per-gene expression comparison (model feature,
  # expression bins, filters) uses it. `baseMean_raw` (raw mean) is kept only for
  # the total-reads summary in 03_Aggregate_Results.R.
  raw_mean <- as.numeric(Matrix::rowMeans(counts))
  res <- data.frame(
    baseMean = raw_mean,
    baseMean_raw = raw_mean,
    median_raw = as.numeric(matrixStats::rowMedians(as.matrix(counts))),
    log2FoldChange = as.numeric(tab[[logfc_col]]),
    padj = as.numeric(tab[[padj_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(pvalue_col) && pvalue_col %in% colnames(tab)) {
    res$pvalue <- as.numeric(tab[[pvalue_col]])
  }
  rownames(res) <- common
  attr(res, "de_method") <- method
  res
}

#' Run DESeq2 differential expression analysis.
run_deseq_worker <- function(agg_result,
                             sex_col = SEX_COL,
                             age_col = AGE_COL,
                             cat_vars = CATEGORICAL_VARS,
                             save_path_results,
                             save_path_log,
                             min_gene_count = MIN_GENE_COUNT,
                             num_vars = if (exists("NUMERIC_VARS")) NUMERIC_VARS else character(),
                             fit_type = if (exists("DE_FIT_TYPE")) DE_FIT_TYPE else "parametric",
                             prepped = NULL) {
  save_de_inputs(agg_result, save_path_log)
  if (is.null(prepped)) {
    prepped <- prepare_de_colData(agg_result, sex_col, age_col, cat_vars, num_vars)
  }

  dds <- tryCatch({
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = prepped$counts,
      colData = prepped$col_data,
      design = prepped$design
    ))
  }, error = function(e) {
    log_warning("Model '%s' failed: %s. Trying sex-only model.",
                prepped$formula_str, e$message)
    prepped <<- sex_only_prepped(prepped)
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = prepped$counts,
      colData = prepped$col_data,
      design = prepped$design
    ))
  })

  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_gene_count, ]

  dds <- tryCatch({
    suppressMessages(DESeq2::DESeq(dds, fitType = fit_type, parallel = FALSE, quiet = TRUE))
  }, error = function(e) {
    log_warning("DESeq2 failed with model '%s': %s. Trying sex-only model.",
                prepped$formula_str, e$message)
    prepped <- sex_only_prepped(prepped)
    DESeq2::design(dds) <- prepped$design
    suppressMessages(DESeq2::DESeq(dds, fitType = fit_type, parallel = FALSE, quiet = TRUE))
  })

  coef_available <- prepped$deseq_coef_name %in% DESeq2::resultsNames(dds)
  if (coef_available) {
    res <- suppressMessages(DESeq2::lfcShrink(dds, coef = prepped$deseq_coef_name,
                                               type = "ashr", quiet = TRUE))
  } else {
    log_warning("Coefficient '%s' not found, using ashr with contrast", prepped$deseq_coef_name)
    res <- suppressMessages(DESeq2::lfcShrink(dds,
                                               contrast = c(prepped$sex_col, MALE_VAL, FEMALE_VAL),
                                               type = "ashr", quiet = TRUE))
  }

  attr(res, "shrinkage_type") <- "ashr"
  attr(res, "de_method") <- "deseq2"

  # Shared cross-method expression ruler, computed on the same min-count-filtered
  # raw counts as canonical_de_table(): `median_raw` (median raw pseudobulk count
  # across donors) for all per-gene comparisons, `baseMean_raw` (raw mean) for the
  # total-reads summary only. DESeq2's normalised `baseMean` is retained unchanged.
  raw_counts <- DESeq2::counts(dds)
  res$baseMean_raw <- as.numeric(Matrix::rowMeans(raw_counts)[rownames(res)])
  res$median_raw <- matrixStats::rowMedians(raw_counts)[match(rownames(res), rownames(raw_counts))]

  dir.create(dirname(save_path_results), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(res, save_path_results)

  if (coef_available) "SUCCESS" else "SUCCESS_FALLBACK"
}

run_edger_worker <- function(agg_result,
                             save_path_results,
                             save_path_log,
                             min_gene_count = MIN_GENE_COUNT,
                             prepped = NULL) {
  save_de_inputs(agg_result, save_path_log)
  if (is.null(prepped)) prepped <- prepare_de_colData(agg_result)
  counts <- filtered_counts(prepped, min_gene_count)

  fit_once <- function(prep) {
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge <- edgeR::estimateDisp(dge, prep$design_matrix)
    fit <- edgeR::glmQLFit(dge, prep$design_matrix, robust = TRUE)
    edgeR::glmQLFTest(fit, coef = prep$sex_coef_name)
  }

  qlf <- tryCatch({
    fit_once(prepped)
  }, error = function(e) {
    log_warning("edgeR failed with model '%s': %s. Trying sex-only model.",
                prepped$formula_str, e$message)
    prepped <<- sex_only_prepped(prepped)
    fit_once(prepped)
  })

  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  # Uniform ashr LFC shrinkage: edgeR has no native SE, so recover it from the
  # single-df QL F statistic (SE = |logFC| / sqrt(F) on the log2 scale).
  tab$logFC <- ashr_shrink_logfc(tab$logFC, se_from_edger_f(tab$logFC, tab$F))
  res <- canonical_de_table(tab, counts, "logFC", "FDR", "PValue", "edger")
  dir.create(dirname(save_path_results), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(res, save_path_results)
  if (identical(prepped$formula_str, paste0("~ ", prepped$sex_col))) "SUCCESS_FALLBACK" else "SUCCESS"
}

run_limma_worker <- function(agg_result,
                             save_path_results,
                             save_path_log,
                             min_gene_count = MIN_GENE_COUNT,
                             prepped = NULL) {
  save_de_inputs(agg_result, save_path_log)
  if (is.null(prepped)) prepped <- prepare_de_colData(agg_result)
  counts <- filtered_counts(prepped, min_gene_count)

  fit_once <- function(prep) {
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    v <- limma::voom(dge, prep$design_matrix, plot = FALSE)
    fit <- limma::lmFit(v, prep$design_matrix)
    limma::eBayes(fit, robust = TRUE)
  }

  fit <- tryCatch({
    fit_once(prepped)
  }, error = function(e) {
    log_warning("limma-voom failed with model '%s': %s. Trying sex-only model.",
                prepped$formula_str, e$message)
    prepped <<- sex_only_prepped(prepped)
    fit_once(prepped)
  })

  tab <- limma::topTable(fit, coef = prepped$sex_coef_name, number = Inf,
                         sort.by = "none")
  # Uniform ashr LFC shrinkage: limma gives an exact posterior SE for the
  # coefficient (unscaled SD times sqrt of the moderated residual variance).
  limma_se <- fit$stdev.unscaled[, prepped$sex_coef_name] * sqrt(fit$s2.post)
  limma_se <- limma_se[rownames(tab)]
  tab$logFC <- ashr_shrink_logfc(tab$logFC, limma_se)
  res <- canonical_de_table(tab, counts, "logFC", "adj.P.Val", "P.Value", "limma")
  dir.create(dirname(save_path_results), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(res, save_path_results)
  if (identical(prepped$formula_str, paste0("~ ", prepped$sex_col))) "SUCCESS_FALLBACK" else "SUCCESS"
}

run_de_worker <- function(agg_result,
                          method,
                          save_path_results,
                          save_path_log,
                          min_gene_count = MIN_GENE_COUNT,
                          prepped = NULL) {
  switch(method,
    deseq2 = run_deseq_worker(
      agg_result = agg_result,
      save_path_results = save_path_results,
      save_path_log = save_path_log,
      min_gene_count = min_gene_count,
      prepped = prepped
    ),
    edger = run_edger_worker(
      agg_result = agg_result,
      save_path_results = save_path_results,
      save_path_log = save_path_log,
      min_gene_count = min_gene_count,
      prepped = prepped
    ),
    limma = run_limma_worker(
      agg_result = agg_result,
      save_path_results = save_path_results,
      save_path_log = save_path_log,
      min_gene_count = min_gene_count,
      prepped = prepped
    ),
    stop("Unknown DE method: ", method)
  )
}

run_all_de_methods <- function(agg_result,
                               methods = DE_METHODS,
                               save_path_results,
                               save_path_log,
                               min_gene_count = MIN_GENE_COUNT,
                               sex_col = SEX_COL,
                               age_col = AGE_COL,
                               cat_vars = CATEGORICAL_VARS,
                               num_vars = if (exists("NUMERIC_VARS")) NUMERIC_VARS else character()) {
  stopifnot(all(methods %in% IMPLEMENTED_DE_METHODS))
  if (is.null(names(save_path_results))) names(save_path_results) <- methods
  if (is.null(names(save_path_log))) names(save_path_log) <- methods

  prepped <- prepare_de_colData(agg_result, sex_col, age_col, cat_vars, num_vars)
  status <- setNames(character(length(methods)), methods)

  for (method in methods) {
    status[[method]] <- tryCatch({
      run_de_worker(
        agg_result = agg_result,
        method = method,
        save_path_results = save_path_results[[method]],
        save_path_log = save_path_log[[method]],
        min_gene_count = min_gene_count,
        prepped = prepped
      )
    }, error = function(e) {
      paste0("ERROR: ", conditionMessage(e))
    })
  }

  status
}

# ==============================================================================
# 4. DONOR METADATA FUNCTIONS
# ==============================================================================

#' Load donor metadata from manifest and shards
#' 
#' @param manifest_path Path to manifest file
#' @param shard_dir Path to shard directory
#' @return Data frame with donor_id, sex, cell_count, total_reads
load_donor_metadata <- function(manifest_path, shard_dir) {
  
  if (!file.exists(manifest_path)) {
    log_error("Manifest not found: %s", manifest_path)
    return(NULL)
  }
  
  manifest_df <- qs::qread(manifest_path)$manifest

  # MINI_RUN: cap the donor universe up front so the (slow) per-shard metadata
  # scan only touches a handful of shards. Take 2x the target so the sex-balanced
  # trim below still has a pool to choose from.
  if (exists("MINI_RUN") && isTRUE(MINI_RUN)) {
    cap <- 2L * MINI_MAX_DONORS
    if (nrow(manifest_df) > cap) {
      manifest_df <- manifest_df[seq_len(cap), , drop = FALSE]
    }
  }

  donor_meta <- do.call(rbind, lapply(manifest_df$file_path, function(f) {
    # Handle relative paths
    if (!file.exists(f)) {
      f <- file.path(shard_dir, basename(f))
    }
    
    if (!file.exists(f)) {
      log_warning("Shard file not found: %s", f)
      return(NULL)
    }
    
    tryCatch({
      tmp <- qs::qread(f)
      cd <- SummarizedExperiment::colData(tmp)
      
      # Compute total reads (UMI counts) for reads-mode eligibility filtering
      # This adds ~1 second per shard but keeps logic simple
      total_reads <- sum(SummarizedExperiment::assay(tmp, "counts"))
      
      data.frame(
        donor_id = as.character(cd[[DONOR_ID_COL]][1]),
        sex = as.character(cd[[SEX_COL]][1]),
        cell_count = ncol(tmp),
        total_reads = total_reads,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      log_error("Failed to read shard %s: %s", f, e$message)
      return(NULL)
    })
  }))
  
  if (is.null(donor_meta) || nrow(donor_meta) == 0) {
    log_error("No valid donor metadata loaded")
    return(NULL)
  }

  # MINI_RUN: trim to a deterministic, sex-balanced universe of MINI_MAX_DONORS.
  # Seeded with SEED_BASE so both 01_Generate_Ground_Truth and
  # 02_Run_Downsampling_Grid resolve the identical universe and GT split.
  if (exists("MINI_RUN") && isTRUE(MINI_RUN) && nrow(donor_meta) > MINI_MAX_DONORS) {
    set.seed(SEED_BASE)
    n_each  <- MINI_MAX_DONORS %/% 2L
    males   <- which(donor_meta$sex == MALE_VAL)
    females <- which(donor_meta$sex == FEMALE_VAL)
    keep <- c(
      sample(males,   min(n_each, length(males))),
      sample(females, min(n_each, length(females)))
    )
    if (length(keep) < MINI_MAX_DONORS) {
      remaining <- setdiff(seq_len(nrow(donor_meta)), keep)
      keep <- c(keep, sample(remaining, min(MINI_MAX_DONORS - length(keep), length(remaining))))
    }
    donor_meta <- donor_meta[sort(keep), , drop = FALSE]
    log_info("MINI_RUN: donor universe trimmed to %d donors (%d male, %d female)",
             nrow(donor_meta),
             sum(donor_meta$sex == MALE_VAL), sum(donor_meta$sex == FEMALE_VAL))
  }

  return(donor_meta)
}

#' Calculate sex proportions from donor metadata
#' 
#' @param donor_meta Data frame from load_donor_metadata()
#' @return List with prop_male, prop_female, n_male, n_female, n_total
get_sex_proportions <- function(donor_meta) {
  n_male <- sum(donor_meta$sex == MALE_VAL)
  n_female <- sum(donor_meta$sex == FEMALE_VAL)
  n_total <- nrow(donor_meta)
  
  list(
    prop_male = n_male / n_total,
    prop_female = n_female / n_total,
    n_male = n_male,
    n_female = n_female,
    n_total = n_total
  )
}

# ==============================================================================
# 5. JOB PLANNING FUNCTIONS
# ==============================================================================

#' Generate job plan for downsampling grid
#' 
#' @param val_df Validation donor data frame (donors not in ground truth)
#' @param donor_steps Vector of donor counts to test
#' @param cell_steps Vector of cell counts to test (NA = max) - used in cells mode
#' @param read_steps Vector of read counts to test (NA = max) - used in reads mode
#' @param mode Downsampling mode: "cells" or "reads"
#' @param n_subsamples Number of iterations per grid point
#' @param prop_male Proportion of males in population
#' @param run_dir Output directory for this run
#' @param shard_dir Path to shard directory
#' @return List of job specifications
generate_job_plan <- function(val_df, 
                               donor_steps = DONOR_STEPS,
                               cell_steps = CELL_STEPS,
                               read_steps = READ_STEPS,
                               mode = DOWNSAMPLE_MODE,
                               n_subsamples = N_SUBSAMPLES,
                               prop_male,
                               run_dir,
                               shard_dir) {
  
  # Select Y-axis steps based on mode
  y_steps <- if (mode == "reads") read_steps else cell_steps
  y_col_name <- if (mode == "reads") "n_reads" else "n_cells"
  
  # Create full grid
  job_grid <- expand.grid(
    n_donors = donor_steps,
    y_value = y_steps,
    iter = seq_len(n_subsamples)
  )
  names(job_grid)[2] <- y_col_name  # Rename column based on mode
  
  # Generate job specifications
  plan <- lapply(seq_len(nrow(job_grid)), function(i) {
    row <- job_grid[i, ]
    
    # Filter eligible donors based on mode
    if (mode == "reads") {
      # Reads mode: filter by total_reads
      if (is.na(row$n_reads)) {
        eligible <- val_df
      } else {
        eligible <- val_df[val_df$total_reads >= row$n_reads, ]
      }
    } else {
      # Cells mode: filter by cell_count
      if (is.na(row$n_cells)) {
        eligible <- val_df
      } else {
        eligible <- val_df[val_df$cell_count >= row$n_cells, ]
      }
    }
    
    # Calculate sex-balanced donor counts
    n_m <- round(row$n_donors * prop_male)
    if (n_m == 0 && row$n_donors > 0 && prop_male > 0) n_m <- 1
    n_f <- row$n_donors - n_m
    
    # Get eligible donors by sex
    m_el <- eligible$donor_id[eligible$sex == MALE_VAL]
    f_el <- eligible$donor_id[eligible$sex == FEMALE_VAL]
    
    # Check if we have enough donors
    if (length(m_el) < n_m || length(f_el) < n_f) {
      return(NULL)
    }
    
    # Generate output file paths (mode-specific naming)
    if (mode == "reads") {
      # Format reads as non-scientific string to avoid "1e+05" in filenames
      y_str <- ifelse(is.na(row$n_reads), "NA", 
                      format(row$n_reads, scientific = FALSE, trim = TRUE))
      out_file <- file.path(run_dir, sprintf("N%d_R%s_I%d.qs", row$n_donors, y_str, row$iter))
      log_file <- file.path(run_dir, sprintf("N%d_R%s_I%d_sampling_metadata.qs", row$n_donors, y_str, row$iter))
    } else {
      y_str <- ifelse(is.na(row$n_cells), "NA", as.character(row$n_cells))
      out_file <- file.path(run_dir, sprintf("N%d_C%s_I%d.qs", row$n_donors, y_str, row$iter))
      log_file <- file.path(run_dir, sprintf("N%d_C%s_I%d_sampling_metadata.qs", row$n_donors, y_str, row$iter))
    }
    
    # Build job spec (include both n_cells and n_reads for compatibility)
    job_spec <- list(
      n_donors = row$n_donors,
      iter = row$iter,
      n_male = n_m,
      n_female = n_f,
      male_pool = m_el,
      female_pool = f_el,
      out = out_file,
      log_out = log_file,
      shard_dir = shard_dir,
      mode = mode
    )
    
    # Add mode-specific field
    if (mode == "reads") {
      job_spec$target_reads <- row$n_reads
      job_spec$n_cells <- NA  # For compatibility
    } else {
      job_spec$n_cells <- row$n_cells
      job_spec$target_reads <- NA  # For compatibility
    }
    
    return(job_spec)
  })
  
  # Remove NULL jobs (insufficient donors)
  plan <- plan[!sapply(plan, is.null)]
  
  return(plan)
}

#' Sample donors for a job
#' 
#' @param job Job specification from generate_job_plan()
#' @param seed Random seed
#' @return Character vector of sampled donor IDs
sample_donors_for_job <- function(job, seed) {
  set.seed(seed)
  c(
    sample(job$male_pool, job$n_male),
    sample(job$female_pool, job$n_female)
  )
}

# ==============================================================================
# 5b. HPC JOB MANIFEST FUNCTIONS
# ==============================================================================

is_hpc_chunk_active <- function(job_chunk = JOB_CHUNK) {
  !is.na(job_chunk)
}

job_y_value <- function(job) {
  if (identical(job$mode, "reads")) job$target_reads else job$n_cells
}

job_y_label <- function(job) {
  y_value <- job_y_value(job)
  if (is.na(y_value)) return("NA")
  format(as.numeric(y_value), scientific = FALSE, trim = TRUE)
}

job_manifest_key <- function(cell_type, run_id, job) {
  paste(
    cell_type,
    run_id,
    job$mode,
    job$n_donors,
    job_y_label(job),
    job$iter,
    sep = "|"
  )
}

estimate_job_cost <- function(job) {
  y_value <- job_y_value(job)
  if (is.na(y_value)) {
    y_steps <- get_y_steps(job$mode)
    y_value <- max(y_steps[!is.na(y_steps)])
  }
  as.numeric(job$n_donors) * as.numeric(y_value)
}

assign_cost_balanced_chunks <- function(manifest, n_chunks = JOB_N_CHUNKS) {
  if (nrow(manifest) == 0) return(manifest)

  order_idx <- order(manifest$est_cost, decreasing = TRUE)
  chunk_loads <- rep(0, n_chunks)
  chunk_ids <- integer(nrow(manifest))

  for (idx in order_idx) {
    chunk_id <- which.min(chunk_loads)
    chunk_ids[idx] <- chunk_id
    chunk_loads[chunk_id] <- chunk_loads[chunk_id] + manifest$est_cost[idx]
  }

  manifest$chunk_id <- chunk_ids
  manifest[order(manifest$chunk_id, -manifest$est_cost, manifest$cell_type,
                 manifest$run_id, manifest$n_donors, manifest$iter), ]
}

build_job_manifest <- function(cell_types = get_cell_types(),
                               n_chunks = JOB_N_CHUNKS,
                               manifest_path = JOB_MANIFEST_PATH) {
  rows <- list()
  row_i <- 0L

  for (cell_type in cell_types) {
    log_section(sprintf("Manifest cell type: %s", cell_type))

    paths_by_method <- get_paths_by_method(cell_type)
    paths <- paths_by_method[[DE_METHODS[1]]]

    donor_meta <- load_donor_metadata(paths$manifest_path, paths$shard_dir)
    if (is.null(donor_meta)) {
      stop(sprintf("Could not load donor metadata for %s", cell_type))
    }
    sex_props <- get_sex_proportions(donor_meta)

    for (run_id in RUN_IDS) {
      run_dir <- file.path(paths$results_dir, paste0("Run_", run_id))
      val_meta_path <- file.path(run_dir, "Validation_Donors.qs")
      if (!file.exists(val_meta_path)) {
        stop(sprintf(
          "Validation donors not found for %s run %d: %s. Run Step 1 first.",
          cell_type, run_id, val_meta_path
        ))
      }

      val_df <- qs::qread(val_meta_path)
      plan <- generate_job_plan(
        val_df = val_df,
        donor_steps = DONOR_STEPS,
        cell_steps = CELL_STEPS,
        read_steps = READ_STEPS,
        mode = DOWNSAMPLE_MODE,
        n_subsamples = N_SUBSAMPLES,
        prop_male = sex_props$prop_male,
        run_dir = run_dir,
        shard_dir = paths$shard_dir
      )

      for (job in plan) {
        row_i <- row_i + 1L
        y_value <- job_y_value(job)
        rows[[row_i]] <- data.frame(
          chunk_id = NA_integer_,
          job_key = job_manifest_key(cell_type, run_id, job),
          cell_type = cell_type,
          run_id = run_id,
          mode = job$mode,
          n_donors = job$n_donors,
          y_value = if (is.na(y_value)) NA_real_ else as.numeric(y_value),
          y_label = job_y_label(job),
          iter = job$iter,
          n_male = job$n_male,
          n_female = job$n_female,
          out_basename = basename(job$out),
          est_cost = estimate_job_cost(job),
          stringsAsFactors = FALSE
        )
      }

      log_info("  run %d: %d manifest jobs", run_id, length(plan))
    }
  }

  manifest <- if (length(rows) > 0) {
    do.call(rbind, rows)
  } else {
    data.frame()
  }

  manifest <- assign_cost_balanced_chunks(manifest, n_chunks)

  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(manifest, manifest_path, row.names = FALSE, na = "NA")

  log_info("Wrote manifest: %s", manifest_path)
  if (nrow(manifest) > 0) {
    chunk_summary <- aggregate(est_cost ~ chunk_id, manifest, sum)
    log_info("Manifest jobs: %d | chunks: %d | cost range: %.0f to %.0f",
             nrow(manifest), n_chunks, min(chunk_summary$est_cost),
             max(chunk_summary$est_cost))
  } else {
    log_warning("Manifest contains zero jobs")
  }

  invisible(manifest)
}

filter_plan_for_hpc_chunk <- function(plan,
                                      cell_type,
                                      run_id,
                                      job_chunk = JOB_CHUNK,
                                      manifest_path = JOB_MANIFEST_PATH) {
  if (!is_hpc_chunk_active(job_chunk)) return(plan)

  if (!file.exists(manifest_path)) {
    stop(sprintf(
      "JOB_CHUNK=%d is set but manifest does not exist: %s. Run Step 1.5 first.",
      job_chunk, manifest_path
    ))
  }

  manifest <- read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    na.strings = character()
  )

  chunk_rows <- manifest[
    manifest$chunk_id == job_chunk &
      manifest$cell_type == cell_type &
      manifest$run_id == run_id &
      manifest$mode == DOWNSAMPLE_MODE,
    ,
    drop = FALSE
  ]

  if (nrow(chunk_rows) == 0) return(plan[FALSE])

  chunk_keys <- chunk_rows$job_key
  plan_keys <- vapply(plan, function(job) {
    job_manifest_key(cell_type, run_id, job)
  }, character(1))

  plan[plan_keys %in% chunk_keys]
}

# ==============================================================================
# 6. VISUALIZATION HELPERS
# ==============================================================================

#' Clean cells factor for plotting
#' 
#' @param x_cells Vector of cell values
#' @param cell_scales Reference scale vector
#' @return Factor with proper levels
clean_cells_factor <- function(x_cells, cell_scales = CELL_STEPS) {
  x_cells <- as.character(x_cells)
  x_cells[x_cells == "NA" | is.na(x_cells)] <- "Max"
  
  # Create levels from cell_scales
  levels_vec <- as.character(cell_scales)
  levels_vec[is.na(cell_scales)] <- "Max"
  
  factor(x_cells, levels = levels_vec)
}

#' Clean donors factor for plotting
#' 
#' @param x_donors Vector of donor values
#' @param donor_scales Reference scale vector
#' @return Factor with proper levels
clean_donors_factor <- function(x_donors, donor_scales = DONOR_STEPS) {
  factor(x_donors, levels = as.character(donor_scales))
}

#' Clean reads factor for plotting (handles scientific notation)
#' 
#' Ported from archived reads-based downsampling script.
#' Ensures read values are formatted as non-scientific strings to prevent
#' R from auto-converting large numbers (e.g., 10000000 -> "1e+07") in factor levels.
#' 
#' @param x_reads Vector of read values
#' @param read_scales Reference scale vector
#' @return Factor with proper levels
clean_reads_factor <- function(x_reads, read_scales = READ_STEPS) {
  x_reads <- as.character(x_reads)
  x_reads[x_reads == "NA" | is.na(x_reads)] <- "Max"
  
  # Format levels strictly as non-scientific strings
  # This prevents R from converting 10000000 to "1e+07" in factor levels
  clean_levels <- sapply(read_scales, function(x) {
    if (is.na(x) || x == "Max") return("Max")
    return(format(as.numeric(x), scientific = FALSE, trim = TRUE))
  })
  
  # Normalize x_reads to match clean_levels format
  # (handles case where input has scientific notation like "7.5e+06")
  x_reads_normalized <- sapply(x_reads, function(x) {
    if (x == "Max") return("Max")
    tryCatch({
      format(as.numeric(x), scientific = FALSE, trim = TRUE)
    }, error = function(e) x)
  })
  
  factor(x_reads_normalized, levels = clean_levels)
}

# ==============================================================================
# 7. CHROMOSOME MAPPING
# ==============================================================================

#' Get chromosome mapping for genes
#' 
#' @param gene_ids Character vector of Ensembl gene IDs
#' @return Named character vector mapping gene IDs to chromosomes
get_chr_map <- function(gene_ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    log_warning("org.Hs.eg.db not available, chromosome mapping disabled")
    return(NULL)
  }
  
  tryCatch({
    map <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = gene_ids,
        column = "CHR",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    return(map)
  }, error = function(e) {
    log_warning("Chromosome mapping failed: %s", e$message)
    return(NULL)
  })
}

# ==============================================================================
# 8. FEASIBILITY DATA
# ==============================================================================

#' Load feasibility map for a cell type
#' 
#' @param cell_type Cell type name
#' @param feasibility_dir Directory containing feasibility CSVs
#' @return Long-format data frame or NULL if not found
load_feasibility_map <- function(cell_type, feasibility_dir = FEASIBILITY_DIR) {
  fname <- file.path(feasibility_dir, paste0("Feasibility_", cell_type, ".csv"))
  
  if (!file.exists(fname)) {
    log_info("No feasibility data found for %s", cell_type)
    return(NULL)
  }
  
  raw <- read.csv(fname)
  
  # Handle row names column if present
  if (ncol(raw) == length(DONOR_STEPS) + 1) {
    raw <- raw[, -1]
  }
  
  # Check if feasibility data matches current DONOR_STEPS
  if (ncol(raw) != length(DONOR_STEPS)) {
    log_warning("Feasibility data for %s has %d columns but DONOR_STEPS has %d elements. Skipping feasibility overlay.", 
                cell_type, ncol(raw), length(DONOR_STEPS))
    return(NULL)
  }
  
  colnames(raw) <- as.character(DONOR_STEPS)
  
  # Create cell labels
  cell_labels <- as.character(CELL_STEPS)
  cell_labels[is.na(CELL_STEPS)] <- "Max"
  raw$Cells <- cell_labels
  
  # Pivot to long format
  long_df <- tidyr::pivot_longer(
    raw,
    cols = -Cells,
    names_to = "Donors",
    values_to = "Feasibility_Count"
  )
  
  long_df$Donors <- clean_donors_factor(long_df$Donors)
  long_df$Cells <- clean_cells_factor(long_df$Cells)
  
  return(long_df)
}

# ==============================================================================
# 9. ITERATION RESULT PARSING
# ==============================================================================

#' Parse iteration result filename to extract grid parameters
#'
#' @param filename Filename, not full path
#' @param mode "cells" or "reads"
#' @return Named list with n_donors, y_value, iteration, or NULL if parse fails
parse_iteration_filename <- function(filename, mode) {
  fp <- get_file_pattern(mode)

  m <- regmatches(filename, regexec(fp$parse_regex, filename))[[1]]
  if (length(m) < 4) return(NULL)

  y_raw <- m[3]
  y_value <- if (y_raw == "NA") NA_real_ else as.numeric(y_raw)

  list(
    n_donors = as.numeric(m[2]),
    y_value = y_value,
    iteration = as.numeric(m[4])
  )
}

#' List all iteration result files in a run directory
#'
#' @param run_dir Path to run directory
#' @param mode "cells" or "reads"
#' @return Data frame with columns: file_path, n_donors, y_value, iteration
list_iteration_files <- function(run_dir, mode) {
  fp <- get_file_pattern(mode)

  all_files <- list.files(run_dir, pattern = fp$pattern, full.names = FALSE)
  all_files <- all_files[!grepl("_sampling_metadata", all_files)]
  if (length(all_files) == 0) return(NULL)

  parsed <- lapply(all_files, function(f) {
    p <- parse_iteration_filename(f, mode)
    if (is.null(p)) return(NULL)
    data.frame(
      file_path = file.path(run_dir, f),
      n_donors = p$n_donors,
      y_value = p$y_value,
      iteration = p$iteration,
      stringsAsFactors = FALSE
    )
  })

  parsed <- parsed[!sapply(parsed, is.null)]
  if (length(parsed) == 0) return(NULL)

  do.call(rbind, parsed)
}

# ==============================================================================
# 10. SAFE DIRECTORY CREATION
# ==============================================================================

#' Create directory if it does not exist
#'
#' @param path Directory path
safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    log_info("Created directory: %s", path)
  }
}

# ==============================================================================
# 11. SCRIPT DIRECTORY DETECTION
# ==============================================================================

#' Get the directory containing the current script
#' 
#' Works in multiple contexts: Rscript, source(), RStudio
#' @return Path to script directory
get_script_dir <- function() {
  # Method 1: Rscript command line
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  
  # Method 2: source() or RStudio
  if (exists("ofile", envir = sys.frame(1))) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  
  # Method 3: RStudio active document
  if (Sys.getenv("RSTUDIO") == "1") {
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      if (rstudioapi::isAvailable()) {
        path <- tryCatch(
          rstudioapi::getActiveDocumentContext()$path,
          error = function(e) ""
        )
        if (nzchar(path)) return(dirname(normalizePath(path)))
      }
    }
  }
  
  # Fallback: working directory
  return(getwd())
}
