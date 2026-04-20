# ==============================================================================
# UTILS: Full Ground Truth Pipeline - Shared Utility Functions
# ==============================================================================
# Functions for pseudobulk aggregation, DESeq2 execution, and cross-dataset
# diagnostic analyses. Source this file after config.R.
# ==============================================================================

# ==============================================================================
# 1. PSEUDOBULK AGGREGATION
# ==============================================================================

#' Aggregate donor shards into pseudobulk counts
#'
#' Reads individual donor .qs shard files, sums counts per donor (pseudobulk),
#' and returns a SingleCellExperiment with donors as columns.
#'
#' @param shard_dir Path to directory containing donor_*.qs files
#' @param design_vars Character vector of metadata columns to extract
#' @param n_cores Number of cores for parallel processing
#' @return SingleCellExperiment with pseudobulk counts and metadata
aggregate_shards_parallel <- function(shard_dir, design_vars, n_cores) {

  # Find shard files
  shard_files <- list.files(shard_dir, pattern = "\\.qs$", full.names = TRUE)
  shard_files <- shard_files[!grepl("MANIFEST", shard_files, ignore.case = TRUE)]
  shard_files <- shard_files[!grepl("donor_NA\\.qs$", shard_files)]

  if (length(shard_files) == 0) {
    stop(sprintf("No valid shard files found in %s", shard_dir))
  }

  log_info("Aggregating %d donor shards using %d cores...", length(shard_files), n_cores)

  local_param <- BiocParallel::MulticoreParam(
    workers = n_cores,
    stop.on.error = FALSE
  )

  pb_list <- BiocParallel::bplapply(shard_files, function(f_path) {
    if (!file.exists(f_path)) return(NULL)

    sce_shard <- tryCatch({
      qs::qread(f_path)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(sce_shard)) return(NULL)

    assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(sce_shard)) {
      "counts"
    } else {
      SummarizedExperiment::assayNames(sce_shard)[1]
    }

    # FIX: Convert to integer to ensure DESeq2 compatibility
    # Matrix::rowSums() returns numeric even when summing integers
    count_mat <- SummarizedExperiment::assay(sce_shard, assay_name)
    counts_vec <- as.integer(Matrix::rowSums(count_mat))
    # Use SCE rownames; fall back to assay matrix rownames if NULL
    gene_names <- rownames(sce_shard)
    if (is.null(gene_names)) gene_names <- rownames(count_mat)
    names(counts_vec) <- gene_names

    full_meta <- as.data.frame(SummarizedExperiment::colData(sce_shard))

    # Ensure all requested vars exist, fill NA if missing
    for (mv in setdiff(design_vars, colnames(full_meta))) {
      full_meta[[mv]] <- NA
    }

    meta_row <- full_meta[1, design_vars, drop = FALSE]

    return(list(counts = counts_vec, meta = meta_row, genes = names(counts_vec)))
  }, BPPARAM = local_param)

  # Remove NULLs
  pb_list <- pb_list[!sapply(pb_list, is.null)]
  if (length(pb_list) == 0) stop("No valid shards loaded.")

  log_info("Successfully loaded %d / %d shards", length(pb_list), length(shard_files))

  # Combine into matrix
  template_genes <- pb_list[[1]]$genes
  combined_counts <- do.call(cbind, lapply(pb_list, function(x) {
    if (identical(x$genes, template_genes)) return(x$counts)
    return(x$counts[match(template_genes, x$genes)])
  }))

  # Ensure integer type for DESeq2
  combined_counts <- as.matrix(combined_counts)
  storage.mode(combined_counts) <- "integer"

  combined_meta <- do.call(rbind, lapply(pb_list, `[[`, "meta"))
  rownames(combined_counts) <- template_genes

  sce_out <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = combined_counts),
    colData = S4Vectors::DataFrame(combined_meta)
  )

  return(sce_out)
}

# ==============================================================================
# 2. DESEQ2 EXECUTION
# ==============================================================================

#' Run DESeq2 differential expression analysis on pseudobulk data
#'
#' Handles formula construction, covariate processing, age scaling,
#' rank-deficiency detection, and apeglm LFC shrinkage.
#'
#' @param pb_sce SingleCellExperiment with pseudobulk counts
#' @param cfg Dataset configuration list from DATASET_CONFIG
#' @param de_param BiocParallel parameter for DESeq2
#' @param min_gene_count Minimum total count to keep a gene
#' @return List with: dds (DESeqDataSet), results (DESeqResults), res_df (data.frame)
run_deseq_for_dataset <- function(pb_sce, cfg, de_param, min_gene_count = 10) {

  meta_df <- SummarizedExperiment::colData(pb_sce)
  sex_col <- cfg$cols$sex
  age_col <- cfg$cols$age
  male_val <- cfg$sex_values$male
  female_val <- cfg$sex_values$female

  # --- Check sex column ---
  if (!sex_col %in% colnames(meta_df)) {
    stop(sprintf("Sex column '%s' not found. Available: %s",
                 sex_col, paste(colnames(meta_df), collapse = ", ")))
  }

  # --- Check replicates ---
  n_m <- sum(meta_df[[sex_col]] == male_val, na.rm = TRUE)
  n_f <- sum(meta_df[[sex_col]] == female_val, na.rm = TRUE)
  log_info("Sample counts - Male (%s): %d, Female (%s): %d",
           male_val, n_m, female_val, n_f)

  if (n_m < 2 || n_f < 2) {
    stop(sprintf("Insufficient replicates (M: %d, F: %d). Need >= 2 each.", n_m, n_f))
  }

  # --- Process AGE (standardized) ---
  if (!age_col %in% colnames(meta_df)) {
    log_warn("Age column '%s' not found - using placeholder", age_col)
    meta_df$internal_age_scaled <- 0
  } else {
    raw_age <- as.character(meta_df[[age_col]])
    num_age <- cfg$age_parser(raw_age)

    if (all(is.na(num_age))) {
      log_warn("All ages are NA. Using 0 for age covariate.")
      meta_df$internal_age_scaled <- 0
    } else {
      if (any(is.na(num_age))) {
        n_na <- sum(is.na(num_age))
        log_warn("%d ages failed parsing. Imputing mean.", n_na)
        num_age[is.na(num_age)] <- mean(num_age, na.rm = TRUE)
      }

      if (sd(num_age, na.rm = TRUE) > 0) {
        meta_df$internal_age_scaled <- as.numeric(scale(num_age))
      } else {
        log_warn("Age has zero variance. Using 0 for age covariate.")
        meta_df$internal_age_scaled <- 0
      }
    }
  }

  # --- Process categorical covariates ---
  design_terms <- c()
  for (cv in cfg$covariates) {
    if (cv %in% colnames(meta_df)) {
      n_levels <- length(unique(na.omit(meta_df[[cv]])))
      if (n_levels > 1) {
        design_terms <- c(design_terms, cv)
        log_info("Added categorical covariate: %s (%d levels)", cv, n_levels)
      } else {
        log_warn("Covariate '%s' has only 1 level. Skipping.", cv)
      }
    } else {
      log_warn("Covariate '%s' not found in metadata. Skipping.", cv)
    }
  }

  # --- Process numeric covariates (standardized) ---
  for (nv in cfg$numeric_covariates) {
    if (nv %in% colnames(meta_df)) {
      raw_val <- as.numeric(as.character(meta_df[[nv]]))

      if (all(is.na(raw_val))) {
        log_warn("Numeric covariate '%s' is all NA. Skipping.", nv)
        next
      }

      if (any(is.na(raw_val))) {
        raw_val[is.na(raw_val)] <- mean(raw_val, na.rm = TRUE)
      }

      if (sd(raw_val, na.rm = TRUE) > 0) {
        scaled_name <- paste0(nv, "_scaled")
        meta_df[[scaled_name]] <- as.numeric(scale(raw_val))
        design_terms <- c(design_terms, scaled_name)
        log_info("Added numeric covariate: %s (standardized)", nv)
      } else {
        log_warn("Numeric covariate '%s' has zero variance. Skipping.", nv)
      }
    } else {
      log_warn("Numeric covariate '%s' not found in metadata. Skipping.", nv)
    }
  }

  # --- Add age if it has variance ---
  if (any(meta_df$internal_age_scaled != 0)) {
    design_terms <- c(design_terms, "internal_age_scaled")
  }

  # --- Always add sex last (main variable of interest) ---
  design_terms <- unique(c(setdiff(design_terms, sex_col), sex_col))
  if (length(design_terms) == 0) design_terms <- sex_col

  formula_str <- paste0("~ ", paste(design_terms, collapse = " + "))
  log_info("Design formula: %s", formula_str)

  # --- Update metadata in object ---
  SummarizedExperiment::colData(pb_sce) <- meta_df

  # --- Create DESeqDataSet ---
  log_info("Creating DESeqDataSet...")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(pb_sce, "counts"),
    colData   = SummarizedExperiment::colData(pb_sce),
    design    = as.formula(formula_str)
  )

  # Filter low count genes
  n_before <- nrow(dds)
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_gene_count, ]
  n_after <- nrow(dds)
  log_info("Gene filter: %d -> %d (removed %d low-count genes)",
           n_before, n_after, n_before - n_after)

  # --- Run DESeq2 ---
  log_info("Running DESeq2 (this may take a while)...")
  dds <- suppressMessages(DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = de_param))

  # --- LFC Shrinkage ---
  res_names <- DESeq2::resultsNames(dds)
  log_info("Available coefficients: %s", paste(res_names, collapse = ", "))

  target_coef <- res_names[grepl(paste0("^", sex_col), res_names)]
  if (length(target_coef) == 0) {
    stop(sprintf("Could not find sex coefficient. Looking for '%s' in: %s",
                 sex_col, paste(res_names, collapse = ", ")))
  }
  target_coef <- target_coef[1]

  log_info("Applying apeglm shrinkage on: %s", target_coef)
  res <- DESeq2::lfcShrink(dds, coef = target_coef, type = "apeglm",
                            parallel = TRUE, BPPARAM = de_param)

  # --- Convert to data.frame ---
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$pvalue), ]

  n_sig <- sum(res_df$padj < FDR_THRESHOLD, na.rm = TRUE)
  log_info("Total genes tested: %d, Significant (padj < %g): %d",
           nrow(res_df), FDR_THRESHOLD, n_sig)

  return(list(dds = dds, results = res, res_df = res_df))
}

# ==============================================================================
# 3. GENE ID HELPERS
# ==============================================================================

#' Strip version numbers from Ensembl IDs
#'
#' @param x Character vector of gene IDs (e.g., "ENSG00000229807.5")
#' @return Character vector with versions removed (e.g., "ENSG00000229807")
clean_gene_ids <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*", "", x)
  return(trimws(x))
}

# ==============================================================================
# 5. CROSS-DATASET DIAGNOSTIC HELPERS
# ==============================================================================

#' Load DE results for all datasets in a cell type directory
#'
#' @param output_dir Path to cell type output directory
#' @param suffix_pattern File suffix pattern for DE results
#' @return Named list of data.frames (one per dataset)
load_all_de_results <- function(output_dir, suffix_pattern = "_DE_results.csv") {
  files <- list.files(
    path = output_dir,
    pattern = paste0(".*", suffix_pattern, "$"),
    full.names = TRUE
  )

  all_results <- list()

  for (f in files) {
    sample_name <- gsub(suffix_pattern, "", basename(f))

    df <- tryCatch({
      read.csv(f, stringsAsFactors = FALSE)
    }, error = function(e) {
      log_warn("Failed to read %s: %s", basename(f), e$message)
      NULL
    })

    if (is.null(df)) next

    required_cols <- c("gene", "padj", "log2FoldChange")
    if (!all(required_cols %in% colnames(df))) {
      log_warn("File %s missing required columns. Skipping.", basename(f))
      next
    }

    # Clean gene IDs and set as row names
    df$gene_id <- clean_gene_ids(df$gene)
    df <- df[!duplicated(df$gene_id), ]
    rownames(df) <- df$gene_id

    all_results[[sample_name]] <- df
    log_info("Loaded %s (%d genes)", sample_name, nrow(df))
  }

  return(all_results)
}

#' Extract significant genes from DE results
#'
#' @param all_results Named list of DE result data.frames
#' @param padj_cutoff Adjusted p-value threshold
#' @param lfc_cutoff Absolute log2 fold change threshold
#' @return Named list of character vectors (significant gene IDs per dataset)
extract_significant_genes <- function(all_results,
                                      padj_cutoff = FDR_THRESHOLD,
                                      lfc_cutoff = LFC_THRESHOLD) {
  sig_list <- list()
  for (name in names(all_results)) {
    df <- all_results[[name]]
    sig_df <- df[!is.na(df$padj) & df$padj < padj_cutoff &
                   abs(df$log2FoldChange) > lfc_cutoff, ]
    sig_list[[name]] <- rownames(sig_df)
  }
  return(sig_list)
}

#' Build pairwise overlap matrix of significant genes
#'
#' @param sig_list Named list of character vectors
#' @return Numeric matrix of pairwise overlaps
build_overlap_matrix <- function(sig_list) {
  n <- length(sig_list)
  mat <- matrix(0, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- names(sig_list)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- length(intersect(sig_list[[i]], sig_list[[j]]))
    }
  }
  return(mat)
}

#' Build LFC matrix for a set of genes across datasets
#'
#' @param gene_ids Character vector of gene IDs
#' @param all_results Named list of DE result data.frames
#' @return Numeric matrix (genes x datasets)
build_lfc_matrix <- function(gene_ids, all_results) {
  mat <- matrix(NA, nrow = length(gene_ids), ncol = length(all_results))
  rownames(mat) <- gene_ids
  colnames(mat) <- names(all_results)

  for (d in names(all_results)) {
    df <- all_results[[d]]
    matched <- intersect(gene_ids, rownames(df))
    mat[matched, d] <- df[matched, "log2FoldChange"]
  }

  return(mat)
}

#' Compute concordance analysis vs reference dataset
#'
#' @param sig_list Named list of significant gene vectors
#' @param reference Name of reference dataset
#' @return data.frame with concordance metrics
compute_concordance <- function(sig_list, reference = REFERENCE_DATASET) {
  if (!reference %in% names(sig_list)) {
    log_warn("Reference dataset '%s' not found in results", reference)
    return(data.frame())
  }

  ref_sig <- sig_list[[reference]]
  results <- data.frame()

  for (dataset in names(sig_list)) {
    if (dataset == reference) next

    ds_sig <- sig_list[[dataset]]
    shared <- intersect(ref_sig, ds_sig)

    results <- rbind(results, data.frame(
      Dataset = dataset,
      Reference_Unique = length(setdiff(ref_sig, ds_sig)),
      Dataset_Unique = length(setdiff(ds_sig, ref_sig)),
      Shared = length(shared),
      Concordance_Rate = if (length(ref_sig) > 0) {
        length(shared) / length(ref_sig)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

#' Build sign agreement matrix
#'
#' @param lfc_matrix Numeric matrix (genes x datasets)
#' @return Numeric matrix (datasets x datasets) of sign agreement proportions
build_sign_agreement <- function(lfc_matrix) {
  sign_mat <- sign(lfc_matrix)
  n <- ncol(sign_mat)
  agreement <- matrix(NA, nrow = n, ncol = n)
  rownames(agreement) <- colnames(agreement) <- colnames(sign_mat)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      valid <- !is.na(sign_mat[, i]) & !is.na(sign_mat[, j])
      if (sum(valid) > 0) {
        agreement[i, j] <- mean(sign_mat[valid, i] == sign_mat[valid, j])
      }
    }
  }

  return(agreement)
}

# ==============================================================================
# 6. EXPANDED LFC AGREEMENT MATRIX HELPERS
# ==============================================================================

#' Compute MI-optimal flip counts for sex label permutation
#'
#' Given n_male true males and n_female true females, computes the number of
#' true males that should retain the "male" label such that the assigned labels
#' are statistically independent of true sex (MI ≈ 0).
#'
#' Formula: a* = round(n_male^2 / N), where N = n_male + n_female.
#' This ensures P(assigned_male | true_male) = P(assigned_male) = n_male/N.
#'
#' @param n_male Integer, number of true males
#' @param n_female Integer, number of true females
#' @return List with: a_star (males keeping "male"), n_flip_males (males -> "female"),
#'         n_flip_females (females -> "male"), mi_deviation (should be near 0)
compute_mi_optimal_flip <- function(n_male, n_female) {
  n_total <- n_male + n_female

  if (n_total == 0) stop("No donors provided")
  if (n_male == 0 || n_female == 0) {
    log_warn("Only one sex present (M=%d, F=%d). Cannot flip labels.", n_male, n_female)
    return(list(a_star = n_male, n_flip_males = 0, n_flip_females = 0,
                mi_deviation = NA_real_))
  }

  # Optimal number of true males to keep as "male"
  a_star <- round(n_male^2 / n_total)

  # Number of males to flip to "female"
  n_flip_males <- n_male - a_star

  # Number of females to flip to "male" (to maintain marginal ratio)
  n_flip_females <- n_male - a_star

  # Deviation from perfect independence (should be near 0)
  p_am_given_tm <- a_star / n_male
  p_am <- n_male / n_total
  mi_deviation <- abs(p_am_given_tm - p_am)

  log_info("MI-optimal flip: %d males keep 'male' (of %d), flip %d males and %d females",
           a_star, n_male, n_flip_males, n_flip_females)
  log_info("  P(assigned_M | true_M) = %.4f, P(assigned_M) = %.4f, |diff| = %.6f",
           p_am_given_tm, p_am, mi_deviation)

  return(list(
    a_star         = a_star,
    n_flip_males   = n_flip_males,
    n_flip_females = n_flip_females,
    mi_deviation   = mi_deviation
  ))
}

#' Apply MI-minimised sex label flip to donor metadata
#'
#' Flips sex labels for a subset of donors such that the assigned labels are
#' approximately independent of true sex (MI ≈ 0). The marginal sex ratio is
#' preserved.
#'
#' @param meta_df Data.frame with at least a sex column
#' @param sex_col Name of the sex column
#' @param male_val Value representing male in the sex column
#' @param female_val Value representing female in the sex column
#' @param seed Random seed for selecting which donors to flip
#' @return Data.frame with modified sex column (original sex saved as 'true_sex')
flip_sex_labels_mi_minimised <- function(meta_df, sex_col, male_val, female_val, seed) {

  # Save original sex
  meta_df$true_sex <- meta_df[[sex_col]]

  male_idx   <- which(meta_df[[sex_col]] == male_val)
  female_idx <- which(meta_df[[sex_col]] == female_val)

  n_male   <- length(male_idx)
  n_female <- length(female_idx)

  # Compute optimal flip
  flip_info <- compute_mi_optimal_flip(n_male, n_female)

  if (flip_info$n_flip_males == 0 && flip_info$n_flip_females == 0) {
    log_warn("No flipping needed (single sex). Returning unchanged.")
    return(meta_df)
  }

  # Randomly select which males to flip and which females to flip
  set.seed(seed)
  males_to_flip   <- sample(male_idx, flip_info$n_flip_males)
  females_to_flip <- sample(female_idx, flip_info$n_flip_females)

  # Apply flips
  meta_df[[sex_col]][males_to_flip]   <- female_val
  meta_df[[sex_col]][females_to_flip] <- male_val

  # Verify marginals preserved
  new_n_male   <- sum(meta_df[[sex_col]] == male_val)
  new_n_female <- sum(meta_df[[sex_col]] == female_val)
  log_info("After flip: %d assigned male, %d assigned female (was %d M, %d F)",
           new_n_male, new_n_female, n_male, n_female)

  if (new_n_male != n_male || new_n_female != n_female) {
    log_error("Marginal sex ratio NOT preserved! Something went wrong.")
  }

  return(meta_df)
}

#' Run DESeq2 on a subset of Roussos donors (for upper/lower bound)
#'
#' Aggregates shards for specified donor IDs, optionally applies a label flip
#' function, then runs DESeq2 using the Roussos dataset configuration.
#'
#' @param donor_ids Character vector of donor IDs to include
#' @param cell_type Cell type name (e.g., "Astro")
#' @param label_flip_fn Optional function(meta_df) -> meta_df that modifies sex labels
#' @param de_param BiocParallel parameter for DESeq2
#' @param cache_path Optional path to cache/load the DESeq2 results (.qs)
#' @param force_recalc Force recalculation even if cache exists
#' @return Data.frame of DE results (gene, log2FoldChange, padj, etc.)
run_deseq_for_roussos_subset <- function(donor_ids, cell_type,
                                          label_flip_fn = NULL,
                                          de_param = NULL,
                                          cache_path = NULL,
                                          force_recalc = FALSE) {

  # Check cache
  if (!is.null(cache_path) && file.exists(cache_path) && !force_recalc) {
    log_info("[CACHE] Loading cached DE results: %s", basename(cache_path))
    cached <- qs::qread(cache_path)
    if (is.data.frame(cached) && "log2FoldChange" %in% colnames(cached)) {
      return(cached)
    }
    log_warn("Cached file invalid, regenerating...")
  }

  cfg <- DATASET_CONFIG[["Roussos"]]
  shard_dir <- get_shard_dir(cell_type, "Roussos")

  if (!dir.exists(shard_dir)) {
    log_error("Shard directory not found: %s", shard_dir)
    return(NULL)
  }

  # Aggregate shards for the specified donors
  vars_requested <- unique(c(
    cfg$cols$donor, cfg$cols$sex, cfg$cols$age,
    cfg$covariates, cfg$numeric_covariates
  ))

  log_info("Aggregating %d Roussos donors for %s...", length(donor_ids), cell_type)

  # Filter shard files to only include requested donors
  all_shard_files <- list.files(shard_dir, pattern = "\\.qs$", full.names = TRUE)
  all_shard_files <- all_shard_files[!grepl("MANIFEST", all_shard_files, ignore.case = TRUE)]
  all_shard_files <- all_shard_files[!grepl("donor_NA\\.qs$", all_shard_files)]

  # Extract donor IDs from filenames and filter
  shard_donor_ids <- gsub("^donor_", "", gsub("\\.qs$", "", basename(all_shard_files)))
  keep_mask <- shard_donor_ids %in% donor_ids
  target_files <- all_shard_files[keep_mask]

  if (length(target_files) == 0) {
    log_error("No matching shard files found for %d donor IDs", length(donor_ids))
    return(NULL)
  }

  log_info("Found %d / %d shard files", length(target_files), length(donor_ids))

  # Create a temporary directory with symlinks to only the target shards
  tmp_shard_dir <- tempfile(pattern = "gt_subset_")
  dir.create(tmp_shard_dir)
  on.exit(unlink(tmp_shard_dir, recursive = TRUE), add = TRUE)

  for (f in target_files) {
    file.symlink(f, file.path(tmp_shard_dir, basename(f)))
  }

  # Aggregate
  pb_sce <- aggregate_shards_parallel(
    shard_dir   = tmp_shard_dir,
    design_vars = vars_requested,
    n_cores     = min(N_CORES_AGG, 6)
  )

  # Apply label flip if provided
  if (!is.null(label_flip_fn)) {
    log_info("Applying label flip function...")
    modified_meta <- label_flip_fn(as.data.frame(SummarizedExperiment::colData(pb_sce)))
    SummarizedExperiment::colData(pb_sce) <- S4Vectors::DataFrame(modified_meta)
  }

  # Setup DE param if not provided
  if (is.null(de_param)) {
    de_param <- setup_parallel(min(N_CORES_DE, 6))
  }

  # Run DESeq2
  log_info("Running DESeq2 on %d donors, %d genes...", ncol(pb_sce), nrow(pb_sce))
  de_result <- run_deseq_for_dataset(
    pb_sce         = pb_sce,
    cfg            = cfg,
    de_param       = de_param,
    min_gene_count = MIN_GENE_COUNT
  )

  res_df <- de_result$res_df

  # Clean gene IDs and set as row names (consistent with load_all_de_results)
  res_df$gene_id <- clean_gene_ids(res_df$gene)
  res_df <- res_df[!duplicated(res_df$gene_id), ]
  rownames(res_df) <- res_df$gene_id

  # Cache if path provided
  if (!is.null(cache_path)) {
    cache_dir <- dirname(cache_path)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    qs::qsave(res_df, cache_path)
    log_info("[CACHE] Saved DE results: %s", basename(cache_path))
  }

  # Clean up
  rm(pb_sce, de_result)
  gc(verbose = FALSE)

  return(res_df)
}

#' Get donor IDs from cartesian downsampling run metadata
#'
#' Reads the Ground_Truth_sampling_metadata.qs and Validation_Donors.qs from
#' a cartesian downsampling run to get the GT and validation donor IDs.
#'
#' @param cell_type Cell type name
#' @param run_id Run number (1-6)
#' @return List with gt_donors and val_donors character vectors
get_cartesian_split_donors <- function(cell_type, run_id) {
  run_dir <- file.path(CARTESIAN_BASE, cell_type, CARTESIAN_DATASET,
                        CARTESIAN_MODE, paste0("Run_", run_id))

  if (!dir.exists(run_dir)) {
    log_error("Cartesian run directory not found: %s", run_dir)
    return(NULL)
  }

  # Load validation donors
  val_path <- file.path(run_dir, "Validation_Donors.qs")
  if (!file.exists(val_path)) {
    log_error("Validation donors file not found: %s", val_path)
    return(NULL)
  }

  val_meta <- qs::qread(val_path)
  val_donors <- as.character(val_meta$donor_id)

  # Load GT sampling metadata to get GT donor IDs
  gt_log_path <- file.path(run_dir, "Ground_Truth_sampling_metadata.qs")
  if (!file.exists(gt_log_path)) {
    log_error("GT sampling metadata not found: %s", gt_log_path)
    return(NULL)
  }

  gt_log <- qs::qread(gt_log_path)
  gt_donors <- names(gt_log)

  log_info("Run %d: %d GT donors, %d validation donors",
           run_id, length(gt_donors), length(val_donors))

  return(list(gt_donors = gt_donors, val_donors = val_donors))
}

#' Load DE results from a cartesian Ground_Truth.qs file
#'
#' Reads the DESeqResults object from a cartesian run's Ground_Truth.qs and
#' converts it to a data.frame consistent with load_all_de_results() format.
#'
#' @param cell_type Cell type name
#' @param run_id Run number
#' @return Data.frame with gene_id, log2FoldChange, padj, etc.
load_cartesian_gt_results <- function(cell_type, run_id) {
  
  # --- Defensive Check 1: Validate CARTESIAN_BASE exists in environment ---
  if (!exists("CARTESIAN_BASE")) {
    stop(sprintf(
      "[ERROR] CARTESIAN_BASE not found in environment. Did you source config.R?\n  Expected config variable: CARTESIAN_BASE\n  Call: load_cartesian_gt_results(cell_type = '%s', run_id = %d)",
      cell_type, run_id
    ))
  }
  
  # --- Defensive Check 2: Validate CARTESIAN_BASE directory exists ---
  if (!dir.exists(CARTESIAN_BASE)) {
    stop(sprintf(
      "[ERROR] CARTESIAN_BASE directory does not exist on disk:\n  Path: %s\n  Check config.R and verify the cartesian downsampling pipeline has run.\n  Call: load_cartesian_gt_results(cell_type = '%s', run_id = %d)",
      CARTESIAN_BASE, cell_type, run_id
    ))
  }
  
  # --- Construct file path ---
  gt_path <- file.path(CARTESIAN_BASE, cell_type, CARTESIAN_DATASET,
                        CARTESIAN_MODE, paste0("Run_", run_id), "Ground_Truth.qs")

  # --- Defensive Check 3: Validate file exists ---
  if (!file.exists(gt_path)) {
    log_error("Ground truth file not found: %s", gt_path)
    return(NULL)
  }

  # --- Defensive Check 4: Try reading with explicit error handling ---
  res <- tryCatch({
    qs::qread(gt_path)
  }, error = function(e) {
    stop(sprintf(
      "[ERROR] Failed to read Ground_Truth.qs file:\n  Path: %s\n  Error: %s\n  Possible causes:\n    - File corruption (try regenerating with cartesian pipeline Step 1)\n    - qs package version mismatch (check qs::qread() version)\n    - File permissions (check read access)\n  Call: load_cartesian_gt_results(cell_type = '%s', run_id = %d)",
      gt_path, e$message, cell_type, run_id
    ))
  })
  
  # --- Defensive Check 5: Validate object class ---
  if (!"DESeqResults" %in% class(res)) {
    stop(sprintf(
      "[ERROR] Ground_Truth.qs did not contain a DESeqResults object:\n  Path: %s\n  Actual class: %s\n  Expected: DESeqResults\n  This suggests the file was created by a different pipeline or is corrupted.\n  Call: load_cartesian_gt_results(cell_type = '%s', run_id = %d)",
      gt_path, paste(class(res), collapse = ", "), cell_type, run_id
    ))
  }

  # Convert DESeqResults to data.frame
  res_df <- as.data.frame(res)
  
  # --- Defensive Check 6: Validate required columns ---
  required_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
  missing_cols <- setdiff(required_cols, colnames(res_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "[ERROR] Ground_Truth DESeqResults is missing required columns:\n  Path: %s\n  Missing: %s\n  Present: %s\n  Call: load_cartesian_gt_results(cell_type = '%s', run_id = %d)",
      gt_path, paste(missing_cols, collapse = ", "),
      paste(colnames(res_df), collapse = ", "), cell_type, run_id
    ))
  }
  
  res_df$gene <- rownames(res_df)
  res_df$gene_id <- clean_gene_ids(res_df$gene)
  res_df <- res_df[!duplicated(res_df$gene_id), ]
  rownames(res_df) <- res_df$gene_id
  
  # --- Defensive Check 7: Validate non-zero genes ---
  if (nrow(res_df) == 0) {
    log_warn("Ground truth file loaded but contains 0 genes after processing: %s", gt_path)
    return(NULL)
  }

  log_info("Loaded cartesian GT Run %d: %d genes", run_id, nrow(res_df))
  return(res_df)
}

#' Get all Roussos donor IDs for a cell type
#'
#' Fast path: reads donor IDs from the Roussos pseudobulk checkpoint (single
#' file, already computed in Step 1). Falls back to listing shard filenames if
#' the checkpoint does not exist.
#'
#' @param cell_type Cell type name
#' @return Character vector of donor IDs
get_all_roussos_donors <- function(cell_type) {

  # --- Fast path: read from pseudobulk checkpoint ---
  pb_path <- file.path(OUTPUT_BASE, cell_type, "Roussos_pb_checkpoint.qs")
  if (file.exists(pb_path)) {
    tryCatch({
      pb_sce <- qs::qread(pb_path)
      cfg <- DATASET_CONFIG[["Roussos"]]
      ids <- as.character(
        SummarizedExperiment::colData(pb_sce)[[cfg$cols$donor]]
      )
      log_info("Found %d Roussos donors for %s (from pb checkpoint)",
               length(ids), cell_type)
      return(ids)
    }, error = function(e) {
      log_warn("Failed to read pb checkpoint, falling back to shard listing: %s",
               e$message)
    })
  }

  # --- Fallback: list shard filenames ---
  shard_dir <- get_shard_dir(cell_type, "Roussos")

  if (!dir.exists(shard_dir)) {
    log_error("Roussos shard directory not found: %s", shard_dir)
    return(NULL)
  }

  shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = FALSE)
  shard_files <- shard_files[!grepl("MANIFEST", shard_files, ignore.case = TRUE)]
  shard_files <- shard_files[!grepl("donor_NA\\.qs$", shard_files)]

  donor_ids <- gsub("^donor_", "", gsub("\\.qs$", "", shard_files))
  log_info("Found %d Roussos donors for %s (from shard listing)",
           length(donor_ids), cell_type)

  return(donor_ids)
}

#' Get Roussos donor sex from metadata
#'
#' Fast path: reads donor_id and sex from the Roussos pseudobulk checkpoint
#' (single file, already computed in Step 1). Falls back to reading individual
#' shard files if the checkpoint does not exist.
#'
#' @param cell_type Cell type name
#' @param donor_ids Character vector of donor IDs (if NULL, uses all)
#' @return Data.frame with donor_id and sex columns
get_roussos_donor_sex <- function(cell_type, donor_ids = NULL) {

  # --- Fast path: read from pseudobulk checkpoint ---
  pb_path <- file.path(OUTPUT_BASE, cell_type, "Roussos_pb_checkpoint.qs")
  if (file.exists(pb_path)) {
    tryCatch({
      pb_sce <- qs::qread(pb_path)
      cfg <- DATASET_CONFIG[["Roussos"]]
      cd <- as.data.frame(SummarizedExperiment::colData(pb_sce))
      sex_df <- data.frame(
        donor_id = as.character(cd[[cfg$cols$donor]]),
        sex      = as.character(cd[[cfg$cols$sex]]),
        stringsAsFactors = FALSE
      )
      if (!is.null(donor_ids)) {
        sex_df <- sex_df[sex_df$donor_id %in% donor_ids, ]
      }
      log_info("Got sex for %d Roussos donors in %s (from pb checkpoint)",
               nrow(sex_df), cell_type)
      return(sex_df)
    }, error = function(e) {
      log_warn("Failed to read pb checkpoint, falling back to shard reading: %s",
               e$message)
    })
  }

  # --- Fallback: read individual shard files ---
  shard_dir <- get_shard_dir(cell_type, "Roussos")

  if (is.null(donor_ids)) {
    donor_ids <- get_all_roussos_donors(cell_type)
  }

  sex_df <- do.call(rbind, lapply(donor_ids, function(id) {
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))

    if (!file.exists(shard_path)) return(NULL)

    tryCatch({
      sce <- qs::qread(shard_path)
      cd <- SummarizedExperiment::colData(sce)
      data.frame(
        donor_id = id,
        sex = as.character(cd[["sex"]][1]),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      log_warn("Failed to read sex for donor %s: %s", id, e$message)
      NULL
    })
  }))

  return(sex_df)
}

# ==============================================================================
# 7. CHECKPOINT / CACHE HELPERS
# ==============================================================================

#' Validate a cached pseudobulk SCE object
#'
#' @param pb_sce SingleCellExperiment object
#' @param cfg Dataset configuration
#' @return TRUE if valid, FALSE if should regenerate
validate_cached_pseudobulk <- function(pb_sce, cfg) {
  counts <- SummarizedExperiment::assay(pb_sce, "counts")

  # Check 1: Counts are integers
  if (!all(counts == floor(counts), na.rm = TRUE)) {
    non_int_pct <- mean(counts != floor(counts), na.rm = TRUE) * 100
    log_warn("Cache invalid: %.2f%% non-integer counts", non_int_pct)
    return(FALSE)
  }

  # Check 2: Non-empty
  if (ncol(pb_sce) == 0 || nrow(pb_sce) == 0) {
    log_warn("Cache invalid: empty pseudobulk object")
    return(FALSE)
  }

  # Check 3: Required columns exist
  required <- c(cfg$cols$donor, cfg$cols$sex)
  if (!all(required %in% colnames(SummarizedExperiment::colData(pb_sce)))) {
    log_warn("Cache invalid: missing required metadata columns")
    return(FALSE)
  }

  return(TRUE)
}

#' Save results in both .qs and .csv formats
#'
#' @param res_df Data.frame of DE results
#' @param dds DESeqDataSet object
#' @param res DESeqResults object
#' @param pb_sce Pseudobulk SCE object
#' @param output_dir Output directory
#' @param dataset_name Dataset name prefix
#' @param cell_type Cell type name
save_de_outputs <- function(res_df, dds, res, pb_sce, output_dir,
                            dataset_name, cell_type) {
  prefix <- paste0(dataset_name)

  # CSV results
  csv_path <- file.path(output_dir, paste0(prefix, "_DE_results.csv"))
  write.csv(res_df, csv_path, row.names = FALSE, quote = TRUE)
  log_info("Saved DE results: %s", basename(csv_path))

  # Full DESeq2 results object
  qs::qsave(res, file.path(output_dir, paste0(prefix, "_Full_DESeq2_Results.qs")))

  # DDS checkpoint
  qs::qsave(dds, file.path(output_dir, paste0(prefix, "_dds_checkpoint.qs")))

  # Pseudobulk checkpoint
  qs::qsave(pb_sce, file.path(output_dir, paste0(prefix, "_pb_checkpoint.qs")))

  log_info("Saved all outputs for %s/%s", dataset_name, cell_type)
}

# ==============================================================================
# 8. GENE AUDIT MATRIX HELPERS (Step 04)
# ==============================================================================

#' Extract gene IDs from all shards in a directory (union)
#'
#' Loads every donor shard .qs file in the directory and takes the union of
#' all rownames (gene IDs). This handles edge cases where shards might have
#' slightly different gene sets.
#'
#' @param shard_dir Path to directory containing donor_*.qs files
#' @return Character vector of unique gene IDs, or NULL if no shards found
extract_shard_gene_ids <- function(shard_dir) {

  if (!dir.exists(shard_dir)) {
    log_warn("Shard directory not found: %s", shard_dir)
    return(NULL)
  }

  shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  shard_files <- shard_files[!grepl("MANIFEST", shard_files, ignore.case = TRUE)]
  shard_files <- shard_files[!grepl("donor_NA\\.qs$", shard_files)]

  if (length(shard_files) == 0) {
    log_warn("No shard files found in: %s", shard_dir)
    return(NULL)
  }

  all_genes <- character(0)

  for (f in shard_files) {
    tryCatch({
      sce <- qs::qread(f)
      # Try SCE rownames first; fall back to assay matrix rownames if NULL
      gene_ids_local <- rownames(sce)
      if (is.null(gene_ids_local) || length(gene_ids_local) == 0) {
        gene_ids_local <- rownames(SummarizedExperiment::assay(sce, 1))
      }
      all_genes <- union(all_genes, gene_ids_local)
      rm(sce)
    }, error = function(e) {
      log_warn("Failed to read shard %s: %s", basename(f), e$message)
    })
  }

  gc(verbose = FALSE)

  if (length(all_genes) == 0) {
    log_warn("No gene IDs extracted from shards in: %s", shard_dir)
    return(NULL)
  }

  log_info("  Extracted %d unique genes from %d shards in %s",
           length(all_genes), length(shard_files), basename(shard_dir))
  return(all_genes)
}

#' Query BioMart for gene annotations including chromosome band
#'
#' Fetches ensembl_gene_id, chromosome_name, gene_biotype, and band from
#' BioMart. Results are cached to a CSV file for subsequent runs.
#'
#' @param gene_ids Character vector of Ensembl gene IDs (no version suffix)
#' @param cache_path Path to cache CSV file
#' @return data.frame with columns: ensembl_gene_id, chromosome_name,
#'         gene_biotype, band
get_audit_annotations <- function(gene_ids, cache_path) {

  gene_ids <- unique(gene_ids)

  # Check cache first
  if (file.exists(cache_path)) {
    log_info("Loading audit annotation cache from: %s", basename(cache_path))
    cached <- read.csv(cache_path, stringsAsFactors = FALSE)

    # Check coverage
    missing <- setdiff(gene_ids, cached$ensembl_gene_id)
    if (length(missing) == 0) {
      log_info("  Cache covers all %d genes", length(gene_ids))
      return(cached[cached$ensembl_gene_id %in% gene_ids, ])
    }
    log_info("  Cache covers %d/%d genes, querying %d missing",
             length(gene_ids) - length(missing), length(gene_ids), length(missing))
    query_ids <- missing
  } else {
    log_info("No audit cache found, querying BioMart for %d genes", length(gene_ids))
    cached <- NULL
    query_ids <- gene_ids
  }

  # Query BioMart
  bm_result <- query_biomart_audit(query_ids)

  if (!is.null(bm_result) && nrow(bm_result) > 0) {
    # Merge with cache
    if (!is.null(cached)) {
      all_annotations <- rbind(cached, bm_result)
      all_annotations <- all_annotations[!duplicated(all_annotations$ensembl_gene_id), ]
    } else {
      all_annotations <- bm_result
    }

    # Save updated cache
    write.csv(all_annotations, cache_path, row.names = FALSE)
    log_info("  Updated audit cache: %d total genes", nrow(all_annotations))

    return(all_annotations[all_annotations$ensembl_gene_id %in% gene_ids, ])
  } else {
    log_warn("BioMart audit query failed or returned no results")
    if (!is.null(cached)) {
      return(cached[cached$ensembl_gene_id %in% gene_ids, ])
    }
    return(NULL)
  }
}

#' Query BioMart for chromosome, biotype, and band
#'
#' @param gene_ids Character vector of Ensembl gene IDs
#' @return data.frame with ensembl_gene_id, chromosome_name, gene_biotype, band
query_biomart_audit <- function(gene_ids) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    log_warn("biomaRt package not available")
    return(NULL)
  }

  # Try primary host
  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = AUDIT_BIOMART_DATASET,
                     host = AUDIT_BIOMART_HOST)
  }, error = function(e) {
    log_warn("BioMart primary host failed: %s", e$message)
    tryCatch({
      biomaRt::useMart("ensembl", dataset = AUDIT_BIOMART_DATASET,
                       host = AUDIT_BIOMART_HOST_FALLBACK)
    }, error = function(e2) {
      log_warn("BioMart fallback host also failed: %s", e2$message)
      return(NULL)
    })
  })

  if (is.null(mart)) return(NULL)

  tryCatch({
    result <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "chromosome_name", "gene_biotype", "band"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )

    # Deduplicate: some genes map to multiple bands (e.g., spanning regions)
    # Keep first occurrence per gene
    result <- result[!duplicated(result$ensembl_gene_id), ]

    log_info("  BioMart returned annotations for %d/%d genes (with band)",
             nrow(result), length(gene_ids))
    return(result)
  }, error = function(e) {
    log_warn("BioMart audit query failed: %s", e$message)
    return(NULL)
  })
}

#' Classify genes into overlapping categories
#'
#' Takes a BioMart annotation data.frame and returns a data.frame with boolean
#' columns for each gene category. Categories are NOT mutually exclusive — a
#' gene can be e.g., both X-chromosome and pseudogene.
#'
#' @param annotation_df data.frame from get_audit_annotations() with columns:
#'        ensembl_gene_id, chromosome_name, gene_biotype, band
#' @return data.frame with ensembl_gene_id and boolean category columns
classify_gene_categories <- function(annotation_df) {

  if (is.null(annotation_df) || nrow(annotation_df) == 0) {
    log_warn("No annotations to classify")
    return(NULL)
  }

  result <- data.frame(
    ensembl_gene_id = annotation_df$ensembl_gene_id,
    stringsAsFactors = FALSE
  )

  chr <- annotation_df$chromosome_name
  biotype <- annotation_df$gene_biotype
  band <- annotation_df$band

  # Chromosome-based categories

  result$is_X_chr <- chr == "X"
  result$is_Y_chr <- chr == "Y"
  result$is_autosomal <- chr %in% AUTOSOMAL_CHROMOSOMES
  result$is_mitochondrial <- chr == "MT"

  # PAR identification via band annotation

  # PAR1: Xp22.33 or Yp11.32
  # PAR2: Xq28 or Yq12
  # Note: band matching uses startsWith to handle sub-bands
  is_par1 <- (chr == "X" & !is.na(band) & startsWith(band, PAR1_BANDS_X)) |
             (chr == "Y" & !is.na(band) & startsWith(band, PAR1_BANDS_Y))
  is_par2 <- (chr == "X" & !is.na(band) & startsWith(band, PAR2_BANDS_X)) |
             (chr == "Y" & !is.na(band) & startsWith(band, PAR2_BANDS_Y))
  result$is_PAR <- is_par1 | is_par2

  # Biotype-based categories
  result$is_protein_coding <- !is.na(biotype) & biotype == "protein_coding"
  result$is_pseudogene <- !is.na(biotype) & biotype %in% PSEUDOGENE_BIOTYPES
  result$is_lncRNA <- !is.na(biotype) & biotype == "lncRNA"

  return(result)
}

#' Build audit matrix from gene lists and category classifications
#'
#' Takes a named list of gene ID vectors (one per dataset) and the category
#' classification, returns a data.frame (rows = datasets, columns = categories).
#'
#' @param gene_lists Named list of character vectors (gene IDs per dataset)
#' @param category_df data.frame from classify_gene_categories()
#' @return data.frame with Dataset, Total, and one column per category
build_audit_matrix <- function(gene_lists, category_df) {

  if (is.null(category_df) || nrow(category_df) == 0) {
    log_warn("No category data available for audit matrix")
    return(NULL)
  }

  # Category columns (all boolean columns except ensembl_gene_id)
  cat_cols <- setdiff(colnames(category_df), "ensembl_gene_id")

  # Clean column names for output (remove is_ prefix)
  clean_names <- gsub("^is_", "", cat_cols)

  rows <- list()

  for (ds_name in names(gene_lists)) {
    ds_genes <- gene_lists[[ds_name]]

    # Match to category data
    matched <- category_df[category_df$ensembl_gene_id %in% ds_genes, ]

    # Count per category
    counts <- sapply(cat_cols, function(col) {
      sum(matched[[col]], na.rm = TRUE)
    })
    names(counts) <- clean_names

    # Also count genes NOT in BioMart (unannotated)
    n_unannotated <- length(setdiff(ds_genes, category_df$ensembl_gene_id))

    row_data <- c(
      Dataset = ds_name,
      Total = length(ds_genes),
      counts,
      Unannotated = n_unannotated
    )

    rows[[ds_name]] <- row_data
  }

  result <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)

  # Convert numeric columns from character
  numeric_cols <- setdiff(colnames(result), "Dataset")
  for (col in numeric_cols) {
    result[[col]] <- as.integer(result[[col]])
  }

  rownames(result) <- NULL
  return(result)
}

# ==============================================================================
# 7. INTERSECTION OVER UNION (IoU) ANALYSIS
# ==============================================================================

#' Compute Intersection over Union (IoU) for two gene sets
#'
#' IoU = |Intersection| / |Union| = |A ∩ B| / |A ∪ B|
#' This is equivalent to the Jaccard index for binary sets.
#'
#' @param set1_sig Character vector of significant genes in set 1
#' @param set2_sig Character vector of significant genes in set 2
#' @return IoU value in [0, 1], or NA if union is empty
compute_iou <- function(set1_sig, set2_sig) {
  intersect_size <- length(intersect(set1_sig, set2_sig))
  union_size <- length(union(set1_sig, set2_sig))
  
  if (union_size == 0) {
    return(NA_real_)
  }
  
  return(intersect_size / union_size)
}

#' Build IoU matrix for multiple result sets
#'
#' Computes pairwise IoU for all combinations of result sets,
#' using the first result set as the universe for gene selection.
#'
#' @param combined_results Named list of DE result data frames
#' @param ref_name Name of reference dataset (for universe selection)
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param lfc_threshold LFC threshold for significance (default: 1.0)
#' @return Named list containing:
#'   - iou_matrix: Matrix of IoU values
#'   - summary_df: Data frame with detailed statistics per pair
build_iou_matrix <- function(combined_results, ref_name = REFERENCE_DATASET,
                              fdr_threshold = 0.05, lfc_threshold = 1.0) {
  
  n_sets <- length(combined_results)
  set_names <- names(combined_results)
  
  if (n_sets == 0) {
    return(list(iou_matrix = NULL, summary_df = NULL))
  }
  
  # Extract significant genes for each set
  sig_genes_list <- lapply(combined_results, function(df) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    
    # Get gene IDs from rownames if not in a column
    if ("gene" %in% colnames(df)) {
      genes <- df$gene
    } else if ("gene_id" %in% colnames(df)) {
      genes <- df$gene_id
    } else {
      genes <- rownames(df)
    }
    
    # Filter by significance thresholds
    sig_idx <- which(df$padj < fdr_threshold & abs(df$log2FoldChange) > lfc_threshold)
    return(genes[sig_idx])
  })
  names(sig_genes_list) <- set_names
  
  # Compute pairwise IoU
  iou_mat <- matrix(NA, nrow = n_sets, ncol = n_sets)
  rownames(iou_mat) <- set_names
  colnames(iou_mat) <- set_names
  
  summary_rows <- list()
  
  for (i in seq_len(n_sets)) {
    for (j in seq_len(n_sets)) {
      set1 <- sig_genes_list[[i]]
      set2 <- sig_genes_list[[j]]
      
      iou_val <- compute_iou(set1, set2)
      iou_mat[i, j] <- iou_val
      
      # Build summary row
      if (!is.na(iou_val)) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Set1 = set_names[i],
          Set2 = set_names[j],
          Set1_Sig_Count = length(set1),
          Set2_Sig_Count = length(set2),
          Intersection_Size = length(intersect(set1, set2)),
          Union_Size = length(union(set1, set2)),
          IoU = iou_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  summary_df <- do.call(rbind, summary_rows)
  
  return(list(iou_matrix = iou_mat, summary_df = summary_df))
}

#' Add group annotations to IoU summary
#'
#' Adds Group column to IoU summary based on set names (Upper Bound, Cross-Dataset, Lower Bound)
#'
#' @param summary_df Data frame from build_iou_matrix
#' @param upper_bound_names Character vector of upper bound set names
#' @param null_names Character vector of null set names
#' @return Summary data frame with added Group column
annotate_iou_groups <- function(summary_df, upper_bound_names = NULL, null_names = NULL) {
  summary_df$Group <- "Cross-Dataset"
  
  if (!is.null(upper_bound_names)) {
    summary_df$Group[summary_df$Set1 %in% upper_bound_names | 
                       summary_df$Set2 %in% upper_bound_names] <- "Upper Bound"
  }
  
  if (!is.null(null_names)) {
    summary_df$Group[summary_df$Set1 %in% null_names | 
                       summary_df$Set2 %in% null_names] <- "Lower Bound"
  }
  
  summary_df$Group <- factor(summary_df$Group, 
                            levels = c("Upper Bound",
                                       "Cross-Dataset",
                                       "Lower Bound"))
  return(summary_df)
}

#' Build IoU matrix with universe control (count matrix based)
#'
#' Computes pairwise IoU restricted to genes present in BOTH count matrices.
#' This "universe-controlled" version uses the full count matrix gene set as the
#' universe - genes that DESeq2 filtered out (NA padj) are counted as "not detected"
#' rather than excluded. This penalizes datasets with shallow sequencing that cannot
#' detect genes that deeper datasets can.
#'
#' @param combined_results Named list of DE result data frames
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param lfc_threshold LFC threshold for significance (default: 1.0)
#' @return Named list containing:
#'   - iou_matrix: Matrix of universe-controlled IoU values
#'   - summary_df: Data frame with detailed statistics per pair
build_iou_matrix_universe_controlled <- function(combined_results,
                                                  fdr_threshold = 0.05, 
                                                  lfc_threshold = 1.0) {
  
  n_sets <- length(combined_results)
  set_names <- names(combined_results)
  
  if (n_sets == 0) {
    return(list(iou_matrix = NULL, summary_df = NULL))
  }
  
  # For each result set, identify the full universe (all genes in count matrix)
  # and significant genes within that universe
  universes <- lapply(combined_results, function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(list(universe = character(0), sig = character(0)))
    }
    
    # Get gene IDs
    if ("gene" %in% colnames(df)) {
      genes <- df$gene
    } else if ("gene_id" %in% colnames(df)) {
      genes <- df$gene_id
    } else {
      genes <- rownames(df)
    }
    
    # Universe: all genes present in the count matrix (results dataframe)
    # This includes genes that DESeq2 filtered out (NA padj) - they count as "not detected"
    universe_idx <- seq_len(nrow(df))
    universe_genes <- genes[universe_idx]
    
    # Significant genes (padj < threshold, |LFC| > threshold)
    # Genes with NA padj (DESeq2 filtered) will NOT be in sig_idx - counted as "not detected"
    sig_idx <- which(df$padj < fdr_threshold & abs(df$log2FoldChange) > lfc_threshold)
    sig_genes <- genes[sig_idx]
    
    return(list(universe = universe_genes, sig = sig_genes))
  })
  names(universes) <- set_names
  
  # Compute pairwise universe-controlled IoU
  iou_mat <- matrix(NA, nrow = n_sets, ncol = n_sets)
  rownames(iou_mat) <- set_names
  colnames(iou_mat) <- set_names
  
  summary_rows <- list()
  
  for (i in seq_len(n_sets)) {
    for (j in seq_len(n_sets)) {
      # Get full universes (count matrix genes) for both sets
      universe_i <- universes[[i]]$universe
      universe_j <- universes[[j]]$universe
      
      # Shared universe: genes present in BOTH count matrices
      shared_universe <- intersect(universe_i, universe_j)
      
      # Significant genes within the shared universe
      sig_i <- intersect(universes[[i]]$sig, shared_universe)
      sig_j <- intersect(universes[[j]]$sig, shared_universe)
      
      # Compute IoU on shared universe
      iou_val <- compute_iou(sig_i, sig_j)
      iou_mat[i, j] <- iou_val
      
      # Build summary row
      if (!is.na(iou_val)) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Set1 = set_names[i],
          Set2 = set_names[j],
          Set1_Universe_Size = length(universe_i),
          Set2_Universe_Size = length(universe_j),
          Shared_Universe_Size = length(shared_universe),
          Set1_Sig_Count = length(sig_i),
          Set2_Sig_Count = length(sig_j),
          Intersection_Size = length(intersect(sig_i, sig_j)),
          Union_Size = length(union(sig_i, sig_j)),
          IoU_Universe_Controlled = iou_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  summary_df <- do.call(rbind, summary_rows)
  
  return(list(iou_matrix = iou_mat, summary_df = summary_df))
}
