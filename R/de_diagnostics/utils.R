# ==============================================================================
# UTILS: DE Diagnostics Pipeline - Shared Utility Functions
# ==============================================================================
# Functions for loading DE results, BioMart annotation, pseudobulk aggregation,
# age parsing, and plotting helpers. Source this file after config.R.
# ==============================================================================

# ==============================================================================
# 1. GENE ID CLEANING
# ==============================================================================

#' Strip Ensembl version suffixes from gene IDs
#' @param x Character vector of gene IDs (e.g., "ENSG00000229807.11")
#' @return Character vector with versions removed (e.g., "ENSG00000229807")
clean_gene_ids <- function(x) {
  sub("\\..*", "", as.character(x))
}

# ==============================================================================
# 2. DE RESULT LOADING
# ==============================================================================

#' Load DE results CSV for a cell type
#'
#' Reads the DESeq2 results CSV, cleans gene IDs, and filters to significant
#' genes based on padj and LFC thresholds.
#'
#' @param cell_type Cell type name
#' @param padj_cutoff FDR threshold
#' @param lfc_min Minimum absolute LFC
#' @return data.frame with filtered DE results, or NULL on failure
load_de_results <- function(cell_type, padj_cutoff = PADJ_CUTOFF,
                            lfc_min = LFC_MIN_HIST) {

  csv_path <- get_de_results_path(cell_type)

  if (!file.exists(csv_path)) {
    log_warn("DE results not found: %s", csv_path)
    return(NULL)
  }

  df <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Validate required columns
  req_cols <- c("padj", "log2FoldChange", "lfcSE", "baseMean", "gene")
  missing <- setdiff(req_cols, colnames(df))
  if (length(missing) > 0) {
    log_warn("Missing columns in %s: %s", basename(csv_path), paste(missing, collapse = ", "))
    return(NULL)
  }

  # Clean gene IDs (from 'gene' column, not rownames)
  df$gene_id <- clean_gene_ids(df$gene)
  rownames(df) <- df$gene_id

  # Filter to significant genes
  df_filt <- df[!is.na(df$padj) & !is.na(df$log2FoldChange) & !is.na(df$lfcSE), ]
  df_filt <- df_filt[df_filt$padj < padj_cutoff & abs(df_filt$log2FoldChange) > lfc_min, ]

  if (nrow(df_filt) == 0) {
    log_warn("No significant genes for %s (padj < %g, |LFC| > %g)",
             cell_type, padj_cutoff, lfc_min)
    return(NULL)
  }

  df_filt$Cell_Type <- cell_type
  log_info("  %s: %d significant genes (of %d tested)", cell_type, nrow(df_filt), nrow(df))

  return(df_filt)
}

#' Load dispersion values from dds checkpoint
#'
#' @param cell_type Cell type name
#' @return Named numeric vector of dispersions (gene_id -> dispersion), or NULL
load_dispersions <- function(cell_type) {

  dds_path <- get_dds_checkpoint_path(cell_type)

  if (!file.exists(dds_path)) {
    log_warn("DDS checkpoint not found: %s", dds_path)
    return(NULL)
  }

  tryCatch({
    dds <- suppressMessages(qs::qread(dds_path))
    disp_vec <- DESeq2::dispersions(dds)
    names(disp_vec) <- clean_gene_ids(rownames(dds))
    log_info("  %s: loaded dispersions for %d genes", cell_type, length(disp_vec))
    return(disp_vec)
  }, error = function(e) {
    log_warn("Failed to load dispersions for %s: %s", cell_type, e$message)
    return(NULL)
  })
}

# ==============================================================================
# 3. BIOMART ANNOTATION
# ==============================================================================

#' Get gene biotype annotations, with CSV caching
#'
#' Queries BioMart for gene biotype (protein_coding vs other), caching results
#' to a CSV file for subsequent runs.
#'
#' @param gene_ids Character vector of Ensembl gene IDs (no version suffix)
#' @param cache_path Path to cache CSV file
#' @return data.frame with columns: ensembl_gene_id, gene_biotype, Biotype_Class
get_biotype_annotations <- function(gene_ids, cache_path) {

  # Check cache first
  if (file.exists(cache_path)) {
    log_info("Loading biotype cache from: %s", basename(cache_path))
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
    log_info("No biotype cache found, querying BioMart for %d genes", length(gene_ids))
    cached <- NULL
    query_ids <- gene_ids
  }

  # Query BioMart
  bm_result <- query_biomart(query_ids)

  if (!is.null(bm_result) && nrow(bm_result) > 0) {
    bm_result$Biotype_Class <- ifelse(
      !is.na(bm_result$gene_biotype) & bm_result$gene_biotype == "protein_coding",
      "Coding", "Non-Coding"
    )

    # Merge with cache
    if (!is.null(cached)) {
      all_annotations <- rbind(cached, bm_result)
      all_annotations <- all_annotations[!duplicated(all_annotations$ensembl_gene_id), ]
    } else {
      all_annotations <- bm_result
    }

    # Save updated cache
    write.csv(all_annotations, cache_path, row.names = FALSE)
    log_info("  Updated biotype cache: %d total genes", nrow(all_annotations))

    return(all_annotations[all_annotations$ensembl_gene_id %in% gene_ids, ])
  } else {
    log_warn("BioMart query failed or returned no results")
    if (!is.null(cached)) {
      return(cached[cached$ensembl_gene_id %in% gene_ids, ])
    }
    return(NULL)
  }
}

#' Query BioMart for gene biotype
#' @param gene_ids Character vector of Ensembl gene IDs
#' @return data.frame with ensembl_gene_id and gene_biotype columns, or NULL
query_biomart <- function(gene_ids) {

  # BIOMART_OFFLINE=true: use the shipped annotation cache only (cluster
  # compute nodes have no internet; the reproduction pipeline sets this).
  if (tolower(Sys.getenv("BIOMART_OFFLINE", "false")) %in% c("true", "1", "yes")) {
    log_info("  BIOMART_OFFLINE=true: not querying BioMart for %d genes", length(gene_ids))
    return(NULL)
  }


  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    log_warn("biomaRt package not available")
    return(NULL)
  }

  # Try primary host
  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = BIOMART_DATASET, host = BIOMART_HOST)
  }, error = function(e) {
    log_warn("BioMart primary host failed: %s", e$message)
    # Try fallback
    tryCatch({
      biomaRt::useMart("ensembl", dataset = BIOMART_DATASET, host = BIOMART_HOST_FALLBACK)
    }, error = function(e2) {
      log_warn("BioMart fallback host also failed: %s", e2$message)
      return(NULL)
    })
  })

  if (is.null(mart)) return(NULL)

  tryCatch({
    result <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "gene_biotype"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )
    log_info("  BioMart returned annotations for %d/%d genes",
             nrow(result), length(gene_ids))
    return(result)
  }, error = function(e) {
    log_warn("BioMart query failed: %s", e$message)
    return(NULL)
  })
}

# ==============================================================================
# 4. DENSITY PLOT HELPERS
# ==============================================================================

#' Compute auto-scaled x-axis limits from data
#'
#' Returns a range padded by a fraction of the data span on each side.
#' Optionally enforces a fixed lower or upper bound.
#'
#' @param df data.frame containing the metric column
#' @param col Name of the metric column
#' @param padding Fraction of data range to pad on each side (default 0.05)
#' @param lim_lower Optional fixed lower bound (overrides auto lower)
#' @param lim_upper Optional fixed upper bound (overrides auto upper)
#' @return Numeric vector of length 2, or NULL if no valid data
compute_auto_limits <- function(df, col, padding = 0.05,
                                lim_lower = NULL, lim_upper = NULL) {
  vals <- df[[col]]
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) return(NULL)

  rng <- range(vals)
  pad <- diff(rng) * padding

  lo <- if (!is.null(lim_lower)) lim_lower else rng[1] - pad
  hi <- if (!is.null(lim_upper)) lim_upper else rng[2] + pad

  c(lo, hi)
}

#' Generate median stats label for density plot subtitle
#'
#' @param df data.frame with metric column and grouping column
#' @param val_col Name of the metric column
#' @param group_col Name of the grouping column (e.g., "Cell_Type")
#' @return Character string with median values per group
get_stats_label <- function(df, val_col, group_col) {
  stats <- df %>%
    dplyr::group_by(!!rlang::sym(group_col)) %>%
    dplyr::summarise(Med = median(!!rlang::sym(val_col), na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(Label = sprintf("%s: %.3f", .data[[group_col]], Med))
  paste(stats$Label, collapse = " | ")
}

#' Generate aggregated density plot (cell types overlaid)
#'
#' @param df data.frame with metric and Cell_Type columns
#' @param x_col Name of the metric column
#' @param title_suffix Metric name for title
#' @param x_lab X-axis label
#' @param xlim_val Numeric vector of length 2 for x-axis limits, or NULL for auto
#' @return ggplot object or NULL
plot_density_aggregated <- function(df, x_col, title_suffix, x_lab, xlim_val) {
  df_c <- df[!is.na(df[[x_col]]), ]
  if (nrow(df_c) == 0) return(NULL)

  stats_str <- get_stats_label(df_c, x_col, "Cell_Type")

  p <- ggplot2::ggplot(df_c, ggplot2::aes(x = .data[[x_col]],
                                      fill = Cell_Type,
                                      color = Cell_Type)) +
    ggplot2::geom_density(alpha = 0.3, position = "identity", linewidth = 0.8) +
    ggplot2::scale_fill_brewer(palette = DENSITY_PALETTE) +
    ggplot2::scale_color_brewer(palette = DENSITY_PALETTE) +
    ggplot2::labs(
      title = sprintf("PsychAD: %s by Cell Type", title_suffix),
      subtitle = paste("Median:", stats_str),
      x = x_lab,
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())

  # Apply x-axis limits only if explicitly provided
  if (!is.null(xlim_val)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim_val)
  }

  p
}

#' Generate stratified density plot (faceted by coding/non-coding)
#'
#' @param df data.frame with metric, Cell_Type, and Biotype_Class columns
#' @param x_col Name of the metric column
#' @param title_suffix Metric name for title
#' @param x_lab X-axis label
#' @param xlim_val Numeric vector of length 2 for x-axis limits, or NULL for auto
#' @return ggplot object or NULL
plot_density_stratified <- function(df, x_col, title_suffix, x_lab, xlim_val) {
  df_c <- df[!is.na(df[[x_col]]) & !is.na(df$Biotype_Class), ]
  if (nrow(df_c) == 0) return(NULL)

  p <- ggplot2::ggplot(df_c, ggplot2::aes(x = .data[[x_col]],
                                      fill = Cell_Type,
                                      color = Cell_Type)) +
    ggplot2::geom_density(alpha = 0.3, position = "identity", linewidth = 0.8) +
    ggplot2::facet_wrap(~Biotype_Class, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_brewer(palette = DENSITY_PALETTE) +
    ggplot2::scale_color_brewer(palette = DENSITY_PALETTE) +
    ggplot2::labs(
      title = sprintf("PsychAD: %s (Stratified by Biotype)", title_suffix),
      x = x_lab,
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())

  # Apply x-axis limits only if explicitly provided
  if (!is.null(xlim_val)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim_val)
  }

  p
}

# ==============================================================================
# 5. AGE PARSING (reused from full_ground_truth)
# ==============================================================================

#' Parse verbose age strings to numeric years
#'
#' Handles formats like "80 year-old stage", "6 month-old stage", etc.
#' Returns NA for "unknown", "not available", etc.
#'
#' @param x Character vector of age strings
#' @return Numeric vector of ages in years
parse_age_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
  is_unknown <- grepl("^(unknown|not available|na|n/a|nan)$", s)
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))

  is_year  <- grepl("year[ -]?old.*stage", s)
  is_month <- grepl("month[ -]?old.*stage", s)
  is_week  <- grepl("week[ -]?old.*stage", s)
  is_day   <- grepl("day[ -]?old.*stage", s)

  y <- rep(NA_real_, length(s))
  y[is_year]    <- n[is_year]
  y[is_month]   <- n[is_month] / 12
  y[is_week]    <- n[is_week]  / 52.1429
  y[is_day]     <- n[is_day]   / 365.25
  y[is_unknown] <- NA_real_

  # Fallback: if no specific unit found but it is a number, assume years
  if (all(is.na(y)) && !any(is_unknown) && any(!is.na(n))) y <- n
  return(y)
}

# ==============================================================================
# 6. PSEUDOBULK AGGREGATION (for PCA)
# ==============================================================================

#' Aggregate donor shards into pseudobulk for PCA
#'
#' Loads individual donor .qs shard files, sums counts per donor, and returns
#' a SingleCellExperiment with donors as columns.
#'
#' @param donor_ids Character vector of donor IDs to load
#' @param shard_dir Path to shard directory
#' @param n_cores Number of cores for parallel loading
#' @return SingleCellExperiment with pseudobulk counts and metadata, or NULL
aggregate_pseudobulk <- function(donor_ids, shard_dir, n_cores = N_CORES) {

  local_param <- BiocParallel::MulticoreParam(
    workers = n_cores,
    stop.on.error = FALSE
  )

  pb_list <- BiocParallel::bplapply(donor_ids, function(donor_id) {
    safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))

    if (!file.exists(shard_path)) return(NULL)

    sce_shard <- tryCatch({
      qs::qread(shard_path)
    }, error = function(e) return(NULL))

    if (is.null(sce_shard)) return(NULL)

    assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(sce_shard)) {
      "counts"
    } else {
      SummarizedExperiment::assayNames(sce_shard)[1]
    }

    counts_vec <- as.integer(Matrix::rowSums(
      SummarizedExperiment::assay(sce_shard, assay_name)
    ))
    names(counts_vec) <- rownames(sce_shard)

    meta_row <- as.data.frame(SummarizedExperiment::colData(sce_shard))[1, , drop = FALSE]

    return(list(counts = counts_vec, meta = meta_row, genes = names(counts_vec)))
  }, BPPARAM = local_param)

  # Remove NULLs
  pb_list <- pb_list[!sapply(pb_list, is.null)]
  if (length(pb_list) == 0) return(NULL)

  # Intersect genes
  all_gene_sets <- lapply(pb_list, `[[`, "genes")
  common_genes <- Reduce(intersect, all_gene_sets)

  if (length(common_genes) == 0) return(NULL)

  # Build count matrix
  count_mat <- do.call(cbind, lapply(pb_list, function(x) {
    x$counts[match(common_genes, x$genes)]
  }))
  rownames(count_mat) <- common_genes

  # Build metadata
  meta_df <- do.call(rbind, lapply(pb_list, `[[`, "meta"))

  # Use donor_id as column names
  if (DONOR_ID_COL %in% colnames(meta_df)) {
    col_ids <- as.character(meta_df[[DONOR_ID_COL]])
  } else {
    col_ids <- paste0("donor_", seq_len(nrow(meta_df)))
  }
  colnames(count_mat) <- col_ids
  rownames(meta_df) <- col_ids

  # Create SCE
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = count_mat),
    colData = S4Vectors::DataFrame(meta_df)
  )

  return(sce)
}

#' Extract GT donor IDs from cartesian downsampling metadata
#'
#' @param metadata_path Path to Ground_Truth_sampling_metadata.qs
#' @return Character vector of GT donor IDs, or NULL
extract_gt_donor_ids <- function(metadata_path) {
  if (!file.exists(metadata_path)) {
    log_warn("GT metadata not found: %s", metadata_path)
    return(NULL)
  }

  tryCatch({
    gt_obj <- qs::qread(metadata_path)
    # GT metadata is a named list — donor IDs are the names
    gt_ids <- names(gt_obj)
    if (is.null(gt_ids)) {
      # Fallback: try as vector
      if (is.vector(gt_obj) && is.atomic(gt_obj)) {
        gt_ids <- as.character(gt_obj)
      } else {
        gt_ids <- unlist(gt_obj)
      }
    }
    return(as.character(gt_ids))
  }, error = function(e) {
    log_warn("Failed to read GT metadata: %s", e$message)
    return(NULL)
  })
}

#' Get all donor IDs from manifest
#'
#' @param shard_dir Path to shard directory containing 00_MASTER_MANIFEST.qs
#' @return Character vector of all donor IDs, or NULL
get_all_donor_ids <- function(shard_dir) {
  manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
  if (!file.exists(manifest_path)) {
    log_warn("Manifest not found: %s", manifest_path)
    return(NULL)
  }

  tryCatch({
    manifest <- qs::qread(manifest_path)$manifest
    return(as.character(manifest$donor_id))
  }, error = function(e) {
    log_warn("Failed to read manifest: %s", e$message)
    return(NULL)
  })
}

# ==============================================================================
# 7. SPLIT GT DE RESULT LOADING (for Step 3)
# ==============================================================================

#' Load cartesian downsampling split GT DE results
#'
#' Reads a Ground_Truth.qs file (DESeqResults object), converts to data.frame,
#' and cleans gene IDs.
#'
#' @param gt_path Path to Ground_Truth.qs
#' @return data.frame with DE results, or NULL
load_split_gt_results <- function(gt_path) {
  if (!file.exists(gt_path)) {
    log_warn("Split GT file not found: %s", gt_path)
    return(NULL)
  }

  tryCatch({
    res <- qs::qread(gt_path)
    df <- as.data.frame(res)
    df$gene_id <- clean_gene_ids(rownames(df))
    rownames(df) <- df$gene_id
    return(df)
  }, error = function(e) {
    log_warn("Failed to load split GT: %s", e$message)
    return(NULL)
  })
}

# ==============================================================================
# 8. CHROMOSOME ANNOTATION (for Step 3)
# ==============================================================================

#' Get chromosome annotations, extending the BioMart cache
#'
#' Queries BioMart for chromosome_name in addition to gene_biotype.
#' Extends the existing cache file if chromosome info is missing.
#'
#' @param gene_ids Character vector of Ensembl gene IDs
#' @param cache_path Path to cache CSV file
#' @return data.frame with ensembl_gene_id, gene_biotype, Biotype_Class,
#'         chromosome_name, Chr_Type columns
get_chromosome_annotations <- function(gene_ids, cache_path) {

  # Check if cache already has chromosome info
  if (file.exists(cache_path)) {
    cached <- read.csv(cache_path, stringsAsFactors = FALSE)
    if ("chromosome_name" %in% colnames(cached)) {
      missing <- setdiff(gene_ids, cached$ensembl_gene_id)
      if (length(missing) == 0) {
        log_info("  Cache covers all %d genes with chromosome info", length(gene_ids))
        result <- cached[cached$ensembl_gene_id %in% gene_ids, ]
        result$Chr_Type <- classify_chromosome(result$chromosome_name)
        return(result)
      }
      log_info("  Cache missing chromosome info for %d genes", length(missing))
      query_ids <- missing
    } else {
      log_info("  Cache exists but lacks chromosome_name — re-querying all %d genes",
               length(gene_ids))
      query_ids <- gene_ids
    }
  } else {
    log_info("  No cache found — querying BioMart for %d genes", length(gene_ids))
    cached <- NULL
    query_ids <- gene_ids
  }

  # Query BioMart with chromosome info
  bm_result <- query_biomart_with_chr(query_ids)

  if (!is.null(bm_result) && nrow(bm_result) > 0) {
    bm_result$Biotype_Class <- ifelse(
      !is.na(bm_result$gene_biotype) & bm_result$gene_biotype == "protein_coding",
      "Coding", "Non-Coding"
    )

    # Merge with existing cache
    if (!is.null(cached) && "chromosome_name" %in% colnames(cached)) {
      all_annotations <- rbind(cached, bm_result)
      all_annotations <- all_annotations[!duplicated(all_annotations$ensembl_gene_id), ]
    } else if (!is.null(cached)) {
      # Cache exists but without chromosome — replace with full data
      # Keep cached entries that aren't in the new query
      old_only <- cached[!(cached$ensembl_gene_id %in% bm_result$ensembl_gene_id), ]
      if (nrow(old_only) > 0) {
        old_only$chromosome_name <- NA_character_
      }
      all_annotations <- rbind(bm_result, old_only)
    } else {
      all_annotations <- bm_result
    }

    # Save updated cache
    write.csv(all_annotations, cache_path, row.names = FALSE)
    log_info("  Updated cache with chromosome info: %d total genes", nrow(all_annotations))

    result <- all_annotations[all_annotations$ensembl_gene_id %in% gene_ids, ]
    result$Chr_Type <- classify_chromosome(result$chromosome_name)
    return(result)
  } else {
    log_warn("BioMart chromosome query failed")
    if (!is.null(cached)) {
      result <- cached[cached$ensembl_gene_id %in% gene_ids, ]
      if (!"chromosome_name" %in% colnames(result)) {
        result$chromosome_name <- NA_character_
      }
      result$Chr_Type <- classify_chromosome(result$chromosome_name)
      return(result)
    }
    return(NULL)
  }
}

#' Query BioMart for gene biotype AND chromosome name
#' @param gene_ids Character vector of Ensembl gene IDs
#' @return data.frame with ensembl_gene_id, gene_biotype, chromosome_name
query_biomart_with_chr <- function(gene_ids) {

  # BIOMART_OFFLINE=true: use the shipped annotation cache only (cluster
  # compute nodes have no internet; the reproduction pipeline sets this).
  if (tolower(Sys.getenv("BIOMART_OFFLINE", "false")) %in% c("true", "1", "yes")) {
    log_info("  BIOMART_OFFLINE=true: not querying BioMart for %d genes", length(gene_ids))
    return(NULL)
  }


  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    log_warn("biomaRt package not available")
    return(NULL)
  }

  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = BIOMART_DATASET, host = BIOMART_HOST)
  }, error = function(e) {
    log_warn("BioMart primary host failed: %s", e$message)
    tryCatch({
      biomaRt::useMart("ensembl", dataset = BIOMART_DATASET, host = BIOMART_HOST_FALLBACK)
    }, error = function(e2) {
      log_warn("BioMart fallback host also failed: %s", e2$message)
      return(NULL)
    })
  })

  if (is.null(mart)) return(NULL)

  tryCatch({
    result <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "gene_biotype", "chromosome_name"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )
    log_info("  BioMart returned %d annotations (with chromosome)", nrow(result))
    return(result)
  }, error = function(e) {
    log_warn("BioMart query failed: %s", e$message)
    return(NULL)
  })
}

#' Classify chromosome names into Autosomal or SexChr
#'
#' @param chr_names Character vector of chromosome names
#' @return Character vector: "Autosomal", "SexChr", or "Other"
classify_chromosome <- function(chr_names) {
  ifelse(chr_names %in% CHR_AUTOSOMAL, "Autosomal",
         ifelse(chr_names %in% CHR_SEX, "SexChr", "Other"))
}

# ==============================================================================
# 9. OVERLAP COMPUTATION (for Step 3)
# ==============================================================================

#' Compute bidirectional overlap between two gene sets
#'
#' @param full_gt_genes Character vector of significant gene IDs from full GT
#' @param split_genes Character vector of significant gene IDs from split GT
#' @return List with overlap stats
compute_overlap_stats <- function(full_gt_genes, split_genes) {
  overlap <- intersect(full_gt_genes, split_genes)
  full_only <- setdiff(full_gt_genes, split_genes)
  split_only <- setdiff(split_genes, full_gt_genes)

  n_full <- length(full_gt_genes)
  n_split <- length(split_genes)
  n_overlap <- length(overlap)

  recall <- if (n_full > 0) n_overlap / n_full else NA_real_
  precision <- if (n_split > 0) n_overlap / n_split else NA_real_

  list(
    N_FullGT = n_full,
    N_Split = n_split,
    N_Overlap = n_overlap,
    N_FullOnly = length(full_only),
    N_SplitOnly = length(split_only),
    Recall = recall,
    Precision = precision,
    overlap_genes = overlap,
    full_only_genes = full_only,
    split_only_genes = split_only
  )
}

#' Assign LFC magnitude bin to genes
#'
#' @param abs_lfc Numeric vector of absolute LFC values
#' @param bins Named list of c(lower, upper) bounds
#' @return Character vector of bin labels
assign_lfc_bin <- function(abs_lfc, bins = RECOVERY_LFC_BINS) {
  result <- rep(NA_character_, length(abs_lfc))
  for (bin_name in names(bins)) {
    bounds <- bins[[bin_name]]
    in_bin <- abs_lfc >= bounds[1] & abs_lfc < bounds[2]
    result[in_bin] <- bin_name
  }
  return(result)
}

# ==============================================================================
# 10. SPLIT VIOLIN GEOM (for Step 6)
# ==============================================================================

#' Custom split violin geom for ggplot2
#'
#' Draws half-violins: one side for each level of the `split` aesthetic.
#' The first level is drawn on the left (negative x offset), the second on
#' the right (positive x offset). Requires exactly 2 levels in the split
#' aesthetic per x position.
#'
#' Based on the well-known community pattern (jan-glx / StackOverflow).
#'
#' @usage geom_split_violin(mapping, data, stat, position, ...)
#' @examples
#' ggplot(df, aes(x = group, y = value, fill = split_var)) +
#'   geom_split_violin(alpha = 0.7)

GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    # Compute the density-scaled left and right edges of the violin.
    # violinwidth varies per row (proportional to density at each y value),
    # while x, xmin, xmax are constants (the x-position and bounding box).
    data <- transform(data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )

    # Determine which half this group is (left or right).
    # The group aesthetic alternates: odd = left, even = right.
    grp <- data$group[1]

    # Store the centre x-position before overwriting x with density edges.
    centre_x <- data$x[1]

    if (grp %% 2 == 1) {
      # Left half: density curve on the left, flat centre edge on the right.
      # Set x to the density-scaled left edge and order y ascending.
      newdata <- transform(data, x = xminv)
      newdata <- newdata[order(newdata$y), ]
    } else {
      # Right half: flat centre edge on the left, density curve on the right.
      # Set x to the density-scaled right edge and order y descending.
      newdata <- transform(data, x = xmaxv)
      newdata <- newdata[order(-newdata$y), ]
    }

    # Close the polygon by adding boundary points at the centre x-position.
    # The first and last rows are duplicated, then their x is set to the
    # centre — this creates the flat vertical edge of the half-violin.
    newdata <- rbind(
      newdata[1, , drop = FALSE],
      newdata,
      newdata[nrow(newdata), , drop = FALSE]
    )
    newdata$x[c(1, nrow(newdata))] <- centre_x

    # Draw using GeomPolygon (filled polygon from x, y coordinates)
    ggplot2::GeomPolygon$draw_panel(newdata, ...)
  }
)

#' Split violin geom layer function
#'
#' @param mapping Aesthetic mapping (requires fill or colour for split)
#' @param data Data frame
#' @param stat Statistical transformation (default: "ydensity")
#' @param position Position adjustment (default: "identity")
#' @param ... Additional arguments passed to layer
#' @param draw_quantiles Quantiles to draw (default: NULL)
#' @param trim Trim violins to data range (default: TRUE)
#' @param scale Scale violins (default: "area")
#' @param na.rm Remove NAs (default: FALSE)
#' @param show.legend Show legend (default: NA)
#' @param inherit.aes Inherit aesthetics (default: TRUE)
#' @return A ggplot2 layer
geom_split_violin <- function(mapping = NULL, data = NULL,
                               stat = "ydensity", position = "identity",
                               ...,
                               draw_quantiles = NULL,
                               trim = TRUE,
                               scale = "area",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}

# ==============================================================================
# 11. GENE SYMBOL LOOKUP (for Step 6)
# ==============================================================================

#' Get gene symbols from Ensembl IDs via BioMart with caching
#'
#' @param gene_ids Character vector of Ensembl gene IDs (no version suffix)
#' @param cache_path Path to cache CSV file (will be created/extended)
#' @return Named character vector: ensembl_id -> gene_symbol
get_gene_symbols <- function(gene_ids, cache_path) {

  # Check cache first
  if (file.exists(cache_path)) {
    cached <- read.csv(cache_path, stringsAsFactors = FALSE)
    if (all(c("ensembl_gene_id", "external_gene_name") %in% colnames(cached))) {
      missing <- setdiff(gene_ids, cached$ensembl_gene_id)
      if (length(missing) == 0) {
        log_info("Gene symbol cache covers all %d genes", length(gene_ids))
        result <- cached$external_gene_name[match(gene_ids, cached$ensembl_gene_id)]
        names(result) <- gene_ids
        # Replace empty strings with Ensembl ID
        result[is.na(result) | result == ""] <- gene_ids[is.na(result) | result == ""]
        return(result)
      }
      log_info("Gene symbol cache missing %d/%d genes, querying BioMart",
               length(missing), length(gene_ids))
      query_ids <- missing
    } else {
      query_ids <- gene_ids
      cached <- NULL
    }
  } else {
    query_ids <- gene_ids
    cached <- NULL
  }

  # Query BioMart
  bm_result <- tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      log_warn("biomaRt not available for gene symbol lookup")
      return(NULL)
    }

    mart <- tryCatch({
      biomaRt::useMart("ensembl", dataset = BIOMART_DATASET, host = BIOMART_HOST)
    }, error = function(e) {
      tryCatch({
        biomaRt::useMart("ensembl", dataset = BIOMART_DATASET,
                          host = BIOMART_HOST_FALLBACK)
      }, error = function(e2) NULL)
    })

    if (is.null(mart)) return(NULL)

    biomaRt::getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = query_ids,
      mart = mart
    )
  }, error = function(e) {
    log_warn("BioMart gene symbol query failed: %s", e$message)
    NULL
  })

  if (!is.null(bm_result) && nrow(bm_result) > 0) {
    # Merge with cache
    if (!is.null(cached)) {
      # Add external_gene_name to cached if missing
      if (!"external_gene_name" %in% colnames(cached)) {
        cached$external_gene_name <- NA_character_
      }
      # Ensure bm_result has all columns from cached
      for (col in setdiff(colnames(cached), colnames(bm_result))) {
        bm_result[[col]] <- NA
      }
      for (col in setdiff(colnames(bm_result), colnames(cached))) {
        cached[[col]] <- NA
      }
      all_data <- rbind(cached, bm_result[, colnames(cached)])
      all_data <- all_data[!duplicated(all_data$ensembl_gene_id), ]
    } else {
      all_data <- bm_result
    }

    # Save updated cache
    write.csv(all_data, cache_path, row.names = FALSE)
    log_info("Updated gene symbol cache: %d total genes", nrow(all_data))

    result <- all_data$external_gene_name[match(gene_ids, all_data$ensembl_gene_id)]
    names(result) <- gene_ids
    result[is.na(result) | result == ""] <- gene_ids[is.na(result) | result == ""]
    return(result)
  }

  # Fallback: return Ensembl IDs as symbols
  log_warn("Gene symbol lookup failed — using Ensembl IDs as labels")
  result <- gene_ids
  names(result) <- gene_ids
  return(result)
}
