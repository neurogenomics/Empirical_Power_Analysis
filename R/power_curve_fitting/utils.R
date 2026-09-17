# ==============================================================================
# UTILS: Power Curve Fitting Pipeline v4.0.0 - Shared Utility Functions
# ==============================================================================
# This file contains all shared utility functions for the power curve fitting
# pipeline. Source this file after config.R.
#
# v4.0.0: Switched from beta regression to logistic regression (glm, family=binomial).
#          Removed squeeze functions and safe_betareg. Added safe_glm_binomial.
# v3.1.0: Added metadata parsing, shard cache, unified formula builder.
# ==============================================================================

# ==============================================================================
# 1. AGE PARSING (reused from cartesian_downsampling)
# ==============================================================================

#' Parse age stage strings to numeric years
#'
#' @param x Character vector of age strings
#' @return Numeric vector of ages in years
parse_age_stage_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))

  is_year    <- grepl("year[ -]?old", s)
  is_month   <- grepl("month[ -]?old", s)
  is_week    <- grepl("week[ -]?old", s)
  is_day     <- grepl("day[ -]?old", s)
  is_unknown <- grepl("^(unknown|not available|na|n/a|nan|)$", s)

  y <- rep(NA_real_, length(s))
  y[is_year]    <- n[is_year]
  y[is_month]   <- n[is_month] / 12
  y[is_week]    <- n[is_week] / 52.1429
  y[is_day]     <- n[is_day] / 365.25
  y[is_unknown] <- NA_real_

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
  }
  as.numeric(scale(age_vec, center = TRUE, scale = TRUE))
}

# ==============================================================================
# 2. SHARD AGGREGATION (reused from cartesian_downsampling)
# ==============================================================================

#' Aggregate counts from donor shards into pseudobulk
#'
#' @param donor_ids Character vector of donor IDs to aggregate
#' @param shard_dir Path to directory containing donor shards
#' @param cell_indices_map Named list of cell indices per donor (from sampling metadata)
#' @param seed Random seed for reproducibility
#' @param verbose Whether to print progress messages
#' @return List with 'sce' (SingleCellExperiment) and 'sampling_log'
aggregate_from_shards <- function(donor_ids, shard_dir,
                                   cell_indices_map = NULL,
                                   seed = 42, verbose = FALSE) {

  results_list <- lapply(donor_ids, function(id) {
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))

    if (!file.exists(shard_path)) {
      log_error("Missing shard: %s", shard_path)
      return(NULL)
    }

    sce_shard <- tryCatch({
      qs::qread(shard_path)
    }, error = function(e) {
      log_error("Failed to load shard %s: %s", shard_path, e$message)
      return(NULL)
    })

    if (is.null(sce_shard)) return(NULL)

    # If cell indices are provided, subset to those cells
    if (!is.null(cell_indices_map) && id %in% names(cell_indices_map)) {
      idx <- cell_indices_map[[id]]
      if (max(idx) <= ncol(sce_shard)) {
        sce_shard <- sce_shard[, idx]
      }
    }

    # Aggregate to pseudobulk (sum counts across cells)
    counts_vec <- Matrix::rowSums(SummarizedExperiment::assay(sce_shard, "counts"))
    meta_row <- as.data.frame(SummarizedExperiment::colData(sce_shard)[1, , drop = FALSE])

    return(list(
      counts = counts_vec,
      meta = meta_row,
      genes = names(counts_vec),
      donor_id = id
    ))
  })

  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    log_warning("No valid donors after aggregation")
    return(NULL)
  }

  all_genes <- results_list[[1]]$genes

  combined_counts <- do.call(cbind, lapply(results_list, function(x) {
    if (!identical(x$genes, all_genes)) {
      return(x$counts[match(all_genes, x$genes)])
    }
    return(x$counts)
  }))

  combined_meta <- do.call(rbind, lapply(results_list, `[[`, "meta"))

  rownames(combined_counts) <- all_genes
  colnames(combined_counts) <- as.character(combined_meta[[DONOR_ID_COL]])
  rownames(combined_meta) <- colnames(combined_counts)

  sce_out <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = combined_counts),
    colData = S4Vectors::DataFrame(combined_meta)
  )

  return(list(sce = sce_out))
}

# ==============================================================================
# 3. DESEQ2 WORKER (adapted from cartesian_downsampling)
# ==============================================================================

#' Run DESeq2 and return the DDS object (not just results)
#'
#' Unlike the cartesian_downsampling version, this saves the full DDS object
#' so we can extract dispersion estimates.
#'
#' @param agg_result Result from aggregate_from_shards()
#' @param save_path_dds Path to save DESeqDataSet object
#' @param save_path_results Path to save DESeqResults (optional, NULL to skip)
#' @param min_gene_count Minimum total count for gene filtering
#' @return "SUCCESS" or error message
run_deseq_for_dds <- function(agg_result,
                               save_path_dds,
                               save_path_results = NULL,
                               min_gene_count = MIN_GENE_COUNT) {

  pb_sce <- agg_result$sce

  # Prepare colData
  cts <- as.data.frame(SummarizedExperiment::colData(pb_sce))

  # Parse and standardize age
  raw_age <- as.character(cts[[AGE_COL]])
  num_age <- parse_age_stage_to_years(raw_age)
  age_has_variance <- length(unique(num_age[!is.na(num_age)])) > 1

  if (age_has_variance) {
    cts$age_scaled <- standardize_age(num_age, impute_na = TRUE)
  }

  # Process numeric covariates (e.g. pmi). Each is mean-imputed for NAs,
  # z-scaled, and entered as <name>_scaled in the design.
  num_vars <- if (exists("NUMERIC_VARS")) NUMERIC_VARS else character()
  numeric_terms <- character()
  for (nv in num_vars) {
    if (!nv %in% colnames(cts)) next
    raw_val <- suppressWarnings(as.numeric(as.character(cts[[nv]])))
    if (all(is.na(raw_val))) next
    if (any(is.na(raw_val))) raw_val[is.na(raw_val)] <- mean(raw_val, na.rm = TRUE)
    if (sd(raw_val, na.rm = TRUE) == 0) next
    scaled_name <- paste0(nv, "_scaled")
    cts[[scaled_name]] <- as.numeric(scale(raw_val))
    numeric_terms <- c(numeric_terms, scaled_name)
  }

  # Identify valid categorical covariates
  valid_vars <- c()
  for (v in CATEGORICAL_VARS) {
    if (v %in% colnames(cts) && length(unique(cts[[v]])) > 1) {
      sex_var_table <- table(cts[[SEX_COL]], cts[[v]])
      if (all(rowSums(sex_var_table > 0) > 1) && all(colSums(sex_var_table > 0) > 1)) {
        valid_vars <- c(valid_vars, v)
      } else {
        log_warning("Covariate '%s' appears collinear with sex, excluding", v)
      }
    }
  }

  # Build formula
  formula_parts <- c(SEX_COL)
  if (age_has_variance) formula_parts <- c(formula_parts, "age_scaled")
  formula_parts <- c(formula_parts, numeric_terms, valid_vars)
  formula_str <- paste0("~ ", paste(formula_parts, collapse = " + "))

  # Create DESeqDataSet
  dds <- tryCatch({
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(pb_sce, "counts"),
      colData = cts,
      design = as.formula(formula_str)
    ))
  }, error = function(e) {
    log_warning("Model '%s' failed: %s. Trying sex-only model.", formula_str, e$message)
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(pb_sce, "counts"),
      colData = cts,
      design = as.formula(paste0("~ ", SEX_COL))
    ))
  })

  # Filter low-count genes
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_gene_count, ]

  # Run DESeq2 with full parallelization (24 cores)
  dds <- tryCatch({
    suppressMessages(DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, quiet = TRUE))
  }, error = function(e) {
    log_warning("DESeq2 failed: %s. Trying sex-only model.", e$message)
    DESeq2::design(dds) <- as.formula(paste0("~ ", SEX_COL))
    suppressMessages(DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, quiet = TRUE))
  })

  # Save DDS checkpoint
  qs::qsave(dds, save_path_dds)
  log_info("Saved DDS checkpoint: %s", basename(save_path_dds))

  # Optionally save shrunk results
  coef_available <- TRUE
  if (!is.null(save_path_results)) {
    coef_name <- paste0(SEX_COL, "_", MALE_VAL, "_vs_", FEMALE_VAL)
    coef_available <- coef_name %in% DESeq2::resultsNames(dds)
    if (coef_available) {
      res <- suppressMessages(DESeq2::lfcShrink(dds, coef = coef_name,
                                                 type = "ashr", quiet = TRUE))
    } else {
      log_warning("Coefficient '%s' not found, using ashr with contrast", coef_name)
      res <- suppressMessages(DESeq2::lfcShrink(dds,
                                                 contrast = c(SEX_COL, MALE_VAL, FEMALE_VAL),
                                                 type = "ashr", quiet = TRUE))
    }
    attr(res, "shrinkage_type") <- "ashr"
    qs::qsave(res, save_path_results)
  }

  return(if (coef_available) "SUCCESS" else "SUCCESS_FALLBACK")
}

# ==============================================================================
# 4. GENE ANNOTATION
# ==============================================================================

#' Load BioMart gene annotation cache
#'
#' @param cache_path Path to BioMart cache CSV
#' @return Data frame with ensembl_gene_id, chromosome_name, gene_biotype, band
load_gene_annotations <- function(cache_path = BIOMART_CACHE) {
  if (!file.exists(cache_path)) {
    log_warning("BioMart cache not found: %s", cache_path)
    return(NULL)
  }

  annot <- read.csv(cache_path, stringsAsFactors = FALSE)
  log_info("Loaded gene annotations: %d genes", nrow(annot))
  return(annot)
}

#' Classify chromosome as autosomal or sex chromosome
#'
#' @param chr_name Character vector of chromosome names
#' @return Character vector: "Autosomal", "SexChr", or NA
classify_chromosome <- function(chr_name) {
  result <- rep(NA_character_, length(chr_name))
  result[chr_name %in% as.character(1:22)] <- "Autosomal"
  result[chr_name %in% c("X", "Y")] <- "SexChr"
  return(result)
}

# ==============================================================================
# 5. ITERATION RESULT PARSING
# ==============================================================================

#' Parse iteration result filename to extract grid parameters
#'
#' @param filename Filename (not full path)
#' @param mode "cells" or "reads"
#' @return Named list with n_donors, y_value, iteration, or NULL if parse fails
parse_iteration_filename <- function(filename, mode) {
  fp <- get_file_pattern(mode)

  m <- regmatches(filename, regexec(fp$parse_regex, filename))[[1]]
  if (length(m) < 4) return(NULL)

  n_donors <- as.numeric(m[2])
  y_raw <- m[3]
  iteration <- as.numeric(m[4])

  # Parse y_value
  if (y_raw == "NA") {
    y_value <- NA_real_
  } else {
    y_value <- as.numeric(y_raw)
  }

  list(
    n_donors = n_donors,
    y_value = y_value,
    iteration = iteration
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

  # Exclude sampling metadata files
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
# 6. SHARD READ COUNT CACHE (v3.1.0)
# ==============================================================================

#' Pre-compute per-cell read counts for a donor shard
#'
#' Computes colSums(counts) for each cell in the shard and saves as a named
#' numeric vector. This avoids loading the full shard repeatedly during
#' metadata parsing.
#'
#' @param shard_path Path to donor shard .qs file
#' @param cache_path Path to save cached read counts
#' @return Named numeric vector of per-cell read counts, or NULL on failure
cache_shard_read_counts <- function(shard_path, cache_path) {
  sce <- tryCatch({
    qs::qread(shard_path)
  }, error = function(e) {
    log_error("Failed to load shard %s: %s", shard_path, e$message)
    return(NULL)
  })

  if (is.null(sce)) return(NULL)

  # Compute per-cell total read counts
  reads_per_cell <- Matrix::colSums(SummarizedExperiment::assay(sce, "counts"))

  # Save cache
  qs::qsave(reads_per_cell, cache_path)

  return(reads_per_cell)
}

#' Load cached per-cell read counts for a donor
#'
#' @param cache_path Path to cached read counts
#' @return Named numeric vector of per-cell read counts, or NULL if not found
load_shard_cache <- function(cache_path) {
  if (!file.exists(cache_path)) return(NULL)
  tryCatch({
    qs::qread(cache_path)
  }, error = function(e) {
    log_error("Failed to load cache %s: %s", cache_path, e$message)
    return(NULL)
  })
}

#' Get the cache path for a donor shard
#'
#' @param donor_id Donor ID string
#' @param cache_dir Directory for shard caches
#' @return Path to cache file
get_shard_cache_path <- function(donor_id, cache_dir) {
  safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
  file.path(cache_dir, paste0("donor_", safe_id, "_reads.qs"))
}

# ==============================================================================
# 7. SAMPLING METADATA PARSING (v3.1.0)
# ==============================================================================

#' Parse a sampling metadata filename to extract grid parameters
#'
#' @param filename Filename (not full path)
#' @param mode "cells" or "reads"
#' @return Named list with n_donors, y_value, iteration, or NULL if parse fails
parse_metadata_filename <- function(filename, mode) {
  fp <- get_metadata_file_pattern(mode)

  m <- regmatches(filename, regexec(fp$parse_regex, filename))[[1]]
  if (length(m) < 4) return(NULL)

  n_donors <- as.numeric(m[2])
  y_raw <- m[3]
  iteration <- as.numeric(m[4])

  if (y_raw == "NA") {
    y_value <- NA_real_
  } else {
    y_value <- as.numeric(y_raw)
  }

  list(
    n_donors = n_donors,
    y_value = y_value,
    iteration = iteration
  )
}

#' List all sampling metadata files in a run directory
#'
#' @param run_dir Path to run directory
#' @param mode "cells" or "reads"
#' @return Data frame with columns: file_path, n_donors, y_value, iteration
list_metadata_files <- function(run_dir, mode) {
  fp <- get_metadata_file_pattern(mode)

  all_files <- list.files(run_dir, pattern = fp$pattern, full.names = FALSE)

  # Exclude the Ground_Truth_sampling_metadata.qs
  all_files <- all_files[!grepl("^Ground_Truth", all_files)]

  if (length(all_files) == 0) return(NULL)

  parsed <- lapply(all_files, function(f) {
    p <- parse_metadata_filename(f, mode)
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

#' Compute actual n_cells and mean reads_per_cell for a single iteration
#'
#' Loads the sampling metadata (donor_id -> cell indices), looks up per-cell
#' read counts from the shard cache, and computes:
#'   - n_cells: mean number of cells across donors in this iteration
#'   - reads_per_cell: mean reads per cell across all cells in this iteration
#'
#' @param metadata_path Path to sampling metadata .qs file
#' @param shard_cache Named list of per-cell read count vectors (donor_id -> numeric vector)
#' @return Named list with n_cells_mean, reads_per_cell_mean, n_donors_actual,
#'         total_cells, total_reads, or NULL on failure
compute_iteration_depth <- function(metadata_path, shard_cache) {
  metadata <- tryCatch({
    qs::qread(metadata_path)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(metadata) || length(metadata) == 0) return(NULL)

  donor_n_cells <- numeric(length(metadata))
  donor_rpc <- numeric(length(metadata))
  valid_donors <- 0

  for (i in seq_along(metadata)) {
    donor_id <- names(metadata)[i]
    cell_indices <- metadata[[i]]

    # Look up cached read counts for this donor
    donor_reads <- shard_cache[[donor_id]]
    if (is.null(donor_reads)) {
      # Try with safe ID
      safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
      donor_reads <- shard_cache[[safe_id]]
    }

    if (is.null(donor_reads)) next

    # Get reads for the selected cells
    # Cell indices from sampling metadata are 1-based column indices
    valid_idx <- cell_indices[cell_indices <= length(donor_reads)]
    if (length(valid_idx) == 0) next

    selected_reads <- donor_reads[valid_idx]

    valid_donors <- valid_donors + 1
    donor_n_cells[valid_donors] <- length(valid_idx)
    donor_rpc[valid_donors] <- mean(selected_reads)
  }

  if (valid_donors == 0) return(NULL)

  donor_n_cells <- donor_n_cells[1:valid_donors]
  donor_rpc <- donor_rpc[1:valid_donors]

  list(
    n_cells_mean = mean(donor_n_cells),
    reads_per_cell_mean = mean(donor_rpc),
    n_donors_actual = valid_donors,
    total_cells = sum(donor_n_cells),
    total_reads = sum(donor_n_cells * donor_rpc)
  )
}

# ==============================================================================
# 8. LOGISTIC REGRESSION HELPERS (v4.0.0)
# ==============================================================================

#' Fit logistic regression with error handling
#'
#' v4.0.0: Replaced beta regression with logistic regression because the
#' response is binary per-iteration detection (0/1). Uses the same logit link
#' and produces identical coefficient interpretation.
#'
#' @param formula Formula for logistic regression (response should be 0/1)
#' @param data Data frame
#' @return glm model object or NULL on failure
safe_glm_binomial <- function(formula, data) {
  tryCatch({
    glm(formula, data = data, family = binomial(link = "logit"))
  }, error = function(e) {
    log_error("Logistic regression failed: %s", e$message)
    return(NULL)
  })
}

#' Compute McFadden's pseudo-R² for a GLM
#'
#' McFadden's R² = 1 - logLik(model) / logLik(null_model)
#'
#' @param model A fitted glm object
#' @return Numeric pseudo-R² value
mcfadden_r2 <- function(model) {
  if (is.null(model)) return(NA_real_)
  1 - model$deviance / model$null.deviance
}

# ==============================================================================
# 9. MODEL FORMULA CONSTRUCTION (v4.0.0)
# ==============================================================================

#' Build the v4.0.0 unified logistic regression formula
#'
#' The formula decomposes depth into breadth (n_cells) and
#' depth-per-cell (reads_per_cell):
#'
#'   detection ~ log_donors + log_n_cells + log_reads_per_cell +
#'               log_abs_gt_lfc + log_gt_basemean + log_gt_dispersion
#'
#' Note on log(|LFC|): |LFC| is already on a log2 scale (from DESeq2), so
#' this is technically a double-log transform. This is intentional: it
#' implements diminishing marginal returns, where doubling |LFC| produces a
#' constant change in log-odds of detection. This is more realistic than
#' raw |LFC|, where the jump from |LFC|=5 to 10 would have the same impact
#' as the jump from 1.5 to 3.0, despite detection being near ceiling at
#' high |LFC| values.
#'
#' The combined model has NO cell_type term — it tests whether these 6
#' predictors are sufficient without knowing cell type.
#'
#' @param include_cell_type Whether to include cell_type as a predictor
#'        (for LRT comparison only — not the primary model)
#' @param interactions Whether to include cell_type interactions
#' @return Formula object
build_model_formula <- function(include_cell_type = FALSE, interactions = FALSE) {

  predictors <- c("log_donors", "log_n_cells", "log_reads_per_cell",
                   "log_abs_gt_lfc", "log_gt_basemean", "log_gt_dispersion")

  if (include_cell_type && interactions) {
    # Full interaction model: cell_type * (all predictors)
    pred_str <- paste0("cell_type * (", paste(predictors, collapse = " + "), ")")
  } else if (include_cell_type) {
    # Additive model: cell_type + all predictors
    pred_str <- paste0("cell_type + ", paste(predictors, collapse = " + "))
  } else {
    # Per-cell-type or combined model: just predictors
    pred_str <- paste(predictors, collapse = " + ")
  }

  as.formula(paste("mean_detection ~", pred_str))
}

# ==============================================================================
# 10. SHARED GENE SET IDENTIFICATION
# ==============================================================================

#' Identify genes present in the truth set of ALL cell types
#'
#' Takes the fitting data for each cell type and returns the intersection
#' of gene_id values. This ensures the shared set reflects genes that
#' survived the truth set criteria (padj < threshold, |LFC| > threshold)
#' in every cell type's ground truth.
#'
#' @param fitting_data_list Named list of fitting data frames (one per cell type)
#' @return Character vector of shared gene IDs
identify_shared_genes <- function(fitting_data_list) {
  gene_sets <- lapply(fitting_data_list, function(df) {
    unique(df$gene_id)
  })

  shared <- Reduce(intersect, gene_sets)
  log_info("Shared gene set: %d genes (from %s)",
           length(shared),
           paste(sprintf("%s=%d", names(gene_sets), sapply(gene_sets, length)),
                 collapse = ", "))
  return(shared)
}

# ==============================================================================
# 11. PREDICTION HELPERS
# ==============================================================================

#' Generate prediction grid for v3.1.0 visualisation
#'
#' Creates a grid over n_donors × n_cells (or n_donors × reads_per_cell)
#' with fixed gene covariates.
#'
#' @param vary_axis "n_cells" or "reads_per_cell" — which depth axis to vary
#' @param fixed_depth Value for the fixed depth axis
#' @param lfc_value Representative |LFC| value
#' @param basemean_value Representative baseMean value
#' @param dispersion_value Representative dispersion value
#' @param n_points Number of points per axis for smooth curves
#' @return Data frame with prediction grid
generate_prediction_grid <- function(vary_axis = "n_cells",
                                      fixed_depth = NULL,
                                      lfc_value = DEFAULT_REPRESENTATIVE_LFC,
                                      basemean_value = DEFAULT_REPRESENTATIVE_BASEMEAN,
                                      dispersion_value = DEFAULT_REPRESENTATIVE_DISPERSION,
                                      n_points = 50) {

  donor_seq <- exp(seq(log(min(DONOR_STEPS)), log(max(DONOR_STEPS)), length.out = n_points))

  if (vary_axis == "n_cells") {
    cell_steps <- CELL_STEPS[!is.na(CELL_STEPS)]
    vary_seq <- exp(seq(log(min(cell_steps)), log(max(cell_steps)), length.out = n_points))
    default_rpc <- 6000  # Reasonable default
    fixed_val <- if (!is.null(fixed_depth)) fixed_depth else default_rpc
  } else {
    # vary reads_per_cell
    vary_seq <- exp(seq(log(2000), log(15000), length.out = n_points))
    default_cells <- 250  # Reasonable default
    fixed_val <- if (!is.null(fixed_depth)) fixed_depth else default_cells
  }

  grid <- expand.grid(
    n_donors = donor_seq,
    vary_value = vary_seq
  )

  grid$log_donors <- log(grid$n_donors)
  grid$log_abs_gt_lfc <- log(lfc_value)
  grid$log_gt_basemean <- log(basemean_value)
  grid$log_gt_dispersion <- log(dispersion_value)

  if (vary_axis == "n_cells") {
    grid$n_cells <- grid$vary_value
    grid$log_n_cells <- log(grid$n_cells)
    grid$reads_per_cell <- fixed_val
    grid$log_reads_per_cell <- log(fixed_val)
  } else {
    grid$reads_per_cell <- grid$vary_value
    grid$log_reads_per_cell <- log(grid$reads_per_cell)
    grid$n_cells <- fixed_val
    grid$log_n_cells <- log(fixed_val)
  }

  grid$vary_value <- NULL
  return(grid)
}

# ==============================================================================
# 12. RESULT FORMATTING
# ==============================================================================

#' Format model coefficients as a clean data frame
#'
#' v4.0.0: Updated for glm objects (was betareg).
#'
#' @param model glm model object
#' @param model_name Name for identification
#' @return Data frame with coefficient details
format_coefficients <- function(model, model_name) {
  if (is.null(model)) return(NULL)

  s <- summary(model)
  coef_df <- as.data.frame(s$coefficients)
  coef_df$term <- rownames(coef_df)
  coef_df$model <- model_name
  rownames(coef_df) <- NULL

  # Rename columns for clarity
  names(coef_df)[names(coef_df) == "Estimate"] <- "estimate"
  names(coef_df)[names(coef_df) == "Std. Error"] <- "std_error"
  names(coef_df)[names(coef_df) == "z value"] <- "z_value"
  names(coef_df)[names(coef_df) == "Pr(>|z|)"] <- "p_value"

  return(coef_df[, c("model", "term", "estimate", "std_error", "z_value", "p_value")])
}

# ==============================================================================
# 13. PROGRESS TRACKING
# ==============================================================================

#' Create a progress tracker closure
#'
#' @param total Total number of items
#' @param label Label for progress messages
#' @return Function that increments and prints progress
create_progress <- function(total, label = "Processing") {
  count <- 0
  last_pct <- -1
  function() {
    count <<- count + 1
    pct <- floor(count / total * 100)
    if (pct %% 10 == 0 && pct != last_pct) {
      log_info("%s: %d/%d (%d%%)", label, count, total, pct)
      last_pct <<- pct
    }
  }
}

# ==============================================================================
# 14. SAFE DIRECTORY CREATION
# ==============================================================================

#' Create directory if it doesn't exist
#'
#' @param path Directory path
safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    log_info("Created directory: %s", path)
  }
}

# ==============================================================================
# 15. DONOR METADATA LOADING
# ==============================================================================

#' Load donor metadata from manifest
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

  # MINI_RUN: cap the donor universe up front so the per-shard scan only touches
  # a handful of shards. Take 2x the target for the sex-balanced trim below.
  if (exists("MINI_RUN") && isTRUE(MINI_RUN)) {
    cap <- 2L * MINI_MAX_DONORS
    if (nrow(manifest_df) > cap) {
      manifest_df <- manifest_df[seq_len(cap), , drop = FALSE]
    }
  }

  donor_meta <- do.call(rbind, lapply(manifest_df$file_path, function(f) {
    if (!file.exists(f)) {
      f <- file.path(shard_dir, basename(f))
    }
    if (!file.exists(f)) return(NULL)

    tryCatch({
      tmp <- qs::qread(f)
      cd <- SummarizedExperiment::colData(tmp)
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

  # MINI_RUN: trim to a deterministic, sex-balanced universe of MINI_MAX_DONORS,
  # seeded with SEED_BASE so the GT split matches the cartesian mini run.
  if (exists("MINI_RUN") && isTRUE(MINI_RUN) &&
      !is.null(donor_meta) && nrow(donor_meta) > MINI_MAX_DONORS) {
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
  }

  return(donor_meta)
}
