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
# 3. DESEQ2 WORKER FUNCTION
# ==============================================================================

#' Run DESeq2 differential expression analysis
#' 
#' @param agg_result Result from aggregate_from_shards()
#' @param sex_col Column name for sex variable
#' @param age_col Column name for age variable
#' @param cat_vars Character vector of categorical covariate column names
#' @param save_path_results Path to save DESeq2 results
#' @param save_path_log Path to save sampling log
#' @param min_gene_count Minimum total count for gene filtering
#' @return "SUCCESS" or error message
run_deseq_worker <- function(agg_result, 
                              sex_col = SEX_COL, 
                              age_col = AGE_COL,
                              cat_vars = CATEGORICAL_VARS,
                              save_path_results, 
                              save_path_log,
                              min_gene_count = MIN_GENE_COUNT) {
  
  # Unpack aggregation result
  pb_sce <- agg_result$sce
  sampling_log <- agg_result$sampling_log
  
  # Save sampling log
  qs::qsave(sampling_log, save_path_log)
  
  # Prepare colData
  cts <- as.data.frame(SummarizedExperiment::colData(pb_sce))
  
  # Parse and standardize age
  raw_age <- as.character(cts[[age_col]])
  num_age <- parse_age_stage_to_years(raw_age)
  
  # Check if age has variance (need at least 2 unique values)
  age_has_variance <- length(unique(num_age[!is.na(num_age)])) > 1
  
  if (age_has_variance) {
    cts$age_scaled <- standardize_age(num_age, impute_na = TRUE)
  } else {
    # If all ages are the same, don't include age in model
    log_warning("Age has no variance, excluding from model")
  }
  
  # Identify valid categorical covariates (present and have >1 level)
  # Also check that adding the covariate doesn't make the model rank-deficient
  valid_vars <- c()
  for (v in cat_vars) {
    if (v %in% colnames(cts) && length(unique(cts[[v]])) > 1) {
      # Check if this variable is not perfectly collinear with sex
      sex_var_table <- table(cts[[sex_col]], cts[[v]])
      # If any row or column is all zeros (except diagonal), it's collinear
      if (all(rowSums(sex_var_table > 0) > 1) && all(colSums(sex_var_table > 0) > 1)) {
        valid_vars <- c(valid_vars, v)
      } else {
        log_warning("Covariate '%s' appears collinear with sex, excluding", v)
      }
    }
  }
  
  # Build formula - start with just sex
  formula_parts <- c(sex_col)
  
  # Add age if it has variance
  if (age_has_variance) {
    formula_parts <- c(formula_parts, "age_scaled")
  }
  
  # Add valid categorical covariates
  formula_parts <- c(formula_parts, valid_vars)
  
  formula_str <- paste0("~ ", paste(formula_parts, collapse = " + "))
  
  # Create DESeqDataSet with error handling for rank-deficient models
  dds <- tryCatch({
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(pb_sce, "counts"),
      colData = cts,
      design = as.formula(formula_str)
    ))
  }, error = function(e) {
    # If model is rank-deficient, try simpler model with just sex
    log_warning("Model '%s' failed: %s. Trying sex-only model.", formula_str, e$message)
    suppressMessages(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(pb_sce, "counts"),
      colData = cts,
      design = as.formula(paste0("~ ", sex_col))
    ))
  })
  
  # Filter low-count genes
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_gene_count, ]
  
  # Run DESeq2 with error handling
  dds <- tryCatch({
    suppressMessages(DESeq2::DESeq(dds, parallel = FALSE, quiet = TRUE))
  }, error = function(e) {
    # If DESeq fails, try with simpler model
    log_warning("DESeq2 failed with current model: %s. Trying sex-only model.", e$message)
    DESeq2::design(dds) <- as.formula(paste0("~ ", sex_col))
    suppressMessages(DESeq2::DESeq(dds, parallel = FALSE, quiet = TRUE))
  })
  
  # Get coefficient name for sex contrast
  coef_name <- paste0(sex_col, "_", MALE_VAL, "_vs_", FEMALE_VAL)
  
  # LFC shrinkage (lfcShrink is in DESeq2, apeglm is just the method)
  if (coef_name %in% DESeq2::resultsNames(dds)) {
    res <- suppressMessages(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE))
  } else {
    # Fallback to normal shrinkage if coefficient not found
    log_warning("Coefficient '%s' not found, using normal shrinkage", coef_name)
    res <- suppressMessages(DESeq2::lfcShrink(dds, contrast = c(sex_col, MALE_VAL, FEMALE_VAL), 
                                               type = "normal", quiet = TRUE))
  }
  
  # Save results
  qs::qsave(res, save_path_results)
  
  return("SUCCESS")
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
# 9. SCRIPT DIRECTORY DETECTION
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
