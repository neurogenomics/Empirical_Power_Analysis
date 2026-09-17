# ==============================================================================
# SCRIPT: PCA Diagnostics - Cohort Analysis
# PURPOSE: Generate PCA plots colored by metadata variables to assess batch
#          effects, biological variation, and cohort structure
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: PCA plots in Datalink/Exploratory_Analysis/{CellType}/{Dataset}/PCA/
# ==============================================================================

suppressPackageStartupMessages({

  library(SingleCellExperiment)
  library(qs)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(viridis)
  library(parallel)
})

# Load shared utilities for prettified labels
# Use hardcoded fallback path for Rscript compatibility
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Load config - skip if SCRIPT_DIR already set (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({ dirname(sys.frame(1)$ofile) }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "eda_roussos")
  })
  source(file.path(script_path, "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### STEP 01: PCA DIAGNOSTICS ###\n")

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

process_cell_type_pca <- function(cell_type, current_dataset, current_label) {
  
  message(sprintf("\n=== Processing: %s [%s] ===", cell_type, current_dataset))
  
  # Setup paths
  shard_dir <- file.path(SHARD_BASE, cell_type, current_dataset)
  output_paths <- get_output_paths(cell_type, current_dataset)
  create_output_dirs(output_paths)
  
  # Check if shard directory exists
  if (!dir.exists(shard_dir)) {
    warning(sprintf("Shard directory not found: %s", shard_dir))
    return(NULL)
  }
  
  # Load donor IDs from manifest
  tryCatch({
    all_donors <- load_donor_ids_from_manifest(shard_dir, DONOR_ID_COL)
    message(sprintf("  Found %d donors", length(all_donors)))
  }, error = function(e) {
    warning(sprintf("Failed to load manifest for %s: %s", cell_type, e$message))
    return(NULL)
  })
  
  # Build pseudobulk object
  message("  Aggregating pseudobulk...")
  sce_full <- get_pb_from_list(all_donors, "Full_Cohort", shard_dir, DONOR_ID_COL)
  
  if (is.null(sce_full)) {
    warning(sprintf("Could not aggregate pseudobulk for %s", cell_type))
    return(NULL)
  }
  
  message(sprintf("  Pseudobulk: %d genes × %d donors", 
                  nrow(sce_full), ncol(sce_full)))
  
  # Process age if available
  HAS_AGE <- FALSE
  if (AGE_COL %in% colnames(colData(sce_full))) {
    message("  Parsing age strings...")
    raw_age <- colData(sce_full)[[AGE_COL]]
    if (looks_like_age_stage(raw_age)) {
      colData(sce_full)$Age_Years <- parse_age_stage_to_years(raw_age)
      n_na <- sum(is.na(colData(sce_full)$Age_Years))
      if (n_na > 0) {
        message(sprintf("    %d donors have NA age after parsing", n_na))
      }
      HAS_AGE <- TRUE
    } else {
      message("    Age column doesn't match expected format, skipping")
    }
  } else {
    message(sprintf("  Column '%s' not found, skipping age", AGE_COL))
  }
  
  # Add log2 cell count
  colData(sce_full)$Log2_Cells <- log2(colData(sce_full)$n_cells)
  
  # Filter genes
  message("  Filtering genes...")
  keep <- rowSums(assay(sce_full, "counts")) >= MIN_GENE_COUNTS
  sce_full <- sce_full[keep, ]
  message(sprintf("  Retained %d genes after filtering", nrow(sce_full)))
  
  if (nrow(sce_full) < 100) {
    warning(sprintf("Too few genes remaining for %s", cell_type))
    return(NULL)
  }
  
  # VST normalization
  message("  Running DESeq2 VST normalization...")
  dds <- DESeqDataSetFromMatrix(
    countData = assay(sce_full, "counts"),
    colData = colData(sce_full),
    design = ~ 1
  )
  
  vst_out <- tryCatch({
    vst(dds, blind = TRUE)
  }, error = function(e) {
    warning(sprintf("VST failed for %s: %s", cell_type, e$message))
    return(NULL)
  })
  
  if (is.null(vst_out)) return(NULL)
  
  # Determine which metadata variables are available
  vars_to_check <- PCA_METADATA_COLS
  valid_vars <- vars_to_check[vars_to_check %in% colnames(colData(vst_out))]
  
  if (!HAS_AGE) {
    valid_vars <- setdiff(valid_vars, "Age_Years")
  }
  
  message(sprintf("  PCA will use %d metadata variables: %s", 
                  length(valid_vars), paste(valid_vars, collapse = ", ")))
  
  # Run PCA
  message("  Computing PCA...")
  pca_data <- plotPCA(vst_out, intgroup = valid_vars, 
                      returnData = TRUE, ntop = N_PCA_GENES)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Generate plots for each variable
  message("  Generating PCA plots...")
  
  for (var_name in valid_vars) {

    # Skip if all NA
    if (all(is.na(pca_data[[var_name]]))) {
      message(sprintf("    Skipping %s (all NA)", var_name))
      next
    }

    # Prettify variable name for display
    pretty_var <- prettify_label(var_name)
    pretty_ct <- prettify_label(cell_type)

    # Prettify data values for categorical variables
    plot_data <- pca_data
    if (!is.numeric(plot_data[[var_name]])) {
      plot_data[[var_name]] <- prettify_metadata_values(
        as.character(plot_data[[var_name]]), column = var_name
      )
    }

    # Reorder ethnicity so European (majority) is plotted first, minorities on top
    if (var_name == "self_reported_ethnicity") {
      vals <- plot_data[[var_name]]
      euro_rows <- which(vals == "European")
      non_euro_rows <- which(vals != "European")
      plot_data <- plot_data[c(euro_rows, non_euro_rows), , drop = FALSE]
    }

    # Base plot
    p_base <- ggplot(plot_data, aes(x = PC1, y = PC2, color = !!sym(var_name))) +
      geom_point(size = 1.5, alpha = 0.5) +
      theme_publication(base_size = 16) +
      labs(
        title = sprintf("%s: PCA", pretty_ct),
        x = paste0("PC1: ", percentVar[1], "%"),
        y = paste0("PC2: ", percentVar[2], "%"),
        color = pretty_var
      )

    # Choose color scale based on variable type
    if (is.numeric(plot_data[[var_name]])) {
      p_base <- p_base + scale_color_viridis_c(option = COLOR_PALETTE_CONTINUOUS,
                                                end = 0.9, na.value = "grey80")
    } else {
      n_levels <- length(unique(plot_data[[var_name]]))
      if (n_levels > 20) {
        # Too many levels, hide legend
        p_base <- p_base + theme(legend.position = "none") +
          labs(title = sprintf("%s: PCA (%d levels)", pretty_ct, n_levels))
      } else if (n_levels > 8) {
        # Many levels: use Dark2 for better distinguishability
        p_base <- p_base + scale_color_brewer(palette = "Dark2")
      } else {
        p_base <- p_base + scale_color_brewer(palette = COLOR_PALETTE_CATEGORICAL)
      }
    }
    
    # File naming
    safe_var_name <- gsub("[^[:alnum:]]", "_", var_name)
    safe_label <- gsub("[^[:alnum:]]", "", current_label)
    
    # Save FREE version (default scaling)
    p_free <- p_base
    output_file_free <- file.path(output_paths$pca, 
                                   paste0("PCA_", safe_label, "_", safe_var_name, "_free.png"))
    ggsave(output_file_free, p_free, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    message(sprintf("    Saved: %s", basename(output_file_free)))
    
    # Save SCALED version (coord_fixed for equal axis scaling)
    p_scaled <- p_base + coord_fixed(ratio = 1)
    output_file_scaled <- file.path(output_paths$pca, 
                                     paste0("PCA_", safe_label, "_", safe_var_name, "_scaled.png"))
    ggsave(output_file_scaled, p_scaled, width = PLOT_WIDTH + 2, height = PLOT_HEIGHT)
    message(sprintf("    Saved: %s", basename(output_file_scaled)))
  }
  
  message(sprintf("  ✓ Completed PCA for %s", cell_type))
  return(TRUE)
}

# ==============================================================================
# 3. RUN FOR ALL DATASETS × CELL TYPES
# ==============================================================================

for (current_dataset in SHARD_DATASETS) {
  
  current_label <- COHORT_LABELS[current_dataset]
  message(sprintf("\n--- Dataset: %s (%s) ---", current_dataset, current_label))
  
  # Skip-if-exists: check whether PCA PNGs already exist for ALL cell types
  all_exist <- all(sapply(CELL_TYPES, function(ct) {
    pca_dir <- file.path(OUTPUT_BASE, ct, current_dataset, "PCA")
    length(list.files(pca_dir, pattern = "\\.png$", full.names = FALSE)) > 0
  }))
  
  if (all_exist && !FORCE_RECALC) {
    message(sprintf("[SKIP] PCA outputs already exist for %s (set FORCE_RECALC=TRUE to regenerate)",
                    current_dataset))
    next
  }
  
  message(sprintf("Processing %d cell types: %s", 
                  length(CELL_TYPES), paste(CELL_TYPES, collapse = ", ")))
  
  # Process in parallel
  results <- mclapply(CELL_TYPES, function(ct) {
    tryCatch({
      process_cell_type_pca(ct, current_dataset, current_label)
    }, error = function(e) {
      warning(sprintf("Critical failure for %s: %s", ct, e$message))
      return(NULL)
    })
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  
  # Summary
  n_success <- sum(sapply(results, function(x) !is.null(x) && x == TRUE))
  message(sprintf("\n  %s: %d/%d cell types processed successfully", 
                  current_dataset, n_success, length(CELL_TYPES)))
}

message("\n### STEP 01 COMPLETE ###")
