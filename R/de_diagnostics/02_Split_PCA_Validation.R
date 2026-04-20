# ==============================================================================
# SCRIPT: GT-vs-Validation PCA Diagnostics
# PURPOSE: Validate cartesian downsampling 50:50 donor splits by generating
#          PCA plots coloured by split membership and key covariates.
#          Verifies no systematic bias was introduced by the random split.
# INPUT: Donor shards from Datalink/Shards/{CellType}/Roussos/
#        GT metadata from Datalink/Cartesian_Downsampling/.../Ground_Truth_sampling_metadata.qs
# OUTPUT: PCA plots in Datalink/DE_Diagnostics/PCA/{mode}/{CellType}/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(qs)
  library(Matrix)
  library(BiocParallel)
  library(viridis)
  library(RColorBrewer)
})

# Load shared utilities for prettified labels
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Load config - skip if already loaded
if (!exists("SCRIPT_DIR") || !grepl("de_diagnostics", SCRIPT_DIR)) {
  possible_config_paths <- c(
    "R/de_diagnostics/config.R",
    "config.R",
    "../config.R"
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
    stop("Could not find config.R. Run from project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 2: GT-vs-VALIDATION PCA DIAGNOSTICS")
log_info("Modes: %s", paste(CARTESIAN_MODES, collapse = ", "))
log_info("Runs: %d | Cell types: %s", N_ROBUSTNESS_RUNS, paste(CELL_TYPES, collapse = ", "))
log_info("PCA covariates: %s", paste(PCA_COVARIATES, collapse = ", "))
log_info("Top variable genes: %d", N_PCA_GENES)

# ==============================================================================
# 2. MAIN PROCESSING LOOP
# ==============================================================================

total_plots <- 0
total_skipped <- 0
total_failed <- 0

for (cell_type in CELL_TYPES) {

  log_section(sprintf("Cell Type: %s", cell_type))

  # Get shard directory and all donor IDs
  shard_dir <- get_shard_dir(cell_type)

  if (!dir.exists(shard_dir)) {
    log_warn("Shard directory not found: %s — skipping %s", shard_dir, cell_type)
    next
  }

  all_donor_ids <- get_all_donor_ids(shard_dir)
  if (is.null(all_donor_ids) || length(all_donor_ids) == 0) {
    log_warn("No donors found for %s — skipping", cell_type)
    next
  }

  log_info("Total donors available: %d", length(all_donor_ids))

  # ============================================================================
  # 2a. AGGREGATE PSEUDOBULK ONCE PER CELL TYPE
  # ============================================================================
  # We aggregate ALL donors once, then label GT/Val per run.
  # This avoids re-aggregating for every run.

  log_info("Aggregating pseudobulk for all %d donors...", length(all_donor_ids))
  pb_sce <- aggregate_pseudobulk(all_donor_ids, shard_dir, n_cores = N_CORES)

  if (is.null(pb_sce)) {
    log_warn("Pseudobulk aggregation failed for %s — skipping", cell_type)
    next
  }

  log_info("Pseudobulk: %d genes x %d donors",
           nrow(pb_sce), ncol(pb_sce))

  # Parse age
  if (AGE_COL %in% colnames(colData(pb_sce))) {
    colData(pb_sce)$Age_Years <- parse_age_to_years(colData(pb_sce)[[AGE_COL]])
    n_na_age <- sum(is.na(colData(pb_sce)$Age_Years))
    if (n_na_age > 0) {
      log_info("  %d donors have NA age after parsing", n_na_age)
    }
  } else {
    log_warn("  Age column '%s' not found — Age PCA will be skipped", AGE_COL)
    colData(pb_sce)$Age_Years <- NA_real_
  }

  # ============================================================================
  # 2b. VST NORMALISATION (once per cell type)
  # ============================================================================

  log_info("Running VST normalisation...")
  dds <- tryCatch({
    DESeqDataSetFromMatrix(
      countData = assay(pb_sce, "counts"),
      colData   = colData(pb_sce),
      design    = ~ 1
    )
  }, error = function(e) {
    log_error("DESeqDataSet creation failed: %s", e$message)
    return(NULL)
  })

  if (is.null(dds)) next

  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) >= MIN_GENE_COUNTS, ]
  log_info("  Genes after filtering: %d", nrow(dds))

  vst_out <- tryCatch({
    vst(dds, blind = TRUE)
  }, error = function(e) {
    log_error("VST failed: %s", e$message)
    return(NULL)
  })

  if (is.null(vst_out)) next

  # ============================================================================
  # 2c. COMPUTE PCA (once per cell type)
  # ============================================================================

  log_info("Computing PCA (top %d genes)...", N_PCA_GENES)

  # Build intgroup from available covariates
  available_vars <- intersect(
    c(SEX_COL, "Age_Years", SOURCE_COL),
    colnames(colData(vst_out))
  )
  # Filter to vars that actually have non-NA values
  available_vars <- Filter(function(v) {
    !all(is.na(colData(vst_out)[[v]]))
  }, available_vars)

  if (length(available_vars) == 0) {
    log_warn("No valid metadata variables for PCA — using sex only")
    available_vars <- SEX_COL
  }

  pca_data <- tryCatch({
    plotPCA(vst_out, intgroup = available_vars,
            returnData = TRUE, ntop = N_PCA_GENES)
  }, error = function(e) {
    log_error("PCA computation failed: %s", e$message)
    return(NULL)
  })

  if (is.null(pca_data)) next

  percentVar <- round(100 * attr(pca_data, "percentVar"))
  log_info("  PC1: %d%%, PC2: %d%%", percentVar[1], percentVar[2])

  # ============================================================================
  # 2d. GENERATE PCA PLOTS FOR EACH MODE x RUN
  # ============================================================================

  for (mode in CARTESIAN_MODES) {
    log_info("--- Mode: %s ---", mode)

    for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {

      # Get GT donor IDs for this run
      gt_meta_path <- get_gt_metadata_path(cell_type, mode, run_id)
      gt_ids <- extract_gt_donor_ids(gt_meta_path)

      if (is.null(gt_ids) || length(gt_ids) == 0) {
        log_warn("  Run %d: no GT metadata found — skipping", run_id)
        total_skipped <- total_skipped + 1
        next
      }

      # Create output directory
      pca_out_dir <- get_pca_output_dir(mode, cell_type)
      if (!dir.exists(pca_out_dir)) {
        dir.create(pca_out_dir, recursive = TRUE)
      }

      # Check skip-if-exists
      if (!FORCE_RECALC) {
        # Check for new-style filenames with _free/_scaled suffixes
        expected_files <- c(
          paste0("Run_", run_id, "_PCA_Split_free.png"),
          paste0("Run_", run_id, "_PCA_Split_scaled.png"),
          paste0("Run_", run_id, "_PCA_Sex_free.png"),
          paste0("Run_", run_id, "_PCA_Sex_scaled.png"),
          paste0("Run_", run_id, "_PCA_Age_free.png"),
          paste0("Run_", run_id, "_PCA_Age_scaled.png"),
          paste0("Run_", run_id, "_PCA_Source_free.png"),
          paste0("Run_", run_id, "_PCA_Source_scaled.png")
        )
        existing <- file.exists(file.path(pca_out_dir, expected_files))
        if (all(existing)) {
          log_info("  Run %d: plots already exist — skipping", run_id)
          total_skipped <- total_skipped + length(expected_files)
          next
        }
      }

      # Label donors as GT or Validation
      pca_run <- pca_data

      # Robust matching: strip non-alphanumeric chars for comparison
      pca_names <- as.character(pca_run$name)
      clean_pca <- gsub("[^[:alnum:]]", "", pca_names)
      clean_gt <- gsub("[^[:alnum:]]", "", gt_ids)

      # Also try stripping "donor" prefix for numeric-only matching
      clean_pca_nums <- gsub("(?i)donor", "", clean_pca)
      clean_gt_nums <- gsub("(?i)donor", "", clean_gt)

      is_gt <- clean_pca %in% clean_gt | clean_pca_nums %in% clean_gt_nums
      pca_run$Split <- ifelse(is_gt, "Ground_Truth", "Validation")

      n_gt <- sum(is_gt)
      n_val <- sum(!is_gt)

      if (n_gt == 0) {
        log_warn("  Run %d: 0 GT matches (GT example: '%s', PCA example: '%s')",
                 run_id, gt_ids[1], pca_names[1])
        total_failed <- total_failed + 1
        next
      }

      log_info("  Run %d: GT=%d, Val=%d", run_id, n_gt, n_val)

      # --- PLOT 1: SPLIT (with density contours) ---
      tryCatch({
        pretty_ct <- prettify_label(cell_type)
        
        p_split_base <- ggplot(pca_run, aes(x = PC1, y = PC2, color = Split)) +
          geom_point(size = 2, alpha = 0.6) +
          geom_density_2d(linewidth = 0.6, bins = 5) +
          scale_color_manual(values = SPLIT_COLOURS,
                             labels = c("Ground_Truth" = "Ground Truth", 
                                        "Validation" = "Validation")) +
          theme_publication(base_size = 16) +
          labs(
            title = sprintf("%s — Run %d: GT vs Validation", pretty_ct, run_id),
            subtitle = sprintf("GT (n = %d) | Val (n = %d)", n_gt, n_val),
            x = sprintf("PC1: %d%%", percentVar[1]),
            y = sprintf("PC2: %d%%", percentVar[2])
          ) +
          theme(legend.position = "bottom")
        
        # Save FREE version
        ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Split_free.png", run_id)),
               p_split_base, width = PCA_PLOT_WIDTH, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
        
        # Save SCALED version (coord_fixed)
        p_split_scaled <- p_split_base + coord_fixed(ratio = 1)
        ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Split_scaled.png", run_id)),
               p_split_scaled, width = PCA_PLOT_WIDTH + 2, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
        
        total_plots <- total_plots + 2
      }, error = function(e) {
        log_warn("  Run %d Split plot failed: %s", run_id, e$message)
        total_failed <<- total_failed + 1
      })

      # --- PLOT 2: SEX ---
      if (SEX_COL %in% colnames(pca_run)) {
        tryCatch({
          p_sex_base <- ggplot(pca_run, aes(x = PC1, y = PC2,
                                        color = .data[[SEX_COL]],
                                        shape = Split)) +
            geom_point(size = 2.5, alpha = 0.7) +
            scale_color_brewer(palette = SOURCE_PALETTE) +
            scale_shape_manual(values = c("Ground_Truth" = 16, "Validation" = 17),
                               labels = c("Ground_Truth" = "Ground Truth", 
                                          "Validation" = "Validation")) +
            theme_publication(base_size = 16) +
            labs(
              title = sprintf("%s — Run %d: Sex", pretty_ct, run_id),
              x = sprintf("PC1: %d%%", percentVar[1]),
              y = sprintf("PC2: %d%%", percentVar[2])
            ) +
            theme(legend.position = "bottom")
          
          # Save FREE version
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Sex_free.png", run_id)),
                 p_sex_base, width = PCA_PLOT_WIDTH, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          # Save SCALED version
          p_sex_scaled <- p_sex_base + coord_fixed(ratio = 1)
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Sex_scaled.png", run_id)),
                 p_sex_scaled, width = PCA_PLOT_WIDTH + 2, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          total_plots <- total_plots + 2
        }, error = function(e) {
          log_warn("  Run %d Sex plot failed: %s", run_id, e$message)
          total_failed <<- total_failed + 1
        })
      }

      # --- PLOT 3: AGE ---
      if ("Age_Years" %in% colnames(pca_run) && !all(is.na(pca_run$Age_Years))) {
        tryCatch({
          p_age_base <- ggplot(pca_run, aes(x = PC1, y = PC2,
                                        color = Age_Years,
                                        shape = Split)) +
            geom_point(size = 2.5, alpha = 0.7) +
            scale_color_viridis_c(option = "magma", end = 0.9, na.value = "grey80") +
            scale_shape_manual(values = c("Ground_Truth" = 16, "Validation" = 17),
                               labels = c("Ground_Truth" = "Ground Truth", 
                                          "Validation" = "Validation")) +
            theme_publication(base_size = 16) +
            labs(
              title = sprintf("%s — Run %d: Age", pretty_ct, run_id),
              color = "Age (Years)",
              x = sprintf("PC1: %d%%", percentVar[1]),
              y = sprintf("PC2: %d%%", percentVar[2])
            ) +
            theme(legend.position = "bottom")
          
          # Save FREE version
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Age_free.png", run_id)),
                 p_age_base, width = PCA_PLOT_WIDTH, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          # Save SCALED version
          p_age_scaled <- p_age_base + coord_fixed(ratio = 1)
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Age_scaled.png", run_id)),
                 p_age_scaled, width = PCA_PLOT_WIDTH + 2, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          total_plots <- total_plots + 2
        }, error = function(e) {
          log_warn("  Run %d Age plot failed: %s", run_id, e$message)
          total_failed <<- total_failed + 1
        })
      }

      # --- PLOT 4: SOURCE ---
      if (SOURCE_COL %in% colnames(pca_run)) {
        tryCatch({
          p_source_base <- ggplot(pca_run, aes(x = PC1, y = PC2,
                                           color = .data[[SOURCE_COL]],
                                           shape = Split)) +
            geom_point(size = 2.5, alpha = 0.7) +
            scale_color_brewer(palette = SOURCE_PALETTE) +
            scale_shape_manual(values = c("Ground_Truth" = 16, "Validation" = 17),
                               labels = c("Ground_Truth" = "Ground Truth", 
                                          "Validation" = "Validation")) +
            theme_publication(base_size = 16) +
            labs(
              title = sprintf("%s — Run %d: Source", pretty_ct, run_id),
              x = sprintf("PC1: %d%%", percentVar[1]),
              y = sprintf("PC2: %d%%", percentVar[2])
            ) +
            theme(legend.position = "bottom")
          
          # Save FREE version
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Source_free.png", run_id)),
                 p_source_base, width = PCA_PLOT_WIDTH, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          # Save SCALED version
          p_source_scaled <- p_source_base + coord_fixed(ratio = 1)
          ggsave(file.path(pca_out_dir, sprintf("Run_%d_PCA_Source_scaled.png", run_id)),
                 p_source_scaled, width = PCA_PLOT_WIDTH + 2, height = PCA_PLOT_HEIGHT, dpi = PLOT_DPI)
          
          total_plots <- total_plots + 2
        }, error = function(e) {
          log_warn("  Run %d Source plot failed: %s", run_id, e$message)
          total_failed <<- total_failed + 1
        })
      }

    } # end run loop
  } # end mode loop

  # Clean up large objects per cell type
  rm(pb_sce, dds, vst_out, pca_data)
  gc()

} # end cell type loop

# ==============================================================================
# 3. COMPLETION REPORT
# ==============================================================================

log_section("STEP 2 COMPLETION REPORT")
log_info("Total PCA plots generated: %d", total_plots)
log_info("Total skipped (already exist): %d", total_skipped)
log_info("Total failed: %d", total_failed)

expected_total <- length(CELL_TYPES) * length(CARTESIAN_MODES) *
  N_ROBUSTNESS_RUNS * length(PCA_COVARIATES)
log_info("Expected total: %d (3 cell types x %d modes x %d runs x %d covariates)",
         expected_total, length(CARTESIAN_MODES), N_ROBUSTNESS_RUNS,
         length(PCA_COVARIATES))

if (total_failed > 0) {
  log_warn("%d plots failed — check warnings above", total_failed)
}

log_info("Step 2 complete")
