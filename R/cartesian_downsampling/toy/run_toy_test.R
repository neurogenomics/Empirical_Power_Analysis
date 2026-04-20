#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Toy Integration Test
# PURPOSE: Validate end-to-end pipeline with minimal data
# STAGE: 2 (Integration Test)
# ==============================================================================
#
# This test runs the complete pipeline with minimal parameters:
# - 1 cell type (Micro)
# - 1 robustness run
# - 3 donor steps (5, 10, 25)
# - 3 cell steps (50, 100, NA)
# - 2 iterations per grid point
#
# Expected outputs:
# - 1 ground truth file
# - 18 subsample results (3 x 3 x 2)
# - Aggregated summary
# - Heatmaps for visual inspection
#
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 2 TOY INTEGRATION TEST\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# ==============================================================================

# ==============================================================================
# 1. SETUP - OVERRIDE CONFIG FOR TOY TEST
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(stringr)
  library(parallel)
})

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",  # From project root
    "../config.R",                         # From toy/ directory
    "config.R"                             # From script directory
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
    stop("Could not find config.R. Ensure working directory is project root or toy/ directory.")
  }
}

# Source utilities
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# 2. TOY TEST PARAMETERS (Override defaults)
# ==============================================================================

cat("--- Toy Test Configuration ---\n\n")

# Override for toy test
TOY_CELL_TYPES <- c("Micro")  # Smallest dataset
TOY_N_ROBUSTNESS_RUNS <- 1
TOY_DONOR_STEPS <- c(5, 10, 25)
TOY_CELL_STEPS <- c(50, 100, NA)
TOY_N_SUBSAMPLES <- 2
TOY_N_CORES <- min(4, N_CORES)  # Limit cores for toy test

# Toy test outputs go to {OUTPUT_BASE}/TOY_{CellType}_{Dataset}/
# This keeps them separate from production runs

cat(sprintf("Cell types:      %s\n", paste(TOY_CELL_TYPES, collapse = ", ")))
cat(sprintf("Robustness runs: %d\n", TOY_N_ROBUSTNESS_RUNS))
cat(sprintf("Donor steps:     %s\n", paste(TOY_DONOR_STEPS, collapse = ", ")))
cat(sprintf("Cell steps:      %s\n", paste(TOY_CELL_STEPS, collapse = ", ")))
cat(sprintf("Subsamples:      %d\n", TOY_N_SUBSAMPLES))
cat(sprintf("Cores:           %d\n", TOY_N_CORES))
cat(sprintf("Grid points:     %d\n", length(TOY_DONOR_STEPS) * length(TOY_CELL_STEPS)))
cat(sprintf("Total jobs:      %d\n", length(TOY_DONOR_STEPS) * length(TOY_CELL_STEPS) * TOY_N_SUBSAMPLES))
cat("\n")

# Update BPPARAM for toy test
BPPARAM <- MulticoreParam(workers = TOY_N_CORES, progressbar = TRUE, stop.on.error = FALSE)

# ==============================================================================
# 3. VALIDATION TRACKING
# ==============================================================================

validation_results <- list()

validate <- function(name, condition, details = "") {
  status <- if (condition) "PASS" else "FAIL"
  validation_results[[name]] <<- list(status = status, details = details)
  cat(sprintf("[%s] %s", status, name))
  if (nzchar(details)) cat(sprintf(" - %s", details))
  cat("\n")
  return(condition)
}

# ==============================================================================
# 4. RUN PIPELINE STEPS
# ==============================================================================

start_time <- Sys.time()

for (cell_type in TOY_CELL_TYPES) {
  
  cat("\n")
  log_section(sprintf("TOY TEST: %s", cell_type))
  
  # Toy test paths - separate from production
  toy_subdir <- paste0("TOY_", cell_type, "_", SHARD_DATASET)
  paths <- list(
    shard_dir = file.path(SHARD_BASE, cell_type, SHARD_DATASET),
    manifest_path = file.path(SHARD_BASE, cell_type, SHARD_DATASET, "00_MASTER_MANIFEST.qs"),
    results_dir = file.path(OUTPUT_BASE, toy_subdir),
    plots_dir = file.path(OUTPUT_BASE, toy_subdir, "Plots")
  )
  
  # Validate input exists
  if (!validate("Manifest exists", file.exists(paths$manifest_path), paths$manifest_path)) {
    next
  }
  
  # Create output directories
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 1: Generate Ground Truth
  # --------------------------------------------------------------------------
  
  log_info("STEP 1: Generate Ground Truth")
  
  donor_meta <- load_donor_metadata(paths$manifest_path, paths$shard_dir)
  validate("Donor metadata loaded", !is.null(donor_meta) && nrow(donor_meta) > 0,
           sprintf("%d donors", nrow(donor_meta)))
  
  sex_props <- get_sex_proportions(donor_meta)
  validate("Both sexes present", sex_props$n_male > 0 && sex_props$n_female > 0,
           sprintf("%d M, %d F", sex_props$n_male, sex_props$n_female))
  
  m_pool <- donor_meta$donor_id[donor_meta$sex == MALE_VAL]
  f_pool <- donor_meta$donor_id[donor_meta$sex == FEMALE_VAL]
  
  for (run_id in seq_len(TOY_N_ROBUSTNESS_RUNS)) {
    
    run_dir <- file.path(paths$results_dir, paste0("Run_", run_id))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    
    gt_path <- file.path(run_dir, "Ground_Truth.qs")
    gt_log_path <- file.path(run_dir, "Ground_Truth_sampling_metadata.qs")
    
    run_seed <- SEED_BASE + run_id
    set.seed(run_seed)
    
    n_gt_male <- floor(length(m_pool) * GT_SPLIT)
    n_gt_female <- floor(length(f_pool) * GT_SPLIT)
    gt_ids <- c(sample(m_pool, n_gt_male), sample(f_pool, n_gt_female))
    
    log_info("Aggregating ground truth (%d donors)...", length(gt_ids))
    
    agg_result <- aggregate_from_shards(gt_ids, paths$shard_dir, NA, run_seed, TRUE)
    validate("Ground truth aggregation", !is.null(agg_result),
             sprintf("%d genes x %d samples", nrow(agg_result$sce), ncol(agg_result$sce)))
    
    log_info("Running DESeq2 for ground truth...")
    
    status <- tryCatch({
      run_deseq_worker(agg_result, SEX_COL, AGE_COL, CATEGORICAL_VARS, gt_path, gt_log_path, MIN_GENE_COUNT)
    }, error = function(e) e$message)
    
    validate("Ground truth DESeq2", status == "SUCCESS", status)
    validate("Ground truth file created", file.exists(gt_path))
    
    # Save validation donors
    val_ids <- donor_meta$donor_id[!(donor_meta$donor_id %in% gt_ids)]
    val_meta <- donor_meta[donor_meta$donor_id %in% val_ids, ]
    qs::qsave(val_meta, file.path(run_dir, "Validation_Donors.qs"))
    
    # --------------------------------------------------------------------------
    # STEP 2: Run Downsampling Grid
    # --------------------------------------------------------------------------
    
    log_info("STEP 2: Run Downsampling Grid")
    
    plan <- generate_job_plan(
      val_df = val_meta,
      donor_steps = TOY_DONOR_STEPS,
      cell_steps = TOY_CELL_STEPS,
      n_subsamples = TOY_N_SUBSAMPLES,
      prop_male = sex_props$prop_male,
      run_dir = run_dir,
      shard_dir = paths$shard_dir
    )
    
    validate("Job plan generated", length(plan) > 0, sprintf("%d jobs", length(plan)))
    
    log_info("Executing %d jobs in parallel...", length(plan))
    
    results <- bplapply(plan, function(p) {
      job_seed <- run_seed + p$n_donors * 1000 + ifelse(is.na(p$n_cells), 9999, p$n_cells) + p$iter
      
      tryCatch({
        suppressMessages({
          donors <- sample_donors_for_job(p, seed = job_seed)
          agg <- aggregate_from_shards(donors, p$shard_dir, p$n_cells, job_seed, FALSE)
          if (is.null(agg)) return("FAILED_SAMPLING")
          run_deseq_worker(agg, SEX_COL, AGE_COL, CATEGORICAL_VARS, p$out, p$log_out, MIN_GENE_COUNT)
        })
      }, error = function(e) sprintf("ERROR: %s", e$message))
    }, BPPARAM = BPPARAM)
    
    res_vec <- unlist(results)
    n_success <- sum(res_vec == "SUCCESS")
    n_fail <- sum(res_vec != "SUCCESS")
    
    validate("Grid execution", n_success > 0, sprintf("%d success, %d failed", n_success, n_fail))
    
    # Count output files
    result_files <- list.files(run_dir, pattern = "^N\\d+_C.*_I\\d+\\.qs$")
    validate("Result files created", length(result_files) > 0, sprintf("%d files", length(result_files)))
  }
  
  # --------------------------------------------------------------------------
  # STEP 3: Aggregate Results
  # --------------------------------------------------------------------------
  
  log_info("STEP 3: Aggregate Results")
  
  checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
  
  run_folders <- list.dirs(paths$results_dir, recursive = FALSE, full.names = TRUE)
  run_folders <- run_folders[grepl("Run_", basename(run_folders))]
  
  # Simplified aggregation for toy test
  all_results <- list()
  all_xist <- list()
  all_reads <- list()
  
  for (rf in run_folders) {
    run_id <- basename(rf)
    gt_path <- file.path(rf, "Ground_Truth.qs")
    
    if (!file.exists(gt_path)) next
    
    gt_res <- as.data.frame(qs::qread(gt_path))
    gt_all_sig <- rownames(gt_res)[!is.na(gt_res$padj) & gt_res$padj < FDR_THRESHOLD]
    
    sub_files <- list.files(rf, pattern = "^N\\d+_C.*_I\\d+\\.qs$", full.names = TRUE)
    
    for (sf in sub_files) {
      fname <- basename(sf)
      meta <- str_match(fname, "N(\\d+)_C(\\d+|NA)_I(\\d+)")
      n_donors <- as.numeric(meta[2])
      n_cells <- meta[3]
      iter <- as.numeric(meta[4])
      
      sub_res <- tryCatch(as.data.frame(qs::qread(sf)), error = function(e) NULL)
      if (is.null(sub_res)) next
      
      sub_all_sig <- rownames(sub_res)[!is.na(sub_res$padj) & sub_res$padj < FDR_THRESHOLD]
      
      # Simple recall calculation (Total only for toy test)
      for (bin_name in names(LFC_INTERVALS)) {
        bounds <- LFC_INTERVALS[[bin_name]]
        
        gt_in_bin <- rownames(gt_res)[
          !is.na(gt_res$padj) & gt_res$padj < FDR_THRESHOLD &
          abs(gt_res$log2FoldChange) >= bounds[1] & abs(gt_res$log2FoldChange) < bounds[2]
        ]
        
        sub_in_bin <- rownames(sub_res)[
          !is.na(sub_res$padj) & sub_res$padj < FDR_THRESHOLD &
          abs(sub_res$log2FoldChange) >= bounds[1] & abs(sub_res$log2FoldChange) < bounds[2]
        ]
        
        tp <- length(intersect(gt_in_bin, sub_all_sig))
        fn <- length(setdiff(gt_in_bin, sub_all_sig))
        recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
        
        tp_p <- length(intersect(sub_in_bin, gt_all_sig))
        fp_p <- length(setdiff(sub_in_bin, gt_all_sig))
        precision <- if ((tp_p + fp_p) > 0) tp_p / (tp_p + fp_p) else NA_real_
        
        all_results[[length(all_results) + 1]] <- data.frame(
          Run = run_id, Donors = n_donors, Cells = n_cells, Iteration = iter,
          LFC_Bin = bin_name, Gene_Category = "Total",
          Precision = precision, Recall = recall, GT_Count = length(gt_in_bin)
        )
      }
    }
  }
  
  df <- do.call(rbind, all_results)
  df$Cells <- ifelse(df$Cells == "NA", "Max", df$Cells)
  
  validate("Results aggregated", nrow(df) > 0, sprintf("%d rows", nrow(df)))
  
  df_summary <- df %>%
    group_by(Donors, Cells, LFC_Bin, Gene_Category) %>%
    summarise(
      Mean_Precision = mean(Precision, na.rm = TRUE),
      Mean_Recall = mean(Recall, na.rm = TRUE),
      Ref_GT_Count = mean(GT_Count, na.rm = TRUE),
      .groups = "drop"
    )
  
  qs::qsave(list(df_summary = df_summary, df_raw = df), checkpoint_file)
  validate("Checkpoint saved", file.exists(checkpoint_file))
  
  # --------------------------------------------------------------------------
  # STEP 4: Generate Visualization
  # --------------------------------------------------------------------------
  
  log_info("STEP 4: Generate Visualization")
  
  df_summary$Donors <- clean_donors_factor(df_summary$Donors, TOY_DONOR_STEPS)
  df_summary$Cells <- clean_cells_factor(df_summary$Cells, TOY_CELL_STEPS)
  
  # Simple heatmap
  p <- ggplot(df_summary, aes(x = Donors, y = Cells, fill = Mean_Recall)) +
    geom_tile(color = "white") +
    facet_wrap(~LFC_Bin) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50") +
    geom_text(aes(label = sprintf("%.2f", Mean_Recall)), color = "white", size = 3) +
    labs(title = paste0("TOY TEST: Mean Recall - ", cell_type),
         subtitle = "Stage 2 Integration Test",
         x = "Donors", y = "Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_path <- file.path(paths$plots_dir, "TOY_TEST_Recall_Heatmap.png")
  ggsave(plot_path, p, width = 10, height = 6)
  
  validate("Heatmap generated", file.exists(plot_path), plot_path)
  
  # --------------------------------------------------------------------------
  # STEP 5: Cross-Cell-Type Aggregation (Synthetic)
  # --------------------------------------------------------------------------
  # The toy test only runs 1 cell type (Micro). To validate Step 5, we create
  # a synthetic second cell type by duplicating the Micro results with modified
  # values, then run the aggregation logic.
  
  log_info("STEP 5: Cross-Cell-Type Aggregation (Synthetic)")
  
  # Load the checkpoint we just saved
  toy_checkpoint <- qs::qread(checkpoint_file)
  df_raw_micro <- toy_checkpoint$df_raw
  
  if (!is.null(df_raw_micro) && nrow(df_raw_micro) > 0) {
    
    # Tag original as Micro
    df_raw_micro$Cell_Type <- "Micro"
    
    # Create synthetic "Astro" by perturbing Micro values
    set.seed(42)
    df_raw_astro <- df_raw_micro
    df_raw_astro$Cell_Type <- "Astro"
    df_raw_astro$Recall <- pmin(1, pmax(0, df_raw_astro$Recall + runif(nrow(df_raw_astro), -0.1, 0.1)))
    df_raw_astro$Precision <- pmin(1, pmax(0, df_raw_astro$Precision + runif(nrow(df_raw_astro), -0.1, 0.1)))
    df_raw_astro$GT_Count <- pmax(0, df_raw_astro$GT_Count + sample(-3:3, nrow(df_raw_astro), replace = TRUE))
    
    # Ensure Recall_Harm and GT_Harm_Count exist
    if (!"Recall_Harm" %in% colnames(df_raw_micro)) {
      df_raw_micro$Recall_Harm <- df_raw_micro$Recall * 0.9
      df_raw_astro$Recall_Harm <- df_raw_astro$Recall * 0.9
    }
    if (!"GT_Harm_Count" %in% colnames(df_raw_micro)) {
      df_raw_micro$GT_Harm_Count <- round(df_raw_micro$GT_Count * 0.8)
      df_raw_astro$GT_Harm_Count <- round(df_raw_astro$GT_Count * 0.8)
    }
    
    # Combine
    df_combined <- rbind(df_raw_micro, df_raw_astro)
    
    validate("Synthetic combined data created", nrow(df_combined) == 2 * nrow(df_raw_micro),
             sprintf("%d rows (2 x %d)", nrow(df_combined), nrow(df_raw_micro)))
    
    # Source the aggregation function from Step 5
    # (Re-define here to avoid sourcing the full script which would execute)
    compute_pooled_aggregation_toy <- function(df_combined) {
      df_combined %>%
        mutate(
          TP_Recall = Recall * GT_Count,
          Called_Sig = ifelse(!is.na(Precision) & Precision > 0, (Recall * GT_Count) / Precision, NA_real_),
          TP_Precision = ifelse(!is.na(Precision), Recall * GT_Count, NA_real_),
          TP_Recall_Harm = Recall_Harm * GT_Harm_Count,
          GT_Harm_Count_Safe = ifelse(is.na(GT_Harm_Count), 0, GT_Harm_Count)
        ) %>%
        group_by(Donors, Cells, LFC_Bin, Gene_Category) %>%
        summarise(
          Pooled_Recall = ifelse(sum(GT_Count, na.rm = TRUE) > 0,
                                  sum(TP_Recall, na.rm = TRUE) / sum(GT_Count, na.rm = TRUE), 0),
          Pooled_Precision = ifelse(sum(Called_Sig, na.rm = TRUE) > 0,
                                     sum(TP_Precision, na.rm = TRUE) / sum(Called_Sig, na.rm = TRUE), 0),
          Pooled_Recall_Harm = ifelse(sum(GT_Harm_Count_Safe, na.rm = TRUE) > 0,
                                       sum(TP_Recall_Harm, na.rm = TRUE) / sum(GT_Harm_Count_Safe, na.rm = TRUE), 0),
          Total_GT_Observations = sum(GT_Count, na.rm = TRUE),
          N_Cell_Types = n_distinct(Cell_Type),
          N_Observations = n(),
          Cell_Types_Present = paste(sort(unique(Cell_Type)), collapse = "+"),
          .groups = "drop"
        )
    }
    
    # Filter to Total category
    df_total <- df_combined %>% filter(Gene_Category == "Total")
    pooled <- compute_pooled_aggregation_toy(df_total)
    
    validate("Pooled aggregation computed", nrow(pooled) > 0,
             sprintf("%d grid points", nrow(pooled)))
    validate("Pooled recall bounded 0-1", 
             all(pooled$Pooled_Recall >= 0 & pooled$Pooled_Recall <= 1))
    validate("Pooled precision bounded 0-1 or NA",
             all(is.na(pooled$Pooled_Precision) | 
                 (pooled$Pooled_Precision >= 0 & pooled$Pooled_Precision <= 1)))
    validate("Two cell types detected", all(pooled$N_Cell_Types == 2))
    validate("Cell types label correct", all(pooled$Cell_Types_Present == "Astro+Micro"))
    
    # Generate combined heatmap
    toy_global_dir <- file.path(OUTPUT_BASE, paste0("TOY_Global_", SHARD_DATASET), "Plots")
    dir.create(toy_global_dir, recursive = TRUE, showWarnings = FALSE)
    
    pooled$Donors <- clean_donors_factor(pooled$Donors, TOY_DONOR_STEPS)
    pooled$Cells <- clean_cells_factor(pooled$Cells, TOY_CELL_STEPS)
    
    p_agg <- ggplot(pooled, aes(x = Donors, y = Cells, fill = Pooled_Recall)) +
      geom_tile(color = "white") +
      facet_wrap(~LFC_Bin) +
      scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey50") +
      geom_text(aes(label = sprintf("%.2f", Pooled_Recall)), color = "white", size = 3) +
      labs(title = "TOY TEST: Pooled Recall (Synthetic Astro + Micro)",
           subtitle = "Stage 2 Integration Test — Cross-Cell-Type Aggregation",
           x = "Donors", y = "Cells") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    agg_plot_path <- file.path(toy_global_dir, "TOY_TEST_Pooled_Recall_Heatmap.png")
    ggsave(agg_plot_path, p_agg, width = 10, height = 6)
    
    validate("Aggregated heatmap generated", file.exists(agg_plot_path), agg_plot_path)
    
    log_info("Step 5 toy test complete")
    
  } else {
    log_warning("No raw data available for Step 5 toy test (df_raw missing from checkpoint)")
    validate("Step 5 skipped", FALSE, "No raw data in checkpoint")
  }
}

# ==============================================================================
# 5. READS MODE INTEGRATION TEST
# ==============================================================================

log_section("READS MODE INTEGRATION TEST")
log_info("Testing read-depth downsampling with parallel execution")

{
  # Configure reads mode
  DOWNSAMPLE_MODE <- "reads"
  TOY_READ_STEPS <- c(5000, 15000, 40000, NA)
  TOY_N_ITERATIONS <- 2
  
  log_info(sprintf("Mode: %s", DOWNSAMPLE_MODE))
  log_info(sprintf("Read steps: %s", paste(TOY_READ_STEPS, collapse = ", ")))
  log_info(sprintf("Iterations: %d", TOY_N_ITERATIONS))
  
  # Create output directory
  toy_reads_dir <- file.path(TOY_OUTPUT_DIR, sprintf("TOY_%s_%s_reads", TOY_CELLTYPE, TOY_DATASET))
  dir.create(toy_reads_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ==============================================================================
  # Step 1: Generate Ground Truth (reads mode)
  # ==============================================================================
  
  log_info("Step 1: Generating ground truth for reads mode")
  
  gt_reads_path <- file.path(toy_reads_dir, "ground_truth.qs")
  
  if (!file.exists(gt_reads_path)) {
    source(file.path(PIPELINE_DIR, "01_Generate_Ground_Truth.R"))
    
    # Generate ground truth
    gt_reads <- generate_ground_truth(
      shard_dir = TOY_SHARD_DIR,
      output_path = gt_reads_path,
      seed = TOY_SEED
    )
    
    validate("Ground truth generated (reads)", file.exists(gt_reads_path), gt_reads_path)
    log_info(sprintf("Ground truth saved: %d genes, %d donors", 
                     nrow(gt_reads$counts), ncol(gt_reads$counts)))
  } else {
    log_info("Ground truth already exists, loading...")
    gt_reads <- qs::qread(gt_reads_path)
    validate("Ground truth loaded (reads)", !is.null(gt_reads), gt_reads_path)
  }
  
  # ==============================================================================
  # Step 2: Run Downsampling Grid (reads mode with parallel execution)
  # ==============================================================================
  
  log_info("Step 2: Running downsampling grid for reads mode")
  
  source(file.path(PIPELINE_DIR, "02_Run_Downsampling_Grid.R"))
  
  # Build execution plan
  plan_reads <- expand.grid(
    n_cells = NA,
    target_reads = TOY_READ_STEPS,
    iteration = seq_len(TOY_N_ITERATIONS),
    stringsAsFactors = FALSE
  )
  
  plan_reads$mode <- DOWNSAMPLE_MODE
  plan_reads$shard_dir <- TOY_SHARD_DIR
  plan_reads$output_dir <- toy_reads_dir
  plan_reads$donors <- list(TOY_DONORS)
  
  log_info(sprintf("Execution plan: %d jobs", nrow(plan_reads)))
  
  # Execute in parallel (same as production)
  bp_param <- BiocParallel::MulticoreParam(workers = 2, RNGseed = TOY_SEED)
  
  results_reads <- BiocParallel::bplapply(seq_len(nrow(plan_reads)), function(i) {
    p <- plan_reads[i, ]
    
    # Generate job seed
    job_seed <- TOY_SEED + i
    
    # Call aggregate_from_shards with named arguments
    agg_result <- aggregate_from_shards(
      donor_ids = p$donors[[1]],
      shard_dir = p$shard_dir,
      n_cells_to_sample = p$n_cells,
      target_reads = p$target_reads,
      mode = p$mode,
      seed = job_seed,
      verbose = FALSE
    )
    
    # Save output
    output_filename <- sprintf("N%s_R%s_I%d.qs",
                               ifelse(is.na(p$n_cells), "NA", p$n_cells),
                               ifelse(is.na(p$target_reads), "NA", p$target_reads),
                               p$iteration)
    output_path <- file.path(p$output_dir, output_filename)
    
    qs::qsave(agg_result, output_path)
    
    return(list(
      job_id = i,
      output_path = output_path,
      target_reads = p$target_reads,
      actual_reads = sum(agg_result$counts),
      n_genes = nrow(agg_result$counts),
      n_donors = ncol(agg_result$counts)
    ))
  }, BPPARAM = bp_param)
  
  validate("Downsampling grid completed (reads)", length(results_reads) == nrow(plan_reads), 
           sprintf("%d jobs", length(results_reads)))
  
  log_info(sprintf("Generated %d output files", length(results_reads)))
  
  # ==============================================================================
  # Validate Read Depth Differentiation
  # ==============================================================================
  
  log_info("Validating read depth differentiation...")
  
  # Group results by target_reads
  results_by_target <- split(results_reads, sapply(results_reads, function(x) x$target_reads))
  
  # Check file size differentiation
  file_sizes <- sapply(results_reads, function(x) file.size(x$output_path))
  unique_sizes <- length(unique(file_sizes))
  
  validate("File sizes differ across targets", 
           unique_sizes > 1,
           sprintf("%d unique file sizes", unique_sizes))
  
  # Check MD5 checksums
  md5_checksums <- sapply(results_reads, function(x) {
    tools::md5sum(x$output_path)
  })
  unique_checksums <- length(unique(md5_checksums))
  
  validate("MD5 checksums differ across targets",
           unique_checksums == length(results_reads),
           sprintf("%d unique checksums out of %d files", unique_checksums, length(results_reads)))
  
  # Check actual read counts
  actual_reads_by_target <- sapply(results_by_target, function(jobs) {
    mean(sapply(jobs, function(x) x$actual_reads))
  })
  
  log_info("Mean actual reads by target:")
  for (target in names(actual_reads_by_target)) {
    log_info(sprintf("  %s: %.0f reads", target, actual_reads_by_target[target]))
  }
  
  # Verify monotonic relationship
  finite_targets <- names(actual_reads_by_target)[names(actual_reads_by_target) != "NA"]
  finite_targets_num <- as.numeric(finite_targets)
  finite_actual <- actual_reads_by_target[finite_targets]
  
  is_monotonic <- all(diff(finite_actual[order(finite_targets_num)]) > 0)
  
  validate("Actual reads increase with target",
           is_monotonic,
           sprintf("Targets: %s", paste(finite_targets_num[order(finite_targets_num)], collapse=" < ")))
  
  # Check coefficient of variation across targets
  cv_actual <- sd(finite_actual) / mean(finite_actual)
  
  validate("Substantial variation in actual reads",
           cv_actual > 0.3,
           sprintf("CV = %.2f", cv_actual))
  
  # ==============================================================================
  # Step 3: Run DESeq2 (reads mode)
  # ==============================================================================
  
  log_info("Step 3: Running DESeq2 for reads mode")
  
  source(file.path(PIPELINE_DIR, "03_Run_DESeq2.R"))
  
  # Run DESeq2 on all downsampled files
  deseq_results_reads <- BiocParallel::bplapply(results_reads, function(job) {
    agg_data <- qs::qread(job$output_path)
    
    # Run DESeq2
    dds <- run_deseq2(
      counts = agg_data$counts,
      metadata = agg_data$metadata,
      design_formula = ~ Sex
    )
    
    # Extract results
    res <- DESeq2::results(dds, contrast = c("Sex", "M", "F"))
    
    return(list(
      job_id = job$job_id,
      target_reads = job$target_reads,
      n_sig = sum(res$padj < 0.05, na.rm = TRUE),
      n_tested = sum(!is.na(res$padj))
    ))
  }, BPPARAM = bp_param)
  
  validate("DESeq2 completed (reads)", length(deseq_results_reads) == length(results_reads),
           sprintf("%d DESeq2 runs", length(deseq_results_reads)))
  
  # Check that different targets produce different DE results
  n_sig_by_target <- sapply(split(deseq_results_reads, sapply(deseq_results_reads, function(x) x$target_reads)),
                             function(jobs) mean(sapply(jobs, function(x) x$n_sig)))
  
  log_info("Mean significant genes by target:")
  for (target in names(n_sig_by_target)) {
    log_info(sprintf("  %s: %.1f genes", target, n_sig_by_target[target]))
  }
  
  cv_sig <- sd(n_sig_by_target) / mean(n_sig_by_target)
  
  validate("DE results vary across targets",
           cv_sig > 0.1 || length(unique(n_sig_by_target)) > 1,
           sprintf("CV = %.2f", cv_sig))
  
  # ==============================================================================
  # Step 4: Aggregate Results (reads mode)
  # ==============================================================================
  
  log_info("Step 4: Aggregating results for reads mode")
  
  source(file.path(PIPELINE_DIR, "04_Aggregate_Results.R"))
  
  # Create aggregation input
  agg_input_reads <- data.frame(
    target_reads = sapply(deseq_results_reads, function(x) x$target_reads),
    iteration = rep(seq_len(TOY_N_ITERATIONS), each = length(TOY_READ_STEPS)),
    n_sig = sapply(deseq_results_reads, function(x) x$n_sig),
    n_tested = sapply(deseq_results_reads, function(x) x$n_tested)
  )
  
  # Calculate recall (requires ground truth DE genes)
  gt_sig_genes <- rownames(gt_reads$deseq_results)[gt_reads$deseq_results$padj < 0.05]
  
  log_info(sprintf("Ground truth: %d significant genes", length(gt_sig_genes)))
  
  # Compute recall for each job
  recall_reads <- sapply(seq_along(results_reads), function(i) {
    job <- results_reads[[i]]
    agg_data <- qs::qread(job$output_path)
    
    # Run DESeq2
    dds <- run_deseq2(agg_data$counts, agg_data$metadata, ~ Sex)
    res <- DESeq2::results(dds, contrast = c("Sex", "M", "F"))
    
    # Find significant genes
    sig_genes <- rownames(res)[which(res$padj < 0.05)]
    
    # Calculate recall
    tp <- length(intersect(sig_genes, gt_sig_genes))
    recall <- tp / length(gt_sig_genes)
    
    return(recall)
  })
  
  agg_input_reads$recall <- recall_reads
  
  # Aggregate by target
  agg_summary_reads <- agg_input_reads %>%
    group_by(target_reads) %>%
    summarise(
      mean_recall = mean(recall, na.rm = TRUE),
      sd_recall = sd(recall, na.rm = TRUE),
      mean_n_sig = mean(n_sig, na.rm = TRUE),
      .groups = "drop"
    )
  
  log_info("Aggregated recall by target:")
  print(agg_summary_reads)
  
  # Validate recall increases with target reads
  finite_recall <- agg_summary_reads$mean_recall[!is.na(agg_summary_reads$target_reads)]
  finite_targets_recall <- agg_summary_reads$target_reads[!is.na(agg_summary_reads$target_reads)]
  
  is_recall_monotonic <- all(diff(finite_recall[order(finite_targets_recall)]) >= 0)
  
  validate("Recall increases with target reads",
           is_recall_monotonic,
           sprintf("Recall range: %.3f - %.3f", min(finite_recall), max(finite_recall)))
  
  log_info("Reads mode integration test complete")
}

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("================================================================================\n")
cat("  TOY TEST SUMMARY\n")
cat("================================================================================\n")
cat("\n")

n_pass <- sum(sapply(validation_results, function(x) x$status == "PASS"))
n_fail <- sum(sapply(validation_results, function(x) x$status == "FAIL"))

cat(sprintf("  Duration:    %.1f minutes\n", as.numeric(duration)))
cat(sprintf("  Validations: %d passed, %d failed\n", n_pass, n_fail))
cat("\n")

cat("  Validation Details:\n")
for (name in names(validation_results)) {
  r <- validation_results[[name]]
  icon <- if (r$status == "PASS") "[OK]  " else "[FAIL]"
  cat(sprintf("    %s %s", icon, name))
  if (nzchar(r$details)) cat(sprintf(" - %s", r$details))
  cat("\n")
}

cat("\n")

if (n_fail == 0) {
  cat("  STATUS: PASSED\n")
  cat("\n")
  cat("  Outputs for visual inspection:\n")
  toy_output_dir <- paste0("TOY_", TOY_CELL_TYPES[1], "_", SHARD_DATASET)
  cat(sprintf("    - Results: %s\n", file.path(OUTPUT_BASE, toy_output_dir)))
  cat(sprintf("    - Plots:   %s\n", file.path(OUTPUT_BASE, toy_output_dir, "Plots")))
  cat("\n")
  cat("  AWAITING USER APPROVAL before proceeding to Stage 3 (Production)\n")
  toy_status <- "SUCCESS"
} else {
  cat("  STATUS: FAILED\n")
  cat("\n")
  cat("  Fix the above failures before proceeding.\n")
  toy_status <- "FAILED"
}

cat("\n")
cat("================================================================================\n")
cat("\n")

# Finalize logging
if (exists("finalize_logging") && is_logging_active()) {
  finalize_logging(status = toy_status)
}
