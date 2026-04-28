# ==============================================================================
# SCRIPT: Generate Ground Truth DE Results
# PURPOSE: Create ground truth differential expression results using 50% of donors
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Ground_Truth.qs and sampling metadata per run
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(Matrix)
  library(dplyr)
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

log_section("STEP 1: Generate Ground Truth")

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {
  
  log_section(sprintf("Cell Type: %s", cell_type))
  
  # Get paths
  paths <- get_paths(cell_type)

  # Check manifest exists
  if (!file.exists(paths$manifest_path)) {
    log_warning("Manifest not found for %s, skipping", cell_type)
    next
  }

  # Create output directory
  if (!dir.exists(paths$results_dir)) {
    dir.create(paths$results_dir, recursive = TRUE)
    log_info("Created output directory: %s", paths$results_dir)
  }

  # Fast-skip via sentinel: if all GT runs already complete, skip the
  # expensive donor-metadata scan (~1,300 qs::qread calls per cell type).
  gt_sentinel <- file.path(paths$results_dir, ".ground_truth_complete")
  if (file.exists(gt_sentinel) && !FORCE_RECALC) {
    log_info("Sentinel found at %s -- skipping all GT runs for %s",
             gt_sentinel, cell_type)
    next
  }

  # Load donor metadata
  log_info("Loading donor metadata...")
  donor_meta <- load_donor_metadata(paths$manifest_path, paths$shard_dir)
  
  if (is.null(donor_meta) || nrow(donor_meta) == 0) {
    log_warning("No valid donors found for %s, skipping", cell_type)
    next
  }
  
  # Calculate sex proportions
  sex_props <- get_sex_proportions(donor_meta)
  log_info("Donors: %d total (%d M, %d F) | Proportions: %.1f%% M, %.1f%% F",
           sex_props$n_total, sex_props$n_male, sex_props$n_female,
           sex_props$prop_male * 100, sex_props$prop_female * 100)
  
  # Get donor pools by sex
  m_pool <- donor_meta$donor_id[donor_meta$sex == MALE_VAL]
  f_pool <- donor_meta$donor_id[donor_meta$sex == FEMALE_VAL]
  
  # ==============================================================================
  # 3. GENERATE GROUND TRUTH FOR EACH RUN
  # ==============================================================================
  
  for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {
    
    log_info("--- Run %d of %d ---", run_id, N_ROBUSTNESS_RUNS)
    
    # Create run directory
    run_dir <- file.path(paths$results_dir, paste0("Run_", run_id))
    if (!dir.exists(run_dir)) {
      dir.create(run_dir, recursive = TRUE)
    }
    
    # Output paths
    gt_path <- file.path(run_dir, "Ground_Truth.qs")
    gt_log_path <- file.path(run_dir, "Ground_Truth_sampling_metadata.qs")
    
    # Skip if already exists and not forcing recalc
    if (file.exists(gt_path) && !FORCE_RECALC) {
      log_info("Ground truth already exists, skipping")
      next
    }
    
    # Set seed for this run
    run_seed <- SEED_BASE + run_id
    set.seed(run_seed)
    
    # Sample ground truth donors (GT_SPLIT proportion)
    n_gt_male <- floor(length(m_pool) * GT_SPLIT)
    n_gt_female <- floor(length(f_pool) * GT_SPLIT)
    
    gt_males <- sample(m_pool, n_gt_male)
    gt_females <- sample(f_pool, n_gt_female)
    gt_ids <- c(gt_males, gt_females)
    
    log_info("Ground truth: %d donors (%d M, %d F)", 
             length(gt_ids), length(gt_males), length(gt_females))
    
    # Aggregate from shards (use all cells for ground truth)
    log_info("Aggregating pseudobulk counts...")
    agg_result <- aggregate_from_shards(
      donor_ids = gt_ids,
      shard_dir = paths$shard_dir,
      n_cells_to_sample = NA,  # Use all cells
      seed = run_seed,
      verbose = TRUE
    )
    
    if (is.null(agg_result)) {
      log_error("Failed to aggregate ground truth for run %d", run_id)
      next
    }
    
    # Run DESeq2
    log_info("Running DESeq2...")
    status <- tryCatch({
      run_deseq_worker(
        agg_result = agg_result,
        sex_col = SEX_COL,
        age_col = AGE_COL,
        cat_vars = CATEGORICAL_VARS,
        save_path_results = gt_path,
        save_path_log = gt_log_path,
        min_gene_count = MIN_GENE_COUNT
      )
    }, error = function(e) {
      log_error("DESeq2 failed: %s", e$message)
      return(e$message)
    })
    
    log_info("Ground truth status: %s", status)
    
    # Save validation donor list for use by step 2
    val_ids <- donor_meta$donor_id[!(donor_meta$donor_id %in% gt_ids)]
    val_meta <- donor_meta[donor_meta$donor_id %in% val_ids, ]
    qs::qsave(val_meta, file.path(run_dir, "Validation_Donors.qs"))
    
    log_info("Saved %d validation donors for downsampling", nrow(val_meta))
    
    # Clean up
    rm(agg_result)
    gc()
  }
  
  # Write completion sentinel if every run produced both required outputs.
  all_runs_ok <- all(vapply(seq_len(N_ROBUSTNESS_RUNS), function(r) {
    rd <- file.path(paths$results_dir, paste0("Run_", r))
    file.exists(file.path(rd, "Ground_Truth.qs")) &&
      file.exists(file.path(rd, "Validation_Donors.qs"))
  }, logical(1)))

  if (all_runs_ok && !file.exists(gt_sentinel)) {
    writeLines(c(
      sprintf("completed_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      sprintf("cell_type: %s",   cell_type),
      sprintf("n_runs: %d",      N_ROBUSTNESS_RUNS),
      sprintf("gt_split: %s",    GT_SPLIT),
      sprintf("seed_base: %s",   SEED_BASE)
    ), gt_sentinel)
    log_info("Wrote ground-truth sentinel: %s", gt_sentinel)
  }

  log_info("Completed ground truth generation for %s", cell_type)
}

log_section("STEP 1 COMPLETE")
log_info("Ground truth files generated for all cell types and runs")
