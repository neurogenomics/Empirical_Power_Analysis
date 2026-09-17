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
  library(edgeR)
  library(limma)
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
  paths_by_method <- get_paths_by_method(cell_type)
  paths <- paths_by_method[[DE_METHODS[1]]]

  # Check manifest exists
  if (!file.exists(paths$manifest_path)) {
    log_warning("Manifest not found for %s, skipping", cell_type)
    next
  }

  # Create output directories
  for (method in DE_METHODS) {
    if (!dir.exists(paths_by_method[[method]]$results_dir)) {
      dir.create(paths_by_method[[method]]$results_dir, recursive = TRUE)
      log_info("Created output directory: %s", paths_by_method[[method]]$results_dir)
    }
  }

  # Fast-skip the expensive donor-metadata scan (~1,300 qs::qread calls per cell
  # type) only when every REQUESTED run (RUN_IDS) already has GT + validation
  # outputs for every method. Range-aware: a completion sentinel from an earlier
  # 1..5 run must NOT block an incremental 6..10 extension, so we check per-run
  # outputs for the active range rather than the (range-agnostic) sentinel.
  runs_present <- all(vapply(RUN_IDS, function(r) {
    all(vapply(paths_by_method, function(p) {
      rd <- file.path(p$results_dir, paste0("Run_", r))
      file.exists(file.path(rd, "Ground_Truth.qs")) &&
        file.exists(file.path(rd, "Validation_Donors.qs"))
    }, logical(1)))
  }, logical(1)))
  if (runs_present && !FORCE_RECALC) {
    log_info("All requested runs (%d-%d) already complete for %s -- skipping GT scan",
             RUN_START, RUN_END, cell_type)
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
  
  for (run_id in RUN_IDS) {

    log_info("--- Run %d (range %d-%d of %d total) ---",
             run_id, RUN_START, RUN_END, N_ROBUSTNESS_RUNS)
    
    # Create run directories
    run_dirs <- vapply(paths_by_method, function(p) {
      file.path(p$results_dir, paste0("Run_", run_id))
    }, character(1))
    for (run_dir in run_dirs) {
      if (!dir.exists(run_dir)) {
        dir.create(run_dir, recursive = TRUE)
      }
    }
    
    # Output paths (file.path drops names, so re-apply method names for
    # downstream method-keyed indexing in run_all_de_methods)
    gt_paths <- file.path(run_dirs, "Ground_Truth.qs")
    gt_log_paths <- file.path(run_dirs, "Ground_Truth_sampling_metadata.qs")
    names(gt_paths) <- names(run_dirs)
    names(gt_log_paths) <- names(run_dirs)
    
    # Skip if already exists and not forcing recalc
    methods_to_run <- DE_METHODS
    if (!FORCE_RECALC) {
      methods_to_run <- methods_to_run[!file.exists(gt_paths[methods_to_run])]
    }
    if (length(methods_to_run) == 0) {
      log_info("Ground truth already exists for all methods, skipping")
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
    
    # Run all enabled DE methods on the shared pseudobulk matrix
    log_info("Running DE methods: %s", paste(methods_to_run, collapse = ", "))
    status <- tryCatch({
      run_all_de_methods(
        agg_result = agg_result,
        methods = methods_to_run,
        sex_col = SEX_COL,
        age_col = AGE_COL,
        cat_vars = CATEGORICAL_VARS,
        save_path_results = gt_paths[methods_to_run],
        save_path_log = gt_log_paths[methods_to_run],
        min_gene_count = MIN_GENE_COUNT
      )
    }, error = function(e) {
      log_error("DE methods failed: %s", e$message)
      return(e$message)
    })
    
    log_info("Ground truth status: %s",
             paste(sprintf("%s=%s", names(status), status), collapse = " | "))
    
    # Save validation donor list for use by step 2
    val_ids <- donor_meta$donor_id[!(donor_meta$donor_id %in% gt_ids)]
    val_meta <- donor_meta[donor_meta$donor_id %in% val_ids, ]
    for (method in DE_METHODS) {
      qs::qsave(val_meta, file.path(run_dirs[[method]], "Validation_Donors.qs"))
    }
    
    log_info("Saved %d validation donors for downsampling", nrow(val_meta))
    
    # Clean up
    rm(agg_result)
    gc()
  }
  
  # Write completion sentinel if every run produced both required outputs.
  for (method in DE_METHODS) {
    method_paths <- paths_by_method[[method]]
    gt_sentinel <- file.path(method_paths$results_dir, ".ground_truth_complete")
    all_runs_ok <- all(vapply(seq_len(N_ROBUSTNESS_RUNS), function(r) {
      rd <- file.path(method_paths$results_dir, paste0("Run_", r))
      file.exists(file.path(rd, "Ground_Truth.qs")) &&
        file.exists(file.path(rd, "Validation_Donors.qs"))
    }, logical(1)))

    if (all_runs_ok && !file.exists(gt_sentinel)) {
      writeLines(c(
        sprintf("completed_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
        sprintf("cell_type: %s",   cell_type),
        sprintf("method: %s",      method),
        sprintf("n_runs: %d",      N_ROBUSTNESS_RUNS),
        sprintf("gt_split: %s",    GT_SPLIT),
        sprintf("seed_base: %s",   SEED_BASE)
      ), gt_sentinel)
      log_info("Wrote ground-truth sentinel: %s", gt_sentinel)
    }
  }

  log_info("Completed ground truth generation for %s", cell_type)
}

log_section("STEP 1 COMPLETE")
log_info("Ground truth files generated for all cell types and runs")
