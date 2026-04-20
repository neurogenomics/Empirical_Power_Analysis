# ==============================================================================
# SCRIPT: Run Downsampling Grid
# PURPOSE: Execute Cartesian grid of donor x cell downsampling experiments
# INPUT: Ground truth from Step 1, donor shards
# OUTPUT: DESeq2 results for each grid point
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

log_section("STEP 2: Run Downsampling Grid")

# ==============================================================================
# 2. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {
  
  log_section(sprintf("Cell Type: %s", cell_type))
  
  # Get paths
  paths <- get_paths(cell_type)
  
  # Check results directory exists (should be created by Step 1)
  if (!dir.exists(paths$results_dir)) {
    log_warning("Results directory not found for %s, run Step 1 first", cell_type)
    next
  }
  
  # Load donor metadata for sex proportions
  donor_meta <- load_donor_metadata(paths$manifest_path, paths$shard_dir)
  if (is.null(donor_meta)) {
    log_warning("Could not load donor metadata for %s", cell_type)
    next
  }
  
  sex_props <- get_sex_proportions(donor_meta)
  
  # ==============================================================================
  # 3. PROCESS EACH RUN
  # ==============================================================================
  
  for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {
    
    log_info("--- Run %d of %d ---", run_id, N_ROBUSTNESS_RUNS)
    
    run_dir <- file.path(paths$results_dir, paste0("Run_", run_id))
    
    # Check ground truth exists
    gt_path <- file.path(run_dir, "Ground_Truth.qs")
    if (!file.exists(gt_path)) {
      log_warning("Ground truth not found for run %d, skipping", run_id)
      next
    }
    
    # Load validation donors (saved by Step 1)
    val_meta_path <- file.path(run_dir, "Validation_Donors.qs")
    if (!file.exists(val_meta_path)) {
      log_warning("Validation donors not found for run %d, skipping", run_id)
      next
    }
    
    val_df <- qs::qread(val_meta_path)
    log_info("Loaded %d validation donors", nrow(val_df))
    
    # ==============================================================================
    # 4. GENERATE JOB PLAN
    # ==============================================================================
    
    log_info("Building job plan...")
    
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
    
    y_steps <- get_y_steps()
    total_possible <- length(DONOR_STEPS) * length(y_steps) * N_SUBSAMPLES
    valid_jobs <- length(plan)
    skipped_insufficient <- total_possible - valid_jobs
    
    # Filter out already completed jobs
    if (!FORCE_RECALC) {
      plan <- plan[!sapply(plan, function(p) file.exists(p$out))]
    }
    
    jobs_to_run <- length(plan)
    already_done <- valid_jobs - jobs_to_run
    
    log_info("Job Stats: Total: %d | Skipped (insufficient donors): %d | Already done: %d | To run: %d",
             total_possible, skipped_insufficient, already_done, jobs_to_run)
    
    # ==============================================================================
    # 4b. INTERLEAVE JOBS FOR UNIFORM MEMORY USAGE
    # ==============================================================================
    # Without reordering, expand.grid produces all small jobs first and all large
    # jobs last. When N_CORES workers all hit the largest grid points at once,
    # aggregate memory can spike and crash. This interleave sorts by estimated
    # cost (n_donors * y_value), then deals jobs out in a zigzag pattern so that
    # big and small jobs are evenly mixed throughout the queue.
    
    if (jobs_to_run > 1) {
      # Estimate relative cost per job: n_donors * y_value (higher = more memory)
      job_costs <- sapply(plan, function(p) {
        y_val <- if (p$mode == "reads") p$target_reads else p$n_cells
        if (is.na(y_val)) y_val <- 1e7  # "Max" is the most expensive
        p$n_donors * y_val
      })
      
      # Sort by cost ascending
      cost_order <- order(job_costs)
      plan_sorted <- plan[cost_order]
      
      # Zigzag interleave: alternate picking from the cheap end and expensive end
      n <- length(plan_sorted)
      interleaved_idx <- integer(n)
      lo <- 1L
      hi <- n
      for (k in seq_len(n)) {
        if (k %% 2 == 1) {
          interleaved_idx[k] <- lo
          lo <- lo + 1L
        } else {
          interleaved_idx[k] <- hi
          hi <- hi - 1L
        }
      }
      plan <- plan_sorted[interleaved_idx]
      
      log_info("Jobs interleaved for uniform memory usage (cost range: %.0f to %.0f)",
               min(job_costs), max(job_costs))
    }
    
    # ==============================================================================
    # 5. EXECUTE JOBS IN PARALLEL
    # ==============================================================================
    
    if (jobs_to_run > 0) {
      
      log_info("Starting parallel execution on %d cores...", N_CORES)
      
      # Set run seed for reproducible sampling
      run_seed <- SEED_BASE + run_id
      
      results <- bplapply(plan, function(p) {
        
        # Generate unique seed for this job
        # Generate unique seed: use mode-appropriate Y value
        y_val <- if (p$mode == "reads") p$target_reads else p$n_cells
        job_seed <- run_seed + p$n_donors * 1000 + 
                    ifelse(is.na(y_val), 9999, y_val) + 
                    p$iter
        
        result <- tryCatch({
          suppressMessages({
            suppressWarnings({
              
              # Sample donors
              donors <- sample_donors_for_job(p, seed = job_seed)
              
              # Aggregate from shards (mode-aware)
              agg_result <- aggregate_from_shards(
                donor_ids = donors,
                shard_dir = p$shard_dir,
                n_cells_to_sample = p$n_cells,
                target_reads = p$target_reads,
                mode = p$mode,
                seed = job_seed,
                verbose = FALSE
              )
              
              if (is.null(agg_result)) {
                return("FAILED_SAMPLING")
              }
              
              # Run DESeq2
              status <- run_deseq_worker(
                agg_result = agg_result,
                save_path_results = p$out,
                save_path_log = p$log_out
              )
              
              # Clean up
              rm(agg_result)
              gc()
              
              return(status)
            })
          })
        }, error = function(e) {
          return(sprintf("ERROR: %s", e$message))
        })
        
        return(result)
        
      }, BPPARAM = BPPARAM)
      
      # ==============================================================================
      # 6. SUMMARIZE RESULTS
      # ==============================================================================
      
      res_vec <- unlist(results)
      n_success <- sum(res_vec == "SUCCESS")
      n_fail <- sum(res_vec != "SUCCESS")
      
      log_info("Batch complete. Success: %d | Failed: %d", n_success, n_fail)
      
      if (n_fail > 0) {
        # Log unique error types
        error_table <- table(res_vec[res_vec != "SUCCESS"])
        log_warning("Error summary:")
        for (err in names(error_table)) {
          log_warning("  %s: %d occurrences", err, error_table[err])
        }
      }
      
    } else {
      log_info("No new jobs to run")
    }
    
  } # End run loop
  
  log_info("Completed downsampling grid for %s", cell_type)
  
} # End cell type loop

log_section("STEP 2 COMPLETE")
log_info("Downsampling grid executed for all cell types and runs")
