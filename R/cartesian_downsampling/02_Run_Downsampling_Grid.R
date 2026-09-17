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
  library(edgeR)
  library(limma)
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
  paths_by_method <- get_paths_by_method(cell_type)
  paths <- paths_by_method[[DE_METHODS[1]]]
  
  # Check results directory exists (should be created by Step 1)
  if (!dir.exists(paths$results_dir)) {
    log_warning("Results directory not found for %s, run Step 1 first", cell_type)
    next
  }

  # Fast-skip via sentinel: if all downsampling jobs already complete, skip the
  # expensive donor-metadata scan (~1,300 qs::qread calls per cell type).
  grid_sentinels <- vapply(paths_by_method, function(p) {
    file.path(p$results_dir, ".downsampling_grid_complete")
  }, character(1))
  if (all(file.exists(grid_sentinels)) && !FORCE_RECALC) {
    log_info("Downsampling sentinels found for all methods -- skipping all runs for %s",
             cell_type)
    next
  }

  # Track whether any work was actually performed in this cell type. If no jobs
  # remain across all runs, we can write the completion sentinel.
  any_jobs_remaining <- FALSE

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
  
  for (run_id in RUN_IDS) {

    log_info("--- Run %d (range %d-%d of %d total) ---",
             run_id, RUN_START, RUN_END, N_ROBUSTNESS_RUNS)
    
    run_dirs <- vapply(paths_by_method, function(p) {
      file.path(p$results_dir, paste0("Run_", run_id))
    }, character(1))
    run_dir <- run_dirs[[DE_METHODS[1]]]
    
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
    
    if (is_hpc_chunk_active()) {
      before_chunk <- length(plan)
      plan <- filter_plan_for_hpc_chunk(
        plan = plan,
        cell_type = cell_type,
        run_id = run_id
      )
      log_info("HPC chunk %d/%d selected %d of %d valid jobs from %s",
               JOB_CHUNK, JOB_N_CHUNKS, length(plan), before_chunk,
               JOB_MANIFEST_PATH)
    }

    # Filter out already completed jobs
    method_output_paths <- function(p, methods = DE_METHODS) {
      out <- file.path(run_dirs[methods], basename(p$out))
      names(out) <- methods
      out
    }
    method_log_paths <- function(p, methods = DE_METHODS) {
      out <- file.path(run_dirs[methods], basename(p$log_out))
      names(out) <- methods
      out
    }

    chunk_valid_jobs <- length(plan)

    if (!FORCE_RECALC) {
      plan <- plan[!sapply(plan, function(p) {
        all(file.exists(method_output_paths(p)))
      })]
    }

    jobs_to_run <- length(plan)
    already_done <- chunk_valid_jobs - jobs_to_run
    
    log_info("Job Stats: Total: %d | Valid: %d | Skipped (insufficient donors): %d | In chunk: %d | Already done: %d | To run: %d",
             total_possible, valid_jobs, skipped_insufficient, chunk_valid_jobs,
             already_done, jobs_to_run)
    
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
      any_jobs_remaining <- TRUE

      log_info("Starting parallel execution on %d cores...", N_CORES)
      
      # Set run seed for reproducible sampling
      run_seed <- SEED_BASE + run_id
      
      results <- bplapply(plan, function(p) {
        
        # Generate unique seed for this job via hash of a composite key.
        # Replaces the previous additive form (run_seed + n_donors*1000 + y_val + iter),
        # which could alias different (donors, y, iter) tuples once target_reads
        # exceeded 1000 in reads mode. Hashing is deterministic across sessions
        # and collision-resistant in practice.
        y_val <- if (p$mode == "reads") p$target_reads else p$n_cells
        job_key <- paste(run_seed, p$n_donors, p$mode,
                         ifelse(is.na(y_val), "NA", y_val), p$iter, sep = "|")
        job_seed <- digest::digest2int(job_key)
        
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
                return(setNames(rep("FAILED_SAMPLING", length(DE_METHODS)), DE_METHODS))
              }
              
              methods_to_run <- DE_METHODS
              out_paths <- method_output_paths(p)
              log_paths <- method_log_paths(p)
              if (!FORCE_RECALC) {
                methods_to_run <- methods_to_run[!file.exists(out_paths[methods_to_run])]
              }
              if (length(methods_to_run) == 0) {
                return(setNames(rep("SUCCESS", length(DE_METHODS)), DE_METHODS))
              }

              # Run all enabled DE methods on the shared pseudobulk matrix
              status <- run_all_de_methods(
                agg_result = agg_result,
                methods = methods_to_run,
                save_path_results = out_paths[methods_to_run],
                save_path_log = log_paths[methods_to_run]
              )
              
              # Clean up
              rm(agg_result)
              gc()
              
              return(status)
            })
          })
        }, error = function(e) {
          return(setNames(rep(sprintf("ERROR: %s", e$message), length(DE_METHODS)), DE_METHODS))
        })
        
        return(result)
        
      }, BPPARAM = BPPARAM)
      
      # ==============================================================================
      # 6. SUMMARIZE RESULTS
      # ==============================================================================
      
      res_vec <- unlist(results)
      success_codes <- c("SUCCESS", "SUCCESS_FALLBACK")
      n_fallback <- sum(res_vec == "SUCCESS_FALLBACK")
      n_success  <- sum(res_vec %in% success_codes)
      n_fail     <- sum(!res_vec %in% success_codes)

      log_info("Batch complete. Success: %d | Failed: %d", n_success, n_fail)
      method_table <- table(
        method = sub("\\..*$", "", names(res_vec)),
        status = ifelse(res_vec %in% success_codes, "success", "failed")
      )
      log_info("Method status summary:")
      for (method in rownames(method_table)) {
        n_method_success <- if ("success" %in% colnames(method_table)) {
          method_table[method, "success"]
        } else {
          0L
        }
        n_method_failed <- if ("failed" %in% colnames(method_table)) {
          method_table[method, "failed"]
        } else {
          0L
        }
        log_info("  %s: success=%d | failed=%d",
                 method,
                 n_method_success,
                 n_method_failed)
      }
      if (n_success > 0) {
        log_info("  fallback status: %d / %d successful method jobs (%.3f%%)",
                 n_fallback, n_success, 100 * n_fallback / n_success)
      }

      if (n_fail > 0) {
        # Log unique error types (SUCCESS_FALLBACK is a success, not an error)
        error_table <- table(res_vec[!res_vec %in% success_codes])
        log_warning("Error summary:")
        for (err in names(error_table)) {
          log_warning("  %s: %d occurrences", err, error_table[err])
        }
      }
      
    } else {
      log_info("No new jobs to run")
    }
    
  } # End run loop

  # Write completion sentinel if no further jobs were needed in this pass and
  # no errors occurred. Errors set any_jobs_remaining via the bplapply branch
  # above (we ran jobs); successful idempotent re-runs leave it FALSE.
  if (!any_jobs_remaining && !is_hpc_chunk_active()) {
    for (method in DE_METHODS) {
      grid_sentinel <- file.path(paths_by_method[[method]]$results_dir,
                                 ".downsampling_grid_complete")
      if (!file.exists(grid_sentinel)) {
        writeLines(c(
          sprintf("completed_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
          sprintf("cell_type: %s",   cell_type),
          sprintf("method: %s",      method),
          sprintf("n_runs: %d",      N_ROBUSTNESS_RUNS),
          sprintf("mode: %s",        DOWNSAMPLE_MODE),
          sprintf("seed_base: %s",   SEED_BASE)
        ), grid_sentinel)
        log_info("Wrote downsampling-grid sentinel: %s", grid_sentinel)
      }
    }
  } else if (!any_jobs_remaining && is_hpc_chunk_active()) {
    log_info("HPC chunk run complete for %s; not writing global grid sentinel",
             cell_type)
  }

  log_info("Completed downsampling grid for %s", cell_type)

} # End cell type loop

log_section("STEP 2 COMPLETE")
log_info("Downsampling grid executed for all cell types and runs")
