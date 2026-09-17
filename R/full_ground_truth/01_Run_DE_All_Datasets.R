# ==============================================================================
# SCRIPT:  01_Run_DE_All_Datasets.R
# PURPOSE: Run DESeq2 on ALL 7 datasets x 3 cell types using ALL donors
#          (no splitting/sampling) to establish full ground truth DE results.
# INPUT:   Donor shards from Datalink/Shards/{CellType}/{Dataset}/donor_*.qs
# OUTPUT:  Per dataset: _DE_results.csv, _Full_DESeq2_Results.qs,
#          _dds_checkpoint.qs, _pb_checkpoint.qs in Full_Ground_Truth/{CellType}/
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
  library(BiocParallel)
})

# ==============================================================================
# 1. SOURCE CONFIGURATION
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/full_ground_truth/config.R",
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

if (file.exists(SANITY_CHECKS_SCRIPT)) {
  source(SANITY_CHECKS_SCRIPT)
} else {
  log_warn("Sanity checks script not found: %s", SANITY_CHECKS_SCRIPT)
}

# ==============================================================================
# 2. SETUP
# ==============================================================================

log_section("STEP 1: Run DE on All Datasets (Full Ground Truth)")
log_info("Output base: %s", OUTPUT_BASE)
log_info("Shard base:  %s", SHARD_BASE)
log_info("Cell types:  %s", paste(CELL_TYPES, collapse = ", "))
log_info("Datasets:    %s", paste(RUN_DATASETS, collapse = ", "))
log_info("Force recalc: %s", FORCE_RECALC)

de_param <- setup_parallel(N_CORES_DE)

# ==============================================================================
# 3. MAIN LOOP: Cell Type -> Dataset
# ==============================================================================

# Track results across all cell types
global_results <- list()

for (cell_type in CELL_TYPES) {

  log_section(sprintf("CELL TYPE: %s", cell_type))

  paths <- get_paths(cell_type)
  output_dir <- paths$output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_info("Created output directory: %s", output_dir)
  }

  ct_successes <- character(0)
  ct_failures  <- character(0)

  for (d_name in RUN_DATASETS) {

    if (!d_name %in% names(DATASET_CONFIG)) {
      log_warn("Dataset '%s' not in DATASET_CONFIG. Skipping.", d_name)
      next
    }

    cfg <- DATASET_CONFIG[[d_name]]
    shard_dir <- get_shard_dir(cell_type, d_name)

    # ------------------------------------------------------------------
    # 3a. Check shard directory exists
    # ------------------------------------------------------------------
    if (!dir.exists(shard_dir)) {
      log_info("[SKIP] %s/%s - shard directory not found", d_name, cell_type)
      next
    }

    log_section(sprintf("Processing: %s / %s", d_name, cell_type))

    # ------------------------------------------------------------------
    # 3b. Check if DE results already exist
    # ------------------------------------------------------------------
    de_results_file <- file.path(output_dir, paste0(d_name, "_DE_results.csv"))

    if (file.exists(de_results_file) && !FORCE_RECALC) {
      log_info("[SKIP] DE results already exist: %s", basename(de_results_file))
      ct_successes <- c(ct_successes, d_name)
      next
    }

    # ------------------------------------------------------------------
    # 3c-3h. Aggregate, validate, run DESeq2, save
    # ------------------------------------------------------------------
    tryCatch({

      # --- Pseudobulk aggregation or cache load ---
      pb_checkpoint_file <- file.path(output_dir, paste0(d_name, "_pb_checkpoint.qs"))
      pb_sce <- NULL

      # Check for cached pseudobulk
      if (file.exists(pb_checkpoint_file) && !FORCE_RECALC) {
        log_info("[CHECKPOINT] Loading cached pseudobulk...")
        pb_sce <- qs::qread(pb_checkpoint_file)

        if (!validate_cached_pseudobulk(pb_sce, cfg)) {
          log_warn("Cached pseudobulk invalid. Regenerating...")
          pb_sce <- NULL
        } else {
          log_info("[CACHE VALID] Using cached pseudobulk (%d donors, %d genes)",
                   ncol(pb_sce), nrow(pb_sce))
        }
      }

      # Aggregate from shards if no valid cache
      if (is.null(pb_sce)) {
        log_info("Aggregating shards from: %s", shard_dir)

        vars_requested <- unique(c(
          cfg$cols$donor, cfg$cols$sex, cfg$cols$age,
          cfg$covariates, cfg$numeric_covariates
        ))

        pb_sce <- aggregate_shards_parallel(
          shard_dir   = shard_dir,
          design_vars = vars_requested,
          n_cores     = N_CORES_AGG
        )

        # Sanity check aggregated pseudobulk
        if (exists("sanity_check") && exists("pseudobulk_checks")) {
          pb_valid <- sanity_check(
            pb_sce,
            sprintf("Pseudobulk for %s/%s", d_name, cell_type),
            pseudobulk_checks()
          )
          if (!pb_valid) {
            stop("Pseudobulk failed sanity checks")
          }
        }

        # Save pseudobulk checkpoint
        qs::qsave(pb_sce, pb_checkpoint_file)
        log_info("[CHECKPOINT] Saved pseudobulk (%d donors, %d genes)",
                 ncol(pb_sce), nrow(pb_sce))
      }

      # --- Run DESeq2 ---
      log_info("Running DESeq2 for %s/%s...", d_name, cell_type)

      de_result <- run_deseq_for_dataset(
        pb_sce         = pb_sce,
        cfg            = cfg,
        de_param       = de_param,
        min_gene_count = MIN_GENE_COUNT
      )

      # --- Save all outputs ---
      save_de_outputs(
        res_df       = de_result$res_df,
        dds          = de_result$dds,
        res          = de_result$results,
        pb_sce       = pb_sce,
        output_dir   = output_dir,
        dataset_name = d_name,
        cell_type    = cell_type
      )

      ct_successes <- c(ct_successes, d_name)
      log_info("[SUCCESS] %s/%s complete", d_name, cell_type)

    }, error = function(e) {
      log_error("[FAILED] %s/%s: %s", d_name, cell_type, e$message)
      ct_failures <<- c(ct_failures, d_name)
    })

    # ------------------------------------------------------------------
    # 3j. Clean up memory
    # ------------------------------------------------------------------
    if (exists("pb_sce", inherits = FALSE))   rm(pb_sce)
    if (exists("de_result", inherits = FALSE)) rm(de_result)
    gc(verbose = FALSE)

  } # end dataset loop

  # Store per-cell-type summary
  global_results[[cell_type]] <- list(
    successes = ct_successes,
    failures  = ct_failures
  )

} # end cell type loop

# ==============================================================================
# 4. COMPLETION REPORT
# ==============================================================================

log_section("PIPELINE COMPLETION REPORT")

total_success <- 0
total_failure <- 0

for (ct in names(global_results)) {
  res <- global_results[[ct]]
  n_s <- length(res$successes)
  n_f <- length(res$failures)
  total_success <- total_success + n_s
  total_failure <- total_failure + n_f

  log_info("%s: %d succeeded, %d failed", ct, n_s, n_f)
  if (n_s > 0) log_info("  Successes: %s", paste(res$successes, collapse = ", "))
  if (n_f > 0) log_error("  Failures:  %s", paste(res$failures, collapse = ", "))
}

log_info("Total: %d succeeded, %d failed across %d cell types",
         total_success, total_failure, length(CELL_TYPES))
log_info("Output location: %s", OUTPUT_BASE)
log_section("STEP 1 COMPLETE")
