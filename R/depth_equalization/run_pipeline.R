#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Depth Equalization Pipeline
# PURPOSE: Multinomially downsample all cells in Astro/Roussos shards to match
#          Microglia's mean reads-per-cell (~4,426), producing depth-equalized
#          shards for sensitivity analysis of sequencing depth confounding.
# INPUT:   Datalink/Shards/Astro/Roussos/ (1,293 donor shards + manifest)
# OUTPUT:  Datalink/Shards/Astro/Roussos_DepthEq/ (1,293 downsampled shards + manifest)
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================

SCRIPT_DIR <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback for Rscript invocation
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    file.path(PROJECT_ROOT, "R", "depth_equalization")
  }
)

# Source shared logging if available
if (file.exists(LOGGING_PATH)) source(LOGGING_PATH)

# Source pipeline config and utils
source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(qs)
  library(parallel)
})

log_section(sprintf("DEPTH EQUALIZATION PIPELINE v%s", PIPELINE_VERSION))
log_info("Start time: %s", Sys.time())

# ==============================================================================
# 1. VALIDATE SOURCE SHARDS
# ==============================================================================

log_section("STEP 1: Validate Source Shards")

if (!dir.exists(SOURCE_SHARD_DIR)) {
  stop(sprintf("Source shard directory not found: %s", SOURCE_SHARD_DIR))
}

source_manifest_path <- file.path(SOURCE_SHARD_DIR, "00_MASTER_MANIFEST.qs")
if (!file.exists(source_manifest_path)) {
  stop(sprintf("Source manifest not found: %s", source_manifest_path))
}

source_manifest <- qs::qread(source_manifest_path)
source_shard_files <- list.files(SOURCE_SHARD_DIR, pattern = "^donor_.*\\.qs$", 
                                  full.names = TRUE)

n_shards <- length(source_shard_files)
log_info("Found %d source shards in %s", n_shards, SOURCE_SHARD_DIR)
log_info("Manifest reports %d donors, %d total cells", 
         source_manifest$config$n_donors, source_manifest$config$total_cells)

if (n_shards != nrow(source_manifest$manifest)) {
  log_warning("Shard count (%d) != manifest rows (%d). Using file system count.",
              n_shards, nrow(source_manifest$manifest))
}

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

log_section("STEP 2: Create Output Directory")

if (!dir.exists(TARGET_SHARD_DIR)) {
  dir.create(TARGET_SHARD_DIR, recursive = TRUE)
  log_info("Created output directory: %s", TARGET_SHARD_DIR)
} else {
  existing_files <- list.files(TARGET_SHARD_DIR, pattern = "^donor_.*\\.qs$")
  if (length(existing_files) > 0) {
    log_warning("Output directory already contains %d shard files. Will overwrite.", 
                length(existing_files))
  }
}

# ==============================================================================
# 3. DOWNSAMPLE ALL SHARDS
# ==============================================================================

log_section(sprintf("STEP 3: Downsample Shards (fraction = %.4f)", DOWNSAMPLE_FRACTION))

# Sort for deterministic ordering
source_shard_files <- sort(source_shard_files)

# Track progress and statistics
t_start <- Sys.time()

process_one_shard <- function(i) {
  shard_path <- source_shard_files[i]
  shard_basename <- basename(shard_path)
  output_path <- file.path(TARGET_SHARD_DIR, shard_basename)
  
  # Deterministic seed per shard: SEED_BASE + index
  shard_seed <- SEED_BASE + i
  
  tryCatch({
    # Downsample
    sce_ds <- downsample_shard(shard_path, DOWNSAMPLE_FRACTION, shard_seed)
    
    # Compute post-downsampling stats
    reads_per_cell <- Matrix::colSums(assay(sce_ds, "counts"))
    stats <- list(
      donor_id    = as.character(unique(colData(sce_ds)[["donor_id"]])),
      n_cells     = ncol(sce_ds),
      mean_rpc    = mean(reads_per_cell),
      total_reads = sum(reads_per_cell)
    )
    
    # Save
    qs::qsave(sce_ds, output_path)
    
    return(list(success = TRUE, stats = stats, file = shard_basename))
    
  }, error = function(e) {
    return(list(success = FALSE, error = e$message, file = shard_basename))
  })
}

# Run in parallel
log_info("Processing %d shards using %d cores...", n_shards, N_CORES)

results <- mclapply(seq_len(n_shards), process_one_shard, mc.cores = N_CORES)

# Summarise results
n_success <- sum(sapply(results, function(r) r$success))
n_failed <- sum(sapply(results, function(r) !r$success))

log_info("Completed: %d succeeded, %d failed", n_success, n_failed)

if (n_failed > 0) {
  failed_files <- sapply(results[!sapply(results, function(r) r$success)], 
                          function(r) r$file)
  failed_errors <- sapply(results[!sapply(results, function(r) r$success)], 
                           function(r) r$error)
  for (k in seq_along(failed_files)) {
    log_error("  %s: %s", failed_files[k], failed_errors[k])
  }
  stop(sprintf("%d shards failed to process. See errors above.", n_failed))
}

# Collect statistics
all_stats <- lapply(results, function(r) r$stats)
all_mean_rpc <- sapply(all_stats, function(s) s$mean_rpc)
all_total_reads <- sapply(all_stats, function(s) s$total_reads)
all_n_cells <- sapply(all_stats, function(s) s$n_cells)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

log_info("Downsampling completed in %.1f minutes", as.numeric(t_elapsed))
log_info("Post-downsampling statistics:")
log_info("  Mean reads/cell across all shards: %.1f (target: %d)", 
         mean(all_mean_rpc), TARGET_MEAN_RPC)
log_info("  Median reads/cell: %.1f", median(all_mean_rpc))
log_info("  Range of per-donor mean RPC: [%.0f, %.0f]", 
         min(all_mean_rpc), max(all_mean_rpc))
log_info("  Total cells: %d (unchanged from source)", sum(all_n_cells))
log_info("  Total reads: %s (was %s in source)", 
         format(sum(all_total_reads), big.mark = ","),
         format(source_manifest$config$total_cells * SOURCE_MEAN_RPC, big.mark = ","))

# ==============================================================================
# 4. GENERATE MANIFEST
# ==============================================================================

log_section("STEP 4: Generate Manifest")

source_config_for_manifest <- list(
  cell_type           = SOURCE_CELL_TYPE,
  min_age             = source_manifest$config$min_age,
  min_cells           = source_manifest$config$min_cells,
  source_dataset      = SOURCE_DATASET,
  target_mean_rpc     = TARGET_MEAN_RPC,
  downsample_fraction = DOWNSAMPLE_FRACTION
)

manifest <- generate_manifest(TARGET_SHARD_DIR, source_config_for_manifest)

manifest_path <- file.path(TARGET_SHARD_DIR, "00_MASTER_MANIFEST.qs")
qs::qsave(manifest, manifest_path)

log_info("Manifest saved: %s", manifest_path)
log_info("  Donors: %d", manifest$config$n_donors)
log_info("  Total cells: %d", manifest$config$total_cells)

# ==============================================================================
# 5. VALIDATION SPOT-CHECK
# ==============================================================================

log_section("STEP 5: Validation Spot-Check")

# Check 5 random shards
set.seed(42)
check_indices <- sort(sample(seq_len(n_shards), min(5, n_shards)))

for (idx in check_indices) {
  shard_basename <- basename(source_shard_files[idx])
  
  # Load source
  source_stats <- summarise_shard_depth(source_shard_files[idx])
  
  # Load target
  target_path <- file.path(TARGET_SHARD_DIR, shard_basename)
  target_stats <- summarise_shard_depth(target_path)
  
  actual_fraction <- target_stats$mean_rpc / source_stats$mean_rpc
  
  log_info("  %s:", shard_basename)
  log_info("    Source: %d cells, mean RPC = %.0f", 
           source_stats$n_cells, source_stats$mean_rpc)
  log_info("    Target: %d cells, mean RPC = %.0f", 
           target_stats$n_cells, target_stats$mean_rpc)
  log_info("    Actual fraction: %.4f (expected: %.4f)", 
           actual_fraction, DOWNSAMPLE_FRACTION)
  log_info("    Cell count preserved: %s", 
           ifelse(source_stats$n_cells == target_stats$n_cells, "YES", "NO"))
  
  # Sanity checks
  if (abs(actual_fraction - DOWNSAMPLE_FRACTION) > 0.05) {
    log_warning("    Fraction deviation > 5%% for %s!", shard_basename)
  }
  if (source_stats$n_cells != target_stats$n_cells) {
    log_error("    Cell count mismatch for %s!", shard_basename)
  }
}

# ==============================================================================
# 6. SAVE SUMMARY
# ==============================================================================

log_section("STEP 6: Save Summary")

summary_df <- data.frame(
  donor_id    = sapply(all_stats, function(s) s$donor_id),
  n_cells     = all_n_cells,
  mean_rpc    = all_mean_rpc,
  total_reads = all_total_reads,
  stringsAsFactors = FALSE
)

summary_path <- file.path(TARGET_SHARD_DIR, "depth_equalization_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
log_info("Summary CSV saved: %s", summary_path)

# Final summary
log_section("PIPELINE COMPLETE")
log_info("Source: %s (%d donors, mean RPC = %d)", 
         SOURCE_SHARD_DIR, n_shards, SOURCE_MEAN_RPC)
log_info("Target: %s (%d donors, mean RPC ≈ %.0f)", 
         TARGET_SHARD_DIR, n_shards, mean(all_mean_rpc))
log_info("Downsampling fraction: %.4f", DOWNSAMPLE_FRACTION)
log_info("Total processing time: %.1f minutes", as.numeric(t_elapsed))
log_info("Next steps:")
log_info("  1. Run cartesian downsampling with SHARD_DATASET = '%s'", TARGET_DATASET)
log_info("  2. Reuse original GT splits from Roussos runs")
log_info("  3. Refit power curve models and compare coefficients")

# ==============================================================================
# 7. CELL COUNT MATCHING (Astro_Shallow_Fewer Extension)
# ==============================================================================

if (FEWER_EXTENSION$enabled) {
  log_section("STEP 7: Cell Count Matching (Astro_Shallow_Fewer Extension)")
  
  # 7a. Load target (Micro) cell count distribution
  target_manifest_path <- file.path(SHARD_BASE, FEWER_EXTENSION$target_cell_type,
                                    FEWER_EXTENSION$target_dataset,
                                    "00_MASTER_MANIFEST.qs")
  
  if (!file.exists(target_manifest_path)) {
    log_error("Target manifest not found: %s", target_manifest_path)
    log_error("Skipping cell count matching.")
  } else {
    target_summary <- load_cell_count_summary(target_manifest_path)
    log_info("Loaded target distribution: %s (%s) - %d donors, mean %.1f cells",
             FEWER_EXTENSION$target_cell_type, FEWER_EXTENSION$target_dataset,
             nrow(target_summary), mean(target_summary$n_cells))
    
    # 7b. Source distribution (from Step 6 summary)
    source_summary <- summary_df  # Already computed in Step 6
    log_info("Source distribution: %d donors, mean %.1f cells",
             nrow(source_summary), mean(source_summary$n_cells))
    
    # 7c. Compute scaling factor
    matching_stats <- compute_scaling_factor(source_summary, target_summary)
    log_info("Matching statistics:")
    log_info("  Target (Micro) mean cells/donor: %.1f (median: %.1f)",
             matching_stats$target_mean, matching_stats$target_median)
    log_info("  Source (Astro_DepthEq) mean cells/donor: %.1f (median: %.1f)",
             matching_stats$source_mean, matching_stats$source_median)
    log_info("  Scaling factor: %.4f", matching_stats$scaling_factor)
    
    # 7d. Apply cell count matching
    source_shard_dir <- TARGET_SHARD_DIR  # Roussos_DepthEq from Steps 1-6
    target_fewer_dir <- file.path(SHARD_BASE, SOURCE_CELL_TYPE, 
                                  FEWER_EXTENSION$output_dataset)
    
    log_info("Applying cell count matching...")
    log_info("  Source: %s", source_shard_dir)
    log_info("  Target: %s", target_fewer_dir)
    log_info("  Min cells threshold: %d", FEWER_EXTENSION$min_cells_threshold)
    
    t_start_matching <- Sys.time()
    
    match_results <- match_cell_counts(
      source_shard_dir = source_shard_dir,
      target_shard_dir = target_fewer_dir,
      matching_stats = matching_stats,
      min_cells_threshold = FEWER_EXTENSION$min_cells_threshold,
      seed_base = FEWER_EXTENSION$seed_base
    )
    
    t_elapsed_matching <- difftime(Sys.time(), t_start_matching, units = "mins")
    
    # 7e. Generate manifest for Roussos_ShallowFewer
    log_section("STEP 7b: Generate Manifest for Matched Shards")
    
    source_config_for_manifest <- list(
      cell_type = SOURCE_CELL_TYPE,
      min_age = source_manifest$config$min_age,
      min_cells = source_manifest$config$min_cells,
      source_dataset = FEWER_EXTENSION$source_dataset,
      target_dataset = FEWER_EXTENSION$target_dataset,
      target_mean_rpc = TARGET_MEAN_RPC,
      downsample_fraction = DOWNSAMPLE_FRACTION
    )
    
    fewer_manifest <- generate_matched_manifest(
      shard_dir = target_fewer_dir,
      source_config = source_config_for_manifest,
      matching_stats = matching_stats,
      match_results = match_results
    )
    
    fewer_manifest_path <- file.path(target_fewer_dir, "00_MASTER_MANIFEST.qs")
    qs::qsave(fewer_manifest, fewer_manifest_path)
    log_info("Manifest saved: %s", fewer_manifest_path)
    log_info("  Donors retained: %d", fewer_manifest$config$n_donors)
    log_info("  Donors dropped: %d", fewer_manifest$config$donors_dropped)
    log_info("  Total cells: %d", fewer_manifest$config$total_cells)
    
    # 7f. Save matching summary
    summary_match_path <- file.path(target_fewer_dir, 
                                    "cell_count_matching_summary.csv")
    write.csv(match_results, summary_match_path, row.names = FALSE)
    log_info("Matching summary saved: %s", summary_match_path)
    
    # 7g. Report statistics
    log_section("CELL COUNT MATCHING SUMMARY")
    log_info("Matching completed in %.1f minutes", as.numeric(t_elapsed_matching))
    log_info("Donors: %d retained, %d dropped (%.1f%% retained)",
             sum(match_results$success),
             sum(!match_results$success),
             100 * mean(match_results$success))
    
    if (sum(match_results$success) > 0) {
      retained_results <- match_results[match_results$success, ]
      log_info("Cell counts (retained donors):")
      log_info("  Before: mean = %.1f, median = %.1f, range = [%d, %d]",
               mean(retained_results$n_cells_before),
               median(retained_results$n_cells_before),
               min(retained_results$n_cells_before),
               max(retained_results$n_cells_before))
      log_info("  After:  mean = %.1f, median = %.1f, range = [%d, %d]",
               mean(retained_results$n_cells_after),
               median(retained_results$n_cells_after),
               min(retained_results$n_cells_after),
               max(retained_results$n_cells_after))
      log_info("  Target: mean = %.1f, median = %.1f",
               matching_stats$target_mean, matching_stats$target_median)
    }
    
    # 7h. Validation spot-check
    log_section("STEP 7c: Validation Spot-Check (Matched Shards)")
    
    set.seed(42)
    successful_donors <- match_results$donor_id[match_results$success]
    if (length(successful_donors) >= 5) {
      check_donors <- sample(successful_donors, 5)
      
      for (donor in check_donors) {
        shard_basename <- paste0("donor_", donor, ".qs")
        
        source_path <- file.path(source_shard_dir, shard_basename)
        target_path <- file.path(target_fewer_dir, shard_basename)
        
        if (file.exists(source_path) && file.exists(target_path)) {
          source_stats <- summarise_shard_depth(source_path)
          target_stats <- summarise_shard_depth(target_path)
          
          match_info <- match_results[match_results$donor_id == donor, ]
          
          log_info("  %s:", donor)
          log_info("    Source cells: %d, Target cells: %d (expected: %d)",
                   source_stats$n_cells, target_stats$n_cells, 
                   match_info$n_cells_after)
          log_info("    Cell count reduction: %.1f%%",
                   100 * (1 - target_stats$n_cells / source_stats$n_cells))
          log_info("    Read depth preserved: %s",
                   ifelse(abs(target_stats$mean_rpc - source_stats$mean_rpc) < 100,
                            "YES", "WARNING"))
        }
      }
    }
    
    # Final summary with next steps
    log_section("ASTRO_SHALLOW_FEWER EXTENSION COMPLETE")
    log_info("Output: %s", target_fewer_dir)
    log_info("Model name: %s", FEWER_EXTENSION$model_name)
    log_info("Next steps:")
    log_info("  1. Run GT-only regeneration: R/cartesian_gt_only/01_Regenerate_GT.R")
    log_info("  2. Run power curve fitting: R/power_curve_fitting/run_pipeline.R")
    log_info("  3. Compare Astro_Shallow_Fewer vs Astro_Shallow coefficients")
  }
}
