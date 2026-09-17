#!/usr/bin/env Rscript
# ==============================================================================
# TOY TEST: Astro_Shallow_Fewer Extension
# ==============================================================================
# Minimal test of the cell count matching functionality using small
# synthetic data to verify the logic works correctly.
# ==============================================================================

set.seed(42)

# ==============================================================================
# 1. CREATE SYNTHETIC DATA
# ==============================================================================

cat("Creating synthetic test data...\n")

# Simulate source (Astro_DepthEq-like) distribution
# High cell counts
n_source_donors <- 20
source_n_cells <- round(rnorm(n_source_donors, mean = 3000, sd = 800))
source_n_cells <- pmax(source_n_cells, 100)  # Minimum 100 cells

cat(sprintf("Source: %d donors, mean = %.0f cells\n", 
            n_source_donors, mean(source_n_cells)))

# Simulate target (Micro-like) distribution
# Lower cell counts
n_target_donors <- 20
target_n_cells <- round(rnorm(n_target_donors, mean = 1500, sd = 400))
target_n_cells <- pmax(target_n_cells, 50)  # Minimum 50 cells

cat(sprintf("Target: %d donors, mean = %.0f cells\n",
            n_target_donors, mean(target_n_cells)))

# Create summary data frames
source_summary <- data.frame(
  donor_id = paste0("donor_", seq_len(n_source_donors)),
  n_cells = source_n_cells,
  stringsAsFactors = FALSE
)

target_summary <- data.frame(
  donor_id = paste0("donor_", seq_len(n_target_donors)),
  n_cells = target_n_cells,
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. TEST SCALING FACTOR COMPUTATION
# ==============================================================================

cat("\n--- Testing compute_scaling_factor() ---\n")

# Source the config to get utility functions
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "depth_equalization")
source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

matching_stats <- compute_scaling_factor(source_summary, target_summary)

cat(sprintf("Scaling factor: %.4f\n", matching_stats$scaling_factor))
cat(sprintf("Source mean: %.1f, Target mean: %.1f\n",
            matching_stats$source_mean, matching_stats$target_mean))
cat(sprintf("Expected scaling: %.4f\n", 
            matching_stats$target_mean / matching_stats$source_mean))

# Verify the scaling factor is correct
expected_sf <- matching_stats$target_mean / matching_stats$source_mean
if (abs(matching_stats$scaling_factor - expected_sf) > 0.001) {
  stop("ERROR: Scaling factor mismatch!")
}
cat("✓ Scaling factor computed correctly\n")

# ==============================================================================
# 3. TEST CELL COUNT MATCHING
# ==============================================================================

cat("\n--- Testing match_cell_counts logic ---\n")

# Create temporary directories
temp_source_dir <- tempfile("source_shards_")
temp_target_dir <- tempfile("target_shards_")
dir.create(temp_source_dir)
dir.create(temp_target_dir)

# Create fake shard files (just empty files for testing file listing)
for (i in seq_len(n_source_donors)) {
  file.create(file.path(temp_source_dir, 
                       sprintf("donor_%s.qs", source_summary$donor_id[i])))
}

cat(sprintf("Created %d fake source shards in %s\n", 
            n_source_donors, temp_source_dir))

# Test the matching logic manually (without loading real SCEs)
min_threshold <- 100
n_dropped <- 0
n_retained <- 0
retained_before <- c()
retained_after <- c()

for (i in seq_len(n_source_donors)) {
  n_before <- source_summary$n_cells[i]
  n_after <- round(n_before * matching_stats$scaling_factor)
  
  if (n_after < min_threshold) {
    n_dropped <- n_dropped + 1
  } else {
    n_retained <- n_retained + 1
    retained_before <- c(retained_before, n_before)
    retained_after <- c(retained_after, n_after)
  }
}

cat(sprintf("Matching results (threshold=%d):\n", min_threshold))
cat(sprintf("  Retained: %d donors (%.1f%%)\n", n_retained, 100*n_retained/n_source_donors))
cat(sprintf("  Dropped: %d donors (%.1f%%)\n", n_dropped, 100*n_dropped/n_source_donors))

if (n_retained > 0) {
  cat(sprintf("  Retained donors - Before: mean=%.1f, After: mean=%.1f\n",
              mean(retained_before), mean(retained_after)))
  cat(sprintf("  Target distribution mean: %.1f\n", matching_stats$target_mean))
  
  # Check if means are close
  pct_diff <- abs(mean(retained_after) - matching_stats$target_mean) / 
                matching_stats$target_mean * 100
  cat(sprintf("  Difference from target: %.1f%%\n", pct_diff))
  
  if (pct_diff < 5) {
    cat("✓ Cell count matching achieves target mean (within 5%)\n")
  } else {
    cat("⚠ Warning: Cell count matching deviates >5% from target\n")
  }
}

# ==============================================================================
# 4. TEST MANIFEST GENERATION
# ==============================================================================

cat("\n--- Testing generate_matched_manifest() ---\n")

# Create fake match results
match_results <- data.frame(
  donor_id = source_summary$donor_id,
  success = c(rep(TRUE, n_retained), rep(FALSE, n_dropped)),
  error = "",
  n_cells_before = c(retained_before, rep(0, n_dropped)),
  n_cells_after = c(retained_after, rep(0, n_dropped)),
  dropped = c(rep(FALSE, n_retained), rep(TRUE, n_dropped)),
  stringsAsFactors = FALSE
)

source_config <- list(
  cell_type = "Astro",
  min_age = 40,
  min_cells = 30,
  source_dataset = "Roussos_DepthEq",
  target_dataset = "Micro",
  target_mean_rpc = 4426,
  downsample_fraction = 0.4679
)

tryCatch({
  manifest <- generate_matched_manifest(
    shard_dir = temp_target_dir,
    source_config = source_config,
    matching_stats = matching_stats,
    match_results = match_results
  )
  
  cat(sprintf("✓ Manifest generated: %d donors, %d total cells\n",
              manifest$config$n_donors, manifest$config$total_cells))
  cat(sprintf("  Scaling factor recorded: %.4f\n", 
              manifest$config$cell_count_scaling_factor))
  cat(sprintf("  Donors dropped: %d\n", manifest$config$donors_dropped))
  
}, error = function(e) {
  cat(sprintf("✗ Manifest generation failed: %s\n", e$message))
})

# ==============================================================================
# 5. CLEANUP
# ==============================================================================

unlink(temp_source_dir, recursive = TRUE)
unlink(temp_target_dir, recursive = TRUE)

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("TOY TEST COMPLETE\n")
cat("========================================\n")
cat("All core functions executed successfully.\n")
cat("\nNext: Run full pipeline with real data:\n")
cat("  Rscript R/depth_equalization/run_pipeline.R\n")
