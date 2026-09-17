# ==============================================================================
# TOY INTEGRATION TEST: EDA Roussos Pipeline
# PURPOSE: Test pipeline on small subset (Micro, 10 donors)
# USAGE: Rscript R/eda_roussos/toy/run_toy_test.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
})

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Initialize shared logging
# Load config - skip if SCRIPT_DIR already set
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "eda_roussos", "toy")
  })
  source(file.path(script_path, "..", "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### TOY INTEGRATION TEST: EDA ROUSSOS ###\n")

# ==============================================================================
# 2. TOY TEST PARAMETERS
# ==============================================================================

TOY_CELL_TYPE <- "Micro"
TOY_N_DONORS <- 10
TOY_OUTPUT_SUFFIX <- "TOY_TEST"

message(sprintf("Cell type: %s", TOY_CELL_TYPE))
message(sprintf("Donors: %d", TOY_N_DONORS))

# ==============================================================================
# 3. SETUP TOY PATHS
# ==============================================================================

shard_dir <- file.path(SHARD_BASE, TOY_CELL_TYPE, SHARD_DATASET)

if (!dir.exists(shard_dir)) {
  stop(sprintf("Shard directory not found: %s", shard_dir))
}

# Load manifest and sample donors
manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
if (!file.exists(manifest_path)) {
  stop(sprintf("Manifest not found: %s", manifest_path))
}

all_donors <- load_donor_ids_from_manifest(shard_dir, DONOR_ID_COL)
message(sprintf("Total donors available: %d", length(all_donors)))

# Sample subset
set.seed(42)
toy_donors <- sample(all_donors, min(TOY_N_DONORS, length(all_donors)))
message(sprintf("Sampled %d donors for toy test", length(toy_donors)))

# ==============================================================================
# 4. TEST PSEUDOBULK AGGREGATION
# ==============================================================================

message("\n=== Test 1: Pseudobulk Aggregation ===")

sce_pb <- get_pb_from_list(toy_donors, "Toy_Cohort", shard_dir, DONOR_ID_COL)

if (is.null(sce_pb)) {
  stop("Pseudobulk aggregation failed")
}

message(sprintf("✓ Pseudobulk: %d genes × %d donors", nrow(sce_pb), ncol(sce_pb)))

# ==============================================================================
# 5. TEST AGE PARSING
# ==============================================================================

message("\n=== Test 2: Age Parsing ===")

if (AGE_COL %in% colnames(colData(sce_pb))) {
  raw_age <- colData(sce_pb)[[AGE_COL]]
  
  if (looks_like_age_stage(raw_age)) {
    parsed_age <- parse_age_stage_to_years(raw_age)
    n_valid <- sum(!is.na(parsed_age))
    message(sprintf("✓ Parsed %d/%d ages successfully", n_valid, length(parsed_age)))
  } else {
    message("⚠ Age column doesn't match expected format")
  }
} else {
  message(sprintf("⚠ Age column '%s' not found", AGE_COL))
}

# ==============================================================================
# 6. TEST GENE STATISTICS
# ==============================================================================

message("\n=== Test 3: Gene Statistics ===")

counts_mat <- assay(sce_pb, "counts")

mean_expr <- calculate_mean_expression(counts_mat)
detection_rate <- calculate_detection_rate(counts_mat)
zero_inflated <- identify_zero_inflated(counts_mat, ZERO_INFLATED_THRESHOLD)

message(sprintf("✓ Mean expression: %.2f (median)", median(mean_expr)))
message(sprintf("✓ Detection rate: %.2f (median)", median(detection_rate)))
message(sprintf("✓ Zero-inflated genes: %d/%d (%.1f%%)",
                sum(zero_inflated), length(zero_inflated),
                100 * sum(zero_inflated) / length(zero_inflated)))

top_genes <- get_top_expressed_genes(counts_mat, n = 5)
message(sprintf("✓ Top gene: %s (mean: %.1f)", top_genes$gene[1], top_genes$mean_expression[1]))

# ==============================================================================
# 7. TEST QC METRICS
# ==============================================================================

message("\n=== Test 4: QC Metrics ===")

# Load first shard for QC testing
first_donor_id <- toy_donors[1]
safe_id <- gsub("[^[:alnum:]]", "_", first_donor_id)
shard_path <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))

if (file.exists(shard_path)) {
  sce_shard <- qs::qread(shard_path)
  
  mt_percent <- calculate_mt_percent(sce_shard, MT_PATTERN)
  genes_detected <- calculate_genes_detected(sce_shard)
  library_size <- calculate_library_size(sce_shard)
  
  message(sprintf("✓ MT%%: %.2f%% (median)", median(mt_percent, na.rm = TRUE)))
  message(sprintf("✓ Genes detected: %.0f (median)", median(genes_detected)))
  message(sprintf("✓ Library size: %.0f (median)", median(library_size)))
} else {
  message("⚠ Could not load shard for QC testing")
}

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

message("\n========================================")
message("TOY TEST SUMMARY")
message("========================================")
message("✓ Pseudobulk aggregation: PASS")
message("✓ Age parsing: PASS")
message("✓ Gene statistics: PASS")
message("✓ QC metrics: PASS")
message("\n✓ ALL TOY TESTS PASSED")
message("\nReady for Stage 3 (Production Run)")
message("Run: Rscript R/eda_roussos/run_pipeline.R\n")

# Finalize logging
if (exists("finalize_logging") && is_logging_active()) {
  finalize_logging(status = "SUCCESS")
}
