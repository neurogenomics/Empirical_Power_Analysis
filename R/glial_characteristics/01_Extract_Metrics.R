# ==============================================================================
# SCRIPT: Extract Metrics - Compute Per-Dataset Cell/Donor Characteristics
# PURPOSE: Extract UMIs/cell, genes/cell, and cells/donor for each dataset
# INPUT: Donor shards from Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Summary CSVs and raw distributions in Datalink/Glial_Characteristics/Metrics/
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(BiocParallel)
  library(Matrix)
})

# ==============================================================================
# 1. LOAD CONFIGURATION
# ==============================================================================

# Load config - skip if SCRIPT_DIR already set (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR")) {
  script_path <- tryCatch({ dirname(sys.frame(1)$ofile) }, error = function(e) {
    file.path(PROJECT_ROOT, "R", "glial_characteristics")
  })
  source(file.path(script_path, "config.R"))
}
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### STEP 01: EXTRACT METRICS ###\n")

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORIES
# ==============================================================================

create_output_dirs(c(OUTPUT_BASE, OUTPUT_METRICS, OUTPUT_PLOTS, OUTPUT_TESTS, OUTPUT_SUMMARY))

# ==============================================================================
# 3. PROCESS EACH DATASET × CELL TYPE
# ==============================================================================

process_dataset_celltype <- function(dataset, cell_type) {
  
  message(sprintf("\n=== Processing: %s / %s ===", dataset, cell_type))
  
  # Load manifest
  manifest_obj <- load_manifest(cell_type, dataset, SHARD_BASE, DATASETS)
  donor_ids <- get_donor_ids(manifest_obj)
  shard_dir <- manifest_obj$shard_dir
  
  message(sprintf("  Found %d donors in manifest", length(donor_ids)))
  
  # Extract cell-level metrics from all shards
  message(sprintf("  Loading shards and extracting metrics (parallel, %d cores)...", N_CORES))
  
  metrics_result <- extract_all_cell_metrics(donor_ids, shard_dir, n_cores = N_CORES)
  
  cell_metrics <- metrics_result$cell_metrics
  donor_summaries <- metrics_result$donor_summaries
  
  n_cells_total <- nrow(cell_metrics)
  n_donors_loaded <- nrow(donor_summaries)
  
  message(sprintf("  Loaded %s cells from %d donors", 
                  format(n_cells_total, big.mark = ","), n_donors_loaded))
  
  # Add dataset and cell type labels
  cell_metrics$dataset <- dataset
  cell_metrics$cell_type <- cell_type
  donor_summaries$dataset <- dataset
  donor_summaries$cell_type <- cell_type
  
  # ============================================================================
  # 3.1 Compute Summary Statistics
  # ============================================================================
  
  # UMIs per cell
  umi_stats <- data.frame(
    dataset = dataset,
    cell_type = cell_type,
    metric = "umis_per_cell",
    n = n_cells_total,
    mean = mean(cell_metrics$umis),
    median = median(cell_metrics$umis),
    sd = sd(cell_metrics$umis),
    min = min(cell_metrics$umis),
    max = max(cell_metrics$umis),
    q25 = quantile(cell_metrics$umis, 0.25),
    q75 = quantile(cell_metrics$umis, 0.75),
    stringsAsFactors = FALSE
  )
  
  # Genes per cell
  genes_stats <- data.frame(
    dataset = dataset,
    cell_type = cell_type,
    metric = "genes_per_cell",
    n = n_cells_total,
    mean = mean(cell_metrics$genes),
    median = median(cell_metrics$genes),
    sd = sd(cell_metrics$genes),
    min = min(cell_metrics$genes),
    max = max(cell_metrics$genes),
    q25 = quantile(cell_metrics$genes, 0.25),
    q75 = quantile(cell_metrics$genes, 0.75),
    stringsAsFactors = FALSE
  )
  
  # Cells per donor
  cell_counts <- donor_summaries$n_cells
  cells_stats <- data.frame(
    dataset = dataset,
    cell_type = cell_type,
    metric = "cells_per_donor",
    n = length(cell_counts),
    mean = mean(cell_counts),
    median = median(cell_counts),
    sd = sd(cell_counts),
    min = min(cell_counts),
    max = max(cell_counts),
    q25 = quantile(cell_counts, 0.25),
    q75 = quantile(cell_counts, 0.75),
    stringsAsFactors = FALSE
  )
  
  combined_stats <- bind_rows(umi_stats, genes_stats, cells_stats)
  
  message(sprintf("  UMIs/cell:   median = %s, mean = %s", 
                  format(round(umi_stats$median), big.mark = ","),
                  format(round(umi_stats$mean), big.mark = ",")))
  message(sprintf("  Genes/cell:  median = %d, mean = %.1f", 
                  round(genes_stats$median), genes_stats$mean))
  message(sprintf("  Cells/donor: median = %d, mean = %.1f", 
                  round(cells_stats$median), cells_stats$mean))
  
  # ============================================================================
  # 3.2 Save Outputs
  # ============================================================================
  
  # Summary stats CSV
  stats_file <- file.path(OUTPUT_METRICS, 
                          sprintf("Summary_Stats_%s_%s.csv", dataset, cell_type))
  write.csv(combined_stats, stats_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(stats_file)))
  
  # Donor-level summaries
  donor_file <- file.path(OUTPUT_METRICS, 
                          sprintf("Donor_Summaries_%s_%s.csv", dataset, cell_type))
  write.csv(donor_summaries, donor_file, row.names = FALSE)
  message(sprintf("  Saved: %s", basename(donor_file)))
  
  # Cell-level distributions (subsampled for memory)
  cell_metrics_sub <- subsample_df(cell_metrics, MAX_CELLS_PER_CT, seed = SUBSAMPLE_SEED)
  cell_file <- file.path(OUTPUT_METRICS, 
                         sprintf("Cell_Metrics_%s_%s.qs", dataset, cell_type))
  qs::qsave(cell_metrics_sub, cell_file)
  message(sprintf("  Saved: %s (%s cells subsampled to %s)", 
                  basename(cell_file), 
                  format(n_cells_total, big.mark = ","),
                  format(nrow(cell_metrics_sub), big.mark = ",")))
  
  message(sprintf("  ✓ Completed %s / %s", dataset, cell_type))
  
  return(list(
    stats = combined_stats,
    donor_summaries = donor_summaries
  ))
}

# ==============================================================================
# 4. RUN FOR ALL DATASETS × CELL TYPES
# ==============================================================================

all_results <- list()

for (dataset in DATASET_ORDER) {
  message(sprintf("\n--- Dataset: %s (%s) ---", 
                  dataset, DATASETS[[dataset]]$label))
  
  for (cell_type in CELL_TYPES) {
    key <- paste(dataset, cell_type, sep = "_")
    
    result <- tryCatch({
      process_dataset_celltype(dataset, cell_type)
    }, error = function(e) {
      warning(sprintf("Failed for %s / %s: %s", dataset, cell_type, e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      all_results[[key]] <- result
    }
  }
}

# ==============================================================================
# 5. COMBINE ALL RESULTS
# ==============================================================================

message("\n=== Combining Results ===")

# Combine all summary stats
all_stats <- bind_rows(lapply(all_results, `[[`, "stats"))
all_stats_file <- file.path(OUTPUT_SUMMARY, "All_Summary_Stats.csv")
write.csv(all_stats, all_stats_file, row.names = FALSE)
message(sprintf("Saved combined stats: %s", basename(all_stats_file)))

# Combine all donor summaries
all_donors <- bind_rows(lapply(all_results, `[[`, "donor_summaries"))
all_donors_file <- file.path(OUTPUT_SUMMARY, "All_Donor_Summaries.csv")
write.csv(all_donors, all_donors_file, row.names = FALSE)
message(sprintf("Saved combined donor summaries: %s", basename(all_donors_file)))

# ==============================================================================
# 6. COMPUTE INTER-CELL-TYPE RATIOS
# ==============================================================================

message("\n=== Computing Inter-Cell-Type Ratios ===")

compute_ratios_for_dataset <- function(stats_df, dataset_name) {
  
  ds_stats <- stats_df %>% filter(dataset == dataset_name)
  
  ratio_results <- list()
  
  for (metric_name in unique(ds_stats$metric)) {
    metric_stats <- ds_stats %>% filter(metric == metric_name)
    
    # Get medians for each cell type
    medians <- setNames(metric_stats$median, metric_stats$cell_type)
    
    # Reference is Micro
    ref_median <- medians[REFERENCE_CELL_TYPE]
    
    # Compute ratios
    for (ct in setdiff(names(medians), REFERENCE_CELL_TYPE)) {
      ratio_results[[paste(metric_name, ct, sep = "_")]] <- data.frame(
        dataset = dataset_name,
        metric = metric_name,
        numerator_ct = ct,
        denominator_ct = REFERENCE_CELL_TYPE,
        numerator_median = medians[ct],
        denominator_median = ref_median,
        ratio = medians[ct] / ref_median,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(ratio_results)
}

all_ratios <- bind_rows(lapply(DATASET_ORDER, function(ds) {
  compute_ratios_for_dataset(all_stats, ds)
}))

ratios_file <- file.path(OUTPUT_SUMMARY, "Inter_CellType_Ratios.csv")
write.csv(all_ratios, ratios_file, row.names = FALSE)
message(sprintf("Saved ratios: %s", basename(ratios_file)))

# Print ratio summary
message("\n=== Ratio Summary ===")
for (dataset in DATASET_ORDER) {
  ds_ratios <- all_ratios %>% filter(dataset == !!dataset)
  message(sprintf("\n%s:", dataset))
  for (i in seq_len(nrow(ds_ratios))) {
    r <- ds_ratios[i, ]
    message(sprintf("  %s: %s:%s = %.2f×", 
                    METRIC_SHORT[r$metric], r$numerator_ct, r$denominator_ct, r$ratio))
  }
}

# ==============================================================================
# 7. QUICK COMPARISON TABLE
# ==============================================================================

message("\n=== Cross-Dataset Ratio Comparison ===")

# Pivot to compare ratios across datasets
ratio_comparison <- all_ratios %>%
  mutate(ratio_label = sprintf("%s:%s", numerator_ct, denominator_ct)) %>%
  select(metric, ratio_label, dataset, ratio) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = ratio) %>%
  mutate(
    fold_diff = abs(Roussos / DeJager - 1),
    direction = ifelse(Roussos > DeJager, "Roussos higher", "DeJager higher")
  )

comparison_file <- file.path(OUTPUT_SUMMARY, "Ratio_Comparison.csv")
write.csv(ratio_comparison, comparison_file, row.names = FALSE)
message(sprintf("Saved comparison: %s", basename(comparison_file)))

# Print comparison
message("\nRatio Comparison (Roussos vs DeJager):")
print(ratio_comparison, row.names = FALSE)

message("\n### STEP 01 COMPLETE ###")
