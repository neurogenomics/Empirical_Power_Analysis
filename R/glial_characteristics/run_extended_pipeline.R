# ==============================================================================
# SCRIPT: Extended Glial Characteristics Pipeline - All Datasets
# PURPOSE: Extract metrics and generate supplementary table for all 5 datasets
# INPUT: Donor shards from all datasets
# OUTPUT: Comprehensive metrics CSV and Supplementary Table 1
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(BiocParallel)
  library(Matrix)
  library(knitr)
})

# ==============================================================================
# 1. LOAD EXTENDED CONFIGURATION
# ==============================================================================

SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "glial_characteristics")
source(file.path(SCRIPT_DIR, "config_extended.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

message("\n### EXTENDED GLIAL CHARACTERISTICS PIPELINE ###\n")
message("Datasets: ", paste(names(DATASETS), collapse = ", "))
message("Cell Types: ", paste(CELL_TYPES, collapse = ", "))

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORIES
# ==============================================================================

create_output_dirs <- function(dirs) {
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      message("[SETUP] Created: ", dir)
    }
  }
}

create_output_dirs(c(OUTPUT_BASE, OUTPUT_METRICS, OUTPUT_PLOTS, OUTPUT_TESTS, OUTPUT_SUMMARY, OUTPUT_TABLES))

# ==============================================================================
# 3. LOAD MANIFEST AND EXTRACT METRICS
# ==============================================================================

load_manifest_extended <- function(cell_type, dataset_name, shard_base, datasets_config) {
  config <- datasets_config[[dataset_name]]
  suffix <- config$shard_dir_suffix
  
  shard_dir <- file.path(shard_base, cell_type, suffix)
  manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
  
  if (!file.exists(manifest_path)) {
    message("[WARNING] Manifest not found: ", manifest_path)
    return(NULL)
  }
  
  manifest <- qread(manifest_path)
  return(list(
    manifest = manifest,
    shard_dir = shard_dir,
    config = config
  ))
}

extract_metrics_from_shard <- function(shard_path, dataset_config) {
  tryCatch({
    sce <- qread(shard_path)
    counts <- assay(sce, "counts")
    
    # UMIs per cell
    umis_per_cell <- colSums(counts)
    
    # Genes per cell (non-zero counts)
    genes_per_cell <- colSums(counts > 0)
    
    # Cells (already known, but verify)
    n_cells <- ncol(sce)
    
    list(
      umis_per_cell = umis_per_cell,
      genes_per_cell = genes_per_cell,
      n_cells = n_cells,
      donor_id = gsub("\\.qs$", "", basename(shard_path))
    )
  }, error = function(e) {
    message("[ERROR] Failed to process ", shard_path, ": ", e$message)
    return(NULL)
  })
}

# ==============================================================================
# 4. PROCESS ALL DATASET × CELL TYPE COMBINATIONS
# ==============================================================================

all_results <- list()

for (dataset_name in names(DATASETS)) {
  for (cell_type in CELL_TYPES) {
    
    message(sprintf("\n=== Processing: %s / %s ===", dataset_name, cell_type))
    
    # Load manifest
    manifest_obj <- load_manifest_extended(cell_type, dataset_name, SHARD_BASE, DATASETS)
    if (is.null(manifest_obj)) {
      message("[SKIP] No manifest found for ", dataset_name, "/", cell_type)
      next
    }
    
    manifest <- manifest_obj$manifest
    shard_dir <- manifest_obj$shard_dir
    config <- manifest_obj$config
    
    # Get donor IDs
    donor_ids <- manifest$donor_id
    message(sprintf("[INFO] Found %d donors", length(donor_ids)))
    
    # Check if already computed
    summary_file <- file.path(OUTPUT_METRICS, sprintf("Summary_Stats_%s_%s.csv", dataset_name, cell_type))
    if (file.exists(summary_file) && !FORCE_RECALC) {
      message("[SKIP] Summary exists for ", dataset_name, "/", cell_type)
      existing <- read.csv(summary_file)
      all_results[[sprintf("%s_%s", dataset_name, cell_type)]] <- existing
      next
    }
    
    # Process shards
    all_umis <- c()
    all_genes <- c()
    donor_summaries <- list()
    
    for (i in seq_along(donor_ids)) {
      if (i %% 100 == 0) message(sprintf("  Progress: %d/%d donors", i, length(donor_ids)))
      
      donor_id <- donor_ids[i]
      shard_file <- file.path(shard_dir, sprintf("donor_%s.qs", donor_id))
      
      if (!file.exists(shard_file)) {
        # Try alternative naming
        shard_file <- file.path(shard_dir, sprintf("%s.qs", donor_id))
      }
      
      if (!file.exists(shard_file)) {
        # Try without donor_ prefix
        alt_files <- list.files(shard_dir, pattern = sprintf("%s.*\\.qs$", donor_id), full.names = TRUE)
        if (length(alt_files) > 0) {
          shard_file <- alt_files[1]
        }
      }
      
      if (!file.exists(shard_file)) {
        next
      }
      
      result <- extract_metrics_from_shard(shard_file, config)
      if (is.null(result)) next
      
      all_umis <- c(all_umis, result$umis_per_cell)
      all_genes <- c(all_genes, result$genes_per_cell)
      
      donor_summaries[[length(donor_summaries) + 1]] <- data.frame(
        donor_id = result$donor_id,
        n_cells = result$n_cells,
        median_umis = median(result$umis_per_cell),
        mean_umis = mean(result$umis_per_cell),
        median_genes = median(result$genes_per_cell),
        mean_genes = mean(result$genes_per_cell),
        stringsAsFactors = FALSE
      )
    }
    
    # Compute summary statistics
    donor_summary_df <- do.call(rbind, donor_summaries)
    
    summary_stats <- data.frame(
      dataset = dataset_name,
      cell_type = cell_type,
      metric = c("umis_per_cell", "genes_per_cell", "cells_per_donor"),
      n = c(length(all_umis), length(all_genes), nrow(donor_summary_df)),
      mean = c(mean(all_umis), mean(all_genes), mean(donor_summary_df$n_cells)),
      median = c(median(all_umis), median(all_genes), median(donor_summary_df$n_cells)),
      sd = c(sd(all_umis), sd(all_genes), sd(donor_summary_df$n_cells)),
      min = c(min(all_umis), min(all_genes), min(donor_summary_df$n_cells)),
      max = c(max(all_umis), max(all_genes), max(donor_summary_df$n_cells)),
      q25 = c(quantile(all_umis, 0.25), quantile(all_genes, 0.25), quantile(donor_summary_df$n_cells, 0.25)),
      q75 = c(quantile(all_umis, 0.75), quantile(all_genes, 0.75), quantile(donor_summary_df$n_cells, 0.75)),
      stringsAsFactors = FALSE
    )
    
    # Save results
    write.csv(summary_stats, summary_file, row.names = FALSE)
    write.csv(donor_summary_df, file.path(OUTPUT_METRICS, sprintf("Donor_Summaries_%s_%s.csv", dataset_name, cell_type)), row.names = FALSE)
    
    all_results[[sprintf("%s_%s", dataset_name, cell_type)]] <- summary_stats
    
    message(sprintf("[DONE] %s/%s: %d cells, %d donors", dataset_name, cell_type, 
                  length(all_umis), nrow(donor_summary_df)))
  }
}

# ==============================================================================
# 5. COMBINE ALL RESULTS INTO MASTER TABLE
# ==============================================================================

message("\n=== Combining Results ===")

# Combine all summary stats
all_summary <- do.call(rbind, all_results)
write.csv(all_summary, file.path(OUTPUT_SUMMARY, "All_Summary_Stats_Extended.csv"), row.names = FALSE)

# ==============================================================================
# 6. GENERATE SUPPLEMENTARY TABLE 1
# ==============================================================================

message("\n=== Generating Supplementary Table 1 ===")

# Reshape for publication format
supp_table_data <- all_summary %>%
  mutate(
    # Format values with appropriate precision
    mean_formatted = ifelse(metric == "cells_per_donor",
                           sprintf("%.1f", mean),
                           sprintf("%.0f", mean)),
    median_formatted = ifelse(metric == "cells_per_donor",
                             sprintf("%.1f", median),
                             sprintf("%.0f", median)),
    sd_formatted = ifelse(metric == "cells_per_donor",
                         sprintf("%.1f", sd),
                         sprintf("%.0f", sd)),
    range_formatted = sprintf("%.0f-%.0f", min, max),
    # Full stats with SD
    stats_full = sprintf("%s ± %s (median: %s)", mean_formatted, sd_formatted, median_formatted)
  )

# Create wide format table for publication
supp_table_wide <- supp_table_data %>%
  select(dataset, cell_type, metric, n, mean, median, sd, min, max) %>%
  pivot_wider(
    id_cols = c(dataset, cell_type),
    names_from = metric,
    values_from = c(n, mean, median, sd, min, max),
    names_glue = "{metric}_{.value}"
  )

# Add dataset metadata
dataset_metadata <- data.frame(
  dataset = names(DATASETS),
  full_name = sapply(DATASETS, function(d) d$label),
  brain_region = sapply(DATASETS, function(d) d$brain_region),
  stringsAsFactors = FALSE
)

supp_table_with_meta <- supp_table_wide %>%
  left_join(dataset_metadata, by = "dataset")

# Reorder columns for clarity
supp_table_final <- supp_table_with_meta %>%
  select(
    dataset, full_name, brain_region, cell_type,
    # Cells per donor
    cells_per_donor_n, cells_per_donor_mean, cells_per_donor_median, 
    cells_per_donor_sd, cells_per_donor_min, cells_per_donor_max,
    # UMIs per cell  
    umis_per_cell_n, umis_per_cell_mean, umis_per_cell_median,
    umis_per_cell_sd, umis_per_cell_min, umis_per_cell_max,
    # Genes per cell
    genes_per_cell_n, genes_per_cell_mean, genes_per_cell_median,
    genes_per_cell_sd, genes_per_cell_min, genes_per_cell_max
  )

# Save as CSV
write.csv(supp_table_final, file.path(OUTPUT_TABLES, "Supplementary_Table_1_Dataset_Characteristics.csv"), row.names = FALSE)

# Create formatted table for paper skeleton (summary view)
summary_for_paper <- supp_table_data %>%
  select(dataset, cell_type, metric, n, mean, median, sd) %>%
  mutate(
    display = sprintf("%.1f ± %.1f (%.0f)", mean, sd, median)
  ) %>%
  select(-mean, -median, -sd) %>%
  pivot_wider(
    id_cols = c(dataset, cell_type),
    names_from = metric,
    values_from = display
  )

write.csv(summary_for_paper, file.path(OUTPUT_TABLES, "Summary_For_Paper_Methods.csv"), row.names = FALSE)

# ==============================================================================
# 7. PRINT SUMMARY
# ==============================================================================

message("\n=== FINAL SUMMARY ===")
message("\nDataset Characteristics Summary:")
message("================================")

for (metric in c("cells_per_donor", "umis_per_cell", "genes_per_cell")) {
  message(sprintf("\n%s:", gsub("_", " ", metric)))
  metric_data <- supp_table_data %>% filter(metric == !!metric)
  
  for (ct in CELL_TYPES) {
    ct_data <- metric_data %>% filter(cell_type == ct)
    if (nrow(ct_data) > 0) {
      message(sprintf("  %s:", ct))
      for (ds in unique(ct_data$dataset)) {
        ds_row <- ct_data %>% filter(dataset == ds)
        if (nrow(ds_row) > 0) {
          message(sprintf("    %s: %s", ds, ds_row$stats_full[1]))
        }
      }
    }
  }
}

message("\n[COMPLETE] Results saved to:")
message(sprintf("  - %s", OUTPUT_TABLES))
message(sprintf("  - %s", OUTPUT_SUMMARY))
message(sprintf("  - %s", OUTPUT_METRICS))

# ==============================================================================
# 8. GENERATE KNOWLEDGE GRAPH ENTRY (for diary)

message("\n[COMPLETE] Extended Glial Characteristics Pipeline finished!")
