# ==============================================================================
# TOY TEST: Glial Characteristics Pipeline - Quick Validation
# ==============================================================================
# This script validates the pipeline logic with synthetic data without
# loading actual shards. Useful for testing code changes before production run.
# 
# Run with: Rscript toy/run_toy_test.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

message("\n### TOY TEST: Glial Characteristics Pipeline ###\n")

# Load config and utils
SCRIPT_DIR <- dirname(dirname(sys.frame(1)$ofile %||% 
  file.path(PROJECT_ROOT, "R", "glial_characteristics")))
source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# 1. CREATE SYNTHETIC DATA
# ==============================================================================

message("Creating synthetic data...")

set.seed(42)

# Simulate cell-level metrics with realistic Roussos-like distributions
# Based on actual Roussos dataset statistics:
#   Astro: mean UMIs ~9,456, median cells/donor ~466
#   Micro: mean UMIs ~4,426, median cells/donor ~187
#   Oligo: mean UMIs ~5,742, median cells/donor ~1,364

create_synthetic_data <- function(dataset, n_donors = 50) {
  
  cell_data_list <- list()
  donor_data_list <- list()
  
  for (ct in CELL_TYPES) {
    # Set parameters based on cell type
    if (ct == "Astro") {
      umi_mean <- 9500
      umi_sd <- 5000
      cells_median <- 450
      genes_mean <- 3500
    } else if (ct == "Micro") {
      umi_mean <- 4400
      umi_sd <- 2000
      cells_median <- 190
      genes_mean <- 2500
    } else {  # Oligo
      umi_mean <- 5700
      umi_sd <- 2500
      cells_median <- 1350
      genes_mean <- 2800
    }
    
    # Adjust for DeJager (slightly different characteristics)
    if (dataset == "DeJager") {
      umi_mean <- umi_mean * runif(1, 0.85, 1.15)  # 15% variation
      cells_median <- round(cells_median * runif(1, 0.8, 1.2))
    }
    
    for (d in 1:n_donors) {
      n_cells <- max(30, round(rnorm(1, cells_median, cells_median * 0.4)))
      
      cell_data_list[[paste(dataset, ct, d, sep = "_")]] <- data.frame(
        donor_id = sprintf("donor_%s_%s_%d", dataset, ct, d),
        dataset = dataset,
        cell_type = ct,
        umis = pmax(1000, rnorm(n_cells, umi_mean, umi_sd)),
        genes = pmax(500, rnorm(n_cells, genes_mean, 800)),
        stringsAsFactors = FALSE
      )
      
      donor_data_list[[paste(dataset, ct, d, sep = "_")]] <- data.frame(
        donor_id = sprintf("donor_%s_%s_%d", dataset, ct, d),
        dataset = dataset,
        cell_type = ct,
        n_cells = n_cells,
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    cell_metrics = bind_rows(cell_data_list),
    donor_summaries = bind_rows(donor_data_list)
  )
}

roussos_data <- create_synthetic_data("Roussos", n_donors = 50)
dejager_data <- create_synthetic_data("DeJager", n_donors = 40)

all_cell_metrics <- bind_rows(roussos_data$cell_metrics, dejager_data$cell_metrics)
all_donors <- bind_rows(roussos_data$donor_summaries, dejager_data$donor_summaries)

message(sprintf("  Created %s synthetic cells", format(nrow(all_cell_metrics), big.mark = ",")))
message(sprintf("  Created %d synthetic donors", nrow(all_donors)))

# ==============================================================================
# 2. TEST SUMMARY STATISTICS
# ==============================================================================

message("\nTesting summary statistics...")

for (dataset in DATASET_ORDER) {
  for (ct in CELL_TYPES) {
    subset_data <- all_cell_metrics %>% filter(dataset == !!dataset, cell_type == !!ct)
    
    umi_stats <- compute_summary_stats(subset_data$umis, "umis")
    message(sprintf("  %s / %s: UMI median = %.0f", dataset, ct, umi_stats$median))
  }
}

# ==============================================================================
# 3. TEST RATIO COMPUTATION
# ==============================================================================

message("\nTesting ratio computation...")

for (dataset in DATASET_ORDER) {
  for (metric in c("umis", "genes")) {
    astro_median <- all_cell_metrics %>%
      filter(dataset == !!dataset, cell_type == "Astro") %>%
      pull(!!sym(metric)) %>%
      median()
    
    micro_median <- all_cell_metrics %>%
      filter(dataset == !!dataset, cell_type == "Micro") %>%
      pull(!!sym(metric)) %>%
      median()
    
    ratio <- astro_median / micro_median
    message(sprintf("  %s: Astro:Micro %s ratio = %.2f", dataset, metric, ratio))
  }
}

# ==============================================================================
# 4. TEST BOOTSTRAP
# ==============================================================================

message("\nTesting bootstrap functions...")

astro_umis <- all_cell_metrics %>%
  filter(dataset == "Roussos", cell_type == "Astro") %>%
  pull(umis)

micro_umis <- all_cell_metrics %>%
  filter(dataset == "Roussos", cell_type == "Micro") %>%
  pull(umis)

boot_ratios <- bootstrap_median_ratio(astro_umis, micro_umis, n_boot = 1000, seed = 42)
ci <- bootstrap_ci(boot_ratios)

message(sprintf("  Astro:Micro UMI ratio: %.2f [95%% CI: %.2f, %.2f]",
                ci["median"], ci["lower"], ci["upper"]))

# ==============================================================================
# 5. TEST PLOTTING
# ==============================================================================

message("\nTesting plot generation...")

# Create test output directory
toy_output <- file.path(OUTPUT_BASE, "Toy_Test")
dir.create(toy_output, recursive = TRUE, showWarnings = FALSE)

all_cell_metrics$cell_type_pretty <- prettify_label(all_cell_metrics$cell_type)
all_cell_metrics$cell_type_pretty <- factor(all_cell_metrics$cell_type_pretty,
                                            levels = prettify_label(CELL_TYPES))

p_test <- ggplot(all_cell_metrics, aes(x = umis, fill = cell_type_pretty)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  facet_wrap(~dataset) +
  theme_publication(base_size = 14) +
  scale_fill_manual(values = celltype_colours(prettified = TRUE)) +
  labs(title = "TOY TEST: UMIs per Cell", x = "UMIs", y = "Count", fill = "Cell Type")

test_plot_path <- file.path(toy_output, "Toy_Test_Histogram.png")
ggsave(test_plot_path, p_test, width = 10, height = 5, dpi = 150)
message(sprintf("  Saved test plot: %s", test_plot_path))

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

message("\n### TOY TEST COMPLETE ###")
message("All components validated successfully.")
message(sprintf("Test outputs: %s", toy_output))
