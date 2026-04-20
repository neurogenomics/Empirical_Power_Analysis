# ==============================================================================
# SCRIPT: Bootstrap Ratio Tests - Statistical Tests for Ratio Differences
# PURPOSE: Compute bootstrap confidence intervals for inter-cell-type ratios
#          and test whether ratios differ between datasets
# INPUT: Donor summaries from Step 01
# OUTPUT: Bootstrap CIs, forest plots, and ratio difference tests
# ==============================================================================
# NOTE: All bootstrapping is performed at the DONOR level (the independent
# statistical unit). For UMIs/cell and genes/cell, we use per-donor medians
# (computed in Step 01) rather than raw cell-level values. This is both
# statistically correct and computationally efficient.
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
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

message("\n### STEP 03: BOOTSTRAP RATIO TESTS ###\n")

# ==============================================================================
# 2. LOAD DATA FROM STEP 01
# ==============================================================================

message("Loading data from Step 01...")

# Load donor summaries — contains per-donor: n_cells, median_umis, median_genes
all_donors <- read.csv(file.path(OUTPUT_SUMMARY, "All_Donor_Summaries.csv"),
                       stringsAsFactors = FALSE)
message(sprintf("  Loaded donor summaries: %d donors", nrow(all_donors)))

# Quick sanity check
for (ds in DATASET_ORDER) {
  for (ct in CELL_TYPES) {
    n <- sum(all_donors$dataset == ds & all_donors$cell_type == ct)
    message(sprintf("    %s / %s: %d donors", ds, ct, n))
  }
}

# ==============================================================================
# 3. BOOTSTRAP INTER-CELL-TYPE RATIOS (DONOR-LEVEL)
# ==============================================================================

message("\n=== Bootstrapping Inter-Cell-Type Ratios ===")
message(sprintf("  Bootstrap samples: %d, Seed: %d", N_BOOTSTRAP, BOOTSTRAP_SEED))
message("  Unit of resampling: DONOR (independent statistical unit)")

# Define ratio pairs to compute
ratio_pairs <- list(
  c("Astro", "Micro"),
  c("Oligo", "Micro"),
  c("Astro", "Oligo")
)

# Metrics to bootstrap (all at donor level)
# - median_umis: per-donor median UMIs/cell
# - median_genes: per-donor median genes/cell
# - n_cells: cells per donor
metric_configs <- list(
  list(col = "median_umis", label = "UMIs/cell", 
       desc = "Median UMIs per cell (donor-level)"),
  list(col = "median_genes", label = "Genes/cell",
       desc = "Median genes detected per cell (donor-level)"),
  list(col = "n_cells", label = "Cells/donor",
       desc = "Cells per donor")
)

# Function to bootstrap ratio of medians between two cell types within a dataset
bootstrap_ratio_for_metric <- function(donors_df, metric_col, dataset_name, 
                                       numerator_ct, denominator_ct,
                                       n_boot = N_BOOTSTRAP, seed = BOOTSTRAP_SEED) {
  
  # Get donor-level values for each cell type
  num_data <- donors_df %>% 
    filter(dataset == dataset_name, cell_type == numerator_ct) %>% 
    pull(!!sym(metric_col))
  denom_data <- donors_df %>% 
    filter(dataset == dataset_name, cell_type == denominator_ct) %>% 
    pull(!!sym(metric_col))
  
  if (length(num_data) == 0 || length(denom_data) == 0) {
    warning(sprintf("No data for %s or %s in %s", numerator_ct, denominator_ct, dataset_name))
    return(NULL)
  }
  
  # Point estimate
  point_ratio <- median(num_data) / median(denom_data)
  
  # Bootstrap (donor-level resampling — fast: ~450-1300 donors per group)
  boot_ratios <- bootstrap_median_ratio(num_data, denom_data, n_boot = n_boot, seed = seed)
  ci <- bootstrap_ci(boot_ratios, ci_level = CI_LEVEL)
  
  data.frame(
    dataset = dataset_name,
    metric = metric_col,
    numerator_ct = numerator_ct,
    denominator_ct = denominator_ct,
    ratio_label = sprintf("%s:%s", numerator_ct, denominator_ct),
    numerator_median = median(num_data),
    denominator_median = median(denom_data),
    ratio_point = point_ratio,
    ratio_ci_lower = ci["lower"],
    ratio_ci_median = ci["median"],
    ratio_ci_upper = ci["upper"],
    n_numerator = length(num_data),
    n_denominator = length(denom_data),
    stringsAsFactors = FALSE
  )
}

# Run bootstrapping for all combinations
all_bootstrap_results <- list()

for (mc in metric_configs) {
  message(sprintf("\n  Bootstrapping %s ratios...", mc$label))
  
  for (ratio_pair in ratio_pairs) {
    for (dataset in DATASET_ORDER) {
      key <- paste(mc$col, dataset, paste(ratio_pair, collapse = "_"), sep = "_")
      message(sprintf("    %s: %s:%s", dataset, ratio_pair[1], ratio_pair[2]))
      
      result <- bootstrap_ratio_for_metric(
        all_donors, mc$col, dataset, ratio_pair[1], ratio_pair[2]
      )
      if (!is.null(result)) {
        result$metric_label <- mc$label
        all_bootstrap_results[[key]] <- result
      }
    }
  }
}

# Combine all results
bootstrap_df <- bind_rows(all_bootstrap_results)

# Save bootstrap results
bootstrap_file <- file.path(OUTPUT_TESTS, "Bootstrap_Ratio_CIs.csv")
write.csv(bootstrap_df, bootstrap_file, row.names = FALSE)
message(sprintf("\n  Saved: %s", basename(bootstrap_file)))

# ==============================================================================
# 4. COMPARE RATIOS BETWEEN DATASETS
# ==============================================================================

message("\n=== Comparing Ratios Between Datasets ===")

# For each metric × ratio_pair, compare Roussos vs DeJager
ratio_diff_results <- list()

for (metric in unique(bootstrap_df$metric)) {
  for (ratio_pair in ratio_pairs) {
    ratio_label <- sprintf("%s:%s", ratio_pair[1], ratio_pair[2])
    
    roussos_row <- bootstrap_df %>%
      filter(dataset == "Roussos", metric == !!metric, ratio_label == !!ratio_label)
    dejager_row <- bootstrap_df %>%
      filter(dataset == "DeJager", metric == !!metric, ratio_label == !!ratio_label)
    
    if (nrow(roussos_row) == 1 && nrow(dejager_row) == 1) {
      
      # Check CI overlap
      ci_overlap <- !(roussos_row$ratio_ci_upper < dejager_row$ratio_ci_lower ||
                       dejager_row$ratio_ci_upper < roussos_row$ratio_ci_lower)
      
      # Compute ratio of ratios
      ratio_of_ratios <- roussos_row$ratio_point / dejager_row$ratio_point
      
      # Percent difference
      pct_diff <- 100 * (roussos_row$ratio_point / dejager_row$ratio_point - 1)
      
      ratio_diff_results[[paste(metric, ratio_label, sep = "_")]] <- data.frame(
        metric = metric,
        metric_label = roussos_row$metric_label,
        ratio_label = ratio_label,
        roussos_ratio = roussos_row$ratio_point,
        roussos_ci_lower = roussos_row$ratio_ci_lower,
        roussos_ci_upper = roussos_row$ratio_ci_upper,
        dejager_ratio = dejager_row$ratio_point,
        dejager_ci_lower = dejager_row$ratio_ci_lower,
        dejager_ci_upper = dejager_row$ratio_ci_upper,
        ratio_of_ratios = ratio_of_ratios,
        pct_diff = pct_diff,
        ci_overlap = ci_overlap,
        interpretation = ifelse(ci_overlap, "PRESERVED", "DIVERGENT"),
        stringsAsFactors = FALSE
      )
    }
  }
}

ratio_diff_df <- bind_rows(ratio_diff_results)

# Save comparison results
diff_file <- file.path(OUTPUT_TESTS, "Ratio_Difference_Tests.csv")
write.csv(ratio_diff_df, diff_file, row.names = FALSE)
message(sprintf("  Saved: %s", basename(diff_file)))

# ==============================================================================
# 5. GENERATE FOREST PLOT
# ==============================================================================

message("\n=== Generating Forest Plot ===")

# Prepare data for forest plot
forest_data <- bootstrap_df %>%
  mutate(
    dataset = factor(dataset, levels = DATASET_ORDER),
    metric_label = factor(metric_label, levels = c("UMIs/cell", "Genes/cell", "Cells/donor")),
    y_label = sprintf("%s\n%s", ratio_label, metric_label)
  )

# Create forest plot
p_forest <- ggplot(forest_data, 
                   aes(x = ratio_point, y = ratio_label, color = dataset)) +
  geom_errorbarh(aes(xmin = ratio_ci_lower, xmax = ratio_ci_upper),
                 height = 0.2, position = position_dodge(0.4), linewidth = 0.8) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~metric_label, ncol = 1, scales = "free_x") +
  theme_publication(base_size = BASE_FONT_SIZE) +
  scale_color_manual(values = DATASET_COLOURS) +
  labs(
    title = "Inter-Cell-Type Ratios with 95% Bootstrap CIs",
    subtitle = "Comparing Roussos vs De Jager (both DLPFC) | 10,000 bootstrap resamples at donor level",
    x = "Ratio (fold)",
    y = "Cell Type Comparison",
    color = "Dataset"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_TESTS, "Forest_Plot_Ratio_CIs.png"),
       p_forest, width = 10, height = 10, dpi = PLOT_DPI)
message("  Saved: Forest_Plot_Ratio_CIs.png")

# ==============================================================================
# 6. SUMMARY TABLE WITH INTERPRETATION
# ==============================================================================

message("\n=== Final Summary ===")

# Create interpretive summary
summary_interp <- ratio_diff_df %>%
  select(metric_label, ratio_label, 
         roussos_ratio, roussos_ci_lower, roussos_ci_upper,
         dejager_ratio, dejager_ci_lower, dejager_ci_upper,
         pct_diff, interpretation) %>%
  mutate(
    roussos_ci = sprintf("%.2f [%.2f, %.2f]", roussos_ratio, roussos_ci_lower, roussos_ci_upper),
    dejager_ci = sprintf("%.2f [%.2f, %.2f]", dejager_ratio, dejager_ci_lower, dejager_ci_upper),
    pct_diff_fmt = sprintf("%+.1f%%", pct_diff)
  ) %>%
  select(metric_label, ratio_label, roussos_ci, dejager_ci, pct_diff_fmt, interpretation)

summary_file <- file.path(OUTPUT_SUMMARY, "Ratio_Interpretation_Summary.csv")
write.csv(summary_interp, summary_file, row.names = FALSE)
message(sprintf("  Saved: %s", basename(summary_file)))

# Print summary to console
message("\n========================================")
message("INTER-CELL-TYPE RATIO COMPARISON SUMMARY")
message("========================================")

for (i in seq_len(nrow(ratio_diff_df))) {
  r <- ratio_diff_df[i, ]
  
  status_icon <- ifelse(r$interpretation == "PRESERVED", "✓", "✗")
  
  message(sprintf("\n%s %s (%s):",
                  status_icon, r$ratio_label, r$metric_label))
  message(sprintf("  Roussos: %.2f [%.2f, %.2f]",
                  r$roussos_ratio, r$roussos_ci_lower, r$roussos_ci_upper))
  message(sprintf("  DeJager: %.2f [%.2f, %.2f]",
                  r$dejager_ratio, r$dejager_ci_lower, r$dejager_ci_upper))
  message(sprintf("  Difference: %+.1f%% | CIs %s → %s",
                  r$pct_diff,
                  ifelse(r$ci_overlap, "OVERLAP", "DO NOT OVERLAP"),
                  r$interpretation))
}

# Overall summary
n_preserved <- sum(ratio_diff_df$interpretation == "PRESERVED")
n_divergent <- sum(ratio_diff_df$interpretation == "DIVERGENT")
n_total <- nrow(ratio_diff_df)

message("\n========================================")
message(sprintf("OVERALL: %d/%d ratios PRESERVED (%.0f%%)",
                n_preserved, n_total, 100 * n_preserved / n_total))
message(sprintf("         %d/%d ratios DIVERGENT (%.0f%%)",
                n_divergent, n_total, 100 * n_divergent / n_total))
message("========================================")

if (n_preserved == n_total) {
  message("\nCONCLUSION: All inter-cell-type ratios are preserved across datasets.")
  message("           The patterns observed in Roussos are likely BIOLOGICAL INVARIANTS")
  message("           (driven by cell biology and brain composition), not dataset-specific artifacts.")
} else if (n_preserved > n_divergent) {
  message("\nCONCLUSION: Most inter-cell-type ratios are preserved across datasets.")
  message("           The dominant patterns are likely BIOLOGICAL, with some dataset-specific variation.")
} else {
  message("\nCONCLUSION: Many inter-cell-type ratios differ between datasets.")
  message("           The patterns may be DATASET-SPECIFIC (protocol, QC, annotation differences).")
}

message("\n### STEP 03 COMPLETE ###")
