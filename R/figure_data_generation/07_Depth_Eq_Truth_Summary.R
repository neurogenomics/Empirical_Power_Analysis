# ==============================================================================
# STEP 07: Depth-Equalisation Truth Set Summary by LFC Bin
# ==============================================================================
# For each astrocyte model (native, depth-eq, depth+cell-eq), counts truth
# genes per LFC bin, averaged across robustness runs.
#
# Inputs:
#   - Datalink/Power_Curve_Fitting/Fitting_Data/{Model}/unified_fitting_data_per_iteration.qs
#     (produced by power_curve_fitting/01_Build_Fitting_Dataset.R)
#     Models: Astro, Astro_Shallow, Astro_Shallow_Fewer
#
# Outputs:
#   - Datalink/Paper_Source_Data/depth_eq_truth_summary.csv      (mean counts)
#   - Datalink/Paper_Source_Data/depth_eq_truth_summary_sd.csv   (SD across runs)
#
# Format: rows = LFC bins + Total, columns = pretty-printed model names.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(qs)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 07: DEPTH-EQ TRUTH SET SUMMARY BY LFC BIN ###\n")

# ==============================================================================
# Configuration
# ==============================================================================

MODELS <- c("Astro", "Astro_Shallow", "Astro_Shallow_Fewer")

PRETTY_NAMES <- c(
  "Astro"                = "Astrocytes",
  "Astro_Shallow"        = "Astro (Depth-Eq.)",
  "Astro_Shallow_Fewer"  = "Astro (Depth+Cell-Eq.)"
)

LFC_BINS <- list(
  "Tiny (0.2-1.0)"   = c(0.2, 1.0),
  "Small (1.0-2.0)"  = c(1.0, 2.0),
  "Medium (2.0-3.0)"  = c(2.0, 3.0),
  "Large (3.0+)"      = c(3.0, Inf)
)

# Try EPA path first, fall back to parent tree's v2 path (precomputed data)
FITTING_BASE <- file.path(PROJECT_ROOT, "Datalink", "Power_Curve_Fitting",
                           "Fitting_Data")
if (!file.exists(file.path(FITTING_BASE, "Astro", "unified_fitting_data_per_iteration.qs"))) {
  parent_v2 <- file.path(dirname(PROJECT_ROOT), "Datalink",
                          "Power_Curve_Fitting_v2", "Fitting_Data")
  if (dir.exists(parent_v2)) {
    FITTING_BASE <- parent_v2
    message(sprintf("  Using parent tree fitting data: %s", parent_v2))
  }
}

# ==============================================================================
# Helper: bin absolute LFC values
# ==============================================================================

bin_lfc <- function(abs_lfc, bins = LFC_BINS) {
  labels <- names(bins)
  result <- rep(NA_character_, length(abs_lfc))
  for (i in seq_along(bins)) {
    lo <- bins[[i]][1]
    hi <- bins[[i]][2]
    result[abs_lfc >= lo & abs_lfc < hi] <- labels[i]
  }
  # Handle upper boundary of last bin (Inf means >=)
  last <- length(bins)
  if (is.infinite(bins[[last]][2])) {
    result[abs_lfc >= bins[[last]][1]] <- labels[last]
  }
  factor(result, levels = labels)
}

# ==============================================================================
# Extract per-run truth gene counts by LFC bin for each model
# ==============================================================================

all_counts <- list()

for (model in MODELS) {
  data_path <- file.path(FITTING_BASE, model,
                          "unified_fitting_data_per_iteration.qs")

  if (!file.exists(data_path)) {
    message(sprintf("  [SKIP] Fitting data not found: %s", data_path))
    next
  }

  message(sprintf("  Loading %s ...", model))
  fit_data <- qread(data_path)

  # Extract unique truth genes per robustness run
  truth_per_run <- fit_data %>%
    distinct(run_id, gene_id, abs_gt_lfc) %>%
    mutate(lfc_bin = bin_lfc(abs_gt_lfc))

  # Count genes per bin per run
  counts <- truth_per_run %>%
    filter(!is.na(lfc_bin)) %>%
    group_by(run_id, lfc_bin, .drop = FALSE) %>%
    summarise(n = n(), .groups = "drop")

  # Add total row per run
  totals <- truth_per_run %>%
    filter(!is.na(lfc_bin)) %>%
    group_by(run_id) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(lfc_bin = "Total")

  counts <- bind_rows(counts, totals) %>%
    mutate(model = model)

  all_counts[[model]] <- counts
  message(sprintf("  %s: %d runs, mean total = %.0f genes",
                  model,
                  n_distinct(counts$run_id),
                  mean(totals$n)))

  rm(fit_data)
  gc(verbose = FALSE)
}

combined <- bind_rows(all_counts)

# ==============================================================================
# Build wide tables: rows = LFC bin, columns = model (pretty name)
# ==============================================================================

bin_order <- c(names(LFC_BINS), "Total")

# Mean table
mean_wide <- combined %>%
  group_by(model, lfc_bin) %>%
  summarise(mean_n = round(mean(n), 1), .groups = "drop") %>%
  mutate(
    model = PRETTY_NAMES[model],
    lfc_bin = factor(lfc_bin, levels = bin_order)
  ) %>%
  arrange(lfc_bin) %>%
  pivot_wider(names_from = model, values_from = mean_n) %>%
  rename(LFC_Bin = lfc_bin)

# SD table
sd_wide <- combined %>%
  group_by(model, lfc_bin) %>%
  summarise(sd_n = round(sd(n), 2), .groups = "drop") %>%
  mutate(
    model = PRETTY_NAMES[model],
    lfc_bin = factor(lfc_bin, levels = bin_order)
  ) %>%
  arrange(lfc_bin) %>%
  pivot_wider(names_from = model, values_from = sd_n) %>%
  rename(LFC_Bin = lfc_bin)

# ==============================================================================
# Write
# ==============================================================================

out_mean <- file.path(OUTPUT_BASE, "depth_eq_truth_summary.csv")
out_sd   <- file.path(OUTPUT_BASE, "depth_eq_truth_summary_sd.csv")

write.csv(mean_wide, out_mean, row.names = FALSE)
write.csv(sd_wide, out_sd, row.names = FALSE)

message(sprintf("\n[DONE] Wrote mean table to %s", out_mean))
message(sprintf("[DONE] Wrote SD table to %s", out_sd))
message("\nMean truth gene counts by LFC bin:")
print(as.data.frame(mean_wide), row.names = FALSE)
