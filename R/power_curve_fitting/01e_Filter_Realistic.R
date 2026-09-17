# ==============================================================================
# SCRIPT: Filter aggregated fitting data to the realistic design sub-grid
# PURPOSE: The full empirical grid trains on n_cells in {5,10,25,50,...} per
#          donor. Rows with n_cells < 100 collapse to binomial noise and are
#          not representative of any real experimental design. The paper's
#          empirical recall/precision surfaces still show the full grid (that
#          collapse is a finding). The FITTED heuristic power model, however,
#          should be trained only on the sub-grid a design planner might
#          actually target -- n_cells >= 100 per donor.
#
#          This script reads the aggregated table produced by 01b, applies the
#          n_cells >= 100 filter, and writes a parallel realistic .qs that the
#          reweighting sweep (02_Fit_XGBoost.R with XGB_FITTING_DATA=...) reads
#          in place of aggregated_fitting_data_all.qs.
#
#          Named 01e_ to preserve the existing XGBoost step numbering after
#          cross-method scoring moved to
#          R/cartesian_downsampling/10_CrossScore_Methods.R.
#
# INPUT:  Fitting_Data/aggregated_fitting_data_all.qs        (from 01b)
# OUTPUT: Fitting_Data/aggregated_fitting_data_realistic.qs  (filtered)
#         Fitting_Data/filter_summary.csv                    (row & bin diagnostics)
#         Fitting_Data/filter_dropped_by_group.csv           (dropped-row breakdown)
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/power_curve_fitting/config.R",
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
  if (!config_loaded) stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 01e: Filter aggregated fitting data to realistic sub-grid (n_cells >= 100)")

# Threshold (env-overridable for exploration; default locked at 100 per design doc).
MIN_N_CELLS <- as.integer(Sys.getenv("XGB_MIN_N_CELLS", "100"))
log_info("Filter threshold: n_cells >= %d", MIN_N_CELLS)

# ==============================================================================
# 2. LOAD AGGREGATED DATA
# ==============================================================================

fit_dir <- file.path(OUTPUT_BASE, "Fitting_Data")
in_path  <- file.path(fit_dir, "aggregated_fitting_data_all.qs")
if (!file.exists(in_path)) {
  stop(sprintf("Aggregated data not found: %s. Run 01b_Aggregate_Fitting_Data.R first.", in_path))
}

dat_all <- qs::qread(in_path)
log_info("Loaded aggregated table: %d rows, %d genes",
         nrow(dat_all), dplyr::n_distinct(dat_all$gene_id))

if (!"n_cells" %in% colnames(dat_all)) {
  stop("Aggregated table has no n_cells column -- cannot apply n_cells filter.")
}

# ==============================================================================
# 3. APPLY FILTER
# ==============================================================================

n_before <- nrow(dat_all)
dat <- dat_all %>% filter(is.finite(n_cells), n_cells >= MIN_N_CELLS)
n_after  <- nrow(dat)
n_drop   <- n_before - n_after

log_info("Filter kept %d of %d rows (%.1f%%); dropped %d (%.1f%%)",
         n_after, n_before, 100 * n_after / n_before,
         n_drop,  100 * n_drop  / n_before)

# ==============================================================================
# 4. DIAGNOSTICS
# ==============================================================================

# Observed detection-rate bin distribution before / after (see design doc §1).
bin_dist <- function(d, tag) {
  b <- cut(d$detection_rate, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  tab <- as.data.frame(table(b), stringsAsFactors = FALSE)
  names(tab) <- c("obs_bin", "n")
  tab$pct <- 100 * tab$n / sum(tab$n)
  tab$stage <- tag
  tab[, c("stage", "obs_bin", "n", "pct")]
}
bin_summary <- bind_rows(bin_dist(dat_all, "before_filter"),
                         bin_dist(dat,     "after_filter"))

# Top-line row-count summary.
row_summary <- data.frame(
  stage = c("before_filter", "after_filter", "dropped"),
  n_rows = c(n_before, n_after, n_drop),
  pct    = c(100, 100 * n_after / n_before, 100 * n_drop / n_before)
)

# Dropped rows broken down by (cell_type, de_method, n_cells).
group_cols <- intersect(c("cell_type", "de_method", "n_cells"), colnames(dat_all))
if (length(group_cols) > 0) {
  dropped <- dat_all %>%
    filter(!is.finite(n_cells) | n_cells < MIN_N_CELLS) %>%
    count(across(all_of(group_cols)), name = "n_dropped") %>%
    arrange(desc(n_dropped))
} else {
  dropped <- data.frame(n_dropped = n_drop)
}

# ==============================================================================
# 5. WRITE OUTPUTS
# ==============================================================================

out_qs <- file.path(fit_dir, "aggregated_fitting_data_realistic.qs")
qs::qsave(dat, out_qs)
log_info("Wrote filtered .qs -> %s", out_qs)

filter_summary_path <- file.path(fit_dir, "filter_summary.csv")
write.csv(bind_rows(
  data.frame(kind = "rows", stage = row_summary$stage,
             obs_bin = NA_character_, n = row_summary$n_rows,
             pct = row_summary$pct),
  data.frame(kind = "obs_bin", stage = bin_summary$stage,
             obs_bin = bin_summary$obs_bin, n = bin_summary$n,
             pct = bin_summary$pct)
), filter_summary_path, row.names = FALSE)
log_info("Wrote diagnostics -> %s", filter_summary_path)

dropped_path <- file.path(fit_dir, "filter_dropped_by_group.csv")
write.csv(dropped, dropped_path, row.names = FALSE)
log_info("Wrote dropped-row breakdown -> %s", dropped_path)

log_section("STEP 01e COMPLETE")
log_info("Realistic sub-grid: %d rows, %d genes", nrow(dat), dplyr::n_distinct(dat$gene_id))
