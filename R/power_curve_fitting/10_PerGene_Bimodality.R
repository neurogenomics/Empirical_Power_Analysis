# ==============================================================================
# SCRIPT: Per-gene detection bimodality
# PURPOSE: Test the "detection is nearly binary per gene" hypothesis. The
#          global calibration table shows only ~7.8% of ROWS land in the
#          (0.1, 0.9) mid-band. Does this generalise per-gene? Concretely:
#            - Are individual genes mostly "always detected" or "never
#              detected", with only a narrow strip of design points in the
#              transition band?
#            - Or is the low global mid-band fraction driven by many genes
#              each with a wide transition band that averages out?
#
#          Approach:
#            1. Sample N genes (default 1000) stratified by cell type.
#            2. For each gene, aggregate across run_id (mean detection_rate
#               per design point = (cell_type, de_method, mode, n_donors, y_value)).
#            3. Compute per-gene distribution metrics.
#            4. Aggregate cross-gene: % of genes with fraction-in-mid-band
#               below thresholds, gene-archetype breakdown, bimodality score.
#            5. Visualise: ECDF of contested fraction, mean-vs-contested
#               scatter (archetype map), per-gene strip small-multiples for
#               a random 30-gene subsample.
#
#          Runs on iclneuro (no HPC needed). ~2-5 min on the realistic
#          aggregated data.
#
# INPUT:  Fitting_Data/aggregated_fitting_data_realistic.qs   (from 01e)
#         (or aggregated_fitting_data.qs if XGB_BIMODAL_DATA=full)
# OUTPUT: Models/xgboost/realistic/per_gene_bimodality/
#           per_gene_summary.csv         one row per sampled gene
#           threshold_summary.csv        cross-gene % below each threshold
#           archetype_counts.csv         (always/never/mixed) archetype counts
#           ecdf_contested_fraction.png  ECDF across genes
#           mean_vs_contested.png        2D archetype map
#           per_gene_strips.png          strip plot for 30 sampled genes
#           contested_density.png        density of contested fraction
#           README.md                    how to read + caveats
#
# ENV VARS:
#   XGB_MODEL_SUBDIR       (default "realistic") - output dir under Models/xgboost/
#   XGB_BIMODAL_DATA       (default "realistic") - "realistic" or "full"
#   XGB_BIMODAL_N_GENES    (default 1000)         - sample size
#   XGB_BIMODAL_SEED       (default 42)           - sampling seed
#   XGB_BIMODAL_MIN_POINTS (default 20)           - min design points per gene
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
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

log_section("STEP 10: Per-gene detection bimodality")

XGB_MODEL_SUBDIR       <- Sys.getenv("XGB_MODEL_SUBDIR",       "realistic")
XGB_BIMODAL_DATA       <- tolower(Sys.getenv("XGB_BIMODAL_DATA","realistic"))
XGB_BIMODAL_N_GENES    <- as.integer(Sys.getenv("XGB_BIMODAL_N_GENES",    "1000"))
XGB_BIMODAL_SEED       <- as.integer(Sys.getenv("XGB_BIMODAL_SEED",       "42"))
XGB_BIMODAL_MIN_POINTS <- as.integer(Sys.getenv("XGB_BIMODAL_MIN_POINTS", "20"))

data_file <- if (XGB_BIMODAL_DATA == "realistic") {
  "aggregated_fitting_data_realistic.qs"
} else if (XGB_BIMODAL_DATA == "full") {
  "aggregated_fitting_data.qs"
} else {
  stop("XGB_BIMODAL_DATA must be 'realistic' or 'full'")
}

fitting_data_path <- file.path(OUTPUT_BASE, "Fitting_Data", data_file)
if (!file.exists(fitting_data_path)) {
  stop("Fitting data not found: ", fitting_data_path,
       ". Run 01_Build_Fitting_Dataset.R (and 01e_Filter_Realistic.R if using realistic).")
}

out_dir <- file.path(OUTPUT_BASE, "Models", "xgboost", XGB_MODEL_SUBDIR,
                     "per_gene_bimodality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_info("Loading %s", basename(fitting_data_path))
dat <- qs::qread(fitting_data_path)
log_info("  rows        : %d", nrow(dat))
log_info("  cell types  : %s", paste(sort(unique(dat$cell_type)), collapse = ", "))
log_info("  de methods  : %s", paste(sort(unique(dat$de_method)), collapse = ", "))
log_info("  genes total : %d", dplyr::n_distinct(dat$gene_id))

# ==============================================================================
# 2. COLLAPSE ACROSS run_id TO GET ONE detection_rate PER (gene, design point)
# ==============================================================================
# design point = (cell_type, de_method, mode, n_donors, y_value)
# We AVERAGE detection_rate across run_ids. Rationale: different runs are
# random half-donor splits of the same underlying design; the true p is the
# same. Averaging reduces the Bernoulli sampling noise from 20 to up to 20*n_runs.
# We keep gene-level features from the first available run (they are constant
# per (gene, run) but the same gene may have slightly different values across
# runs since GT is estimated per-run; taking the mean smooths).

log_info("Collapsing across run_id...")

if ("mode" %in% colnames(dat)) {
  design_cols <- c("cell_type", "de_method", "mode", "n_donors", "y_value")
} else {
  design_cols <- c("cell_type", "de_method", "n_donors", "y_value")
}

per_point <- dat |>
  group_by(across(all_of(c("gene_id", design_cols)))) |>
  summarise(
    detection_rate = mean(detection_rate, na.rm = TRUE),
    n_runs         = dplyr::n(),
    abs_gt_lfc     = mean(abs_gt_lfc,     na.rm = TRUE),
    gt_basemean    = mean(gt_basemean,    na.rm = TRUE),
    gt_dispersion  = mean(gt_dispersion,  na.rm = TRUE),
    .groups = "drop"
  )

log_info("  design points (gene x setup): %d", nrow(per_point))
log_info("  unique genes                : %d", dplyr::n_distinct(per_point$gene_id))

# ==============================================================================
# 3. SAMPLE GENES (STRATIFIED BY CELL TYPE, MIN POINTS FILTER)
# ==============================================================================

set.seed(XGB_BIMODAL_SEED)

gene_totals <- per_point |>
  group_by(gene_id, cell_type) |>
  summarise(n_points = dplyr::n(), .groups = "drop") |>
  filter(n_points >= XGB_BIMODAL_MIN_POINTS)

log_info("Genes with >= %d design points per cell type: %d",
         XGB_BIMODAL_MIN_POINTS, nrow(gene_totals))

per_ct <- max(1L, floor(XGB_BIMODAL_N_GENES / dplyr::n_distinct(gene_totals$cell_type)))
sampled_keys <- gene_totals |>
  group_by(cell_type) |>
  slice_sample(n = per_ct, replace = FALSE) |>
  ungroup() |>
  select(gene_id, cell_type)

log_info("Sampled %d (gene, cell_type) keys (~%d per cell type)",
         nrow(sampled_keys), per_ct)

sampled <- per_point |> semi_join(sampled_keys, by = c("gene_id", "cell_type"))

# ==============================================================================
# 4. PER-GENE METRICS
# ==============================================================================
# For each (gene_id, cell_type) compute:
#   - n_points                        : number of design points
#   - mean_p                          : mean detection_rate
#   - sd_p                            : sd of detection_rate
#   - frac_in_(0.1,0.9), (0.2,0.8), (0.3,0.7)  : fraction of contested points
#   - frac_below_0.1, frac_above_0.9  : fraction of "definitely not / safe" points
#   - bimodality_coefficient          : sample-corrected BC = (skew^2 + 1)/(kurt + 3(n-1)^2 / ((n-2)(n-3)))
#     Higher BC => more bimodal. Uniform = 5/9 = 0.555; BC > 0.555 => bimodal.
#   - archetype                       : "always_off" (mean < 0.1) / "always_on" (mean > 0.9)
#                                       / "transition" (mean 0.1-0.9)

bimodality_coef <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 8) return(NA_real_)
  m  <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  z  <- (x - m) / s
  skew <- sum(z^3) / n
  kurt <- sum(z^4) / n - 3
  denom <- kurt + 3 * ((n - 1)^2 / ((n - 2) * (n - 3)))
  (skew^2 + 1) / denom
}

log_info("Computing per-gene metrics...")

per_gene <- sampled |>
  group_by(gene_id, cell_type) |>
  summarise(
    n_points          = dplyr::n(),
    mean_p            = mean(detection_rate),
    sd_p              = sd(detection_rate),
    mean_abs_gt_lfc   = mean(abs_gt_lfc,    na.rm = TRUE),
    mean_gt_basemean  = mean(gt_basemean,   na.rm = TRUE),
    mean_gt_dispersion= mean(gt_dispersion, na.rm = TRUE),
    frac_below_0.1    = mean(detection_rate <  0.1),
    frac_above_0.9    = mean(detection_rate >  0.9),
    frac_in_0.1_0.9   = mean(detection_rate >= 0.1 & detection_rate <= 0.9),
    frac_in_0.2_0.8   = mean(detection_rate >= 0.2 & detection_rate <= 0.8),
    frac_in_0.3_0.7   = mean(detection_rate >= 0.3 & detection_rate <= 0.7),
    bimodality_coef   = bimodality_coef(detection_rate),
    .groups = "drop"
  ) |>
  mutate(archetype = case_when(
    mean_p <  0.1                    ~ "always_off",
    mean_p >  0.9                    ~ "always_on",
    TRUE                             ~ "transition"
  ))

log_info("Per-gene table: %d rows", nrow(per_gene))

write.csv(per_gene,
          file.path(out_dir, "per_gene_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# 5. CROSS-GENE THRESHOLD SUMMARY (Chris's asked-for statistic)
# ==============================================================================
# "What proportion of them had a 0.1-0.9 probability less than 5%, 10%, 20%,
# 30% of the time?" -> for each band and threshold, % of genes whose
# fraction-in-band is BELOW the threshold.

thresholds <- c(0.05, 0.10, 0.20, 0.30, 0.50)
bands <- list(
  "0.1-0.9" = "frac_in_0.1_0.9",
  "0.2-0.8" = "frac_in_0.2_0.8",
  "0.3-0.7" = "frac_in_0.3_0.7"
)

threshold_summary <- do.call(rbind, lapply(names(bands), function(bname) {
  col <- bands[[bname]]
  do.call(rbind, lapply(thresholds, function(thr) {
    data.frame(
      band          = bname,
      threshold     = thr,
      n_genes       = nrow(per_gene),
      n_below       = sum(per_gene[[col]] <  thr, na.rm = TRUE),
      pct_below     = 100 * mean(per_gene[[col]] < thr, na.rm = TRUE),
      median_frac   = median(per_gene[[col]], na.rm = TRUE),
      mean_frac     = mean(per_gene[[col]], na.rm = TRUE)
    )
  }))
}))

write.csv(threshold_summary,
          file.path(out_dir, "threshold_summary.csv"),
          row.names = FALSE)

log_info("Threshold summary (%% of genes with fraction-in-band below X):")
print(threshold_summary)

# ==============================================================================
# 6. ARCHETYPE COUNTS
# ==============================================================================

archetype_counts <- per_gene |>
  group_by(cell_type, archetype) |>
  summarise(
    n            = dplyr::n(),
    median_frac_contested = median(frac_in_0.1_0.9, na.rm = TRUE),
    median_bimodality     = median(bimodality_coef, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

write.csv(archetype_counts,
          file.path(out_dir, "archetype_counts.csv"),
          row.names = FALSE)

log_info("Archetype counts by cell type:")
print(archetype_counts)

# ==============================================================================
# 7. PLOTS
# ==============================================================================

# ---- 7a. ECDF of contested fraction across genes -----------------------------

ecdf_data <- per_gene |>
  select(gene_id, cell_type, frac_in_0.1_0.9, frac_in_0.2_0.8, frac_in_0.3_0.7) |>
  pivot_longer(starts_with("frac_in_"),
               names_to = "band",
               values_to = "frac") |>
  mutate(band = recode(band,
                       frac_in_0.1_0.9 = "0.1-0.9",
                       frac_in_0.2_0.8 = "0.2-0.8",
                       frac_in_0.3_0.7 = "0.3-0.7"))

p_ecdf <- ggplot(ecdf_data, aes(x = frac, colour = band)) +
  stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.7) +
  scale_colour_manual(values = c("0.1-0.9" = "#440154",
                                 "0.2-0.8" = "#21918c",
                                 "0.3-0.7" = "#fde725")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Fraction of a gene's design points in the contested band",
       y = "ECDF across sampled genes",
       colour = "Band") +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.25))

ggsave(file.path(out_dir, "ecdf_contested_fraction.png"), p_ecdf,
       width = 7, height = 5, dpi = 150)

# ---- 7b. mean_p vs contested fraction (archetype map) ------------------------

p_map <- ggplot(per_gene, aes(x = mean_p, y = frac_in_0.1_0.9)) +
  geom_bin_2d(bins = 40) +
  scale_fill_viridis_c(trans = "log10", name = "Genes") +
  facet_wrap(~ cell_type) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mean detection rate (per gene, across design points)",
       y = "Fraction of design points in (0.1, 0.9)") +
  theme_minimal()

ggsave(file.path(out_dir, "mean_vs_contested.png"), p_map,
       width = 10, height = 4, dpi = 150)

# ---- 7c. Per-gene strip plot for 30 randomly sampled genes -------------------
# Colours strips by mean_p so you can eyeball "always off" (blue) vs
# "always on" (yellow) vs "transitions" (mid).

set.seed(XGB_BIMODAL_SEED + 1)
strip_keys <- per_gene |>
  slice_sample(n = min(30, nrow(per_gene))) |>
  arrange(mean_p) |>
  mutate(gene_label = sprintf("%s / %s (mean=%.2f)", gene_id, cell_type, mean_p),
         gene_label = factor(gene_label, levels = gene_label))

strip_dat <- sampled |>
  inner_join(strip_keys |> select(gene_id, cell_type, gene_label),
             by = c("gene_id", "cell_type"))

p_strip <- ggplot(strip_dat,
                  aes(x = detection_rate, y = gene_label, colour = detection_rate)) +
  geom_jitter(height = 0.25, size = 0.7, alpha = 0.5) +
  scale_colour_viridis_c(option = "plasma") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "detection_rate across a gene's design points", y = NULL,
       colour = "detection_rate") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right")

ggsave(file.path(out_dir, "per_gene_strips.png"), p_strip,
       width = 9, height = 8, dpi = 150)

# ---- 7d. Density of contested fraction, split by archetype -------------------
# Chris asked: histogram is ugly. Alternative: violin/density of the
# CONTESTED FRACTION statistic, split by archetype so you can see whether
# "always on/off" archetypes have near-zero contested fraction AND
# "transition" archetype has meaningfully non-zero (or is itself narrow).

p_density <- ggplot(per_gene, aes(x = frac_in_0.1_0.9, fill = archetype)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(always_off = "#440154",
                               transition = "#fde725",
                               always_on  = "#21918c")) +
  labs(x = "Fraction of a gene's design points in (0.1, 0.9)",
       y = "Density (across sampled genes)",
       fill = "Archetype") +
  facet_wrap(~ cell_type) +
  theme_minimal()

ggsave(file.path(out_dir, "contested_density.png"), p_density,
       width = 10, height = 4, dpi = 150)

# ==============================================================================
# 8. README
# ==============================================================================

readme <- sprintf(paste(
  "# Per-gene detection bimodality",
  "",
  "**Generated by:** `R/power_curve_fitting/10_PerGene_Bimodality.R`",
  "**Source data:** `%s`",
  "**Sample size:** %d (gene, cell_type) keys stratified across cell types",
  "**Seed:** %d",
  "**Min design points per gene:** %d",
  "",
  "## Question",
  "",
  "The global calibration table shows only ~7.8%% of ROWS land in the (0.1, 0.9) mid-band.",
  "Does that translate to a per-gene picture where detection is nearly binary --",
  "each gene is either 'always detected' or 'never detected' at the design points it appears at,",
  "with only a narrow strip of design points in the transition?",
  "",
  "## Files",
  "",
  "- `per_gene_summary.csv` -- one row per (gene, cell_type). Columns include",
  "  `mean_p` (mean detection across design points), `frac_in_0.1_0.9` (fraction",
  "  of design points in the contested band), `bimodality_coef` (Pearson's sample",
  "  bimodality coefficient; BC > 0.555 suggests bimodality), and `archetype`",
  "  (always_off / always_on / transition based on mean_p thresholds 0.1 and 0.9).",
  "- `threshold_summary.csv` -- direct answer to 'what %% of genes have contested",
  "  fraction below X'. Rows are (band, threshold), value is `pct_below`.",
  "- `archetype_counts.csv` -- breakdown of genes into archetypes by cell type.",
  "",
  "## Plots",
  "",
  "- `ecdf_contested_fraction.png` -- ECDF of contested fraction across genes,",
  "  one line per band definition. A steep ECDF near 0 means most genes have",
  "  almost no contested design points (binary detection holds).",
  "- `mean_vs_contested.png` -- 2D density heatmap. X = mean detection per gene,",
  "  Y = fraction contested. Genes concentrated in corners (0,0) or (1,0) support",
  "  binary detection; a diagonal ridge from (0.5,~1) means transition band is real.",
  "- `contested_density.png` -- density of contested fraction, split by archetype.",
  "  Distinguishes 'mass at 0 because gene is always_off' from 'mass at 0 because",
  "  transition genes have a narrow band'.",
  "- `per_gene_strips.png` -- 30 randomly sampled genes, ordered by mean detection.",
  "  Direct visual check of whether individual genes cluster their design points",
  "  near 0/1 or spread them across the mid-band.",
  "",
  "## Caveats",
  "",
  "- Averaging across `run_id` reduces Bernoulli noise but hides between-run",
  "  variation in the GT-set definition. A gene that appears in only 1 of 10 GT",
  "  runs still contributes to `per_point`, and its detection_rate at any design",
  "  point is computed from that one run's 20 iterations.",
  "- 'Transition' archetype is defined by mean_p in (0.1, 0.9). A gene classed",
  "  as transition may still be bimodal (some design points at ~0, some at ~1)",
  "  with almost NO contested design points -- exactly Chris's hypothesis. That",
  "  case is diagnosed by low `frac_in_0.1_0.9` despite `mean_p` near 0.5.",
  "  Cross-check via `mean_vs_contested.png` archetype map.",
  "- Bimodality coefficient assumes >=8 points and is unreliable at edge cases.",
  sep = "\n"),
  basename(fitting_data_path),
  XGB_BIMODAL_N_GENES,
  XGB_BIMODAL_SEED,
  XGB_BIMODAL_MIN_POINTS
)

writeLines(readme, file.path(out_dir, "README.md"))

log_info("Wrote outputs to %s", out_dir)
log_section("STEP 10 complete")
