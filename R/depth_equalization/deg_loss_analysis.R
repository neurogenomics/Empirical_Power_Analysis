#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  deg_loss_analysis.R
# PURPOSE: Simple, decoupled DEG-loss result for the depth-equalisation section.
#          Runs the 3 DE packages (DESeq2, edgeR, limma) on the ground-truth
#          split (50% of donors, GT_SPLIT, averaged over 10 robustness runs) for
#          each of four astrocyte configurations, and reports how many ground-
#          truth DEGs are lost relative to the native (full) configuration.
#
#          This REPLACES the old derivation (07_Depth_Eq_Truth_Summary.R), which
#          scraped DEG counts out of the heavy power-curve-fitting dataset. Here
#          the number is computed directly and is fully independent of the
#          XGBoost / power-curve machinery.
#
# DESIGN (proper 2x2 of depth x cell-count equalisation):
#   native          Roussos               full cells,   full depth
#   depth_only      Roussos_DepthEq       full cells,   reduced depth (~Micro)
#   cells_only      Roussos_CellsOnly     reduced cells (~Micro), full depth   [BUILT HERE if missing]
#   depth_and_cells Roussos_ShallowFewer  reduced cells (~Micro), reduced depth
#
# GROUND TRUTH per run: set.seed(SEED_BASE + run_id); take GT_SPLIT of each sex;
#   aggregate pseudobulk over ALL cells; run all 3 DE methods. Identical to
#   R/cartesian_downsampling/01_Generate_Ground_Truth.R.
#
# DEG definition: padj < FDR_THRESHOLD & |log2FoldChange| > LFC_CUT, computed at
#   BOTH thresholds (0.2 and 1) from the same cached DE tables. The 0.2 vs 1
#   contrast is what exposes limma-voom's instability at the permissive cutoff.
#
# OUTPUT (Datalink/Depth_Eq_DEG_Loss/):
#   deg_loss_table.csv    ONE display-ready supplementary table. Rows = 4 configs;
#                         6 value cols = {deseq2,edger,limma} x {|LFC|>0.2, |LFC|>1};
#                         each cell = "mean DEG count (% of native count)".
#   deg_loss_counts.csv   tidy numeric backing: (variant, method, lfc_cut) ->
#                         mean_n_deg, sd_n_deg, pct_of_native.
#   deg_loss_per_run.csv  per (variant, method, run) DEG counts at both thresholds.
#   GT_Cache/<variant>/Run_<r>/<method>.qs   cached per-run DE tables (server-only; not synced)
#
# USAGE:  Rscript R/depth_equalization/deg_loss_analysis.R [n_cores]
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(parallel)
})

# ------------------------------------------------------------------------------
# 0. ARGS + ROOT + SHARED MACHINERY
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
N_CORES <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 12L
# Optional arg[2] caps the number of robustness runs (smoke test). Unset = all.
MAX_RUNS <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else NA_integer_

# Resolve project root the same way the pipelines do.
if (!nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  if (requireNamespace("here", quietly = TRUE)) {
    Sys.setenv(POWER_ANALYSIS_ROOT = here::here())
  }
}

# Reuse the cartesian pipeline's config (constants: SEED_BASE, GT_SPLIT,
# MIN_GENE_COUNT, FDR_THRESHOLD, DE_METHODS, SEX/AGE cols, etc.) and utils
# (aggregate_from_shards, run_all_de_methods, load_donor_metadata).
CART_DIR <- file.path(Sys.getenv("POWER_ANALYSIS_ROOT"), "R", "cartesian_downsampling")
source(file.path(CART_DIR, "config.R"))
source(file.path(CART_DIR, "utils.R"))
# Depth-eq utils supply the cell-count matching used to build cells_only.
DE_UTILS_DIR <- file.path(Sys.getenv("POWER_ANALYSIS_ROOT"), "R", "depth_equalization")
source(file.path(DE_UTILS_DIR, "utils.R"))

set.seed(1)  # only affects any incidental RNG; per-run seeds are set explicitly

CELL_TYPE <- "Astro"
METHODS   <- IMPLEMENTED_DE_METHODS          # deseq2, edger, limma
RUN_IDS   <- 1:N_ROBUSTNESS_RUNS             # 1..10
if (!is.na(MAX_RUNS)) RUN_IDS <- 1:min(N_ROBUSTNESS_RUNS, MAX_RUNS)
LFC_CUTS  <- c(0.2, 1)                       # both truth-set thresholds (0.2 exposes limma instability)

OUT_DIR   <- project_path("Datalink", "Depth_Eq_DEG_Loss")
CACHE_DIR <- file.path(OUT_DIR, "GT_Cache")
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

log_section(sprintf("DEG-LOSS ANALYSIS (depth equalisation) | cores=%d | runs=%d",
                    N_CORES, length(RUN_IDS)))

# ------------------------------------------------------------------------------
# 1. BUILD cells_only VARIANT (Roussos_CellsOnly) IF MISSING
#    = native cells randomly subsampled to Micro's mean cell count, FULL depth.
#    Uses the SAME scaling factor and per-donor seeds (seed_base = 9998) as the
#    ShallowFewer build, so cells_only is the exact "cells removed, depth kept"
#    counterpart to depth_and_cells.
# ------------------------------------------------------------------------------
CELLS_ONLY_DATASET <- "Roussos_CellsOnly"
cells_only_dir <- file.path(SHARD_BASE, CELL_TYPE, CELLS_ONLY_DATASET)
cells_only_sentinel <- file.path(cells_only_dir, ".cells_only_complete")

native_dir      <- file.path(SHARD_BASE, CELL_TYPE, "Roussos")
native_manifest <- file.path(native_dir, "00_MASTER_MANIFEST.qs")
micro_manifest  <- file.path(SHARD_BASE, "Micro", "Roussos", "00_MASTER_MANIFEST.qs")

if (!file.exists(cells_only_sentinel)) {
  log_section("STEP 1: Build cells_only variant (Roussos_CellsOnly)")
  if (!file.exists(native_manifest) || !file.exists(micro_manifest)) {
    log_warning("Cannot build cells_only: missing native or micro manifest. Skipping this variant.")
  } else {
    native_summary <- load_cell_count_summary(native_manifest)
    micro_summary  <- load_cell_count_summary(micro_manifest)
    matching_stats <- compute_scaling_factor(native_summary, micro_summary)
    log_info("cells_only scaling factor: %.4f (native mean %.1f cells -> micro mean %.1f cells)",
             matching_stats$scaling_factor, matching_stats$source_mean, matching_stats$target_mean)

    match_results <- match_cell_counts(
      source_shard_dir    = native_dir,
      target_shard_dir    = cells_only_dir,
      matching_stats      = matching_stats,
      min_cells_threshold = 30,      # same threshold as ShallowFewer (FEWER_EXTENSION)
      seed_base           = 9998     # same seed base as ShallowFewer
    )

    src_cfg <- list(
      cell_type = CELL_TYPE, min_age = NA, min_cells = NA,
      source_dataset = "Roussos", target_mean_rpc = NA, downsample_fraction = 1
    )
    # generate_matched_manifest references FEWER_EXTENSION$target_dataset for a
    # config field; depth_equalization/config.R (sourced via its utils? no) —
    # ensure it exists.
    if (!exists("FEWER_EXTENSION")) FEWER_EXTENSION <- list(target_dataset = "Roussos")
    cells_manifest <- generate_matched_manifest(
      shard_dir = cells_only_dir, source_config = src_cfg,
      matching_stats = matching_stats, match_results = match_results
    )
    qs::qsave(cells_manifest, file.path(cells_only_dir, "00_MASTER_MANIFEST.qs"))
    write.csv(match_results, file.path(cells_only_dir, "cell_count_matching_summary.csv"),
              row.names = FALSE)
    writeLines(c(
      sprintf("completed_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      sprintf("source: %s", native_dir),
      sprintf("scaling_factor: %.6f", matching_stats$scaling_factor),
      sprintf("n_donors_retained: %d", sum(match_results$success)),
      sprintf("n_donors_dropped: %d",  sum(!match_results$success))
    ), cells_only_sentinel)
    log_info("Built cells_only: %d donors retained, %d dropped",
             sum(match_results$success), sum(!match_results$success))
  }
} else {
  log_info("cells_only variant already built (sentinel present): %s", cells_only_dir)
}

# ------------------------------------------------------------------------------
# 2. VARIANT DEFINITIONS
# ------------------------------------------------------------------------------
VARIANTS <- c(
  native          = "Roussos",
  depth_only      = "Roussos_DepthEq",
  cells_only      = CELLS_ONLY_DATASET,
  depth_and_cells = "Roussos_ShallowFewer"
)
VARIANTS <- VARIANTS[vapply(VARIANTS, function(ds) {
  file.exists(file.path(SHARD_BASE, CELL_TYPE, ds, "00_MASTER_MANIFEST.qs"))
}, logical(1))]
log_info("Variants available: %s", paste(names(VARIANTS), collapse = ", "))

# ------------------------------------------------------------------------------
# 3. DONOR SEX MAP (scan native shards ONCE, parallel) + per-variant donor pools
#    Sex is not in the manifest, so we must read colData once. Native donor sex
#    applies to every variant (equalisation never changes a donor's sex). Each
#    variant's donor universe/order is taken from its own manifest.
# ------------------------------------------------------------------------------
log_section("STEP 2: Build donor sex map (scan native shards once)")
native_mf <- qs::qread(native_manifest)$manifest
native_files <- native_mf$file_path
native_files[!file.exists(native_files)] <-
  file.path(native_dir, basename(native_files[!file.exists(native_files)]))

sex_list <- mclapply(native_files, function(f) {
  if (!file.exists(f)) return(NULL)
  tryCatch({
    cd <- SummarizedExperiment::colData(qs::qread(f))
    data.frame(donor_id = as.character(cd[[DONOR_ID_COL]][1]),
               sex      = as.character(cd[[SEX_COL]][1]),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}, mc.cores = N_CORES)
sex_map <- do.call(rbind, sex_list)
sex_map <- sex_map[!is.na(sex_map$donor_id), ]
sex_vec <- setNames(sex_map$sex, sex_map$donor_id)
log_info("Sex map built for %d donors (%d M, %d F)",
         nrow(sex_map), sum(sex_map$sex == MALE_VAL), sum(sex_map$sex == FEMALE_VAL))

# Per-variant ordered donor pools (donor order = variant manifest order, matching
# production semantics of load_donor_metadata for that variant).
variant_pools <- lapply(VARIANTS, function(ds) {
  mf <- qs::qread(file.path(SHARD_BASE, CELL_TYPE, ds, "00_MASTER_MANIFEST.qs"))$manifest
  ids <- as.character(mf$donor_id)
  sx  <- sex_vec[ids]
  keep <- !is.na(sx)
  ids <- ids[keep]; sx <- sx[keep]
  list(
    m_pool = ids[sx == MALE_VAL],
    f_pool = ids[sx == FEMALE_VAL],
    shard_dir = file.path(SHARD_BASE, CELL_TYPE, ds)
  )
})
names(variant_pools) <- names(VARIANTS)
for (v in names(variant_pools)) {
  log_info("  %-16s %d donors (%d M, %d F)", v,
           length(variant_pools[[v]]$m_pool) + length(variant_pools[[v]]$f_pool),
           length(variant_pools[[v]]$m_pool), length(variant_pools[[v]]$f_pool))
}

# ------------------------------------------------------------------------------
# 4. RUN GROUND-TRUTH DE FOR EVERY (variant, run) IN PARALLEL
#    Each job: reproduce the GT split, aggregate pseudobulk, run 3 DE methods,
#    cache per-method result tables, return DEG tables in memory.
# ------------------------------------------------------------------------------
log_section("STEP 3: Ground-truth DE across variants x runs x methods")

jobs <- expand.grid(variant = names(variant_pools), run = RUN_IDS,
                    stringsAsFactors = FALSE)
log_info("Dispatching %d (variant x run) jobs across %d cores", nrow(jobs), N_CORES)

read_de_table <- function(path) {
  r  <- qs::qread(path)
  df <- as.data.frame(r)
  gene <- rownames(df)
  bm <- if ("baseMean_raw" %in% names(df)) df$baseMean_raw else df$baseMean
  data.frame(
    gene    = gene,
    baseMean = as.numeric(bm),
    lfc     = as.numeric(df$log2FoldChange),
    padj    = as.numeric(df$padj),
    stringsAsFactors = FALSE
  )
}

run_one_job <- function(j) {
  variant <- jobs$variant[j]
  run_id  <- jobs$run[j]
  pool    <- variant_pools[[variant]]
  run_seed <- SEED_BASE + run_id

  run_cache <- file.path(CACHE_DIR, variant, paste0("Run_", run_id))
  dir.create(run_cache, recursive = TRUE, showWarnings = FALSE)
  res_paths <- setNames(file.path(run_cache, paste0(METHODS, ".qs")), METHODS)
  log_paths <- setNames(file.path(run_cache, paste0(METHODS, "_log.qs")), METHODS)

  # --- reproduce GT donor split (identical to 01_Generate_Ground_Truth.R) ---
  set.seed(run_seed)
  n_gt_m <- floor(length(pool$m_pool) * GT_SPLIT)
  n_gt_f <- floor(length(pool$f_pool) * GT_SPLIT)
  gt_ids <- c(sample(pool$m_pool, n_gt_m), sample(pool$f_pool, n_gt_f))

  # --- aggregate + DE (skip if cached) ---
  need <- !all(file.exists(res_paths))
  if (need) {
    agg <- aggregate_from_shards(
      donor_ids = gt_ids, shard_dir = pool$shard_dir,
      n_cells_to_sample = NA, mode = "cells", seed = run_seed, verbose = FALSE
    )
    if (is.null(agg)) return(NULL)
    run_all_de_methods(
      agg_result = agg, methods = METHODS,
      save_path_results = res_paths, save_path_log = log_paths,
      min_gene_count = MIN_GENE_COUNT
    )
  }

  # --- read back DE tables (thresholding is applied downstream, per LFC cut) ---
  out <- lapply(METHODS, function(m) {
    if (!file.exists(res_paths[[m]])) return(NULL)
    read_de_table(res_paths[[m]])
  })
  names(out) <- METHODS
  list(variant = variant, run = run_id, n_gt = length(gt_ids), tables = out)
}

results <- mclapply(seq_len(nrow(jobs)), run_one_job, mc.cores = N_CORES)
ok <- !vapply(results, is.null, logical(1))
if (any(!ok)) log_warning("%d/%d jobs returned NULL", sum(!ok), length(results))
results <- results[ok]

# ------------------------------------------------------------------------------
# 5. AGGREGATE: net DEG counts at BOTH LFC thresholds (0.2 and 1)
# ------------------------------------------------------------------------------
log_section("STEP 4: Aggregate net DEG counts at |LFC|>0.2 and |LFC|>1")

count_deg <- function(tab, lfc_cut) {
  sum(!is.na(tab$padj) & tab$padj < FDR_THRESHOLD & abs(tab$lfc) > lfc_cut)
}

# --- per (variant, method, run, cut) DEG counts ---
per_run_rows <- list()
for (r in results) {
  for (m in METHODS) {
    tab <- r$tables[[m]]
    if (is.null(tab)) next
    for (cut in LFC_CUTS) {
      per_run_rows[[length(per_run_rows) + 1]] <- data.frame(
        variant = r$variant, method = m, run = r$run,
        n_gt_donors = r$n_gt, lfc_cut = cut,
        n_deg = count_deg(tab, cut), stringsAsFactors = FALSE
      )
    }
  }
}
per_run <- do.call(rbind, per_run_rows)
per_run <- per_run[order(per_run$variant, per_run$method, per_run$lfc_cut, per_run$run), ]
write.csv(per_run, file.path(OUT_DIR, "deg_loss_per_run.csv"), row.names = FALSE)
log_info("Wrote %s (%d rows)", file.path(OUT_DIR, "deg_loss_per_run.csv"), nrow(per_run))

# --- mean over runs per (variant, method, cut) + net count as % of native ---
variant_order <- names(variant_pools)
counts_rows <- list()
for (m in METHODS) {
  for (cut in LFC_CUTS) {
    nat_sub  <- per_run[per_run$variant == "native" & per_run$method == m & per_run$lfc_cut == cut, ]
    nat_mean <- if (nrow(nat_sub) > 0) mean(nat_sub$n_deg) else NA_real_
    for (variant in variant_order) {
      sub <- per_run[per_run$variant == variant & per_run$method == m & per_run$lfc_cut == cut, ]
      if (nrow(sub) == 0) next
      mn <- mean(sub$n_deg)
      counts_rows[[length(counts_rows) + 1]] <- data.frame(
        variant = variant, method = m, lfc_cut = cut, n_runs = nrow(sub),
        mean_n_deg = mn, sd_n_deg = stats::sd(sub$n_deg),
        pct_of_native = if (!is.na(nat_mean) && nat_mean > 0) 100 * mn / nat_mean else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
counts <- do.call(rbind, counts_rows)
write.csv(counts, file.path(OUT_DIR, "deg_loss_counts.csv"), row.names = FALSE)
log_info("Wrote %s (%d rows)", file.path(OUT_DIR, "deg_loss_counts.csv"), nrow(counts))

# ------------------------------------------------------------------------------
# 6. BUILD ONE DISPLAY-READY SUPPLEMENTARY TABLE
#    Rows = 4 configs; 6 value cols = {method} x {|LFC|>0.2, |LFC|>1};
#    each cell = "mean DEG count (% of native count)". Native row = 100%.
# ------------------------------------------------------------------------------
log_section("STEP 5: Build display-ready supplementary table")

variant_labels <- c(
  native          = "Native (full depth & cells)",
  depth_only      = "Depth-equalised",
  cells_only      = "Cell-count-equalised",
  depth_and_cells = "Depth + cell equalised"
)
method_labels <- c(deseq2 = "DESeq2", edger = "edgeR", limma = "limma")

fmt_cell <- function(variant, m, cut) {
  row <- counts[counts$variant == variant & counts$method == m & counts$lfc_cut == cut, ]
  if (nrow(row) == 0 || is.na(row$mean_n_deg)) return(NA_character_)
  sprintf("%d (%.1f%%)", as.integer(round(row$mean_n_deg)), row$pct_of_native)
}

col_specs <- expand.grid(method = METHODS, cut = LFC_CUTS, stringsAsFactors = FALSE)
col_specs <- col_specs[order(match(col_specs$method, METHODS), col_specs$cut), ]

table_df <- data.frame(Configuration = unname(variant_labels[variant_order]),
                       stringsAsFactors = FALSE, check.names = FALSE)
for (i in seq_len(nrow(col_specs))) {
  m <- col_specs$method[i]; cut <- col_specs$cut[i]
  col_name <- sprintf("%s (|LFC| > %s)", method_labels[[m]], format(cut))
  table_df[[col_name]] <- vapply(variant_order, function(v) fmt_cell(v, m, cut), character(1))
}
write.csv(table_df, file.path(OUT_DIR, "deg_loss_table.csv"), row.names = FALSE)
log_info("Wrote %s (%d rows x %d cols)", file.path(OUT_DIR, "deg_loss_table.csv"),
         nrow(table_df), ncol(table_df))

# ------------------------------------------------------------------------------
# 7. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
log_section("HEADLINE: mean DEG count (pct of native) by config x method x cut")
for (variant in variant_order) {
  for (m in METHODS) {
    parts <- vapply(LFC_CUTS, function(cut) {
      row <- counts[counts$variant == variant & counts$method == m & counts$lfc_cut == cut, ]
      if (nrow(row) == 0) return("NA")
      sprintf("|LFC|>%s: %.1f (%.1f%%)", format(cut), row$mean_n_deg, row$pct_of_native)
    }, character(1))
    log_info("  %-16s %-7s  %s", variant, m, paste(parts, collapse = "  "))
  }
}
log_section("DEG-LOSS ANALYSIS COMPLETE")
