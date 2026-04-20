# ==============================================================================
# STEP 06: Verify Paper Claims Against Source Data
# ==============================================================================
# Extracts all UNVERIFIED values from the attribution checklist into a CSV.
# ==============================================================================

cat("=== Step 06: Verify Paper Claims ===\n")

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
})

DATALINK <- file.path(PROJECT_ROOT, "Datalink")
OUTPUT_DIR <- file.path(DATALINK, "Paper_Source_Data")

results <- list()

# ==============================================================================
# 1. TRUTH SET SIZES (Items #36-38, #76, #169-170)
# ==============================================================================
cat("\n--- 1. Truth set sizes per run ---\n")

for (ct in c("Astro", "Micro", "Oligo", "Astro_Shallow", "Astro_Shallow_Fewer")) {
  fitting_file <- file.path(DATALINK, "Power_Curve_Fitting", "Fitting_Data", ct,
                            "unified_fitting_data_per_iteration.qs")
  if (!file.exists(fitting_file)) {
    cat(sprintf("  SKIP: %s\n", ct))
    next
  }

  cat(sprintf("  Loading %s...\n", ct))
  fd <- qread(fitting_file)
  fd <- as.data.table(fd)

  # Unique truth genes per run
  truth_per_run <- fd[, .(n_truth = uniqueN(gene_id)), by = .(run_id)]
  cat(sprintf("  %s: mean=%.0f, min=%d, max=%d\n",
              ct, mean(truth_per_run$n_truth), min(truth_per_run$n_truth), max(truth_per_run$n_truth)))
  results[[paste0("truth_", ct, "_mean")]] <- round(mean(truth_per_run$n_truth), 1)
  results[[paste0("truth_", ct, "_min")]] <- min(truth_per_run$n_truth)
  results[[paste0("truth_", ct, "_max")]] <- max(truth_per_run$n_truth)

  # LFC bin distribution
  gene_lfc <- fd[, .(abs_lfc = abs_gt_lfc[1]), by = .(run_id, gene_id)]
  lfc_bins <- gene_lfc[, .(
    tiny = sum(abs_lfc >= 0.2 & abs_lfc < 1.0),
    small = sum(abs_lfc >= 1.0 & abs_lfc < 2.0),
    medium = sum(abs_lfc >= 2.0 & abs_lfc < 3.0),
    large = sum(abs_lfc >= 3.0),
    total = .N
  ), by = .(run_id)]
  lfc_bins[, pct_tiny := tiny / total * 100]

  cat(sprintf("  %s bins (mean): tiny=%.0f (%.1f%%), small=%.0f, med=%.0f, large=%.0f, total=%.0f\n",
              ct, mean(lfc_bins$tiny), mean(lfc_bins$pct_tiny),
              mean(lfc_bins$small), mean(lfc_bins$medium),
              mean(lfc_bins$large), mean(lfc_bins$total)))

  results[[paste0("bins_", ct, "_pct_tiny")]] <- round(mean(lfc_bins$pct_tiny), 1)
  results[[paste0("bins_", ct, "_small")]] <- round(mean(lfc_bins$small), 1)
  results[[paste0("bins_", ct, "_medium")]] <- round(mean(lfc_bins$medium), 1)
  results[[paste0("bins_", ct, "_large")]] <- round(mean(lfc_bins$large), 1)

  # Genes at |LFC| > 1
  lfc_gt1 <- gene_lfc[, .(n = sum(abs_lfc >= 1.0)), by = run_id]
  cat(sprintf("  %s |LFC|>1: mean=%.1f, range=%d-%d\n",
              ct, mean(lfc_gt1$n), min(lfc_gt1$n), max(lfc_gt1$n)))
  results[[paste0("n_lfc_gt1_", ct)]] <- round(mean(lfc_gt1$n), 1)

  rm(fd, gene_lfc, lfc_bins, lfc_gt1)
  gc(verbose = FALSE)
}

# ==============================================================================
# 2. MEDIAN GENE COVARIATES (Items #68-69)
# ==============================================================================
cat("\n--- 2. Median gene covariates ---\n")

for (ct in c("Astro", "Micro", "Oligo")) {
  fitting_file <- file.path(DATALINK, "Power_Curve_Fitting", "Fitting_Data", ct,
                            "unified_fitting_data_per_iteration.qs")
  fd <- qread(fitting_file)
  fd <- as.data.table(fd)

  # Per-gene characteristics (one row per gene per run)
  gene_chars <- fd[, .(
    abs_lfc = abs_gt_lfc[1],
    basemean = gt_basemean[1],
    dispersion = gt_dispersion[1]
  ), by = .(run_id, gene_id)]

  med_lfc <- median(gene_chars$abs_lfc)
  med_bm <- median(gene_chars$basemean)
  med_disp <- median(gene_chars$dispersion)

  cat(sprintf("  %s: median |LFC|=%.3f, baseMean=%.1f, disp=%.3f\n",
              ct, med_lfc, med_bm, med_disp))

  results[[paste0("median_lfc_", ct)]] <- round(med_lfc, 3)
  results[[paste0("median_basemean_", ct)]] <- round(med_bm, 1)
  results[[paste0("median_disp_", ct)]] <- round(med_disp, 3)

  rm(fd, gene_chars)
  gc(verbose = FALSE)
}

# ==============================================================================
# 3. RECALL AT SPECIFIC GRID POINTS
# ==============================================================================
cat("\n--- 3. Recall at specific grid points ---\n")

for (ct in c("Oligo", "Micro")) {
  fitting_file <- file.path(DATALINK, "Power_Curve_Fitting", "Fitting_Data", ct,
                            "unified_fitting_data_per_iteration.qs")
  cat(sprintf("  Loading %s...\n", ct))
  fd <- qread(fitting_file)
  fd <- as.data.table(fd)

  # Assign LFC bins
  fd[, lfc_bin := fifelse(abs_gt_lfc >= 3.0, "Large",
                  fifelse(abs_gt_lfc >= 2.0, "Medium",
                  fifelse(abs_gt_lfc >= 1.0, "Small", "Tiny")))]

  # Mean recall per grid point per LFC bin (across runs and iterations)
  recall <- fd[, .(recall = mean(detected), n_obs = .N, n_genes = uniqueN(gene_id)),
               by = .(n_donors, n_cells, lfc_bin)]

  if (ct == "Oligo") {
    # Grid points the paper claims values for:
    claims <- list(
      c(15, 750, "Large", "oligo_15d_750c_large"),
      c(50, 250, "Large", "oligo_50d_250c_large"),
      c(100, 250, "Large", "oligo_100d_250c_large"),
      c(75, 500, "Medium", "oligo_75d_500c_medium"),
      c(50, 100, "Large", "oligo_50d_100c_large"),
      c(100, 100, "Large", "oligo_100d_100c_large"),
      c(10, 1500, "Large", "oligo_10d_1500c_large"),
      c(50, 250, "Large", "oligo_50d_250c_large2")
    )
    for (cl in claims) {
      r <- recall[n_donors == as.numeric(cl[1]) & n_cells == as.numeric(cl[2]) & lfc_bin == cl[3]]
      if (nrow(r) > 0) {
        cat(sprintf("  %s %s %sd/%sc: recall=%.4f\n", ct, cl[3], cl[1], cl[2], r$recall))
        results[[paste0("recall_", cl[4])]] <- round(r$recall, 4)
      } else {
        cat(sprintf("  %s %s %sd/%sc: NOT FOUND in grid\n", ct, cl[3], cl[1], cl[2]))
      }
    }

    # Max recall for small-effect genes
    max_small <- recall[lfc_bin == "Small", max(recall)]
    cat(sprintf("  Oligo Small max recall: %.4f\n", max_small))
    results[["oligo_max_recall_small"]] <- round(max_small, 4)

    # Show top 10 Small recall values to understand the max
    cat("  Oligo Small recall - top 10 configs:\n")
    small_top <- recall[lfc_bin == "Small"][order(-recall)][1:10]
    for (i in seq_len(nrow(small_top))) {
      cat(sprintf("    %dd/%sc: %.4f (n_obs=%d, n_genes=%d)\n",
                  small_top$n_donors[i], as.character(small_top$n_cells[i]),
                  small_top$recall[i], small_top$n_obs[i], small_top$n_genes[i]))
    }

    # Gene counts per LFC bin (at |LFC| > 1)
    gene_lfc <- fd[, .(abs_lfc = abs_gt_lfc[1]), by = .(run_id, gene_id)]
    bin_counts <- gene_lfc[abs_lfc >= 1.0, .(
      small = sum(abs_lfc >= 1.0 & abs_lfc < 2.0),
      medium = sum(abs_lfc >= 2.0 & abs_lfc < 3.0),
      large = sum(abs_lfc >= 3.0)
    ), by = run_id]
    cat(sprintf("  Oligo |LFC|>1 bins (mean): small=%.1f, med=%.1f, large=%.1f\n",
                mean(bin_counts$small), mean(bin_counts$medium), mean(bin_counts$large)))
    results[["oligo_lfc_bins_small"]] <- round(mean(bin_counts$small), 1)
    results[["oligo_lfc_bins_medium"]] <- round(mean(bin_counts$medium), 1)
    results[["oligo_lfc_bins_large"]] <- round(mean(bin_counts$large), 1)
  }

  if (ct == "Micro") {
    # Max config observation: micro at 15 donors
    # "Max" cells = the highest n_cells value in the grid
    max_cells_val <- max(recall$n_cells, na.rm = TRUE)
    r_15_max_med <- recall[n_donors == 15 & n_cells == max_cells_val & lfc_bin == "Medium"]
    r_15_500_med <- recall[n_donors == 15 & n_cells == 500 & lfc_bin == "Medium"]

    if (nrow(r_15_max_med) > 0)
      cat(sprintf("  Micro Medium 15d/Max(%d)c: recall=%.4f\n", max_cells_val, r_15_max_med$recall))
    if (nrow(r_15_500_med) > 0)
      cat(sprintf("  Micro Medium 15d/500c: recall=%.4f\n", r_15_500_med$recall))

    results[["micro_15d_max_medium"]] <- if (nrow(r_15_max_med) > 0) round(r_15_max_med$recall, 4) else NA
    results[["micro_15d_500c_medium"]] <- if (nrow(r_15_500_med) > 0) round(r_15_500_med$recall, 4) else NA

    # Also check all cell counts at 15 donors for Medium bin
    cat("  Micro Medium recall at 15 donors, all cell counts:\n")
    r_15_med <- recall[n_donors == 15 & lfc_bin == "Medium"][order(n_cells)]
    for (i in seq_len(nrow(r_15_med))) {
      cat(sprintf("    %s cells: %.4f\n", as.character(r_15_med$n_cells[i]), r_15_med$recall[i]))
    }

    # Also get all recall for Small bin at max configs to check "never exceeded 80%"
    cat("  Micro Small recall at high configs:\n")
    r_small <- recall[lfc_bin == "Small" & n_donors >= 300][order(-recall)]
    for (i in seq_len(min(5, nrow(r_small)))) {
      cat(sprintf("    %dd/%sc: %.4f\n", r_small$n_donors[i], as.character(r_small$n_cells[i]), r_small$recall[i]))
    }
  }

  rm(fd, recall)
  gc(verbose = FALSE)
}

# ==============================================================================
# 4. INDEPENDENT FILTERING (Items #5, #32-35)
# ==============================================================================
cat("\n--- 4. Independent filtering ---\n")

for (ct in c("Astro", "Micro", "Oligo")) {
  gt_dir <- file.path(DATALINK, "Power_Curve_Fitting", "GT_DDS_Checkpoints", ct)
  gt_files <- list.files(gt_dir, pattern = "GT_DDS_Run.*\\.qs$", full.names = TRUE)
  if (length(gt_files) == 0) next

  cat(sprintf("  Loading %s GT DDS (Run 1)...\n", ct))
  suppressPackageStartupMessages(library(DESeq2))
  dds <- qread(gt_files[1])

  n_total <- nrow(dds)
  tryCatch({
    res <- results(dds)
    n_tested <- sum(!is.na(res$padj))
    n_filtered <- n_total - n_tested
    pct_filtered <- round(n_filtered / n_total * 100, 1)
    pct_testable <- round(n_tested / n_total * 100, 1)

    cat(sprintf("  %s: %d total, %d tested (%.1f%% testable), %d filtered (%.1f%%)\n",
                ct, n_total, n_tested, pct_testable, n_filtered, pct_filtered))
    results[[paste0("filter_", ct, "_pct_filtered")]] <- pct_filtered
    results[[paste0("filter_", ct, "_pct_testable")]] <- pct_testable
    results[[paste0("filter_", ct, "_n_total")]] <- n_total
    results[[paste0("filter_", ct, "_n_tested")]] <- n_tested

    # Get significant genes for cross-CT testability check
    sig <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
    results[[paste0("sig_genes_", ct)]] <- paste(sig, collapse = ";")
    cat(sprintf("  %s significant (padj<0.05, |LFC|>1): %d genes\n", ct, length(sig)))
  }, error = function(e) {
    cat(sprintf("  %s: ERROR - %s\n", ct, e$message))
  })

  rm(dds)
  gc(verbose = FALSE)
}

# ==============================================================================
# 5. PRE-FILTER DONOR COUNT
# ==============================================================================
cat("\n--- 5. Pre-filter shard counts ---\n")
for (ct in c("Astro", "Micro", "Oligo")) {
  shard_dir <- file.path(DATALINK, "Shards", ct, "Roussos")
  if (dir.exists(shard_dir)) {
    n <- length(list.files(shard_dir, pattern = "\\.qs$"))
    cat(sprintf("  %s: %d shards\n", ct, n))
    results[[paste0("prefilter_", ct)]] <- n
  }
}

# ==============================================================================
# 6. POWER PREDICTIONS WITH EXPLICIT PARAMETERS
# ==============================================================================
cat("\n--- 6. Power predictions ---\n")
coefs <- fread(file.path(DATALINK, "Paper_Source_Data", "model_coefficients.csv"))

predict_p <- function(ct_name, nd, nc, rpc, lfc, bm, disp) {
  r <- coefs[cell_type == ct_name & model_type == "per_cell_type"]
  if (nrow(r) == 0) return(NA)
  logit <- r$intercept + r$beta_log_donors * log(nd) + r$beta_log_n_cells * log(nc) +
    r$beta_log_reads_per_cell * log(rpc) + r$beta_log_lfc * log(lfc) +
    r$beta_log_basemean * log(bm) + r$beta_log_dispersion * log(disp)
  1 / (1 + exp(-logit))
}

find_d80 <- function(ct_name, nc, rpc, lfc, bm, disp) {
  for (n in 1:1000) {
    if (predict_p(ct_name, n, nc, rpc, lfc, bm, disp) >= 0.80) return(n)
  }
  NA
}

# Moderate: LFC=1.5, baseMean=50, disp=0.5, 250 cells
cat("\n  Moderate (LFC=1.5, bm=50, disp=0.5, 250 cells):\n")
for (rpc in c(5000, 5742, 5827, 6000)) {
  cat(sprintf("    RPC=%d: ", rpc))
  for (ct in c("Oligo", "Astro", "Micro")) {
    d <- find_d80(ct, 250, rpc, 1.5, 50, 0.5)
    cat(sprintf("%s=%d  ", ct, d))
  }
  cat("\n")
}

# Large: LFC=3.0, baseMean=100, disp=0.5, 250 cells
cat("\n  Large (LFC=3.0, bm=100, disp=0.5, 250 cells):\n")
for (rpc in c(5000, 5742, 5827, 6000)) {
  cat(sprintf("    RPC=%d: ", rpc))
  for (ct in c("Oligo", "Astro", "Micro")) {
    d <- find_d80(ct, 250, rpc, 3.0, 100, 0.5)
    cat(sprintf("%s=%d  ", ct, d))
  }
  cat("\n")
}

# Depth-eq astro moderate
cat("\n  Depth-eq Astro moderate (LFC=1.5, bm=50, disp=0.5, 250 cells):\n")
for (rpc in c(4426, 5000, 5742, 6000)) {
  d <- find_d80("Astro_Shallow_Fewer", 250, rpc, 1.5, 50, 0.5)
  cat(sprintf("    Astro_Shallow_Fewer at RPC=%d: %d donors\n", rpc, d))
}

# Median covariates
cat("\n  Median (LFC=0.3, bm=15, disp=0.34, 75 cells, RPC=5827):\n")
for (ct in c("Oligo", "Astro", "Micro")) {
  p <- predict_p(ct, 600, 75, 5827, 0.3, 15, 0.34)
  cat(sprintf("    %s at 600d: %.4f\n", ct, p))
}

# ==============================================================================
# 7. TOTAL CELLS
# ==============================================================================
cat("\n--- 7. Total cells ---\n")
ds <- fread(file.path(DATALINK, "Paper_Source_Data", "dataset_statistics.csv"))
total <- sum(ds$N_Cells)
cat(sprintf("  Total cells: %d (%.2fM)\n", total, total / 1e6))
donors_max <- ds[, max(N_Donors), by = Dataset]
cat(sprintf("  Total donors (max per dataset): %d\n", sum(donors_max$V1)))

# ==============================================================================
# 8. GTEx FILES
# ==============================================================================
cat("\n--- 8. GTEx/Bulk files ---\n")
bulk_dir <- file.path(DATALINK, "Bulk_Analysis")
if (dir.exists(bulk_dir)) {
  all_files <- list.files(bulk_dir, recursive = TRUE)
  cat(sprintf("  %d files found:\n", length(all_files)))
  for (f in all_files) cat(sprintf("    %s\n", f))

  # Read key results
  results_dir <- file.path(bulk_dir, "Bulk_Results")
  if (dir.exists(results_dir)) {
    for (f in list.files(results_dir, pattern = "\\.csv$")) {
      d <- fread(file.path(results_dir, f), nrows = 2)
      cat(sprintf("  %s: %d cols (%s)\n", f, ncol(d), paste(names(d)[1:min(5,ncol(d))], collapse=", ")))
    }
  }
}

# ==============================================================================
# SAVE
# ==============================================================================
cat("\n=== Saving ===\n")
results_df <- data.frame(
  claim = names(results),
  value = as.character(unlist(results)),
  stringsAsFactors = FALSE
)
# Remove very long values (sig gene lists)
results_df$value <- ifelse(nchar(results_df$value) > 200,
                           paste0(substr(results_df$value, 1, 100), "..."),
                           results_df$value)
write.csv(results_df, file.path(OUTPUT_DIR, "paper_claim_verification.csv"), row.names = FALSE)
cat(sprintf("Saved %d results\n", nrow(results_df)))
cat("\n=== Step 06 Complete ===\n")
