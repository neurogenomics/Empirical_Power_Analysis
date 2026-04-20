# Verify truth gene counts for paper verification
# Extracts per-run truth gene counts from fitting data

library(qs)

base_dir <- file.path(PROJECT_ROOT, "Datalink", "Power_Curve_Fitting")

cat("=== Truth Gene Counts (|LFC| > 0.2) ===\n\n")

for (ct in c("Oligo", "Astro", "Micro", "Astro_Shallow", "Astro_Shallow_Fewer")) {
  fitting_file <- file.path(base_dir, "Fitting_Data", ct, "unified_fitting_data_per_iteration.qs")

  if (!file.exists(fitting_file)) {
    cat(sprintf("SKIP: %s (file not found)\n", ct))
    next
  }

  cat(sprintf("--- %s ---\n", ct))
  dat <- qread(fitting_file)

  # Get unique truth genes per run (truth genes are those in the fitting data)
  # Each gene appears multiple times (once per iteration per grid point)
  # The truth set is the same across iterations within a run
  # Get unique genes per run
  runs <- unique(dat$run_id)

  cat(sprintf("Total observations: %s\n", format(nrow(dat), big.mark=",")))
  cat(sprintf("Columns: %s\n", paste(names(dat), collapse=", ")))

  # Get truth genes per run (unique gene_id values)
  for (r in sort(runs)[1:min(2, length(runs))]) {
    run_dat <- dat[dat$run_id == r, ]
    genes <- unique(run_dat$gene_id)
    cat(sprintf("  Run %d: %d unique truth genes\n", r, length(genes)))

    # Get the GT LFC for each gene (from the first occurrence)
    gene_info <- run_dat[!duplicated(run_dat$gene_id), c("gene_id", "abs_lfc")]
    if ("abs_gt_lfc" %in% names(run_dat)) {
      gene_info <- run_dat[!duplicated(run_dat$gene_id), c("gene_id", "abs_gt_lfc")]
      names(gene_info)[2] <- "abs_lfc"
    }

    # Bin by LFC
    tiny <- sum(gene_info$abs_lfc >= 0.2 & gene_info$abs_lfc < 1.0)
    small <- sum(gene_info$abs_lfc >= 1.0 & gene_info$abs_lfc < 2.0)
    medium <- sum(gene_info$abs_lfc >= 2.0 & gene_info$abs_lfc < 3.0)
    large <- sum(gene_info$abs_lfc >= 3.0)
    cat(sprintf("    Tiny (0.2-1): %d, Small (1-2): %d, Medium (2-3): %d, Large (3+): %d = Total: %d\n",
                tiny, small, medium, large, tiny + small + medium + large))
  }

  # Average across all runs
  gene_counts <- sapply(runs, function(r) {
    length(unique(dat$gene_id[dat$run_id == r]))
  })
  cat(sprintf("  Mean truth genes across %d runs: %.1f (range: %d-%d)\n\n",
              length(runs), mean(gene_counts), min(gene_counts), max(gene_counts)))

  rm(dat)
  gc(verbose = FALSE)
}

cat("\n=== Done ===\n")
