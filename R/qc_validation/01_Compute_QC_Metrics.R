# ==============================================================================
# SCRIPT: Compute QC Metrics Across All Roussos Shards
# PURPOSE: For each cell type, compute per-cell QC metrics (UMIs, genes, mito%)
#          for every donor shard and check against Fullard et al. 2025 thresholds.
#          Implements fail-early logic to stop if data is obviously raw/pre-QC.
# INPUT: Donor shard files (.qs) from Datalink/Shards/{CellType}/Roussos_Full/
# OUTPUT: Per-cell-type QC metric files in Datalink/QC_Validation/{CellType}/
# ==============================================================================

if (!exists("PROJECT_ROOT")) source("config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

message("=" |> rep(70) |> paste0(collapse = ""))
message("STEP 01: Compute QC Metrics Across All Roussos Shards")
message("=" |> rep(70) |> paste0(collapse = ""))

# ==============================================================================
# 1. Get mitochondrial gene list
# ==============================================================================

message("\n[STEP 1.1] Fetching mitochondrial gene list from Ensembl...")
mito_genes <- get_mito_genes(cache_file = MITO_GENE_CACHE)

if (is.null(mito_genes) || length(mito_genes) == 0) {
  message("[WARN] Ensembl query failed. Will use ^MT- prefix fallback for each shard.")
}

# ==============================================================================
# 2. Process each cell type
# ==============================================================================

for (cell_type in CELL_TYPES) {
  message(sprintf("\n%s", "=" |> rep(60) |> paste0(collapse = "")))
  message(sprintf("Processing cell type: %s", cell_type))
  message("=" |> rep(60) |> paste0(collapse = ""))

  shard_dir <- file.path(SHARD_BASE, cell_type, SHARD_DATASET)
  output_dir <- file.path(OUTPUT_BASE, cell_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Check for checkpoint
  checkpoint_file <- file.path(output_dir, "qc_all_cells.qs")
  donor_summary_file <- file.path(output_dir, "qc_summary_per_donor.csv")

  if (file.exists(checkpoint_file) && file.exists(donor_summary_file)) {
    message(sprintf("[CHECKPOINT] Found existing QC results for %s. Skipping.", cell_type))
    message(sprintf("  Delete %s to force recomputation.", checkpoint_file))
    next
  }

  # List shard files
  shard_files <- list_donor_shards(shard_dir)
  n_shards <- length(shard_files)
  message(sprintf("  Found %d donor shards in: %s", n_shards, shard_dir))

  if (n_shards == 0) {
    message(sprintf("[WARN] No shards found for %s. Skipping.", cell_type))
    next
  }

  # Accumulators
  all_cell_qc <- list()
  all_donor_summaries <- list()
  n_processed <- 0
  total_cells <- 0
  failed_donors <- character(0)

  for (i in seq_along(shard_files)) {
    shard_path <- shard_files[i]
    donor_id <- extract_donor_id(shard_path)

    # Progress reporting
    if (i %% 100 == 0 || i <= 10 || i == n_shards) {
      message(sprintf("  [%d/%d] Processing donor: %s", i, n_shards, donor_id))
    }

    tryCatch({
      sce <- qs::qread(shard_path)

      # Sanity check on first donor
      if (i == 1) {
        counts_mat <- assay(sce, "counts")
        sanity_ok <- sanity_check_counts(counts_mat,
          sprintf("%s - first donor (%s)", cell_type, donor_id))
        if (!sanity_ok) {
          message("[WARN] Sanity check failed on first donor. Continuing but flagging.")
        }
      }

      # Compute QC metrics
      qc_df <- compute_cell_qc_metrics(sce, mito_genes)
      qc_df <- check_thresholds(qc_df,
        min_umis = QC_MIN_UMIS, max_umis = QC_MAX_UMIS,
        min_genes = QC_MIN_GENES, max_genes = QC_MAX_GENES,
        max_mito_pct = QC_MAX_MITO_PCT
      )
      qc_df$donor_id <- donor_id
      qc_df$cell_type <- cell_type

      # Store
      all_cell_qc[[i]] <- qc_df
      all_donor_summaries[[i]] <- summarise_donor(qc_df, donor_id)

      n_processed <- n_processed + 1
      total_cells <- total_cells + nrow(qc_df)

      # Fail-early check
      if (n_processed %in% c(10, 50)) {
        combined_so_far <- do.call(rbind, all_cell_qc[!sapply(all_cell_qc, is.null)])
        fail_check <- check_fail_early(combined_so_far, n_processed,
          threshold_10 = FAIL_EARLY_10_DONORS_PCT,
          threshold_50 = FAIL_EARLY_50_DONORS_PCT
        )

        if (fail_check$should_stop) {
          message(sprintf("\n🔴 %s", fail_check$reason))
          message(sprintf("   Processed %d donors, %d cells so far for %s",
            n_processed, total_cells, cell_type))

          # Save partial results before stopping
          partial_cells <- do.call(rbind, all_cell_qc[!sapply(all_cell_qc, is.null)])
          partial_donors <- do.call(rbind, all_donor_summaries[!sapply(all_donor_summaries, is.null)])

          qs::qsave(partial_cells, file.path(output_dir, "qc_all_cells_PARTIAL.qs"))
          write.csv(partial_donors, file.path(output_dir, "qc_summary_per_donor_PARTIAL.csv"),
            row.names = FALSE)

          # Write fail-early report
          fail_report <- c(
            sprintf("FAIL-EARLY TRIGGERED for %s", cell_type),
            sprintf("Reason: %s", fail_check$reason),
            sprintf("Donors processed: %d / %d", n_processed, n_shards),
            sprintf("Cells processed: %d", total_cells),
            "",
            "Violation summary:",
            sprintf("  fail_low_umis:  %d (%.2f%%)",
              sum(partial_cells$fail_low_umis), mean(partial_cells$fail_low_umis) * 100),
            sprintf("  fail_high_umis: %d (%.2f%%)",
              sum(partial_cells$fail_high_umis), mean(partial_cells$fail_high_umis) * 100),
            sprintf("  fail_low_genes: %d (%.2f%%)",
              sum(partial_cells$fail_low_genes), mean(partial_cells$fail_low_genes) * 100),
            sprintf("  fail_high_genes: %d (%.2f%%)",
              sum(partial_cells$fail_high_genes), mean(partial_cells$fail_high_genes) * 100),
            sprintf("  fail_high_mito: %d (%.2f%%)",
              sum(partial_cells$fail_high_mito, na.rm = TRUE),
              mean(partial_cells$fail_high_mito, na.rm = TRUE) * 100),
            sprintf("  fail_any:       %d (%.2f%%)",
              sum(partial_cells$fail_any), mean(partial_cells$fail_any) * 100)
          )
          writeLines(fail_report, file.path(output_dir, "FAIL_EARLY_REPORT.txt"))

          # Stop processing this cell type but continue to next
          break
        } else {
          violation_rate <- mean(combined_so_far$fail_any, na.rm = TRUE) * 100
          message(sprintf("  ✅ Fail-early check passed at %d donors (%.2f%% violation rate)",
            n_processed, violation_rate))
        }
      }

    }, error = function(e) {
      message(sprintf("  [ERROR] Failed to process donor %s: %s", donor_id, e$message))
      failed_donors <<- c(failed_donors, donor_id)
    })
  }

  # Check if we stopped early
  if (n_processed < n_shards && file.exists(file.path(output_dir, "FAIL_EARLY_REPORT.txt"))) {
    message(sprintf("\n[%s] Stopped early. See FAIL_EARLY_REPORT.txt", cell_type))
    next
  }

  # Combine results
  message(sprintf("\n  Combining results for %s: %d donors, %d total cells",
    cell_type, n_processed, total_cells))

  all_cells_df <- do.call(rbind, all_cell_qc[!sapply(all_cell_qc, is.null)])
  all_donors_df <- do.call(rbind, all_donor_summaries[!sapply(all_donor_summaries, is.null)])

  # Save
  qs::qsave(all_cells_df, checkpoint_file)
  write.csv(all_donors_df, donor_summary_file, row.names = FALSE)

  message(sprintf("  Saved: %s", checkpoint_file))
  message(sprintf("  Saved: %s", donor_summary_file))

  # Print summary
  violations <- summarise_violations(all_cells_df)
  message(sprintf("\n  === %s QC Summary ===", cell_type))
  message(sprintf("  Total cells: %d", violations$n_total))
  for (vname in names(violations$counts)) {
    message(sprintf("  %s: %d (%.4f%%)", vname,
      violations$counts[[vname]], violations$percentages[[vname]]))
  }

  if (length(failed_donors) > 0) {
    message(sprintf("  Failed donors: %d (%s)", length(failed_donors),
      paste(head(failed_donors, 5), collapse = ", ")))
  }
}

message("\n[DONE] Step 01 complete.")
