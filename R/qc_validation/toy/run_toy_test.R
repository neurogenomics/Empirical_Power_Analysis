# ==============================================================================
# SCRIPT: QC Validation Pipeline - Toy Integration Test
# PURPOSE: Stage 2 validation — run full pipeline on 10 donors per cell type
# RUN: Rscript R/qc_validation/toy/run_toy_test.R
# ==============================================================================

# ==============================================================================
# ==============================================================================
message("=== QC VALIDATION PIPELINE — TOY INTEGRATION TEST ===")
message(sprintf("Started: %s", Sys.time()))
message("Minimal run: 10 donors per cell type, Roussos_Full shards")
message("")

# ==============================================================================
# Source config and utils
# ==============================================================================
source(file.path(PROJECT_ROOT, "R", "qc_validation", "config.R"))
source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# Override output to toy directory
TOY_OUTPUT <- file.path(SCRIPT_DIR, "toy", "output")
dir.create(TOY_OUTPUT, recursive = TRUE, showWarnings = FALSE)

# Override fail-early to be more sensitive for toy test
TOY_N_DONORS <- 10

# ==============================================================================
# Step 1: Get mito genes
# ==============================================================================
log_section("Step 1: Mitochondrial Gene Lookup")

mito_genes <- get_mito_genes(cache_file = MITO_GENE_CACHE)
if (is.null(mito_genes)) {
  log_warn("Ensembl query failed. Using ^MT- prefix fallback.")
} else {
  log_info("Got %d mitochondrial genes from Ensembl", length(mito_genes))
}

# ==============================================================================
# Step 2: Process 10 donors per cell type
# ==============================================================================
log_section("Step 2: Compute QC Metrics (10 donors per cell type)")

for (cell_type in CELL_TYPES) {
  message(sprintf("\n--- %s ---", cell_type))

  shard_dir <- file.path(SHARD_BASE, cell_type, SHARD_DATASET)
  shard_files <- list_donor_shards(shard_dir)

  if (length(shard_files) == 0) {
    log_warn("No shards found for %s at %s", cell_type, shard_dir)
    next
  }

  # Take first 10
  toy_files <- head(shard_files, TOY_N_DONORS)
  log_info("Processing %d donors for %s (out of %d available)",
    length(toy_files), cell_type, length(shard_files))

  all_qc <- list()
  all_summaries <- list()

  for (i in seq_along(toy_files)) {
    donor_id <- extract_donor_id(toy_files[i])
    message(sprintf("  [%d/%d] %s", i, length(toy_files), donor_id))

    tryCatch({
      sce <- qs::qread(toy_files[i])

      # Sanity check on first donor
      if (i == 1) {
        counts_mat <- assay(sce, "counts")
        sanity_check_counts(counts_mat, sprintf("%s first donor", cell_type))

        # Report basic info
        log_info("%s: %d genes x %d cells in first shard",
          cell_type, nrow(sce), ncol(sce))
        log_info("Gene name examples: %s",
          paste(head(rownames(sce), 5), collapse = ", "))

        # Check for mito gene presence
        gene_names <- rownames(sce)
        is_mito <- identify_mito_genes(gene_names, mito_genes)
        log_info("Mito genes found: %d out of %d total genes",
          sum(is_mito), length(gene_names))
        if (sum(is_mito) > 0) {
          log_info("Mito gene names: %s",
            paste(gene_names[is_mito], collapse = ", "))
        }
      }

      qc_df <- compute_cell_qc_metrics(sce, mito_genes)
      qc_df <- check_thresholds(qc_df)
      qc_df$donor_id <- donor_id
      qc_df$cell_type <- cell_type

      all_qc[[i]] <- qc_df
      all_summaries[[i]] <- summarise_donor(qc_df, donor_id)

    }, error = function(e) {
      log_error("Failed on donor %s: %s", donor_id, e$message)
    })
  }

  # Combine
  if (length(all_qc) == 0) {
    log_error("No donors processed for %s", cell_type)
    next
  }

  cells_df <- do.call(rbind, all_qc[!sapply(all_qc, is.null)])
  donors_df <- do.call(rbind, all_summaries[!sapply(all_summaries, is.null)])

  # Save toy results
  toy_ct_dir <- file.path(TOY_OUTPUT, cell_type)
  dir.create(toy_ct_dir, recursive = TRUE, showWarnings = FALSE)
  qs::qsave(cells_df, file.path(toy_ct_dir, "qc_all_cells.qs"))
  write.csv(donors_df, file.path(toy_ct_dir, "qc_summary_per_donor.csv"), row.names = FALSE)

  # Print summary
  violations <- summarise_violations(cells_df)
  log_info("%s: %d cells from %d donors", cell_type, violations$n_total, nrow(donors_df))
  for (vname in names(violations$counts)) {
    log_info("  %s: %d (%.4f%%)", vname,
      violations$counts[[vname]], violations$percentages[[vname]])
  }

  # Quick histogram for visual inspection
  tryCatch({
    p <- ggplot(cells_df, aes(x = n_umis)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = QC_MIN_UMIS, colour = "red", linetype = "dashed") +
      geom_vline(xintercept = QC_MAX_UMIS, colour = "red", linetype = "dashed") +
      scale_x_log10(labels = scales::comma) +
      labs(title = sprintf("UMI Distribution — %s (Toy: %d donors)", cell_type, TOY_N_DONORS),
           x = "UMIs per Cell", y = "Count") +
      theme_minimal()
    ggsave(file.path(toy_ct_dir, "Toy_Histogram_UMIs.png"), p, width = 8, height = 5, dpi = 150)

    p <- ggplot(cells_df, aes(x = n_genes)) +
      geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
      geom_vline(xintercept = QC_MIN_GENES, colour = "red", linetype = "dashed") +
      geom_vline(xintercept = QC_MAX_GENES, colour = "red", linetype = "dashed") +
      labs(title = sprintf("Gene Count Distribution — %s (Toy: %d donors)", cell_type, TOY_N_DONORS),
           x = "Genes Detected", y = "Count") +
      theme_minimal()
    ggsave(file.path(toy_ct_dir, "Toy_Histogram_Genes.png"), p, width = 8, height = 5, dpi = 150)

    if (!all(is.na(cells_df$pct_mito))) {
      p <- ggplot(cells_df, aes(x = pct_mito)) +
        geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
        geom_vline(xintercept = QC_MAX_MITO_PCT, colour = "red", linetype = "dashed") +
        labs(title = sprintf("Mito %% Distribution — %s (Toy: %d donors)", cell_type, TOY_N_DONORS),
             x = "Mitochondrial UMI %", y = "Count") +
        theme_minimal()
      ggsave(file.path(toy_ct_dir, "Toy_Histogram_Mito.png"), p, width = 8, height = 5, dpi = 150)
    }

    log_info("Plots saved to %s", toy_ct_dir)
  }, error = function(e) {
    log_warn("Plot generation failed: %s", e$message)
  })
}

# ==============================================================================
# Summary
# ==============================================================================
log_section("Toy Test Summary")

message("\n=== TOY INTEGRATION TEST COMPLETE ===")
message(sprintf("Output directory: %s", TOY_OUTPUT))
message("Check the following:")
message("  1. qc_summary_per_donor.csv — are violation counts reasonable?")
message("  2. Toy_Histogram_UMIs.png — do cells fall within threshold lines?")
message("  3. Toy_Histogram_Genes.png — do cells fall within threshold lines?")
message("  4. Toy_Histogram_Mito.png — are cells below 5% mito?")
message("")
message("Awaiting user visual confirmation before Stage 3 production run.")

finalize_logging(status = "SUCCESS")
