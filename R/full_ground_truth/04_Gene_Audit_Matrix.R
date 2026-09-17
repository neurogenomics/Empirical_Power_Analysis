# ==============================================================================
# SCRIPT 04: Gene Audit Matrix — Gene Composition Across Datasets
# ==============================================================================
# Audits the gene composition of each dataset by categorising genes into
# overlapping groups: X chromosome, Y chromosome, PAR, protein-coding,
# pseudogene, lncRNA, mitochondrial, autosomal.
#
# Produces TWO audit matrices:
#   A. Full Gene Universe — all genes present in the shard count matrices
#   B. Significant Gene Universe — genes passing padj < 0.05 & |LFC| > 1
#
# Each matrix has rows = datasets, columns = gene categories (+ Total).
# Categories are NOT mutually exclusive (a gene can be X-chr AND pseudogene).
#
# Outputs (per cell type):
#   {OUTPUT_BASE}/{CellType}/Gene_Audit_Full_Universe.csv
#   {OUTPUT_BASE}/{CellType}/Gene_Audit_Significant.csv
#
# Combined outputs:
#   {OUTPUT_BASE}/Gene_Audit_Full_Universe_Combined.csv
#   {OUTPUT_BASE}/Gene_Audit_Significant_Combined.csv
#
# Prerequisites:
#   - Shards exist in Datalink/Shards/{CellType}/{Dataset}/
#   - Step 01 DE results exist in {OUTPUT_BASE}/{CellType}/{Dataset}_DE_results.csv
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/full_ground_truth/config.R",
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

  if (!config_loaded) {
    stop("Could not find config.R. Run from project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(qs)
  library(SingleCellExperiment)
  library(dplyr)
})

log_section("STEP 04: GENE AUDIT MATRIX")

# ==============================================================================
# 1. COLLECT GENE LISTS FROM ALL DATASETS × CELL TYPES
# ==============================================================================

log_section("STEP 4A: FULL GENE UNIVERSE AUDIT")

# Storage for per-cell-type results
all_full_matrices <- list()
all_sig_matrices <- list()

# Collect ALL gene IDs across all datasets and cell types for a single
# BioMart query (much more efficient than per-dataset queries)
all_gene_ids_global <- character(0)

# First pass: extract gene IDs from shards
full_gene_lists <- list()  # nested: full_gene_lists[[cell_type]][[dataset]]

for (cell_type in CELL_TYPES) {
  log_info("--- Extracting gene IDs for %s ---", cell_type)

  full_gene_lists[[cell_type]] <- list()

  for (d_name in names(DATASET_CONFIG)) {
    shard_dir <- get_shard_dir(cell_type, d_name)

    if (!dir.exists(shard_dir)) {
      log_warn("  [SKIP] %s/%s — shard directory not found", cell_type, d_name)
      next
    }

    gene_ids <- extract_shard_gene_ids(shard_dir)

    if (is.null(gene_ids) || length(gene_ids) == 0) {
      log_warn("  [SKIP] %s/%s — no genes extracted", cell_type, d_name)
      next
    }

    # Clean gene IDs (strip version suffixes)
    gene_ids <- clean_gene_ids(gene_ids)
    full_gene_lists[[cell_type]][[d_name]] <- gene_ids
    all_gene_ids_global <- union(all_gene_ids_global, gene_ids)

    log_info("  %s/%s: %d genes", cell_type, d_name, length(gene_ids))
  }
}

log_info("Total unique genes across all datasets and cell types: %d",
         length(all_gene_ids_global))

# ==============================================================================
# 2. QUERY BIOMART FOR ANNOTATIONS (SINGLE BATCH)
# ==============================================================================

log_section("BIOMART ANNOTATION")

cache_path <- file.path(OUTPUT_BASE, AUDIT_BIOMART_CACHE_FILE)
annotations <- get_audit_annotations(all_gene_ids_global, cache_path)

if (is.null(annotations) || nrow(annotations) == 0) {
  stop("[CRITICAL] Failed to obtain BioMart annotations. Cannot proceed.")
}

log_info("Annotations obtained for %d / %d genes (%.1f%%)",
         nrow(annotations), length(all_gene_ids_global),
         100 * nrow(annotations) / length(all_gene_ids_global))

# ==============================================================================
# 3. CLASSIFY GENES INTO CATEGORIES
# ==============================================================================

log_info("Classifying genes into categories...")
category_df <- classify_gene_categories(annotations)

if (is.null(category_df)) {
  stop("[CRITICAL] Gene classification failed. Cannot proceed.")
}

cat_cols <- setdiff(colnames(category_df), "ensembl_gene_id")
log_info("Categories: %s", paste(gsub("^is_", "", cat_cols), collapse = ", "))

# ==============================================================================
# 4A. BUILD FULL GENE UNIVERSE AUDIT MATRICES
# ==============================================================================

log_section("BUILDING FULL GENE UNIVERSE MATRICES")

for (cell_type in CELL_TYPES) {
  log_info("--- %s ---", cell_type)

  gene_lists <- full_gene_lists[[cell_type]]

  if (length(gene_lists) == 0) {
    log_warn("  No datasets available for %s. Skipping.", cell_type)
    next
  }

  audit_matrix <- build_audit_matrix(gene_lists, category_df)

  if (is.null(audit_matrix)) {
    log_warn("  Failed to build audit matrix for %s.", cell_type)
    next
  }

  # Save per-cell-type CSV
  output_dir <- file.path(OUTPUT_BASE, cell_type)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_path <- file.path(output_dir, AUDIT_FULL_UNIVERSE_FILE)
  write.csv(audit_matrix, out_path, row.names = FALSE)
  log_info("  Saved: %s", out_path)

  # Store for combined output
  audit_matrix$CellType <- cell_type
  all_full_matrices[[cell_type]] <- audit_matrix

  # Print summary
  log_info("  Full Universe Audit for %s:", cell_type)
  for (i in seq_len(nrow(audit_matrix))) {
    log_info("    %s: %d total genes", audit_matrix$Dataset[i], audit_matrix$Total[i])
  }
}

# ==============================================================================
# 4B. BUILD SIGNIFICANT GENE UNIVERSE AUDIT MATRICES
# ==============================================================================

log_section("STEP 4B: SIGNIFICANT GENE UNIVERSE AUDIT")

for (cell_type in CELL_TYPES) {
  log_info("--- %s ---", cell_type)

  sig_gene_lists <- list()

  for (d_name in names(DATASET_CONFIG)) {
    # Load DE results CSV
    de_path <- file.path(OUTPUT_BASE, cell_type,
                         paste0(d_name, "_DE_results.csv"))

    if (!file.exists(de_path)) {
      log_warn("  [SKIP] %s/%s — DE results not found: %s",
               cell_type, d_name, basename(de_path))
      next
    }

    de_df <- tryCatch({
      read.csv(de_path, stringsAsFactors = FALSE)
    }, error = function(e) {
      log_warn("  [SKIP] %s/%s — failed to read DE results: %s",
               cell_type, d_name, e$message)
      NULL
    })

    if (is.null(de_df)) next

    # Validate required columns
    if (!all(c("gene", "padj", "log2FoldChange") %in% colnames(de_df))) {
      log_warn("  [SKIP] %s/%s — missing required columns", cell_type, d_name)
      next
    }

    # Filter to significant genes
    sig_df <- de_df[!is.na(de_df$padj) &
                      de_df$padj < AUDIT_FDR_THRESHOLD &
                      abs(de_df$log2FoldChange) > AUDIT_LFC_THRESHOLD, ]

    if (nrow(sig_df) == 0) {
      log_warn("  %s/%s: 0 significant genes (padj < %g, |LFC| > %g)",
               cell_type, d_name, AUDIT_FDR_THRESHOLD, AUDIT_LFC_THRESHOLD)
      next
    }

    sig_genes <- clean_gene_ids(sig_df$gene)
    sig_gene_lists[[d_name]] <- sig_genes
    log_info("  %s/%s: %d significant genes", cell_type, d_name, length(sig_genes))
  }

  if (length(sig_gene_lists) == 0) {
    log_warn("  No significant gene lists for %s. Skipping.", cell_type)
    next
  }

  # Build audit matrix for significant genes
  sig_audit_matrix <- build_audit_matrix(sig_gene_lists, category_df)

  if (is.null(sig_audit_matrix)) {
    log_warn("  Failed to build significant audit matrix for %s.", cell_type)
    next
  }

  # Save per-cell-type CSV
  output_dir <- file.path(OUTPUT_BASE, cell_type)
  out_path <- file.path(output_dir, AUDIT_SIGNIFICANT_FILE)
  write.csv(sig_audit_matrix, out_path, row.names = FALSE)
  log_info("  Saved: %s", out_path)

  # Store for combined output
  sig_audit_matrix$CellType <- cell_type
  all_sig_matrices[[cell_type]] <- sig_audit_matrix

  # Print summary
  log_info("  Significant Gene Audit for %s:", cell_type)
  for (i in seq_len(nrow(sig_audit_matrix))) {
    log_info("    %s: %d significant genes (%d X_Chr, %d Y_Chr, %d PAR)",
             sig_audit_matrix$Dataset[i],
             sig_audit_matrix$Total[i],
             sig_audit_matrix$X_chr[i],
             sig_audit_matrix$Y_chr[i],
             sig_audit_matrix$PAR[i])
  }
}

# ==============================================================================
# 5. COMBINED SUMMARY TABLES
# ==============================================================================

log_section("COMBINED SUMMARY TABLES")

# Full Universe Combined
if (length(all_full_matrices) > 0) {
  combined_full <- do.call(rbind, all_full_matrices)
  rownames(combined_full) <- NULL

  # Reorder columns: CellType first
  combined_full <- combined_full[, c("CellType", setdiff(colnames(combined_full), "CellType"))]

  out_path <- file.path(OUTPUT_BASE, AUDIT_COMBINED_FULL_FILE)
  write.csv(combined_full, out_path, row.names = FALSE)
  log_info("Saved combined full universe audit: %s", out_path)

  # Print summary table
  log_info("=== FULL GENE UNIVERSE AUDIT (Combined) ===")
  for (ct in unique(combined_full$CellType)) {
    ct_data <- combined_full[combined_full$CellType == ct, ]
    log_info("  --- %s ---", ct)
    for (i in seq_len(nrow(ct_data))) {
      log_info("    %-10s: Total=%6d  X=%4d  Y=%3d  PAR=%3d  Coding=%6d  Pseudo=%5d  lncRNA=%5d  MT=%3d  Auto=%6d  Unann=%4d",
               ct_data$Dataset[i], ct_data$Total[i],
               ct_data$X_chr[i], ct_data$Y_chr[i], ct_data$PAR[i],
               ct_data$protein_coding[i], ct_data$pseudogene[i],
               ct_data$lncRNA[i], ct_data$mitochondrial[i],
               ct_data$autosomal[i], ct_data$Unannotated[i])
    }
  }
} else {
  log_warn("No full universe audit matrices generated.")
}

# Significant Combined
if (length(all_sig_matrices) > 0) {
  combined_sig <- do.call(rbind, all_sig_matrices)
  rownames(combined_sig) <- NULL

  # Reorder columns: CellType first
  combined_sig <- combined_sig[, c("CellType", setdiff(colnames(combined_sig), "CellType"))]

  out_path <- file.path(OUTPUT_BASE, AUDIT_COMBINED_SIG_FILE)
  write.csv(combined_sig, out_path, row.names = FALSE)
  log_info("Saved combined significant gene audit: %s", out_path)

  # Print summary table
  log_info("=== SIGNIFICANT GENE AUDIT (Combined) ===")
  for (ct in unique(combined_sig$CellType)) {
    ct_data <- combined_sig[combined_sig$CellType == ct, ]
    log_info("  --- %s ---", ct)
    for (i in seq_len(nrow(ct_data))) {
      log_info("    %-10s: Total=%4d  X=%3d  Y=%3d  PAR=%3d  Coding=%4d  Pseudo=%4d  lncRNA=%4d  MT=%2d  Auto=%4d  Unann=%3d",
               ct_data$Dataset[i], ct_data$Total[i],
               ct_data$X_chr[i], ct_data$Y_chr[i], ct_data$PAR[i],
               ct_data$protein_coding[i], ct_data$pseudogene[i],
               ct_data$lncRNA[i], ct_data$mitochondrial[i],
               ct_data$autosomal[i], ct_data$Unannotated[i])
    }
  }
} else {
  log_warn("No significant gene audit matrices generated.")
}

# ==============================================================================
# 6. CROSS-CELL-TYPE GENE SET COMPARISON
# ==============================================================================

log_section("CROSS-CELL-TYPE GENE SET COMPARISON")

# Check whether gene sets differ across cell types within each dataset
for (d_name in names(DATASET_CONFIG)) {
  ct_gene_sets <- list()

  for (cell_type in CELL_TYPES) {
    if (!is.null(full_gene_lists[[cell_type]][[d_name]])) {
      ct_gene_sets[[cell_type]] <- full_gene_lists[[cell_type]][[d_name]]
    }
  }

  if (length(ct_gene_sets) < 2) next

  # Check if all cell types have identical gene sets
  all_identical <- TRUE
  ref_set <- ct_gene_sets[[1]]

  for (i in 2:length(ct_gene_sets)) {
    if (!setequal(ref_set, ct_gene_sets[[i]])) {
      all_identical <- FALSE
      break
    }
  }

  if (all_identical) {
    log_info("  %s: Gene sets IDENTICAL across %d cell types (%d genes)",
             d_name, length(ct_gene_sets), length(ref_set))
  } else {
    # Report differences
    all_union <- Reduce(union, ct_gene_sets)
    all_intersect <- Reduce(intersect, ct_gene_sets)
    log_info("  %s: Gene sets DIFFER across cell types", d_name)
    log_info("    Union: %d genes, Intersection: %d genes, Difference: %d",
             length(all_union), length(all_intersect),
             length(all_union) - length(all_intersect))

    for (ct_name in names(ct_gene_sets)) {
      unique_to_ct <- setdiff(ct_gene_sets[[ct_name]], Reduce(union, ct_gene_sets[names(ct_gene_sets) != ct_name]))
      if (length(unique_to_ct) > 0) {
        log_info("    %s has %d unique genes", ct_name, length(unique_to_ct))
      }
    }
  }
}

log_section("STEP 04 COMPLETE")
