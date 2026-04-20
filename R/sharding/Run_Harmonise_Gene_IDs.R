# ==============================================================================
# Shared config (resolves PROJECT_ROOT via here::here() or POWER_ANALYSIS_ROOT env var).
# ==============================================================================
if (!exists("PROJECT_ROOT")) {
  .sharding_config <- tryCatch({
    # commandArgs path works under `Rscript R/sharding/<this>.R`
    .args <- commandArgs(trailingOnly = FALSE)
    .file_arg <- grep("--file=", .args, value = TRUE)
    if (length(.file_arg)) {
      file.path(dirname(sub("--file=", "", .file_arg[1])), "config.R")
    } else if (requireNamespace("here", quietly = TRUE)) {
      file.path(here::here(), "R", "sharding", "config.R")
    } else {
      file.path(getwd(), "R", "sharding", "config.R")
    }
  }, error = function(e) file.path(getwd(), "R", "sharding", "config.R"))
  source(.sharding_config)
  rm(.sharding_config)
}

# ==============================================================================
# SCRIPT: Harmonise Gene IDs Across All Datasets
# PURPOSE: Convert gene symbols to Ensembl IDs in DeJager and Johansen shards,
#          delete stale checkpoints/DE results for affected datasets, so the
#          full_ground_truth pipeline can be re-run with consistent Ensembl IDs.
# INPUT: Existing shard files in Datalink/Shards/{CellType}/{Dataset}/
# OUTPUT: Updated shard files with Ensembl ID rownames; stale files deleted
# ==============================================================================
#
# BACKGROUND:
#   - Roussos, Allen: Already use Ensembl IDs (no action needed)
#   - Rexach: Shards have Ensembl IDs but stale PB checkpoint has NULL rownames
#     (delete checkpoint + DE results, re-run pipeline)
#   - DeJager: Shards use gene symbols (e.g., XIST, FAM87B)
#     (convert symbols → Ensembl, delete checkpoint + DE results)
#   - Johansen: Shards use gene symbols (e.g., WASH7P, LOC105379213)
#     (convert symbols → Ensembl, delete checkpoint + DE results)
#
# USAGE:
#   Rscript R/sharding/Run_Harmonise_Gene_IDs.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(SingleCellExperiment)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_SHARD_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")
BASE_GT_DIR    <- file.path(PROJECT_ROOT, "Datalink", "Full_Ground_Truth")
CELL_TYPES     <- c("Astro", "Micro", "Oligo")

# Datasets that need symbol → Ensembl conversion
# NOTE (2026-02-20): Johansen already converted in previous run; only DeJager remains
CONVERT_DATASETS <- c("DeJager")

# Datasets that only need checkpoint/DE deletion (shards are fine)
# NOTE (2026-02-20): Disabled — DE results are correct (used cached pseudobulk with Ensembl IDs)
CHECKPOINT_ONLY_DATASETS <- c()

# All affected datasets
ALL_AFFECTED <- c(CONVERT_DATASETS, CHECKPOINT_ONLY_DATASETS)

# ==============================================================================
# 1. BUILD SYMBOL → ENSEMBL MAPPING
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  GENE ID HARMONISATION\n")
cat("================================================================\n\n")

cat("--- Building symbol → Ensembl mapping from org.Hs.eg.db ---\n")

# Get all symbol → Ensembl mappings
# org.Hs.eg.db maps: SYMBOL → ENTREZID → ENSEMBL
# We go SYMBOL → ENSEMBL directly via the ENSEMBL column

# Method: Use mapIds with SYMBOL as key, ENSEMBL as column
# This handles the many-to-many mapping by taking the first match
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
cat(sprintf("  Total symbols in org.Hs.eg.db: %d\n", length(all_symbols)))

# Build mapping: symbol → ensembl (first match for each symbol)
symbol_to_ensembl <- mapIds(
  org.Hs.eg.db,
  keys = all_symbols,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove NAs
symbol_to_ensembl <- symbol_to_ensembl[!is.na(symbol_to_ensembl)]
cat(sprintf("  Symbols with Ensembl mapping: %d (%.1f%%)\n",
            length(symbol_to_ensembl),
            100 * length(symbol_to_ensembl) / length(all_symbols)))

# Also build alias → ensembl for genes that might use aliases
# (e.g., LOC* identifiers, older gene names)
all_aliases <- tryCatch({
  keys(org.Hs.eg.db, keytype = "ALIAS")
}, error = function(e) character(0))

if (length(all_aliases) > 0) {
  cat(sprintf("  Total aliases in org.Hs.eg.db: %d\n", length(all_aliases)))
  
  alias_to_ensembl <- mapIds(
    org.Hs.eg.db,
    keys = all_aliases,
    column = "ENSEMBL",
    keytype = "ALIAS",
    multiVals = "first"
  )
  alias_to_ensembl <- alias_to_ensembl[!is.na(alias_to_ensembl)]
  cat(sprintf("  Aliases with Ensembl mapping: %d\n", length(alias_to_ensembl)))
} else {
  alias_to_ensembl <- character(0)
}

# Combined lookup function: try SYMBOL first, then ALIAS
lookup_ensembl <- function(gene_symbols) {
  # First pass: direct symbol lookup
  mapped <- symbol_to_ensembl[gene_symbols]
  
  # Second pass: alias lookup for unmapped genes
  unmapped_mask <- is.na(mapped)
  if (any(unmapped_mask)) {
    unmapped_genes <- gene_symbols[unmapped_mask]
    alias_mapped <- alias_to_ensembl[unmapped_genes]
    mapped[unmapped_mask] <- alias_mapped
  }
  
  return(mapped)
}

# ==============================================================================
# 2. CONVERT SHARDS (DeJager, Johansen)
# ==============================================================================

total_converted <- 0
total_dropped <- 0
total_shards <- 0

for (dataset in CONVERT_DATASETS) {
  cat(sprintf("\n────────────────────────────────────────────────────────────────\n"))
  cat(sprintf("  Converting: %s\n", dataset))
  cat(sprintf("────────────────────────────────────────────────────────────────\n\n"))
  
  for (ct in CELL_TYPES) {
    shard_dir <- file.path(BASE_SHARD_DIR, ct, dataset)
    
    if (!dir.exists(shard_dir)) {
      cat(sprintf("  [SKIP] %s/%s - directory not found\n", ct, dataset))
      next
    }
    
    shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
    # Exclude donor_NA.qs
    shard_files <- shard_files[!grepl("donor_NA\\.qs$", shard_files)]
    
    if (length(shard_files) == 0) {
      cat(sprintf("  [SKIP] %s/%s - no shard files\n", ct, dataset))
      next
    }
    
    cat(sprintf("  %s/%s: %d shards\n", ct, dataset, length(shard_files)))
    
    # Load first shard to build the mapping for this dataset/celltype
    # (all shards in same dir should have same gene set)
    first_sce <- qs::qread(shard_files[1])
    original_genes <- rownames(first_sce)
    n_original <- length(original_genes)
    rm(first_sce); gc(verbose = FALSE)
    
    # Map symbols to Ensembl
    ensembl_ids <- lookup_ensembl(original_genes)
    
    # Statistics
    n_mapped <- sum(!is.na(ensembl_ids))
    n_unmapped <- sum(is.na(ensembl_ids))
    n_dup <- sum(duplicated(ensembl_ids[!is.na(ensembl_ids)]))
    
    cat(sprintf("    Genes: %d total, %d mapped (%.1f%%), %d unmapped, %d duplicate Ensembl IDs\n",
                n_original, n_mapped, 100 * n_mapped / n_original, n_unmapped, n_dup))
    
    # Show some unmapped examples
    if (n_unmapped > 0) {
      unmapped_examples <- head(original_genes[is.na(ensembl_ids)], 10)
      cat(sprintf("    Unmapped examples: %s\n", paste(unmapped_examples, collapse = ", ")))
    }
    
    # Build the keep mask: mapped AND not duplicate
    keep_mask <- !is.na(ensembl_ids)
    
    # Handle duplicates: keep first occurrence only
    if (n_dup > 0) {
      mapped_ensembl <- ensembl_ids[keep_mask]
      dup_mask_within_mapped <- duplicated(mapped_ensembl)
      # Convert back to full-length mask
      keep_indices <- which(keep_mask)
      drop_indices <- keep_indices[dup_mask_within_mapped]
      keep_mask[drop_indices] <- FALSE
      cat(sprintf("    After dedup: keeping %d genes\n", sum(keep_mask)))
    }
    
    final_ensembl <- ensembl_ids[keep_mask]
    n_final <- length(final_ensembl)
    
    cat(sprintf("    Final: %d genes (dropped %d = %.1f%%)\n",
                n_final, n_original - n_final,
                100 * (n_original - n_final) / n_original))
    
    # Process each shard
    for (sf in shard_files) {
      sce <- qs::qread(sf)
      
      # Subset to keep genes and rename
      sce <- sce[keep_mask, ]
      rownames(sce) <- final_ensembl
      
      # Clear rowData gene_ids column if it exists (it would still have symbols)
      if ("gene_ids" %in% colnames(rowData(sce))) {
        rowData(sce)$gene_ids <- final_ensembl
      }
      
      # Save back
      qs::qsave(sce, sf)
      total_shards <- total_shards + 1
    }
    
    # Also handle donor_NA.qs if it exists
    na_shard <- file.path(shard_dir, "donor_NA.qs")
    if (file.exists(na_shard)) {
      sce <- qs::qread(na_shard)
      sce <- sce[keep_mask, ]
      rownames(sce) <- final_ensembl
      if ("gene_ids" %in% colnames(rowData(sce))) {
        rowData(sce)$gene_ids <- final_ensembl
      }
      qs::qsave(sce, na_shard)
    }
    
    total_converted <- total_converted + n_final * length(shard_files)
    total_dropped <- total_dropped + (n_original - n_final) * length(shard_files)
    
    cat(sprintf("    [OK] %d shards updated\n", length(shard_files)))
    gc(verbose = FALSE)
  }
}

# ==============================================================================
# 3. DELETE STALE CHECKPOINTS AND DE RESULTS
# ==============================================================================

cat(sprintf("\n────────────────────────────────────────────────────────────────\n"))
cat(sprintf("  Deleting stale checkpoints and DE results\n"))
cat(sprintf("────────────────────────────────────────────────────────────────\n\n"))

files_deleted <- 0

for (dataset in ALL_AFFECTED) {
  for (ct in CELL_TYPES) {
    gt_dir <- file.path(BASE_GT_DIR, ct)
    
    if (!dir.exists(gt_dir)) next
    
    # Files to delete for this dataset
    patterns <- c(
      sprintf("%s_pb_checkpoint.qs", dataset),
      sprintf("%s_dds_checkpoint.qs", dataset),
      sprintf("%s_DE_results.csv", dataset),
      sprintf("%s_Full_DESeq2_Results.qs", dataset),
      sprintf("%s_%s_log.txt", dataset, ct)
    )
    
    for (pat in patterns) {
      f <- file.path(gt_dir, pat)
      if (file.exists(f)) {
        unlink(f)
        cat(sprintf("  [DELETED] %s/%s\n", ct, pat))
        files_deleted <- files_deleted + 1
      }
    }
  }
}

# Also check the old output location (Cross-Dataset-Analysis)
old_gt_dir <- file.path(PROJECT_ROOT, "Datalink", "Cross-Dataset-Analysis", "Full_Set_Outputs")
if (dir.exists(old_gt_dir)) {
  for (dataset in ALL_AFFECTED) {
    for (ct in CELL_TYPES) {
      old_ct_dir <- file.path(old_gt_dir, ct)
      if (!dir.exists(old_ct_dir)) next
      
      patterns <- c(
        sprintf("%s_pb_checkpoint.qs", dataset),
        sprintf("%s_dds_checkpoint.qs", dataset)
      )
      
      for (pat in patterns) {
        f <- file.path(old_ct_dir, pat)
        if (file.exists(f)) {
          unlink(f)
          cat(sprintf("  [DELETED] (old location) %s/%s\n", ct, pat))
          files_deleted <- files_deleted + 1
        }
      }
    }
  }
}

# ==============================================================================
# 4. VERIFICATION
# ==============================================================================

cat(sprintf("\n────────────────────────────────────────────────────────────────\n"))
cat(sprintf("  Verification\n"))
cat(sprintf("────────────────────────────────────────────────────────────────\n\n"))

# Spot-check one shard from each converted dataset
for (dataset in CONVERT_DATASETS) {
  shard_dir <- file.path(BASE_SHARD_DIR, "Astro", dataset)
  shard_files <- list.files(shard_dir, pattern = "^donor_(?!NA).*\\.qs$", 
                            full.names = TRUE, perl = TRUE)
  if (length(shard_files) > 0) {
    sce <- qs::qread(shard_files[1])
    rn <- head(rownames(sce), 5)
    all_ensg <- all(grepl("^ENSG", rownames(sce)))
    cat(sprintf("  %s: %d genes, all ENSG=%s, sample=[%s]\n",
                dataset, nrow(sce), all_ensg, paste(rn, collapse = ", ")))
    rm(sce); gc(verbose = FALSE)
  }
}

# Verify Rexach shards are untouched (should already have ENSG)
rexach_dir <- file.path(BASE_SHARD_DIR, "Astro", "Rexach")
rexach_files <- list.files(rexach_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
if (length(rexach_files) > 0) {
  sce <- qs::qread(rexach_files[1])
  rn <- head(rownames(sce), 5)
  all_ensg <- all(grepl("^ENSG", rownames(sce)))
  cat(sprintf("  Rexach: %d genes, all ENSG=%s, sample=[%s]\n",
              nrow(sce), all_ensg, paste(rn, collapse = ", ")))
  rm(sce); gc(verbose = FALSE)
}

# Verify no stale checkpoints remain
cat("\n  Checking for remaining checkpoints:\n")
for (dataset in ALL_AFFECTED) {
  for (ct in CELL_TYPES) {
    pb_path <- file.path(BASE_GT_DIR, ct, sprintf("%s_pb_checkpoint.qs", dataset))
    de_path <- file.path(BASE_GT_DIR, ct, sprintf("%s_DE_results.csv", dataset))
    cat(sprintf("    %s/%s: pb_checkpoint=%s, DE_results=%s\n",
                ct, dataset, file.exists(pb_path), file.exists(de_path)))
  }
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat(sprintf("\n================================================================\n"))
cat(sprintf("  SUMMARY\n"))
cat(sprintf("================================================================\n\n"))
cat(sprintf("  Shards converted: %d (across %s)\n", total_shards, 
            paste(CONVERT_DATASETS, collapse = ", ")))
cat(sprintf("  Files deleted: %d (checkpoints + DE results)\n", files_deleted))
cat(sprintf("\n  NEXT STEPS:\n"))
cat(sprintf("  1. Re-run DE pipeline: Rscript R/full_ground_truth/run_pipeline.R\n"))
cat(sprintf("     (This will re-aggregate from shards and re-run DESeq2 for %s)\n",
            paste(ALL_AFFECTED, collapse = ", ")))
cat(sprintf("  2. Re-run diagnostics: Step 02 will use the new Ensembl IDs\n"))
cat(sprintf("\n================================================================\n"))
