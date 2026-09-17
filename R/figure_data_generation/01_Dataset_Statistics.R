# ==============================================================================
# STEP 01: Dataset Statistics — Compute from Shard Files
# ==============================================================================
# Reads every donor shard for all 5 datasets x 3 cell types and computes:
#   N_Donors, N_Cells, Mean/SD/Median Cells_Per_Donor, UMIs_Per_Cell, Genes_Per_Cell
#
# Output: Paper_Source_Data/dataset_statistics.csv (15 rows: 5 datasets x 3 CTs)
#
# This replaces the hardcoded source_data_1 data.frame in paper_final_draft.rmd,
# which contained placeholder values for Allen, Johansen, and Rexach.
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(dplyr)
  library(Matrix)
  library(BiocParallel)
})

if (!exists("PROJECT_ROOT")) source("config.R")

# Load shard utilities from glial_characteristics pipeline
source(file.path(PROJECT_ROOT, "R", "glial_characteristics", "utils.R"))

message("\n### STEP 01: DATASET STATISTICS ###\n")

dir.create(OUTPUT_BASE, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Process each dataset x cell type
# ==============================================================================

all_rows <- list()

for (ds_name in names(DATASETS)) {
  ds <- DATASETS[[ds_name]]

  for (ct in CELL_TYPES) {
    message(sprintf("\n--- %s / %s ---", ds_name, ct))

    shard_dir <- file.path(SHARD_BASE, ct, ds$shard_dir_suffix)

    if (!dir.exists(shard_dir)) {
      message(sprintf("  [SKIP] Shard directory not found: %s", shard_dir))
      next
    }

    # Load manifest
    manifest_path <- file.path(shard_dir, "00_MASTER_MANIFEST.qs")
    if (!file.exists(manifest_path)) {
      message(sprintf("  [SKIP] No manifest at: %s", manifest_path))
      next
    }

    manifest_obj <- qs::qread(manifest_path)
    donor_ids <- unique(as.character(manifest_obj$manifest$donor_id))
    n_donors <- length(donor_ids)
    message(sprintf("  Found %d donors in manifest", n_donors))

    # Extract metrics from all shards
    message(sprintf("  Loading shards (%d cores)...", N_CORES))
    BPPARAM <- MulticoreParam(workers = min(N_CORES, n_donors), progressbar = FALSE)

    results <- bplapply(donor_ids, function(did) {
      safe_id <- gsub("[^[:alnum:]]", "_", did)
      sp <- file.path(shard_dir, paste0("donor_", safe_id, ".qs"))
      if (!file.exists(sp)) return(NULL)
      sce <- qs::qread(sp)
      counts_mat <- assay(sce, "counts")
      umis <- as.numeric(Matrix::colSums(counts_mat))
      genes <- as.integer(Matrix::colSums(counts_mat > 0))
      list(
        donor_id = did,
        n_cells = ncol(sce),
        umis = umis,
        genes = genes
      )
    }, BPPARAM = BPPARAM)

    results <- results[!sapply(results, is.null)]
    if (length(results) == 0) {
      message("  [SKIP] No valid shards loaded")
      next
    }

    # Aggregate
    n_donors_loaded <- length(results)
    cells_per_donor <- sapply(results, function(r) r$n_cells)
    all_umis <- unlist(lapply(results, function(r) r$umis))
    all_genes <- unlist(lapply(results, function(r) r$genes))
    n_cells_total <- sum(cells_per_donor)

    message(sprintf("  Loaded %s cells from %d donors",
                    format(n_cells_total, big.mark = ","), n_donors_loaded))

    row <- data.frame(
      Dataset = ds$label,
      Brain_Region = ds$brain_region,
      Cell_Type = ct,
      N_Donors = n_donors_loaded,
      N_Cells = n_cells_total,
      Mean_Cells_Per_Donor = round(mean(cells_per_donor), 1),
      SD_Cells_Per_Donor = round(sd(cells_per_donor), 1),
      Median_Cells_Per_Donor = median(cells_per_donor),
      Mean_UMIs_Per_Cell = round(mean(all_umis), 0),
      SD_UMIs_Per_Cell = round(sd(all_umis), 0),
      Median_UMIs_Per_Cell = median(all_umis),
      Mean_Genes_Per_Cell = round(mean(all_genes), 0),
      SD_Genes_Per_Cell = round(sd(all_genes), 0),
      Median_Genes_Per_Cell = median(all_genes),
      stringsAsFactors = FALSE
    )

    all_rows[[length(all_rows) + 1]] <- row
  }
}

# ==============================================================================
# Combine and write
# ==============================================================================

dataset_stats <- bind_rows(all_rows)

out_path <- file.path(OUTPUT_BASE, "dataset_statistics.csv")
write.csv(dataset_stats, out_path, row.names = FALSE)

message(sprintf("\n[DONE] Wrote %d rows to %s", nrow(dataset_stats), out_path))
message(sprintf("  Datasets: %s", paste(unique(dataset_stats$Dataset), collapse = ", ")))
