# create_biomart_cache.R — Generate BioMart annotation cache from existing data
# Uses the DE Diagnostics gene_biotype_cache.csv (which has chromosome info from
# a previous BioMart query) as the primary source. Falls back to gene name matching
# for any genes not in that cache.

source("config.R")
library(qs)
library(dplyr)

# ==============================================================================
# 1. LOAD EXISTING CHROMOSOME ANNOTATIONS FROM DE DIAGNOSTICS CACHE
# ==============================================================================

DE_DIAG_CACHE <- file.path(dirname(FULL_GT_BASE), "DE_Diagnostics", "gene_biotype_cache.csv")

if (!file.exists(DE_DIAG_CACHE)) {
  stop("DE Diagnostics gene_biotype_cache.csv not found at: ", DE_DIAG_CACHE,
       "\nRun the DE diagnostics pipeline first to generate this cache.")
}

existing_cache <- read.csv(DE_DIAG_CACHE, stringsAsFactors = FALSE)
message(sprintf("Loaded %d gene annotations from DE Diagnostics cache", nrow(existing_cache)))
message(sprintf("  X: %d, Y: %d, Autosomal: %d",
                sum(existing_cache$chromosome_name == "X", na.rm = TRUE),
                sum(existing_cache$chromosome_name == "Y", na.rm = TRUE),
                sum(!existing_cache$chromosome_name %in% c("X", "Y"), na.rm = TRUE)))

# ==============================================================================
# 2. COLLECT ALL GENE IDs FROM DE RESULTS + BULK
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")
DATASETS   <- c("Roussos", "Allen", "DeJager", "Johansen", "Rexach")

all_ids <- character(0)

for (ct in CELL_TYPES) {
  for (ds in DATASETS) {
    f <- file.path(FULL_GT_BASE, ct, paste0(ds, "_Full_DESeq2_Results.qs"))
    if (file.exists(f)) {
      res <- tryCatch(as.data.frame(qread(f)), error = function(e) NULL)
      if (!is.null(res)) {
        ids <- sub("\\..*", "", rownames(res))
        all_ids <- c(all_ids, ids)
        message(sprintf("  Loaded %d genes from %s / %s", length(ids), ct, ds))
      }
    }
  }
}

if (file.exists(BULK_ALL_FILE)) {
  bulk <- qread(BULK_ALL_FILE)
  bulk_ids <- sub("\\..*", "", bulk$gene)
  all_ids <- c(all_ids, bulk_ids)
  message(sprintf("  Added %d genes from bulk", length(bulk_ids)))
}

all_ids <- unique(all_ids)
message(sprintf("Total unique Ensembl IDs: %d", length(all_ids)))

# ==============================================================================
# 3. BUILD CACHE — merge with existing annotations
# ==============================================================================

# Get gene names from the annotated CSV
gene_name_map <- character(0)
if (file.exists(ANNOTATED_FILE)) {
  desc <- read.csv(ANNOTATED_FILE, stringsAsFactors = FALSE)
  if (all(c("gene", "external_gene_name") %in% colnames(desc))) {
    clean_ids <- sub("\\..*", "", desc$gene)
    gene_name_map <- setNames(desc$external_gene_name, clean_ids)
  }
}

# Start with the existing cache for all genes that appear in it
in_cache <- all_ids[all_ids %in% existing_cache$ensembl_gene_id]
not_in_cache <- all_ids[!all_ids %in% existing_cache$ensembl_gene_id]

message(sprintf("  %d genes found in DE Diagnostics cache", length(in_cache)))
message(sprintf("  %d genes not in cache (will use default autosomal)", length(not_in_cache)))

# Build the output: known genes get real annotations, unknown get autosomal default
raw_coords <- existing_cache %>%
  filter(ensembl_gene_id %in% all_ids) %>%
  mutate(external_gene_name = ifelse(
    ensembl_gene_id %in% names(gene_name_map),
    gene_name_map[ensembl_gene_id],
    ""
  ))

# Add missing genes as autosomal
if (length(not_in_cache) > 0) {
  missing_df <- data.frame(
    ensembl_gene_id = not_in_cache,
    gene_biotype = "protein_coding",
    chromosome_name = "1",
    Biotype_Class = "Coding",
    external_gene_name = ifelse(
      not_in_cache %in% names(gene_name_map),
      gene_name_map[not_in_cache],
      ""
    ),
    stringsAsFactors = FALSE
  )
  raw_coords <- bind_rows(raw_coords, missing_df)
}

# Ensure we have the columns the bulk vis script expects
raw_coords <- raw_coords %>%
  select(ensembl_gene_id, chromosome_name, external_gene_name, gene_biotype)

# ==============================================================================
# 4. SAVE CACHE
# ==============================================================================

cache_path <- file.path(OUTPUT_DIR, "biomart_gene_coords_cache.qs")
qsave(raw_coords, cache_path)

message(sprintf(
  "\nSaved BioMart cache: %d genes (%d X, %d Y, %d autosomal) to %s",
  nrow(raw_coords),
  sum(raw_coords$chromosome_name == "X", na.rm = TRUE),
  sum(raw_coords$chromosome_name == "Y", na.rm = TRUE),
  sum(!raw_coords$chromosome_name %in% c("X", "Y"), na.rm = TRUE),
  cache_path
))
