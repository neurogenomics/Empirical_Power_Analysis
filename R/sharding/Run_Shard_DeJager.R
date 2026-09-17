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
# SCRIPT: DeJager H5Seurat Shard Creator (Filtered & Metadata-Enriched)
# PURPOSE: Reads .h5seurat files, joins external CSV metadata (DeJagerMeta.csv),
#          filters donors by Age (>=40) and Cell Count (>=30), injects all CSV
#          metadata into colData, and saves individual donor .qs shards.
# INPUT: DeJager .h5seurat files + DeJagerMeta.csv
# OUTPUT: Datalink/Shards/{CellType}/DeJager/donor_*.qs + manifest
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(rhdf5)
  library(Matrix)
  library(qs)
  library(dplyr)
  library(future)       
  library(future.apply) 
})

# ==============================================================================
# 0. PARALLELIZATION SETUP (RStudio Safe)
# ==============================================================================
N_WORKERS <- 16
if (N_WORKERS < 1) N_WORKERS <- 1

is_rstudio <- Sys.getenv("RSTUDIO") == "1"

if (.Platform$OS.type == "windows" || is_rstudio) {
  message(sprintf("Running in %s. Using 'multisession' plan.", 
                  ifelse(is_rstudio, "RStudio", "Windows")))
  plan(multisession, workers = N_WORKERS)
} else {
  message("Running in Terminal (Linux/Mac). Using 'multicore' plan.")
  plan(multicore, workers = N_WORKERS)
}

options(future.globals.maxSize = 50 * 1024^3)

# ==============================================================================
# HELPER FUNCTION: Read h5seurat manually
# ==============================================================================
read_h5seurat_manual <- function(h5_file) {
  
  # 1. Get Cell Names
  cell_names <- tryCatch({
    h5read(h5_file, "cell.names")
  }, error = function(e) { return(NULL) })
  
  # 2. Get Metadata
  meta_contents <- h5ls(h5_file)
  meta_entries <- subset(meta_contents, group == "/meta.data")
  
  meta_df <- data.frame(row.names = cell_names)
  
  if (nrow(meta_entries) > 0) {
    for (col_name in meta_entries$name) {
      col_path <- paste0("/meta.data/", col_name)
      col_val <- h5read(h5_file, col_path)
      
      if (length(col_val) == length(cell_names)) {
        meta_df[[col_name]] <- col_val
      }
    }
  }
  
  # 3. Read Counts Matrix
  assay_root <- "/assays/RNA"
  
  if (any(grepl("counts", subset(meta_contents, group == assay_root)$name))) {
    matrix_path <- paste0(assay_root, "/counts")
  } else {
    matrix_path <- paste0(assay_root, "/data")
  }
  
  message(sprintf("    Reading matrix from: %s", matrix_path))
  
  data_slot    <- h5read(h5_file, paste0(matrix_path, "/data"))
  indices_slot <- h5read(h5_file, paste0(matrix_path, "/indices"))
  indptr_slot  <- h5read(h5_file, paste0(matrix_path, "/indptr"))
  features     <- h5read(h5_file, paste0(assay_root, "/features"))
  
  counts_mat <- Matrix::sparseMatrix(
    i = indices_slot + 1,
    p = indptr_slot,
    x = as.numeric(data_slot),
    dims = c(length(features), length(cell_names)),
    dimnames = list(features, cell_names)
  )
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = meta_df
  )
  
  return(sce)
}

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

INPUT_FILES <- list(
  "Micro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "DeJager_microglia.h5Seurat"),
  "Oligo" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "DeJager_oligodendroglia.h5Seurat"),
  "Astro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "DeJager_astrocytes.h5Seurat")
)

BASE_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")
DONOR_ID_COL    <- "individualID"

# External metadata CSV (contains age_death, msex, pmi, etc.)
METADATA_FILE   <- file.path(PROJECT_ROOT, "Datalink", "raw_data", "DeJagerMeta.csv")
META_ID_COL     <- "individualID"   # Join column in CSV
AGE_COL         <- "age_death"      # Age column in CSV

# Filtering Thresholds
MIN_AGE         <- 40
MIN_CELLS       <- 30

CELL_STEPS      <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA)

# ==============================================================================
# 2. LOAD EXTERNAL METADATA
# ==============================================================================

message(sprintf(">>> Loading external metadata from: %s", METADATA_FILE))
if (!file.exists(path.expand(METADATA_FILE))) {
  stop("Metadata file not found: ", METADATA_FILE)
}

meta_db <- read.csv(path.expand(METADATA_FILE), stringsAsFactors = FALSE)
meta_db[[META_ID_COL]] <- as.character(meta_db[[META_ID_COL]])

if (!AGE_COL %in% names(meta_db)) {
  stop(sprintf("Age column '%s' not found in metadata CSV. Available: %s",
               AGE_COL, paste(names(meta_db), collapse = ", ")))
}

message(sprintf("    Loaded %d donor records with %d columns", 
                nrow(meta_db), ncol(meta_db)))
message(sprintf("    Age range: %.1f - %.1f (column: %s)",
                min(meta_db[[AGE_COL]], na.rm = TRUE),
                max(meta_db[[AGE_COL]], na.rm = TRUE), AGE_COL))

# Pre-filter metadata: identify donors passing age threshold
meta_db$age_numeric <- suppressWarnings(as.numeric(meta_db[[AGE_COL]]))
age_valid_ids <- meta_db[[META_ID_COL]][
  !is.na(meta_db$age_numeric) & meta_db$age_numeric >= MIN_AGE
]

message(sprintf("    Donors passing age filter (>= %d): %d / %d",
                MIN_AGE, length(age_valid_ids), nrow(meta_db)))

# ==============================================================================
# 3. CORE PROCESSING LOOP
# ==============================================================================

for (target_type in names(INPUT_FILES)) {
  
  input_path <- INPUT_FILES[[target_type]]
  safe_type_name <- gsub("[^[:alnum:]]", "_", target_type)
  type_shard_dir <- file.path(BASE_OUTPUT_DIR, safe_type_name, "DeJager")
  
  message(sprintf("\n>>> Processing Cell Type: %s", target_type))
  message(sprintf("    Input File: %s", input_path))
  
  if (!file.exists(path.expand(input_path))) {
    message(sprintf("    [ERROR] File not found: %s. Skipping...", input_path))
    next
  }
  
  if (!dir.exists(type_shard_dir)) dir.create(type_shard_dir, recursive = TRUE)
  
  # --- STEP A: PEEK METADATA (Efficiency Check) ---
  message("    Peeking at H5 metadata to identify donors...")
  
  all_donors_raw <- tryCatch({
    donor_path <- paste0("/meta.data/", DONOR_ID_COL)
    as.character(h5read(path.expand(input_path), donor_path))
  }, error = function(e) {
    message(sprintf("    [ERROR] Could not read donor column: %s", e$message))
    character(0)
  })
  
  # Remove NA donor IDs
  all_donors_raw <- all_donors_raw[!is.na(all_donors_raw) & all_donors_raw != "NA"]
  unique_donors_raw <- unique(all_donors_raw)
  
  if (length(unique_donors_raw) == 0) {
    message("    [SKIP] No valid donors found or metadata read failed.")
    next
  }
  
  # --- STEP B: APPLY AGE FILTER ---
  # Only keep donors that exist in the CSV AND pass the age threshold
  donors_with_meta <- intersect(unique_donors_raw, meta_db[[META_ID_COL]])
  donors_no_meta   <- setdiff(unique_donors_raw, meta_db[[META_ID_COL]])
  
  if (length(donors_no_meta) > 0) {
    message(sprintf("    [WARNING] %d donors have no metadata match (will be excluded): %s",
                    length(donors_no_meta),
                    paste(head(donors_no_meta, 5), collapse = ", ")))
  }
  
  donors_age_valid <- intersect(donors_with_meta, age_valid_ids)
  donors_age_fail  <- setdiff(donors_with_meta, age_valid_ids)
  
  message(sprintf("    [FILTER] Age >= %d: %d pass, %d fail (of %d with metadata)",
                  MIN_AGE, length(donors_age_valid), length(donors_age_fail),
                  length(donors_with_meta)))
  
  if (length(donors_age_valid) == 0) {
    message("    [SKIP] No donors passed age filter.")
    next
  }
  
  # --- STEP C: LOAD H5SEURAT & APPLY CELL COUNT FILTER ---
  message("    Loading .h5seurat file manually (Heavy Load)...")
  
  sce_full <- tryCatch({
    read_h5seurat_manual(path.expand(input_path))
  }, error = function(e) {
    message(sprintf("    [ERROR] Failed to load file: %s", e$message))
    return(NULL)
  })
  
  if (is.null(sce_full)) next
  
  if (!DONOR_ID_COL %in% names(colData(sce_full))) {
    message(sprintf("    [ERROR] Donor column '%s' missing from loaded object.", DONOR_ID_COL))
    next
  }
  
  # Subset to age-valid donors only
  donor_ids_loaded <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  keep_age <- donor_ids_loaded %in% donors_age_valid
  sce_full <- sce_full[, keep_age]
  donor_ids_loaded <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  
  message(sprintf("    After age filter: %d cells from %d donors",
                  ncol(sce_full), length(unique(donor_ids_loaded))))
  
  # Apply cell count filter
  donor_counts <- table(donor_ids_loaded)
  valid_donors <- names(donor_counts[donor_counts >= MIN_CELLS])
  low_count_donors <- names(donor_counts[donor_counts < MIN_CELLS])
  
  if (length(low_count_donors) > 0) {
    message(sprintf("    [FILTER] Cell count >= %d: %d pass, %d fail",
                    MIN_CELLS, length(valid_donors), length(low_count_donors)))
  }
  
  if (length(valid_donors) == 0) {
    message("    [SKIP] No donors passed cell count filter.")
    rm(sce_full); gc()
    next
  }
  
  # Subset to valid donors
  keep_cells <- donor_ids_loaded %in% valid_donors
  sce_full <- sce_full[, keep_cells]
  donor_ids_loaded <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  
  message(sprintf("    [FILTER] Final: %d cells from %d donors (Age >= %d, Cells >= %d)",
                  ncol(sce_full), length(valid_donors), MIN_AGE, MIN_CELLS))
  
  # --- STEP D: INJECT CSV METADATA INTO COLDATA ---
  message("    Injecting external metadata into colData...")
  
  # For each cell, look up its donor in the CSV and inject all columns
  cell_donor_ids <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  match_indices <- match(cell_donor_ids, meta_db[[META_ID_COL]])
  
  n_matched <- sum(!is.na(match_indices))
  if (n_matched != length(match_indices)) {
    message(sprintf("    [WARNING] %d / %d cells could not be matched to CSV metadata",
                    sum(is.na(match_indices)), length(match_indices)))
  }
  
  # Inject all CSV columns (except the join column, to avoid duplication)
  csv_cols_to_inject <- setdiff(names(meta_db), c(META_ID_COL, "age_numeric"))
  for (col in csv_cols_to_inject) {
    colData(sce_full)[[col]] <- meta_db[[col]][match_indices]
  }
  
  message(sprintf("    Injected %d metadata columns: %s",
                  length(csv_cols_to_inject),
                  paste(csv_cols_to_inject, collapse = ", ")))
  
  # --- STEP E: DETERMINE WHICH SHARDS TO WRITE ---
  existing_files <- list.files(type_shard_dir, pattern = "^donor_.*\\.qs$")
  
  donor_map <- data.frame(
    id = valid_donors,
    safe_id = gsub("[^[:alnum:]]", "_", valid_donors),
    stringsAsFactors = FALSE
  )
  donor_map$filename <- paste0("donor_", donor_map$safe_id, ".qs")
  
  donors_to_process <- donor_map$id[!donor_map$filename %in% existing_files]
  
  message(sprintf("    Total Valid Donors: %d | Already Sharded: %d | Remaining: %d",
                  length(valid_donors),
                  length(valid_donors) - length(donors_to_process),
                  length(donors_to_process)))
  
  # --- STEP F: PARALLEL SHARDING ---
  if (length(donors_to_process) == 0) {
    message("    [SKIP] All valid donor shards already exist.")
  } else {
    message(sprintf("    Processing %d donors in parallel (%d workers)...",
                    length(donors_to_process), N_WORKERS))
    
    future_lapply(donors_to_process, function(id) {
      safe_id <- gsub("[^[:alnum:]]", "_", id)
      donor_file <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
      
      donor_sce <- sce_full[, donor_ids_loaded == id]
      qs::qsave(donor_sce, donor_file, nthreads = 1)
      
      return(TRUE)
    }, future.seed = TRUE)
  }
  
  # Clean up memory
  rm(sce_full, donor_ids_loaded)
  gc()
  
  # --- STEP G: GENERATE MANIFEST ---
  message("    Building/Updating Manifest...")
  
  all_shards <- list.files(type_shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  # Exclude donor_NA.qs if it exists from a previous unfiltered run
  all_shards <- all_shards[!grepl("donor_NA\\.qs$", all_shards)]
  
  if (length(all_shards) == 0) {
    message("    [WARNING] No shard files found. Skipping manifest.")
    next
  }
  
  manifest_list <- future_lapply(all_shards, function(f) {
    tmp <- qs::qread(f, nthreads = 1)
    data.frame(
      donor_id   = as.character(unique(colData(tmp)[[DONOR_ID_COL]])),
      cell_count = ncol(tmp),
      file_path  = f,
      stringsAsFactors = FALSE
    )
  }, future.seed = TRUE)
  
  manifest_df <- do.call(rbind, manifest_list)
  
  if (is.null(manifest_df) || nrow(manifest_df) == 0) {
    message("    [WARNING] Manifest empty.")
    next
  }
  
  # Filter manifest to match current thresholds (in case old shards exist)
  manifest_df <- manifest_df %>% filter(cell_count >= MIN_CELLS)
  
  feasibility_map <- lapply(CELL_STEPS, function(n) {
    if (is.na(n)) return(manifest_df$donor_id)
    manifest_df %>% filter(cell_count >= n) %>% pull(donor_id)
  })
  names(feasibility_map) <- paste0("step_", ifelse(is.na(CELL_STEPS), "Max", CELL_STEPS))
  
  sharding_final <- list(
    manifest    = manifest_df,
    feasibility = feasibility_map,
    config      = list(
      donor_col   = DONOR_ID_COL,
      n_donors    = nrow(manifest_df),
      total_cells = sum(manifest_df$cell_count),
      cell_filter = target_type,
      min_age     = MIN_AGE,
      min_cells   = MIN_CELLS
    )
  )
  
  qs::qsave(sharding_final, file.path(type_shard_dir, "00_MASTER_MANIFEST.qs"), nthreads = 4)
  
  message(sprintf("    [DONE] %d donors, %d total cells\n",
                  sharding_final$config$n_donors,
                  sharding_final$config$total_cells))
}

message("\n================================================================")
message(" ALL PROCESSING COMPLETE ")
message(sprintf(" Filters Applied: Age >= %d, Cells >= %d", MIN_AGE, MIN_CELLS))
message(sprintf(" Metadata Source: %s", METADATA_FILE))
message("================================================================")

plan(sequential)
