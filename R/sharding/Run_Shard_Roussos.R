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
# SCRIPT: Optimized Shard Creator with Integrated Pre-Filtering
# DESC:   Reads H5AD files, filters donors by Age/Cell Count/Type *before* #         sharding, and saves individual donor .qs files.
# ==============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(qs)
  library(dplyr)
  library(Matrix)
})

# --- SETUP PYTHON ---
if (!py_module_available("anndata")) {
  reticulate::py_install('anndata')
}
ad <- import("anndata")

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# PATHS
INPUT_FILES <- c(
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-HBCC.h5ad"), 
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-RADC.h5ad"), 
  file.path(PROJECT_ROOT, "Datalink", "raw_data", "PsychAD-MSSM.h5ad")
)

# Output Directory (Will be created if it doesn't exist)
BASE_OUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Roussos_Analysis")
TARGET_CELL_TYPE <- "endothelial cell" # Change this for other runs (e.g. "Oligo", "Micro")

# Define where shards go based on cell type
SHARD_DIR <- file.path(BASE_OUT_DIR, gsub(" ", "_", TARGET_CELL_TYPE), "Donor_Shards")
PROGRESS_LOG <- file.path(SHARD_DIR, "sharding_progress.rds")

# DATA COLUMNS
DONOR_ID_COL  <- "donor_id"   
CELL_TYPE_COL <- "cell_type" 
AGE_COL       <- "development_stage" # Column to check for age strings

# FILTERING THRESHOLDS
MIN_AGE_THRESHOLD   <- 40  # Minimum years
MIN_CELLS_THRESHOLD <- 30  # Minimum cells per donor (within the specific cell type)

# PROCESSING
CHUNK_SIZE <- 50000        
CELL_STEPS <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA) 

if(!dir.exists(SHARD_DIR)) dir.create(SHARD_DIR, recursive = TRUE)

# ==============================================================================
# 2. HELPER FUNCTIONS (Age Parsing & Resume Logic)
# ==============================================================================

# --- Resume Helpers ---
get_progress <- function() {
  if (file.exists(PROGRESS_LOG)) return(readRDS(PROGRESS_LOG))
  return(data.frame(file = character(), chunk = integer(), stringsAsFactors = FALSE))
}

update_progress <- function(file_name, chunk_idx) {
  progress <- get_progress()
  progress <- rbind(progress, data.frame(file = file_name, chunk = chunk_idx))
  saveRDS(progress, PROGRESS_LOG)
}

# --- Age Parsing Helpers ---
is_unknown_string <- function(x) {
  s <- tolower(trimws(as.character(x)))
  grepl("^(unknown|not available|na|n/a|nan)$", s)
}

looks_like_age_stage <- function(x) {
  if (!(is.character(x) || is.factor(x))) return(FALSE)
  s <- tolower(trimws(as.character(x)))
  matches <- grepl("(year|month|week|day)[ -]?old.*stage", s)
  return(sum(matches, na.rm=TRUE) > (0.5 * length(matches))) 
}

parse_age_stage_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
  # Extract numeric part
  n <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s)))
  
  is_year    <- grepl("year[ -]?old.*stage", s)
  is_month   <- grepl("month[ -]?old.*stage", s)
  is_week    <- grepl("week[ -]?old.*stage", s)
  is_day     <- grepl("day[ -]?old.*stage", s)
  is_unknown <- is_unknown_string(s)
  
  y <- rep(NA_real_, length(s))
  y[is_year]  <- n[is_year]
  y[is_month] <- n[is_month] / 12                     
  y[is_week]  <- n[is_week]  / 52.1429                     
  y[is_day]   <- n[is_day]   / 365.25                     
  y[is_unknown] <- -1.0
  
  # Fallback: if purely numeric, assume years
  is_numeric_str <- !is.na(suppressWarnings(as.numeric(s))) & is.na(y)
  y[is_numeric_str] <- as.numeric(s[is_numeric_str])
  
  return(y)
}

# ==============================================================================
# 3. CORE PROCESSING LOOP
# ==============================================================================

for (f_path in INPUT_FILES) {
  if (!file.exists(f_path)) {
    message(paste0("!!! [SKIP] File not found: ", f_path))
    next
  }
  
  f_name <- basename(f_path)
  message(paste0("\n>>> Opening File: ", f_name))
  
  # Open in backed mode to save RAM
  adata <- ad$read_h5ad(f_path, backed = "r")
  
  # --------------------------------------------------------------------------
  # STEP A: Pre-Scan Metadata to Identify Valid Donors
  # --------------------------------------------------------------------------
  message("   [PRE-SCAN] Loading metadata to determine valid donors...")
  
  # Load only metadata into R
  full_meta_df <- as.data.frame(adata$obs)
  
  # 1. Filter by Cell Type
  type_matches <- full_meta_df[[CELL_TYPE_COL]] == TARGET_CELL_TYPE
  target_meta  <- full_meta_df[type_matches, ]
  
  if (nrow(target_meta) == 0) {
    message(sprintf("   [SKIP] No cells of type '%s' found in this file.", TARGET_CELL_TYPE))
    adata$file$close()
    next
  }
  
  # 2. Parse Ages
  message("   [PRE-SCAN] Parsing ages...")
  raw_ages <- as.character(target_meta[[AGE_COL]])
  parsed_ages <- parse_age_stage_to_years(raw_ages)
  
  # 3. Group by Donor to check counts and age
  #    (We assume 1 donor has 1 age; if mixed, we take the first valid one)
  donor_stats <- target_meta %>%
    mutate(age_years = parsed_ages) %>%
    group_by(!!sym(DONOR_ID_COL)) %>%
    summarise(
      n_cells = n(),
      donor_age = max(age_years, na.rm = TRUE) # Use max to avoid -1 if mixed with valid
    )
  
  # 4. Apply Filters
  valid_donors_df <- donor_stats %>%
    filter(n_cells >= MIN_CELLS_THRESHOLD) %>%
    filter(donor_age >= MIN_AGE_THRESHOLD)
  
  valid_donor_ids <- unique(valid_donors_df[[DONOR_ID_COL]])
  
  message(sprintf("   [FILTER] Initial Donors: %d -> Valid Donors: %d", 
                  nrow(donor_stats), length(valid_donor_ids)))
  message(sprintf("   [FILTER] Criteria: Age >= %s, Cells >= %s", 
                  MIN_AGE_THRESHOLD, MIN_CELLS_THRESHOLD))
  
  if (length(valid_donor_ids) == 0) {
    message("   [SKIP] No donors met the criteria.")
    adata$file$close()
    next
  }
  
  # --------------------------------------------------------------------------
  # STEP B: Chunk Processing (Only reading/saving valid donors)
  # --------------------------------------------------------------------------
  
  n_cells  <- as.numeric(adata$n_obs)
  n_chunks <- ceiling(n_cells / CHUNK_SIZE)
  gene_names <- rownames(adata$var)
  
  completed_chunks <- get_progress() %>% 
    filter(file == f_name) %>% 
    pull(chunk)
  
  for (i in seq_len(n_chunks)) {
    if (i %in% completed_chunks) {
      message(sprintf("      [SKIP] Chunk %d/%d (Already processed)", i, n_chunks))
      next
    }
    
    # Calculate Python indices (0-based)
    py_idx_start <- (i - 1) * CHUNK_SIZE
    py_idx_end   <- min(i * CHUNK_SIZE, n_cells)
    
    # Calculate R indices (1-based) for metadata subsetting
    r_idx_start <- py_idx_start + 1
    r_idx_end   <- r_idx_start + (py_idx_end - py_idx_start) - 1
    
    # Get metadata for this chunk
    chunk_meta_all <- full_meta_df[r_idx_start:r_idx_end, , drop = FALSE]
    
    # ----------------------------------------------------------------------
    # CRITICAL: Filter Indices for this Chunk
    # Keep row IF: (Is Target Cell Type) AND (Is Valid Donor)
    # ----------------------------------------------------------------------
    keep_indices <- which(
      chunk_meta_all[[CELL_TYPE_COL]] == TARGET_CELL_TYPE & 
        chunk_meta_all[[DONOR_ID_COL]] %in% valid_donor_ids
    )
    
    if (length(keep_indices) == 0) {
      # No useful cells in this chunk, mark done and continue
      if(i %% 10 == 0) message(sprintf("      [INFO] Chunk %d/%d: No valid cells.", i, n_chunks))
      update_progress(f_name, i)
      next
    }
    
    message(sprintf("      [PROC] Chunk %d/%d: Processing %d valid cells...", 
                    i, n_chunks, length(keep_indices)))
    
    # Pull data from Python (Only load the chunk)
    chunk_slice <- adata[py_idx_start:(py_idx_end - 1)]
    
    # Transpose to Sparse Matrix (Genes x Cells) and subset
    counts_sparse <- t(chunk_slice$X)
    counts_sparse <- counts_sparse[, keep_indices, drop = FALSE]
    
    # Subset metadata
    chunk_meta_filtered <- chunk_meta_all[keep_indices, , drop = FALSE]
    
    rownames(counts_sparse) <- gene_names
    colnames(counts_sparse) <- rownames(chunk_meta_filtered)
    
    # Create temp SCE
    sce_small <- SingleCellExperiment(
      assays  = list(counts = counts_sparse),
      colData = chunk_meta_filtered
    )
    
    # ----------------------------------------------------------------------
    # Save/Append to Donor Shards
    # ----------------------------------------------------------------------
    donor_ids <- as.character(colData(sce_small)[[DONOR_ID_COL]])
    unique_chunk_donors <- unique(donor_ids)
    
    for (id in unique_chunk_donors) {
      safe_id <- gsub("[^[:alnum:]]", "_", id)
      donor_file <- file.path(SHARD_DIR, paste0("donor_", safe_id, ".qs"))
      
      donor_subset <- sce_small[, donor_ids == id]
      
      if (file.exists(donor_file)) {
        # Append to existing
        old_shard <- qs::qread(donor_file)
        # Ensure gene order matches (intersection safety)
        common_genes <- intersect(rownames(old_shard), rownames(donor_subset))
        new_shard <- cbind(old_shard[common_genes, ], donor_subset[common_genes, ])
        qs::qsave(new_shard, donor_file)
      } else {
        # Create new
        qs::qsave(donor_subset, donor_file)
      }
    }
    
    # Cleanup & Progress
    update_progress(f_name, i)
    rm(sce_small, counts_sparse, chunk_slice, chunk_meta_filtered); gc()
  }
  
  adata$file$close()
  rm(full_meta_df); gc()
}

# ==============================================================================
# 4. FINAL MANIFEST GENERATION
# ==============================================================================
message("\n===> Building Global Master Manifest...")

all_shards <- list.files(SHARD_DIR, pattern = "^donor_.*\\.qs$", full.names = TRUE)

if (length(all_shards) == 0) {
  warning("No cells were saved. Please check your filters (Cell Type / Age / Min Cells).")
} else {
  
  manifest_df <- do.call(rbind, lapply(all_shards, function(f) {
    tmp <- qs::qread(f)
    data.frame(
      donor_id   = as.character(unique(colData(tmp)[[DONOR_ID_COL]])),
      cell_count = ncol(tmp),
      file_path  = f,
      stringsAsFactors = FALSE
    )
  }))
  
  # Create Feasibility Map (Which donors pass which cell count thresholds?)
  # Note: Donors are already filtered by MIN_CELLS_THRESHOLD, but this maps higher steps.
  feasibility_map <- lapply(CELL_STEPS, function(n) {
    if(is.na(n)) return(manifest_df$donor_id)
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
      cell_filter = TARGET_CELL_TYPE,
      min_age     = MIN_AGE_THRESHOLD,
      min_cells   = MIN_CELLS_THRESHOLD
    )
  )
  
  qs::qsave(sharding_final, file.path(SHARD_DIR, "00_MASTER_MANIFEST.qs"))
  
  message("\n================================================================")
  message(" SHARDING & FILTERING COMPLETE ")
  message(sprintf(" Cell Type:           %s", TARGET_CELL_TYPE))
  message(sprintf(" Age Threshold:       >= %s", MIN_AGE_THRESHOLD))
  message(sprintf(" Cell Count Thresh:   >= %s", MIN_CELLS_THRESHOLD))
  message(sprintf(" Total Donors Saved:  %d", sharding_final$config$n_donors))
  message(sprintf(" Total Cells Saved:   %d", sharding_final$config$total_cells))
  message("================================================================")
}