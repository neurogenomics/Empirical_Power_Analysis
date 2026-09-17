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
# SCRIPT: Johansen .h5ad.tar Sharder (With Age & Cell Count Filtering)
# INPUT: Compressed .h5ad.tar files
# ACTION: Untars -> Loads -> Age Filter -> QC -> Min Cell Filter -> Shards
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter) # For reading .h5ad
  library(qs)            # For fast saving/loading
  library(dplyr)
  library(Matrix)        # For sparse matrix conversion
  library(future)        
  library(future.apply) 
  library(tools)
})

# ==============================================================================
# HELPER: Age Parser
# ==============================================================================
parse_age_column <- function(age_vec) {
  age_vec <- as.character(age_vec)
  age_vec[age_vec == "unknown"] <- NA # Handle strict "unknown" string
  clean_nums <- gsub("[^0-9]", "", age_vec) # Remove " yrs" etc
  return(as.integer(clean_nums))
}

# ==============================================================================
# 0. PARALLELIZATION SETUP
# ==============================================================================
N_WORKERS <- 8
if (N_WORKERS < 1) N_WORKERS <- 1

is_rstudio <- Sys.getenv("RSTUDIO") == "1"

if (.Platform$OS.type == "windows" || is_rstudio) {
  message("Running in Windows/RStudio. Using 'multisession'.")
  plan(multisession, workers = N_WORKERS)
} else {
  message("Running in Terminal. Using 'multicore'.")
  plan(multicore, workers = N_WORKERS)
}

options(future.globals.maxSize = 100 * 1024^3)

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

# Input Files
INPUT_FILES <- list(
  "Astro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "Johansen_Astro.h5ad.tar"),
  "Micro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "Johansen_micro.h5ad.tar"),
  "Oligo" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "Johansen_oligo.h5ad.tar")
)

# Output Directory
BASE_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")

# The column used to split files (The Donor ID)
DONOR_ID_COL  <- "donor_name" 
AGE_COL       <- "age"

# Metadata columns you MUST have (Script will warn if missing)
REQUIRED_COLS <- c("age", "sex", "subclass_scANVI", 
                   "medical_conditions", "batch_vendor_name", "library_prep_set")

# Filtering Thresholds
MIN_CELLS_PER_DONOR <- 30
MIN_AGE_YEARS       <- 40

# Temp space for extraction
TEMP_EXTRACT_DIR <- file.path(PROJECT_ROOT, "Datalink", "raw_data", "temp_extraction")
CELL_STEPS    <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA) 

# ==============================================================================
# 2. CORE PROCESSING LOOP
# ==============================================================================

if(!dir.exists(TEMP_EXTRACT_DIR)) dir.create(TEMP_EXTRACT_DIR, recursive = TRUE)

for (target_type in names(INPUT_FILES)) {
  
  input_tar_path <- INPUT_FILES[[target_type]]
  safe_type_name <- gsub("[^[:alnum:]]", "_", target_type)
  type_shard_dir <- file.path(BASE_OUTPUT_DIR, safe_type_name, "Johansen")
  
  message(sprintf("\n>>> Processing Dataset: %s", target_type))
  
  if (!file.exists(input_tar_path)) {
    message(sprintf("    [ERROR] File not found: %s", input_tar_path))
    next
  }
  
  if(!dir.exists(type_shard_dir)) dir.create(type_shard_dir, recursive = TRUE)
  
  # --- STEP A: EXTRACT & LOAD ---
  message("    Extracting & Loading .h5ad...")
  
  specific_temp_dir <- file.path(TEMP_EXTRACT_DIR, safe_type_name)
  if(dir.exists(specific_temp_dir)) unlink(specific_temp_dir, recursive = TRUE)
  dir.create(specific_temp_dir)
  
  tryCatch({
    untar(input_tar_path, exdir = specific_temp_dir)
  }, error = function(e) { message("Untar failed"); NULL })
  
  h5ad_file <- list.files(specific_temp_dir, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE)[1]
  
  sce_full <- tryCatch({
    # use_hdf5 = FALSE loads data directly into memory (avoids DelayedArray issues)
    # Server has 208GB RAM, files are max 8.4GB, so this is safe
    zellkonverter::readH5AD(h5ad_file, use_hdf5 = FALSE)
  }, error = function(e) { 
    message(sprintf("    [ERROR] readH5AD failed: %s", e$message))
    NULL 
  })
  
  if(is.null(sce_full)) {
    message("    [ERROR] Failed to load H5AD.")
    unlink(specific_temp_dir, recursive = TRUE)
    next
  }
  
  # --- STEP B: FIX MATRIX (Load into RAM *BEFORE* Deleting File) ---
  # CRITICAL FIX: Use raw UMIs instead of scvi-normalized slot 1
  message("    Loading matrix into memory (dgCMatrix)...")
  
  # Check available assays
  avail_assays <- assayNames(sce_full)
  message(sprintf("    Available assays: %s", paste(avail_assays, collapse = ", ")))
  
  # Use UMIs layer if available (raw counts), otherwise warn
  if ("UMIs" %in% avail_assays) {
    message("    [INFO] Using 'UMIs' assay (raw counts)")
    # Convert UMIs to dgCMatrix and store as "counts"
    assay(sce_full, "counts") <- as(assay(sce_full, "UMIs"), "dgCMatrix")
  } else if ("counts" %in% avail_assays) {
    message("    [INFO] Using existing 'counts' assay")
    assay(sce_full, "counts") <- as(assay(sce_full, "counts"), "dgCMatrix")
  } else {
    warning("    [WARNING] No 'UMIs' or 'counts' assay found! Using slot 1 which may be normalized data.")
    assay(sce_full, "counts") <- as(assay(sce_full, 1), "dgCMatrix")
  }
  
  # Validate that counts are integers (raw counts, not normalized)
  # Note: dgCMatrix stores values as double, but we check if they represent integers
  counts_mat <- assay(sce_full, "counts")
  if (!all(counts_mat@x == floor(counts_mat@x), na.rm = TRUE)) {
    non_int_rate <- mean(counts_mat@x != floor(counts_mat@x), na.rm = TRUE)
    warning(sprintf("    [WARNING] Non-integer values detected in count matrix (%.2f%%). Data may be normalized, not raw counts!",
                    non_int_rate * 100))
  } else {
    message("    [OK] Count matrix validated - contains integer values (raw counts)")
  }
  # Note: dgCMatrix requires x slot to be double - integer conversion happens 
  # later during pseudobulk aggregation when converting to dense matrix
  rm(counts_mat)
  
  unlink(specific_temp_dir, recursive = TRUE) 
  gc() 
  
  # --- STEP C: AGE FILTER (>= 40) ---
  message(sprintf("    Applying Age Filter (Age >= %d)...", MIN_AGE_YEARS))
  
  if (AGE_COL %in% names(colData(sce_full))) {
    # 1. Parse Ages
    age_clean <- parse_age_column(sce_full[[AGE_COL]])
    sce_full$age_numeric <- age_clean 
    
    # 2. Identify INVALID Donors (NA or < 40)
    # We find IDs where age is missing OR too young
    invalid_mask <- is.na(age_clean) | (age_clean < MIN_AGE_YEARS)
    invalid_ids  <- unique(sce_full[[DONOR_ID_COL]][invalid_mask])
    invalid_ids  <- invalid_ids[!is.na(invalid_ids)]
    
    if (length(invalid_ids) > 0) {
      message(sprintf("    [CLEANUP] Found %d donors failing age check (<%d or unknown).", 
                      length(invalid_ids), MIN_AGE_YEARS))
      
      # 3. DELETE EXISTING FILES for these invalid donors
      files_removed <- 0
      for (bid in invalid_ids) {
        safe_id <- gsub("[^[:alnum:]]", "_", bid)
        f_path  <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
        if (file.exists(f_path)) {
          unlink(f_path)
          files_removed <- files_removed + 1
        }
      }
      if(files_removed > 0) message(sprintf("    [CLEANUP] Deleted %d existing shards for age-excluded donors.", files_removed))
      
      # 4. REMOVE FROM OBJECT
      sce_full <- sce_full[, !invalid_mask]
      message(sprintf("    [FILTER] Dropped %d cells (Age < %d or NA).", sum(invalid_mask), MIN_AGE_YEARS))
      
    } else {
      message("    [OK] All donors pass age filtering.")
    }
  } else {
    message(sprintf("    [WARNING] Age column '%s' not found. Skipping age filter.", AGE_COL))
  }
  
  if (ncol(sce_full) == 0) {
    message("    [STOP] No cells remaining after age filtering.")
    next
  }
  
  # --- STEP D: QC FILTERING ---
  message("    Applying QC Filters...")
  n_start <- ncol(sce_full)
  
  is_pass <- function(col_val) {
    if(is.null(col_val)) return(rep(TRUE, length(col_val)))
    as.character(col_val) == "True"
  }
  
  keep_cells <- is_pass(sce_full$QC2_metacell_class_flag) &
    is_pass(sce_full$QC4_metacell_gene_umi_low_flag) &
    is_pass(sce_full$QC4_metacell_gene_umi_high_flag) &
    is_pass(sce_full$QC5_class_prediction_flag)
  
  if("QC6_atleast_two_flag" %in% names(colData(sce_full))) {
    keep_cells <- keep_cells & (sce_full$QC6_atleast_two_flag == 0)
  }
  
  sce_full <- sce_full[, keep_cells]
  
  message(sprintf("    QC Filter: Dropped %d cells. Remaining: %d cells.", n_start - ncol(sce_full), ncol(sce_full)))
  
  if (ncol(sce_full) == 0) next
  
  # --- STEP E: CELL COUNT FILTER (>= 30) ---
  # Done AFTER QC because QC might drop a donor below threshold
  message(sprintf("    Applying Cell Count Filter (>= %d cells per donor)...", MIN_CELLS_PER_DONOR))
  
  current_donors <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  donor_counts   <- table(current_donors)
  
  # Identify donors with too few cells
  low_count_donors <- names(donor_counts[donor_counts < MIN_CELLS_PER_DONOR])
  
  if (length(low_count_donors) > 0) {
    message(sprintf("    [CLEANUP] Found %d donors with < %d cells.", length(low_count_donors), MIN_CELLS_PER_DONOR))
    
    # Delete existing shards for these donors
    files_removed <- 0
    for (bid in low_count_donors) {
      safe_id <- gsub("[^[:alnum:]]", "_", bid)
      f_path  <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
      if (file.exists(f_path)) {
        unlink(f_path)
        files_removed <- files_removed + 1
      }
    }
    if(files_removed > 0) message(sprintf("    [CLEANUP] Deleted %d existing shards for low-cell donors.", files_removed))
    
    # Filter object
    keep_donor_idx <- !(current_donors %in% low_count_donors)
    sce_full <- sce_full[, keep_donor_idx]
    
    message(sprintf("    [FILTER] Removed %d donors. Final valid donors: %d.", 
                    length(low_count_donors), length(unique(current_donors[keep_donor_idx]))))
  } else {
    message("    [OK] All donors have sufficient cells.")
  }
  
  if (ncol(sce_full) == 0) {
    message("    [STOP] No cells remaining after count filtering.")
    next
  }
  
  # --- STEP F: PROCESS REMAINING DONORS ---
  
  final_donors <- as.character(colData(sce_full)[[DONOR_ID_COL]])
  unique_donors <- unique(final_donors[!is.na(final_donors)])
  
  # Only process if file doesn't exist (Check again in case we deleted it above)
  donors_to_process <- Filter(function(id) {
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    f_path <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
    return(!file.exists(f_path))
  }, unique_donors)
  
  if (length(donors_to_process) == 0) {
    message("    [SKIP] All valid donor shards already exist.")
  } else {
    message(sprintf("    Processing %d missing donors in parallel...", length(donors_to_process)))
    
    future_lapply(donors_to_process, function(id) {
      safe_id <- gsub("[^[:alnum:]]", "_", id)
      donor_file <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
      
      idx <- which(colData(sce_full)[[DONOR_ID_COL]] == id)
      donor_sce <- sce_full[, idx]
      
      qs::qsave(donor_sce, donor_file, nthreads = 1)
      return(TRUE)
    }, future.seed = TRUE)
  }
  
  rm(sce_full); gc()
  
  # --- STEP G: GENERATE MANIFEST ---
  message("    Updating Manifest...")
  
  all_shards <- list.files(type_shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  
  if (length(all_shards) > 0) {
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
    
    feasibility_map <- lapply(CELL_STEPS, function(n) {
      if(is.na(n)) return(manifest_df$donor_id)
      manifest_df %>% filter(cell_count >= n) %>% pull(donor_id)
    })
    names(feasibility_map) <- paste0("step_", ifelse(is.na(CELL_STEPS), "Max", CELL_STEPS))
    
    sharding_final <- list(
      manifest    = manifest_df,
      feasibility = feasibility_map,
      config      = list(
        donor_col     = DONOR_ID_COL, 
        required_cols = REQUIRED_COLS,
        n_donors      = nrow(manifest_df),
        total_cells   = sum(manifest_df$cell_count),
        min_age       = MIN_AGE_YEARS,
        min_cells     = MIN_CELLS_PER_DONOR
      )
    )
    
    qs::qsave(sharding_final, file.path(type_shard_dir, "00_MASTER_MANIFEST.qs"), nthreads = 4)
  }
  message("    [DONE]")
}

message("\n================================================================")
message(" JOB COMPLETE ")
message("================================================================")
plan(sequential)