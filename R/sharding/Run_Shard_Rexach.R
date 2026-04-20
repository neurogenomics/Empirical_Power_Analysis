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
# SCRIPT: Single H5AD to Multi-Folder Shard Creator (Filtered & Robust)
# INPUT: One large .h5ad file containing all cell types
# ACTION: Reads metadata -> Filters (Age/Counts) -> Shards by Donor
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
options(future.globals.maxSize = 100 * 1024^3) # 100GB RAM Limit

# ==============================================================================
# HELPER: H5AD Reader (Unchanged)
# ==============================================================================
read_h5ad_obs <- function(h5_file) {
  message("    Scanning H5AD structure...")
  listing <- h5ls(h5_file)
  
  obs_vars <- subset(listing, (group == "/obs" | group == "obs") & 
                       name != "__categories")
  
  if (nrow(obs_vars) == 0) stop("No columns found in /obs.")
  message(sprintf("    Found %d potential columns. Reading...", nrow(obs_vars)))
  
  raw_list <- list()
  lengths_vec <- c()
  
  for(col in obs_vars$name) {
    path <- paste0("/obs/", col)
    val <- tryCatch({ h5read(h5_file, path) }, error = function(e) NULL)
    if (is.null(val)) next
    
    if (is.list(val) && !is.data.frame(val) && "codes" %in% names(val)) {
      raw_codes <- val$codes
      if ("categories" %in% names(val)) {
        lookup <- val$categories
        val <- lookup[raw_codes + 1] 
      } else {
        val <- raw_codes
      }
    }
    
    cat_path <- paste0("/obs/__categories/", col)
    has_cat <- any(listing$group == "/obs/__categories" & listing$name == col)
    if (has_cat && !is.list(val)) { 
      lookup <- tryCatch({ h5read(h5_file, cat_path) }, error=function(e) NULL)
      if (!is.null(lookup)) val <- lookup[val + 1]
    }
    
    raw_list[[col]] <- val
    lengths_vec <- c(lengths_vec, length(val))
  }
  
  true_n_cells <- max(lengths_vec)
  message(sprintf("    Detected Cell Count: %d", true_n_cells))
  
  final_cols <- list()
  for (col_name in names(raw_list)) {
    if (length(raw_list[[col_name]]) == true_n_cells) {
      final_cols[[col_name]] <- raw_list[[col_name]]
    }
  }
  
  meta_df <- as.data.frame(final_cols, stringsAsFactors = FALSE)
  
  if ("_index" %in% names(meta_df)) {
    if (anyDuplicated(meta_df$`_index`) == 0) {
      rownames(meta_df) <- meta_df$`_index`
      meta_df$`_index` <- NULL
    } else {
      rownames(meta_df) <- paste0("cell_", 1:nrow(meta_df))
    }
  } else {
    rownames(meta_df) <- paste0("cell_", 1:nrow(meta_df))
  }
  return(meta_df)
}

read_h5ad_matrix_subset <- function(h5_file, cell_indices_1based, feature_names) {
  message("    Reading full sparse structure (Fast Scan)...")
  x_ls <- h5ls(h5_file)
  
  # CRITICAL FIX: Check for raw counts first (/raw/X), fall back to /X if needed
  has_raw <- any(grepl("^/raw$", x_ls$group))
  
  if (has_raw) {
    message("    [INFO] Found /raw folder - using raw counts from /raw/X")
    matrix_path <- "/raw/X"
    is_sparse <- any(subset(x_ls, group == "/raw/X")$name == "data")
  } else {
    message("    [WARN] No /raw folder found - using /X (may be normalized)")
    matrix_path <- "/X"
    is_sparse <- any(subset(x_ls, group == "/X")$name == "data")
  }
  
  if (!is_sparse) stop("This script currently supports Sparse (CSR/CSC) H5AD files. Found Dense.")
  
  # Read from the appropriate path
  data <- h5read(h5_file, paste0(matrix_path, "/data"))
  indices <- h5read(h5_file, paste0(matrix_path, "/indices"))
  indptr <- h5read(h5_file, paste0(matrix_path, "/indptr"))
  
  n_cells_total <- length(indptr) - 1
  n_features <- length(feature_names)
  
  message("    Constructing Matrix...")
  mat_t <- Matrix::sparseMatrix(
    p = indptr, j = indices, x = as.numeric(data),
    dims = c(n_cells_total, n_features), index1 = FALSE
  )
  
  mat <- Matrix::t(mat_t)
  rownames(mat) <- feature_names
  
  message(sprintf("    Subsetting for %d target cells...", length(cell_indices_1based)))
  mat_subset <- mat[, cell_indices_1based, drop=FALSE]
  return(mat_subset)
}

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

INPUT_H5AD <- file.path(PROJECT_ROOT, "Datalink", "raw_data", "Rexach.h5ad")

CELL_TYPE_MAPPING <- list(
  "microglial cell" = "Micro",
  "oligodendrocyte" = "Oligo",
  "astrocyte"       = "Astro"
)

# Column Setup
CELL_TYPE_COL <- "cell_type"      
DONOR_ID_COL  <- "donor_id"   
AGE_COL       <- "development_stage"        # <--- NEW: Define Age Column

# Filter Thresholds
MIN_AGE             <- 40     # <--- NEW: Min Age
MIN_CELLS_PER_DONOR <- 30     # <--- NEW: Min Cells per Donor

BASE_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")
CELL_STEPS      <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA) 

# ==============================================================================
# 2. PRE-PROCESS: GLOBAL METADATA LOAD
# ==============================================================================

message(sprintf(">>> Opening H5AD: %s", INPUT_H5AD))
if (!file.exists(INPUT_H5AD)) stop("Input file not found.")

message(">>> Reading Global Metadata (obs)...")
global_meta <- read_h5ad_obs(INPUT_H5AD)
available_cols <- names(global_meta)

# --- 2a. COLUMN VALIDATION & AUTO-CORRECTION ---

# Check Cell Type
if (!CELL_TYPE_COL %in% available_cols) {
  match <- available_cols[tolower(available_cols) == tolower(CELL_TYPE_COL)][1]
  if (!is.na(match)) {
    message(sprintf("!!! Auto-Correction: Using '%s' instead of '%s'", match, CELL_TYPE_COL))
    CELL_TYPE_COL <- match
  } else stop(paste("Cell type column", CELL_TYPE_COL, "not found."))
}

# Check Donor ID
if (!DONOR_ID_COL %in% available_cols) {
  match <- available_cols[tolower(available_cols) == tolower(DONOR_ID_COL)][1]
  if (!is.na(match)) {
    message(sprintf("!!! Auto-Correction: Using '%s' instead of '%s'", match, DONOR_ID_COL))
    DONOR_ID_COL <- match
  } else stop(paste("Donor column", DONOR_ID_COL, "not found."))
}

# Check Age
if (!AGE_COL %in% available_cols) {
  match <- available_cols[tolower(available_cols) == tolower(AGE_COL)][1]
  if (!is.na(match)) {
    message(sprintf("!!! Auto-Correction: Using '%s' instead of '%s'", match, AGE_COL))
    AGE_COL <- match
  } else stop(paste("Age column", AGE_COL, "not found."))
}
# ... (Previous parts of script remain unchanged)

# ==============================================================================
# 2b. CUSTOM AGE PARSING (Rexach Specific)
# ==============================================================================
message(sprintf(">>> Parsing Age Column: '%s'", AGE_COL))

# Function to parse specific Rexach formats
parse_rexach_age <- function(x) {
  # 1. Infinite Interval: "80 year-old and over stage" -> Return 80
  if (grepl("over", x)) {
    return(as.numeric(sub("^(\\d+).*", "\\1", x)))
  }
  
  # 2. Interval: "60-79 year-old stage" -> Return Mean (69.5)
  # Regex looks for: Digit(s) - Digit(s)
  if (grepl("^(\\d+)-(\\d+)", x)) {
    parts <- regmatches(x, regexec("^(\\d+)-(\\d+)", x))[[1]]
    low   <- as.numeric(parts[2])
    high  <- as.numeric(parts[3])
    return((low + high) / 2)
  }
  
  # 3. Single Age: "56-year-old stage" -> Return 56
  # Extract first sequence of digits
  if (grepl("\\d+", x)) {
    return(as.numeric(sub("^(\\d+).*", "\\1", x)))
  }
  
  return(NA) # Fallback
}

# Apply the parser safely
raw_ages <- global_meta[[AGE_COL]]

# We convert to character first to be safe, then apply the function
parsed_ages <- sapply(as.character(raw_ages), parse_rexach_age)

# Check for failures
n_na <- sum(is.na(parsed_ages))
if (n_na > 0) {
  # Optional: Print a few examples of what failed to help debug
  failed_examples <- head(raw_ages[is.na(parsed_ages)], 5)
  warning(sprintf("    WARNING: %d cells could not be parsed for age. Examples: %s", 
                  n_na, paste(failed_examples, collapse=", ")))
  # We set them to NA; they will be filtered out by the >= MIN_AGE check later
}

# Overwrite the column with the numeric values
global_meta[[AGE_COL]] <- parsed_ages

message("    Age parsing complete. Range: ", 
        paste(range(parsed_ages, na.rm=TRUE), collapse=" - "))

# Get Feature Names
message(">>> Reading Feature Names (var)...")
tryCatch({
  global_features <- h5read(INPUT_H5AD, "var/_index")
}, error = function(e) {
  global_features <- h5read(INPUT_H5AD, "var")[[1]]
})

# ==============================================================================
# 3. CORE PROCESSING LOOP
# ==============================================================================

for (target_name_in_file in names(CELL_TYPE_MAPPING)) {
  
  target_folder_name <- CELL_TYPE_MAPPING[[target_name_in_file]]
  safe_type_name <- gsub("[^[:alnum:]]", "_", target_folder_name)
  type_shard_dir <- file.path(BASE_OUTPUT_DIR, safe_type_name, "Rexach")
  
  message(sprintf("\n>>> Processing Mapping: '%s' -> '%s'", target_name_in_file, target_folder_name))
  if(!dir.exists(type_shard_dir)) dir.create(type_shard_dir, recursive = TRUE)
  
  # --- STEP A: SUBSET & FILTER ---
  
  # 1. Base Cell Type Filter
  raw_indices <- which(global_meta[[CELL_TYPE_COL]] == target_name_in_file)
  if (length(raw_indices) == 0) {
    message(sprintf("    [SKIP] No cells found for type '%s'.", target_name_in_file))
    next
  }
  
  # 2. Apply Age Filter
  age_valid_indices <- which(global_meta[[AGE_COL]][raw_indices] >= MIN_AGE)
  filtered_indices <- raw_indices[age_valid_indices]
  
  removed_age <- length(raw_indices) - length(filtered_indices)
  if (length(filtered_indices) == 0) {
    message(sprintf("    [SKIP] All cells removed by age filter (< %d).", MIN_AGE))
    next
  }
  
  # 3. Apply Donor Count Filter
  # We must count donors *after* age filtering
  temp_meta <- global_meta[filtered_indices, ]
  donor_counts <- table(temp_meta[[DONOR_ID_COL]])
  
  valid_donors <- names(donor_counts)[donor_counts >= MIN_CELLS_PER_DONOR]
  final_indices <- filtered_indices[ temp_meta[[DONOR_ID_COL]] %in% valid_donors ]
  
  removed_counts <- length(filtered_indices) - length(final_indices)
  
  # Summary Report
  message(sprintf("    Filtering Summary for %s:", target_name_in_file))
  message(sprintf("      - Initial Cells:     %d", length(raw_indices)))
  message(sprintf("      - Removed (Age <%d): %d", MIN_AGE, removed_age))
  message(sprintf("      - Removed (<%d cells/donor): %d", MIN_CELLS_PER_DONOR, removed_counts))
  message(sprintf("      - FINAL COUNT:       %d", length(final_indices)))
  
  if (length(final_indices) == 0) {
    message("    [SKIP] No cells remaining after all filters.")
    next
  }
  
  # Update Metadata for Processing
  type_meta <- global_meta[final_indices, ]
  unique_donors <- unique(type_meta[[DONOR_ID_COL]])
  
  # --- STEP B: DETERMINE MISSING SHARDS ---
  donor_map <- data.frame(
    id = unique_donors,
    safe_id = gsub("[^[:alnum:]]", "_", unique_donors),
    stringsAsFactors = FALSE
  )
  donor_map$filename <- paste0("donor_", donor_map$safe_id, ".qs")
  existing_files <- list.files(type_shard_dir, pattern = "^donor_.*\\.qs$")
  donors_to_process <- donor_map$id[ !donor_map$filename %in% existing_files ]
  
  message(sprintf("    Processing %d valid donors (%d missing shards).", 
                  length(unique_donors), length(donors_to_process)))
  
  # --- STEP C: PROCESS (Load Matrix & Split) ---
  
  if (length(donors_to_process) == 0) {
    message("    [SKIP] All donor shards already exist.")
  } else {
    message("    Loading Matrix Subset...")
    
    # Load Matrix using FINAL indices
    counts_subset <- read_h5ad_matrix_subset(INPUT_H5AD, final_indices, global_features)
    
    # Validate that counts are integers (raw counts, not normalized)
    # Note: dgCMatrix stores values as double, but we check if they represent integers
    if (!all(counts_subset@x == floor(counts_subset@x), na.rm = TRUE)) {
      non_int_rate <- mean(counts_subset@x != floor(counts_subset@x), na.rm = TRUE)
      warning(sprintf("    [WARNING] Non-integer values detected in count matrix (%.2f%%). Data may be normalized, not raw counts!", 
                      non_int_rate * 100))
    } else {
      message("    [OK] Count matrix validated - contains integer values (raw counts)")
    }
    # Note: dgCMatrix requires x slot to be double - integer conversion happens 
    # later during pseudobulk aggregation when converting to dense matrix
    
    sce_subset <- SingleCellExperiment(
      assays = list(counts = counts_subset),
      colData = type_meta
    )
    
    message(sprintf("    Writing %d missing shards...", length(donors_to_process)))
    donor_vec <- as.character(sce_subset[[DONOR_ID_COL]])
    
    future_lapply(donors_to_process, function(id) {
      safe_id <- gsub("[^[:alnum:]]", "_", id)
      donor_file <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
      
      donor_sce <- sce_subset[, donor_vec == id]
      qs::qsave(donor_sce, donor_file, nthreads = 1)
      return(TRUE)
    }, future.seed = TRUE)
    
    rm(counts_subset, sce_subset, donor_vec)
    gc()
  }
  
  # --- STEP D: GENERATE MANIFEST ---
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
    
    # Filter manifest to ensure it matches current stringent filters 
    # (in case old files exist)
    manifest_df <- manifest_df %>% filter(cell_count >= MIN_CELLS_PER_DONOR)
    
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
        cell_filter = target_folder_name,
        min_age     = MIN_AGE,
        min_cells   = MIN_CELLS_PER_DONOR
      )
    )
    qs::qsave(sharding_final, file.path(type_shard_dir, "00_MASTER_MANIFEST.qs"), nthreads = 4)
  }
  message("    [DONE]")
}

message("================================================================")
message(" ALL PROCESSING COMPLETE ")
message("================================================================")

plan(sequential)