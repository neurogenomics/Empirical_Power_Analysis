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
# SCRIPT: SCE Multi-CellType Shard Store Creator (Standardized 'donor_id')
# INPUT: List of .qs files (SingleCellExperiment)
# ACTION: Loads -> Parses Age -> Filters -> Standardizes 'donor_id' -> Shards
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(dplyr)
  library(Matrix)
})

# ==============================================================================
# HELPER: Interval & Limit Age Parser
# ==============================================================================
parse_verbose_age <- function(x) {
  s <- tolower(trimws(as.character(x)))
  y <- rep(NA_real_, length(s))
  
  # 1. INTERVALS ("65 to 77" -> 71)
  interval_pattern <- "(\\d+)\\s*(?:to|-)\\s*(\\d+)"
  is_interval <- grepl(interval_pattern, s)
  
  if (any(is_interval)) {
    matches <- regmatches(s, regexec(interval_pattern, s))
    y[is_interval] <- sapply(matches[is_interval], function(m) {
      (as.numeric(m[2]) + as.numeric(m[3])) / 2
    })
  }
  
  # 2. LIMITS ("Less than 65" -> 65, "90+" -> 90)
  mask_single <- !is_interval & !grepl("^(unknown|not available|na|n/a|nan)$", s)
  if (any(mask_single)) {
    y[mask_single] <- suppressWarnings(as.numeric(sub("^(.*?)(\\d+(?:\\.\\d+)?)(.*)$", "\\2", s[mask_single])))
  }
  return(y)
}

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

INPUT_FILES <- list(
  "Micro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "allen_microglia.qs"),
  "Astro" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "allen_astrocytes.qs"),
  "Oligo" = file.path(PROJECT_ROOT, "Datalink", "raw_data", "allen_oligo.qs")
)

BASE_OUTPUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")

# Candidates for detection (will be renamed to 'donor_id' in output)
CANDIDATE_DONOR_COLS <- c("Donor.ID", "donor_name", "donor_id", "External.ID", "individualID", "Subject", "projid")
CANDIDATE_AGE_COLS   <- c("Age.at.death","age", "Age", "donor_age", "years", "Age_years")

MIN_CELLS_PER_DONOR <- 30
MIN_AGE_YEARS       <- 40

CELL_STEPS    <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA) 

# ==============================================================================
# 2. CORE PROCESSING LOOP
# ==============================================================================

for (dataset_key in names(INPUT_FILES)) {
  
  input_path <- INPUT_FILES[[dataset_key]]
  type_shard_dir <- file.path(BASE_OUTPUT_DIR, dataset_key, "Allen")
  
  message("\n################################################################")
  message(sprintf(" PROCESSING DATASET: %s", dataset_key))
  message("################################################################")
  
  if (!file.exists(input_path)) {
    message(sprintf("    [WARNING] Input file not found: %s. Skipping...", input_path))
    next
  }
  
  if(!dir.exists(type_shard_dir)) dir.create(type_shard_dir, recursive = TRUE)
  
  message(">>> Loading SingleCellExperiment object...")
  sce_full <- qs::qread(input_path)
  
  # --- DETECT COLUMNS ---
  file_cols <- names(colData(sce_full))
  
  donor_match <- intersect(CANDIDATE_DONOR_COLS, file_cols)
  if (length(donor_match) == 0) { message("    [ERROR] No valid donor column found!"); next }
  ACTIVE_DONOR_COL <- donor_match[1]
  
  age_match <- intersect(CANDIDATE_AGE_COLS, file_cols)
  if (length(age_match) == 0) { message("    [ERROR] No valid Age column found!"); next }
  ACTIVE_AGE_COL <- age_match[1]
  
  # --- PARSE AGE ---
  message(">>> Parsing Ages...")
  raw_ages <- colData(sce_full)[[ACTIVE_AGE_COL]]
  sce_full$age <- parse_verbose_age(raw_ages) # Standardize Age Column
  
  # --- FILTERING ---
  message(">>> Applying Filters...")
  initial_cells <- ncol(sce_full)
  
  # 1. Filter Age
  valid_age_idx <- which(!is.na(sce_full$age) & sce_full$age >= MIN_AGE_YEARS)
  if (length(valid_age_idx) == 0) { message("    [SKIP] No cells remaining after Age filtering."); next }
  sce_full <- sce_full[, valid_age_idx]
  
  # 2. Filter Cell Count
  donor_ids <- as.character(colData(sce_full)[[ACTIVE_DONOR_COL]])
  donor_counts <- table(donor_ids)
  valid_donors <- names(donor_counts[donor_counts >= MIN_CELLS_PER_DONOR])
  
  if (length(valid_donors) == 0) { message("    [SKIP] No valid donors (>= 30 cells)."); next }
  sce_full <- sce_full[, donor_ids %in% valid_donors]
  
  message(sprintf("    [FILTER] Final: %d cells from %d donors.", ncol(sce_full), length(valid_donors)))
  
  # --- SHARDING ---
  final_donor_ids <- as.character(colData(sce_full)[[ACTIVE_DONOR_COL]])
  unique_donors <- unique(final_donor_ids)
  
  message(sprintf(">>> Sharding %d donors...", length(unique_donors)))
  
  for (id in unique_donors) {
    safe_id <- gsub("[^[:alnum:]]", "_", id)
    donor_file <- file.path(type_shard_dir, paste0("donor_", safe_id, ".qs"))
    
    idx <- which(final_donor_ids == id)
    donor_sce <- sce_full[, idx]
    
    # *** STANDARDIZATION STEP ***
    # Ensure the output object has a column named 'donor_id'
    colData(donor_sce)$donor_id <- as.character(colData(donor_sce)[[ACTIVE_DONOR_COL]])
    
    qs::qsave(donor_sce, donor_file)
  }
  
  # --- MANIFEST ---
  message("    Building Manifest...")
  all_shards <- list.files(type_shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  
  if (length(all_shards) > 0) {
    manifest_df <- do.call(rbind, lapply(all_shards, function(f) {
      tmp <- qs::qread(f)
      data.frame(
        # READ FROM STANDARDIZED COLUMN
        donor_id   = as.character(unique(colData(tmp)$donor_id)),
        cell_count = ncol(tmp),
        file_path  = f,
        stringsAsFactors = FALSE
      )
    }))
    
    feasibility_map <- lapply(CELL_STEPS, function(n) {
      if(is.na(n)) return(manifest_df$donor_id)
      manifest_df %>% filter(cell_count >= n) %>% pull(donor_id)
    })
    names(feasibility_map) <- paste0("step_", ifelse(is.na(CELL_STEPS), "Max", CELL_STEPS))
    
    sharding_final <- list(
      manifest    = manifest_df,
      feasibility = feasibility_map,
      config      = list(
        dataset            = dataset_key,
        n_donors           = nrow(manifest_df),
        total_cells        = sum(manifest_df$cell_count),
        min_age_cutoff     = MIN_AGE_YEARS,
        min_cell_cutoff    = MIN_CELLS_PER_DONOR
      )
    )
    qs::qsave(sharding_final, file.path(type_shard_dir, "00_MASTER_MANIFEST.qs"))
  }
  
  message("    [DONE]")
  rm(sce_full); gc()
}