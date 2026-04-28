# ==============================================================================
# Shared config (resolves PROJECT_ROOT via here::here() or POWER_ANALYSIS_ROOT env var).
# ==============================================================================
if (!exists("PROJECT_ROOT")) {
  .sharding_config <- tryCatch({
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

#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Roussos Shard Creator (FILTERED — Astro / Micro / Oligo)
# PURPOSE:
#   Two-step orchestrator that produces filtered donor shards for the
#   Roussos (PsychAD) dataset. Designed so that downstream pipelines
#   (cartesian_downsampling, depth_equalization, power_curve_fitting,
#    full_ground_truth) can read shards from the same path/format as
#   the historical layout: Datalink/Shards/{CellType}/Roussos/donor_*.qs
#
# STEP 1 — UNFILTERED SHARDING (per cell type)
#   For each of Astro / Micro / Oligo, invoke Run_Shard_Roussos_Full.R
#   as a subprocess. That script reads PsychAD-{HBCC,RADC,MSSM}.h5ad in
#   chunks via reticulate/anndata (backed mode) and writes one .qs per
#   donor to Datalink/Shards/{CellType}/Roussos_Full/. It maintains its
#   own per-chunk progress log, so partial runs resume cleanly.
#
# STEP 2 — DONOR-LEVEL FILTERING (per cell type)
#   For each cell type, walk the Roussos_Full shards, parse donor age
#   from colData$development_stage, and copy donors that satisfy
#       cells_in_shard >= MIN_CELLS_THRESHOLD  (default 30)
#       donor_age      >= MIN_AGE_THRESHOLD    (default 40 years)
#   into Datalink/Shards/{CellType}/Roussos/. Then build a master
#   manifest at Datalink/Shards/{CellType}/Roussos/00_MASTER_MANIFEST.qs
#   matching the format the cartesian_downsampling pipeline expects
#   (manifest$manifest data.frame + feasibility map + config block).
#
# RESUME / SKIP BEHAVIOUR:
#   - Step 1 is skipped for a cell type when its Roussos_Full manifest
#     exists. (To force a re-shard, delete the Roussos_Full directory.)
#   - Step 2 file copies are skipped per donor when the destination
#     .qs already exists and is newer than the source.
#   - The Step 2 manifest is always rebuilt so it stays consistent with
#     whatever shards are present in Roussos/.
#
# INPUT:  PsychAD-{HBCC,RADC,MSSM}.h5ad (via subprocess to _Full script)
# OUTPUT: Datalink/Shards/{CellType}/Roussos_Full/donor_*.qs   (unfiltered)
#         Datalink/Shards/{CellType}/Roussos/donor_*.qs         (filtered)
#         Datalink/Shards/{CellType}/Roussos/00_MASTER_MANIFEST.qs
#
# USAGE:
#   Rscript R/sharding/Run_Shard_Roussos.R                 # All three cell types
#   Rscript R/sharding/Run_Shard_Roussos.R Astro           # Just Astro
#   Rscript R/sharding/Run_Shard_Roussos.R Astro,Oligo     # Subset
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(dplyr)
})

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# Cell types this orchestrator manages
DEFAULT_CELL_TYPES <- c("Astro", "Micro", "Oligo")

# Parse CLI: optional comma-separated subset
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) {
  CELL_TYPES <- trimws(strsplit(args[1], ",", fixed = TRUE)[[1]])
  invalid <- setdiff(CELL_TYPES, DEFAULT_CELL_TYPES)
  if (length(invalid)) {
    stop(sprintf("Unknown cell type(s): %s. Valid: %s",
                 paste(invalid, collapse = ", "),
                 paste(DEFAULT_CELL_TYPES, collapse = ", ")))
  }
} else {
  CELL_TYPES <- DEFAULT_CELL_TYPES
}

# Filter thresholds (donor must satisfy BOTH)
MIN_AGE_THRESHOLD   <- 40   # years
MIN_CELLS_THRESHOLD <- 30   # cells per donor (within target cell type)

# Cell-count steps used by the feasibility map (mirrors historical script)
CELL_STEPS <- c(10, 15, 20, 30, 40, 50, 75, 100, 250, 500, 1000, NA)

# Metadata column names (must match what _Full writes into colData)
DONOR_ID_COL <- "donor_id"
AGE_COL      <- "development_stage"

# Paths
SHARDING_DIR  <- file.path(PROJECT_ROOT, "R", "sharding")
FULL_SCRIPT   <- file.path(SHARDING_DIR, "Run_Shard_Roussos_Full.R")
SHARDS_BASE   <- file.path(PROJECT_ROOT, "Datalink", "Shards")

if (!file.exists(FULL_SCRIPT)) {
  stop(sprintf("Run_Shard_Roussos_Full.R not found at: %s", FULL_SCRIPT))
}

message("\n================================================================")
message(" ROUSSOS SHARD CREATOR (Filtered)")
message("================================================================")
message(sprintf(" Cell types:       %s", paste(CELL_TYPES, collapse = ", ")))
message(sprintf(" Min age:          >= %s years", MIN_AGE_THRESHOLD))
message(sprintf(" Min cells/donor:  >= %s",       MIN_CELLS_THRESHOLD))
message(sprintf(" Shards base:      %s", SHARDS_BASE))
message("================================================================\n")

# ==============================================================================
# 2. AGE PARSING (mirrors helpers in Run_Shard_Roussos_Full's input metadata)
# ==============================================================================

is_unknown_string <- function(x) {
  s <- tolower(trimws(as.character(x)))
  grepl("^(unknown|not available|na|n/a|nan)$", s)
}

parse_age_stage_to_years <- function(x) {
  s <- tolower(trimws(as.character(x)))
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

  is_numeric_str <- !is.na(suppressWarnings(as.numeric(s))) & is.na(y)
  y[is_numeric_str] <- as.numeric(s[is_numeric_str])

  y
}

# ==============================================================================
# 3. STEP 1 — UNFILTERED SHARDING (subprocess per cell type)
# ==============================================================================

run_full_sharding <- function(cell_type_short) {
  full_dir      <- file.path(SHARDS_BASE, cell_type_short, "Roussos_Full")
  full_manifest <- file.path(full_dir, "00_MASTER_MANIFEST.qs")

  if (file.exists(full_manifest)) {
    message(sprintf("[STEP 1] [%s] SKIP — Roussos_Full manifest exists (%s).",
                    cell_type_short, full_manifest))
    return(invisible(TRUE))
  }

  message(sprintf("\n[STEP 1] [%s] Running unfiltered sharding via %s",
                  cell_type_short, basename(FULL_SCRIPT)))

  rscript_bin <- file.path(R.home("bin"), "Rscript")
  status <- system2(
    command = rscript_bin,
    args    = c(shQuote(FULL_SCRIPT), cell_type_short),
    stdout  = "",
    stderr  = ""
  )

  if (status != 0) {
    stop(sprintf("[STEP 1] [%s] _Full script exited with status %d",
                 cell_type_short, status))
  }

  if (!file.exists(full_manifest)) {
    stop(sprintf("[STEP 1] [%s] _Full script returned 0 but manifest is missing: %s",
                 cell_type_short, full_manifest))
  }

  message(sprintf("[STEP 1] [%s] DONE — %s", cell_type_short, full_manifest))
  invisible(TRUE)
}

# ==============================================================================
# 4. STEP 2 — DONOR-LEVEL FILTERING (Roussos_Full -> Roussos)
# ==============================================================================

inspect_shard <- function(shard_path) {
  # Load the SCE; we need colData to compute donor age + cell count
  sce <- qs::qread(shard_path)
  cd  <- as.data.frame(SummarizedExperiment::colData(sce))

  donor_id <- as.character(unique(cd[[DONOR_ID_COL]]))
  if (length(donor_id) != 1) {
    warning(sprintf("Shard %s has %d donor IDs; using first.",
                    basename(shard_path), length(donor_id)))
    donor_id <- donor_id[1]
  }

  raw_age   <- as.character(cd[[AGE_COL]])
  age_years <- max(parse_age_stage_to_years(raw_age), na.rm = TRUE)
  if (!is.finite(age_years)) age_years <- NA_real_

  list(
    donor_id   = donor_id,
    cell_count = ncol(sce),
    age_years  = age_years
  )
}

filter_cell_type <- function(cell_type_short) {
  full_dir <- file.path(SHARDS_BASE, cell_type_short, "Roussos_Full")
  filt_dir <- file.path(SHARDS_BASE, cell_type_short, "Roussos")
  sentinel <- file.path(filt_dir, ".filter_complete")

  if (file.exists(sentinel)) {
    message(sprintf("[STEP 2] [%s] SKIP — sentinel exists (%s). Delete it to force re-filter.",
                    cell_type_short, sentinel))
    return(invisible(TRUE))
  }

  if (!dir.exists(full_dir)) {
    stop(sprintf("[STEP 2] [%s] Roussos_Full directory missing: %s",
                 cell_type_short, full_dir))
  }
  if (!dir.exists(filt_dir)) dir.create(filt_dir, recursive = TRUE)

  full_shards <- list.files(full_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  if (length(full_shards) == 0) {
    stop(sprintf("[STEP 2] [%s] No donor shards found in %s",
                 cell_type_short, full_dir))
  }

  message(sprintf("\n[STEP 2] [%s] Inspecting %d unfiltered shards...",
                  cell_type_short, length(full_shards)))

  # Inspect each shard (this is the slow part — one qs::qread per donor)
  donor_meta <- vector("list", length(full_shards))
  for (i in seq_along(full_shards)) {
    if (i %% 50 == 0 || i == length(full_shards)) {
      message(sprintf("  [%s] inspected %d / %d", cell_type_short, i, length(full_shards)))
    }
    info <- tryCatch(inspect_shard(full_shards[i]), error = function(e) {
      warning(sprintf("Failed to inspect %s: %s", full_shards[i], conditionMessage(e)))
      NULL
    })
    if (is.null(info)) next
    info$source_path <- full_shards[i]
    donor_meta[[i]]  <- info
  }
  donor_meta <- Filter(Negate(is.null), donor_meta)

  donor_df <- do.call(rbind, lapply(donor_meta, function(x) {
    data.frame(
      donor_id    = x$donor_id,
      cell_count  = x$cell_count,
      age_years   = x$age_years,
      source_path = x$source_path,
      stringsAsFactors = FALSE
    )
  }))

  # Apply filters
  pass_cells <- donor_df$cell_count >= MIN_CELLS_THRESHOLD
  pass_age   <- !is.na(donor_df$age_years) & donor_df$age_years >= MIN_AGE_THRESHOLD
  keep_idx   <- pass_cells & pass_age
  pass_df    <- donor_df[keep_idx, , drop = FALSE]

  message(sprintf("[STEP 2] [%s] Total: %d  |  pass cells (>=%d): %d  |  pass age (>=%d): %d  |  pass both: %d",
                  cell_type_short,
                  nrow(donor_df),
                  MIN_CELLS_THRESHOLD, sum(pass_cells),
                  MIN_AGE_THRESHOLD,   sum(pass_age),
                  nrow(pass_df)))

  if (nrow(pass_df) == 0) {
    warning(sprintf("[STEP 2] [%s] No donors survived filtering — no shards copied.",
                    cell_type_short))
    return(invisible(NULL))
  }

  # Copy passing shards (skip if dest already up to date)
  dest_paths <- file.path(filt_dir, basename(pass_df$source_path))
  copied <- 0L
  skipped <- 0L

  for (j in seq_len(nrow(pass_df))) {
    src <- pass_df$source_path[j]
    dst <- dest_paths[j]
    if (file.exists(dst) &&
        file.info(dst)$mtime >= file.info(src)$mtime) {
      skipped <- skipped + 1L
    } else {
      ok <- file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
      if (!ok) stop(sprintf("[STEP 2] [%s] file.copy failed: %s -> %s",
                             cell_type_short, src, dst))
      copied <- copied + 1L
    }
  }

  message(sprintf("[STEP 2] [%s] Copied %d shards, skipped %d already up-to-date.",
                  cell_type_short, copied, skipped))

  # Remove any donor .qs in the filtered directory that no longer passes
  # (e.g. if filter thresholds were tightened on a re-run)
  existing_filtered <- list.files(filt_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  stale <- setdiff(existing_filtered, dest_paths)
  if (length(stale)) {
    message(sprintf("[STEP 2] [%s] Removing %d stale shards no longer passing filters.",
                    cell_type_short, length(stale)))
    file.remove(stale)
  }

  # Build manifest matching cartesian_downsampling expectations
  manifest_df <- data.frame(
    donor_id   = pass_df$donor_id,
    cell_count = pass_df$cell_count,
    file_path  = dest_paths,
    stringsAsFactors = FALSE
  )

  feasibility_map <- lapply(CELL_STEPS, function(n) {
    if (is.na(n)) return(manifest_df$donor_id)
    manifest_df$donor_id[manifest_df$cell_count >= n]
  })
  names(feasibility_map) <- paste0("step_", ifelse(is.na(CELL_STEPS), "Max", CELL_STEPS))

  sharding_final <- list(
    manifest    = manifest_df,
    feasibility = feasibility_map,
    config      = list(
      donor_col      = DONOR_ID_COL,
      n_donors       = nrow(manifest_df),
      total_cells    = sum(manifest_df$cell_count),
      cell_filter    = cell_type_short,
      min_age        = MIN_AGE_THRESHOLD,
      min_cells      = MIN_CELLS_THRESHOLD,
      filtered       = TRUE,
      source_dir     = full_dir,
      filtered_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  )

  manifest_path <- file.path(filt_dir, "00_MASTER_MANIFEST.qs")
  qs::qsave(sharding_final, manifest_path)

  # Write sentinel so future runs can skip the slow per-donor inspection.
  writeLines(
    c(sprintf("filtered_at: %s", sharding_final$config$filtered_at),
      sprintf("n_donors: %d",   sharding_final$config$n_donors),
      sprintf("total_cells: %d", sharding_final$config$total_cells),
      sprintf("min_age: %d",    MIN_AGE_THRESHOLD),
      sprintf("min_cells: %d",  MIN_CELLS_THRESHOLD)),
    sentinel
  )

  message(sprintf("[STEP 2] [%s] DONE — %d donors, %d cells -> %s",
                  cell_type_short,
                  sharding_final$config$n_donors,
                  sharding_final$config$total_cells,
                  manifest_path))

  invisible(sharding_final)
}

# ==============================================================================
# 5. ORCHESTRATION
# ==============================================================================

# STEP 1: ensure unfiltered shards exist for every requested cell type
for (ct in CELL_TYPES) {
  run_full_sharding(ct)
}

# STEP 2: filter and build manifest for every requested cell type
for (ct in CELL_TYPES) {
  filter_cell_type(ct)
}

message("\n================================================================")
message(" ALL ROUSSOS CELL TYPES PROCESSED ")
message(sprintf(" Cell types: %s", paste(CELL_TYPES, collapse = ", ")))
message("================================================================\n")
