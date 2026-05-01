#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Annotate Roussos shards with donor PMI (post-mortem interval)
# PURPOSE:
#   Reads PMI (minutes) for each PsychAD donor from
#   Datalink/raw_data/media-1.csv (column: PMI, key: DonorID padded
#   to 4 digits, e.g. "Donor_0001"), and writes that scalar into
#   colData(sce)$pmi for every cell of every Roussos shard at:
#     Datalink/Shards/{Astro,Micro,Oligo}/{Roussos,Roussos_Full,Roussos_DepthEq}/donor_*.qs
#
#   PMI is uniform per donor, so it is uniform per shard. We avoid re-running
#   the slow Run_Shard_Roussos_Full.R sharding by mutating shard colData in
#   place (atomic write via qs::qsave to f.tmp then rename).
#
# IDEMPOTENCY:
#   - A sentinel file is written at
#       Datalink/Shards/{CellType}/{Variant}/.pmi_annotated
#     containing the CSV sha256 + donor count + timestamp.
#   - On rerun, if the sentinel exists and the CSV sha256 still matches, the
#     entire variant directory is skipped.
#   - Per-shard, if colData$pmi already exists and matches the CSV value
#     within tolerance, the shard is skipped.
#
# USAGE:
#   Rscript R/sharding/Annotate_Roussos_PMI.R                   # All Astro/Micro/Oligo, all variants
#   Rscript R/sharding/Annotate_Roussos_PMI.R Astro             # Just Astro
#   Rscript R/sharding/Annotate_Roussos_PMI.R Astro,Oligo       # Subset
#   POWER_ANALYSIS_PMI_FORCE=1 Rscript ...                      # Force rerun even if sentinels match
#   POWER_ANALYSIS_PMI_ALLOW_MISSING=1 Rscript ...              # Don't abort on unmatched donors
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(digest)
})

# --- Resolve project root ---
if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
    Sys.getenv("POWER_ANALYSIS_ROOT")
  } else if (requireNamespace("here", quietly = TRUE)) {
    here::here()
  } else {
    getwd()
  }
}

source(file.path(PROJECT_ROOT, "R", "shared", "donor_id_utils.R"))

# --- Constants ---
SHARDS_BASE <- file.path(PROJECT_ROOT, "Datalink", "Shards")
PMI_CSV     <- file.path(PROJECT_ROOT, "Datalink", "raw_data", "media-1.csv")

DEFAULT_CELL_TYPES <- c("Astro", "Micro", "Oligo")
ROUSSOS_VARIANTS   <- c("Roussos", "Roussos_Full", "Roussos_DepthEq")
DONOR_ID_COL       <- "donor_id"
PMI_COL            <- "pmi"        # canonical lower-case to match DeJager
PMI_TOL            <- 1e-6         # tolerance for "already annotated"

FORCE          <- nzchar(Sys.getenv("POWER_ANALYSIS_PMI_FORCE")) &&
                  tolower(Sys.getenv("POWER_ANALYSIS_PMI_FORCE")) %in% c("1", "true", "yes")
ALLOW_MISSING  <- nzchar(Sys.getenv("POWER_ANALYSIS_PMI_ALLOW_MISSING")) &&
                  tolower(Sys.getenv("POWER_ANALYSIS_PMI_ALLOW_MISSING")) %in% c("1", "true", "yes")

# --- Parse CLI args ---
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

# ==============================================================================
# 1. Load CSV and build PMI lookup
# ==============================================================================

if (!file.exists(PMI_CSV)) {
  stop(sprintf("media-1.csv not found at %s", PMI_CSV))
}

csv_sha256 <- digest::digest(file = PMI_CSV, algo = "sha256")
pmi_df <- read.csv(PMI_CSV, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("DonorID", "PMI")
missing_cols <- setdiff(required, colnames(pmi_df))
if (length(missing_cols)) {
  stop(sprintf("CSV missing required columns: %s",
               paste(missing_cols, collapse = ", ")))
}

pmi_lookup <- setNames(suppressWarnings(as.numeric(pmi_df$PMI)), pmi_df$DonorID)

message(sprintf("Loaded PMI for %d donors from %s (sha256=%s...)",
                length(pmi_lookup), PMI_CSV, substr(csv_sha256, 1, 12)))

# ==============================================================================
# 2. Per-variant annotation
# ==============================================================================

annotate_variant <- function(cell_type, variant) {
  variant_dir <- file.path(SHARDS_BASE, cell_type, variant)
  sentinel    <- file.path(variant_dir, ".pmi_annotated")

  if (!dir.exists(variant_dir)) {
    message(sprintf("[%s/%s] SKIP — directory not present.", cell_type, variant))
    return(invisible(NULL))
  }

  if (file.exists(sentinel) && !FORCE) {
    sentinel_lines <- tryCatch(readLines(sentinel), error = function(e) character())
    sha_line <- grep("^csv_sha256:", sentinel_lines, value = TRUE)
    sha_recorded <- if (length(sha_line)) trimws(sub("^csv_sha256:", "", sha_line[1])) else ""
    if (identical(sha_recorded, csv_sha256)) {
      message(sprintf("[%s/%s] SKIP — sentinel matches CSV sha256.",
                      cell_type, variant))
      return(invisible(NULL))
    } else {
      message(sprintf("[%s/%s] Sentinel exists but CSV sha256 differs — re-annotating.",
                      cell_type, variant))
    }
  }

  shards <- list.files(variant_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  if (length(shards) == 0) {
    message(sprintf("[%s/%s] SKIP — no donor shards.", cell_type, variant))
    return(invisible(NULL))
  }

  message(sprintf("\n[%s/%s] Annotating %d shards", cell_type, variant, length(shards)))

  n_written  <- 0L
  n_skipped  <- 0L
  n_missing  <- 0L
  unmatched  <- character()

  for (i in seq_along(shards)) {
    shard_path <- shards[i]

    sce <- tryCatch(qs::qread(shard_path), error = function(e) {
      warning(sprintf("Failed to read %s: %s", shard_path, conditionMessage(e)))
      NULL
    })
    if (is.null(sce)) next

    cd <- SummarizedExperiment::colData(sce)
    donor_raw <- as.character(cd[[DONOR_ID_COL]][1])
    padded <- roussos_pad_donor(donor_raw)
    if (is.na(padded) || !padded %in% names(pmi_lookup)) {
      n_missing <- n_missing + 1L
      unmatched <- c(unmatched, donor_raw)
      next
    }
    pmi_val <- pmi_lookup[[padded]]
    if (is.na(pmi_val)) {
      n_missing <- n_missing + 1L
      unmatched <- c(unmatched, donor_raw)
      next
    }

    # Idempotency: skip if already annotated with matching value
    if (PMI_COL %in% colnames(cd)) {
      existing_vals <- suppressWarnings(as.numeric(as.character(cd[[PMI_COL]])))
      if (length(existing_vals) > 0 &&
          all(!is.na(existing_vals)) &&
          all(abs(existing_vals - pmi_val) <= PMI_TOL)) {
        n_skipped <- n_skipped + 1L
        next
      }
    }

    # Write PMI as a uniform scalar across all cells
    cd[[PMI_COL]] <- rep(pmi_val, nrow(cd))
    SummarizedExperiment::colData(sce) <- cd

    # Record provenance in metadata
    md <- S4Vectors::metadata(sce)
    md$pmi_units        <- "minutes"
    md$pmi_source       <- basename(PMI_CSV)
    md$pmi_csv_sha256   <- csv_sha256
    md$pmi_annotated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    S4Vectors::metadata(sce) <- md

    # Atomic write
    tmp_path <- paste0(shard_path, ".tmp")
    qs::qsave(sce, tmp_path)
    file.rename(tmp_path, shard_path)
    n_written <- n_written + 1L

    if (i %% 200 == 0 || i == length(shards)) {
      message(sprintf("  [%s/%s] processed %d / %d (written=%d skipped=%d missing=%d)",
                      cell_type, variant, i, length(shards),
                      n_written, n_skipped, n_missing))
    }
  }

  if (n_missing > 0) {
    msg <- sprintf("[%s/%s] %d shards had no/NA PMI in CSV. First 10: %s",
                   cell_type, variant, n_missing,
                   paste(utils::head(unique(unmatched), 10), collapse = ", "))
    if (ALLOW_MISSING) {
      warning(msg)
    } else {
      stop(msg, "\nSet POWER_ANALYSIS_PMI_ALLOW_MISSING=1 to proceed despite missing PMI.")
    }
  }

  # Write sentinel
  writeLines(c(
    sprintf("annotated_at: %s",   format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    sprintf("csv_sha256: %s",     csv_sha256),
    sprintf("csv_path: %s",       PMI_CSV),
    sprintf("n_written: %d",      n_written),
    sprintf("n_skipped: %d",      n_skipped),
    sprintf("n_missing: %d",      n_missing),
    sprintf("n_total_shards: %d", length(shards))
  ), sentinel)

  message(sprintf("[%s/%s] DONE — written=%d skipped=%d missing=%d (sentinel: %s)",
                  cell_type, variant, n_written, n_skipped, n_missing, sentinel))
  invisible(list(written = n_written, skipped = n_skipped, missing = n_missing))
}

# ==============================================================================
# 3. Drive every cell-type x variant combination
# ==============================================================================

message("================================================================")
message(" ROUSSOS PMI ANNOTATION ")
message(sprintf(" Cell types: %s", paste(CELL_TYPES, collapse = ", ")))
message(sprintf(" Variants:   %s", paste(ROUSSOS_VARIANTS, collapse = ", ")))
message(sprintf(" Force:      %s", FORCE))
message("================================================================")

for (ct in CELL_TYPES) {
  for (variant in ROUSSOS_VARIANTS) {
    annotate_variant(ct, variant)
  }
}

message("\nALL DONE.")
