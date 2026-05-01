#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Audit Roussos PMI Join (READ-ONLY)
# PURPOSE:
#   Verifies that every Roussos donor shard can be joined to a PMI value in
#   media-1.csv before we mutate any shards. This is a pure read-only check;
#   it writes a single CSV report and never modifies shard files.
#
# CHECKS:
#   1. For every donor_*.qs shard in
#      Datalink/Shards/{Astro,Micro,Oligo}/{Roussos,Roussos_Full,Roussos_DepthEq},
#      load colData$donor_id, pad to CSV form via roussos_pad_donor(), and
#      look up PMI in media-1.csv.
#   2. Report:
#        - donors with NO match in CSV (missing entries)
#        - donors with NA / blank PMI in CSV
#        - donors whose padded ID could not be parsed
#        - CSV donors not represented in any shard (informational)
#
# OUTPUT:
#   Datalink/Shards/.pmi_audit/donor_id_audit.csv
#   Datalink/Shards/.pmi_audit/summary.txt
#
# USAGE:
#   Rscript R/sharding/Audit_Roussos_PMI.R
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(SingleCellExperiment)
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

SHARDS_BASE <- file.path(PROJECT_ROOT, "Datalink", "Shards")
PMI_CSV     <- file.path(PROJECT_ROOT, "Datalink", "raw_data", "media-1.csv")
AUDIT_DIR   <- file.path(SHARDS_BASE, ".pmi_audit")
AUDIT_CSV   <- file.path(AUDIT_DIR, "donor_id_audit.csv")
AUDIT_TXT   <- file.path(AUDIT_DIR, "summary.txt")

CELL_TYPES        <- c("Astro", "Micro", "Oligo")
ROUSSOS_VARIANTS  <- c("Roussos", "Roussos_Full", "Roussos_DepthEq")
DONOR_ID_COL      <- "donor_id"

# ==============================================================================
# 1. Load CSV
# ==============================================================================

if (!file.exists(PMI_CSV)) {
  stop(sprintf("media-1.csv not found at %s", PMI_CSV))
}

pmi_df <- read.csv(PMI_CSV, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("DonorID", "PMI")
missing_cols <- setdiff(required_cols, colnames(pmi_df))
if (length(missing_cols)) {
  stop(sprintf("CSV missing required columns: %s",
               paste(missing_cols, collapse = ", ")))
}

message(sprintf("Loaded %d rows from %s", nrow(pmi_df), PMI_CSV))

pmi_lookup <- setNames(suppressWarnings(as.numeric(pmi_df$PMI)), pmi_df$DonorID)

# ==============================================================================
# 2. Walk every Roussos shard variant and resolve PMI
# ==============================================================================

inspect_one <- function(shard_path, cell_type, variant) {
  out <- data.frame(
    cell_type     = cell_type,
    variant       = variant,
    shard_file    = basename(shard_path),
    shard_donor   = NA_character_,
    padded_donor  = NA_character_,
    in_csv        = FALSE,
    pmi_minutes   = NA_real_,
    pmi_is_na     = NA,
    issue         = NA_character_,
    stringsAsFactors = FALSE
  )
  sce <- tryCatch(qs::qread(shard_path), error = function(e) NULL)
  if (is.null(sce)) {
    out$issue <- "shard_unreadable"
    return(out)
  }
  cd <- SummarizedExperiment::colData(sce)
  if (!DONOR_ID_COL %in% colnames(cd)) {
    out$issue <- "no_donor_id_col"
    return(out)
  }
  donor_raw <- as.character(cd[[DONOR_ID_COL]][1])
  out$shard_donor <- donor_raw
  padded <- roussos_pad_donor(donor_raw)
  out$padded_donor <- padded
  if (is.na(padded)) {
    out$issue <- "donor_id_unparseable"
    return(out)
  }
  hit <- !is.na(pmi_lookup[padded])
  matched <- padded %in% names(pmi_lookup)
  out$in_csv <- matched
  if (!matched) {
    out$issue <- "no_csv_row"
    return(out)
  }
  pmi_val <- pmi_lookup[[padded]]
  out$pmi_minutes <- pmi_val
  out$pmi_is_na <- is.na(pmi_val)
  if (is.na(pmi_val)) out$issue <- "csv_pmi_blank"
  out
}

audit_rows <- list()
for (ct in CELL_TYPES) {
  for (variant in ROUSSOS_VARIANTS) {
    variant_dir <- file.path(SHARDS_BASE, ct, variant)
    if (!dir.exists(variant_dir)) next
    shards <- list.files(variant_dir, pattern = "^donor_.*\\.qs$",
                         full.names = TRUE)
    if (length(shards) == 0) next
    message(sprintf("[%s/%s] %d shards", ct, variant, length(shards)))
    for (i in seq_along(shards)) {
      audit_rows[[length(audit_rows) + 1L]] <- inspect_one(shards[i], ct, variant)
      if (i %% 200 == 0 || i == length(shards)) {
        message(sprintf("  inspected %d / %d", i, length(shards)))
      }
    }
  }
}

if (length(audit_rows) == 0) {
  stop("No Roussos shards found under any variant. Aborting audit.")
}

audit_df <- do.call(rbind, audit_rows)

# ==============================================================================
# 3. Write report
# ==============================================================================

dir.create(AUDIT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(audit_df, AUDIT_CSV, row.names = FALSE)

n_total       <- nrow(audit_df)
n_unparseable <- sum(audit_df$issue %in% "donor_id_unparseable")
n_no_csv_row  <- sum(audit_df$issue %in% "no_csv_row")
n_blank_pmi   <- sum(audit_df$issue %in% "csv_pmi_blank")
n_unreadable  <- sum(audit_df$issue %in% "shard_unreadable")
n_clean       <- sum(is.na(audit_df$issue))

# Donors in CSV but never seen in any shard (informational)
shard_padded_set <- unique(stats::na.omit(audit_df$padded_donor))
csv_only <- setdiff(pmi_df$DonorID, shard_padded_set)

summary_lines <- c(
  sprintf("PMI audit summary — %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("CSV: %s", PMI_CSV),
  sprintf("Shards base: %s", SHARDS_BASE),
  "",
  sprintf("Total shards inspected: %d", n_total),
  sprintf("  clean (PMI joined): %d", n_clean),
  sprintf("  unreadable shards: %d", n_unreadable),
  sprintf("  unparseable donor_id: %d", n_unparseable),
  sprintf("  no matching CSV row: %d", n_no_csv_row),
  sprintf("  CSV row present but PMI blank/NA: %d", n_blank_pmi),
  "",
  sprintf("CSV rows not represented by any shard: %d / %d",
          length(csv_only), nrow(pmi_df))
)
writeLines(summary_lines, AUDIT_TXT)

cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")

if (n_unparseable + n_no_csv_row + n_unreadable > 0) {
  cat("\nDonors with hard issues (unparseable / no_csv_row / unreadable):\n")
  bad <- audit_df[audit_df$issue %in%
                    c("donor_id_unparseable", "no_csv_row", "shard_unreadable"), ,
                  drop = FALSE]
  print(utils::head(bad, 20))
  if (nrow(bad) > 20) cat(sprintf("... %d more (see %s)\n", nrow(bad) - 20, AUDIT_CSV))
  cat(sprintf("\nFull audit written to: %s\n", AUDIT_CSV))
  quit(status = 1)
}

if (n_blank_pmi > 0) {
  cat(sprintf("\n%d donors have a CSV row but PMI is blank/NA — they will be mean-imputed during DE.\n",
              n_blank_pmi))
}

cat(sprintf("\nFull audit written to: %s\n", AUDIT_CSV))
invisible(NULL)
