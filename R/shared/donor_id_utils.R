# ==============================================================================
# Donor ID utilities — shared helpers used across pipelines
# ==============================================================================

#' Pad Roussos/PsychAD donor IDs to 4-digit zero-padded form
#'
#' Shards store `donor_id` like "Donor_1", "Donor_42", "Donor_1303".
#' The PsychAD supplementary table (media-1.csv) uses 4-digit zero-padded form
#' "Donor_0001", "Donor_0042", "Donor_1303".
#'
#' This helper converts shard-format IDs to CSV-format IDs so we can join
#' supplementary metadata (e.g. PMI) against shard colData. Use it everywhere
#' that joins between these two sources to avoid bespoke per-script padders.
#'
#' @param x Character vector of donor IDs (shard format, e.g. "Donor_1")
#' @return Character vector in CSV format (e.g. "Donor_0001"); NA where the
#'   input did not match the expected pattern.
roussos_pad_donor <- function(x) {
  m <- regmatches(x, regexec("^Donor_(\\d+)$", x))
  vapply(seq_along(x), function(i) {
    if (length(m[[i]]) < 2) return(NA_character_)
    sprintf("Donor_%04d", as.integer(m[[i]][2]))
  }, character(1))
}
