# ==============================================================================
# STEP 04: FDP Summary — From Full Ground Truth Diagnostics
# ==============================================================================
# Reads the existing FDP_Summary.csv from the full ground truth pipeline
# and copies it to Paper_Source_Data/ for canonical reference.
#
# Output: Paper_Source_Data/fdp_summary.csv
#
# This replaces the hardcoded fdp_summary data.frame in paper_final_draft.rmd.
# ==============================================================================

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 04: FDP SUMMARY ###\n")

# ==============================================================================
# Read existing FDP summary
# ==============================================================================

fdp_src <- file.path(GT_BASE, "FDP_Summary.csv")

if (!file.exists(fdp_src)) {
  stop("FDP_Summary.csv not found: ", fdp_src)
}

fdp <- read.csv(fdp_src)
message(sprintf("  Read %d rows from %s", nrow(fdp), fdp_src))

# Rename columns for consistency with paper
# Source has: CellType, Source, Group, FDP, N_Sig, N_In_Universe, N_Disagreed
# Paper expects: Cell_Type, Source, Group, FDP, N_Sig, N_Disagreed
if ("CellType" %in% colnames(fdp)) {
  colnames(fdp)[colnames(fdp) == "CellType"] <- "Cell_Type"
}

# Write canonical copy
out_path <- file.path(OUTPUT_BASE, "fdp_summary.csv")
write.csv(fdp, out_path, row.names = FALSE)

message(sprintf("\n[DONE] Wrote %d rows to %s", nrow(fdp), out_path))
message(sprintf("  Cell types: %s", paste(unique(fdp$Cell_Type), collapse = ", ")))
message(sprintf("  Groups: %s", paste(unique(fdp$Group), collapse = ", ")))
