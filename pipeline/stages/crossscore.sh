# shellcheck shell=bash
# Stage crossscore: score every method's downsampled calls against every
# method's ground truth (R/cartesian_downsampling/10_CrossScore_Methods.R), for
# both modes and all robustness runs. Single-threaded, roughly 1 h per run.
# Output: Cartesian_Downsampling/Cross_Score/.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_env.sh"

stage_main() {
  epa_init_ws
  epa_in_array grid
  epa_cartesian_env
  epa_env XSCORE_METHODS "$EPA_DE_METHODS"
  epa_env XSCORE_CELLTYPES "Astro,Micro,Oligo"
  epa_env XSCORE_MODES "cells,reads"
  epa_env XSCORE_RUNS "$(seq -s, 1 "$EPA_N_RUNS")"
  epa_step "cross-method scoring"
  epa_R R/cartesian_downsampling/10_CrossScore_Methods.R
  epa_require Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv
  # 10_CrossScore_Methods.R aggregates across runs (sums tp/gt counts) before
  # writing cross_score_matrix.csv, so the final table has no run column.
  # Per-run skips are logged by the R script; here we only check that both
  # modes are present in the output.
  epa_exec Rscript -e '
    x <- data.table::fread("Datalink/Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv")
    mc <- intersect(c("mode", "Mode"), names(x))[1]
    if (is.na(mc)) stop("cross_score_matrix.csv: no mode column")
    if (!setequal(unique(x[[mc]]), c("cells", "reads"))) stop("cross-score does not cover both modes")
    cat("cross-score covers both modes,", nrow(x), "rows\n")'
  epa_out crossscore Cartesian_Downsampling/Cross_Score
  epa_complete
}
