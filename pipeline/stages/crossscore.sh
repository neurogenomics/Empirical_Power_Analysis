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
  # 10_CrossScore_Methods.R skips (with a log line) any run whose truth is
  # missing; a reproduction must cover every run and both modes, so check.
  epa_exec Rscript -e '
    x <- data.table::fread("Datalink/Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv")
    rc <- intersect(c("run", "run_id", "Run"), names(x))[1]
    mc <- intersect(c("mode", "Mode"), names(x))[1]
    if (is.na(rc) || is.na(mc)) stop("cross_score_matrix.csv: no run/mode column")
    n <- as.integer(Sys.getenv("N_ROBUSTNESS_RUNS"))
    if (!setequal(unique(x[[rc]]), seq_len(n))) stop("cross-score does not cover runs 1..", n)
    if (!setequal(unique(x[[mc]]), c("cells", "reads"))) stop("cross-score does not cover both modes")
    cat("cross-score covers runs 1 ..", n, "and both modes,", nrow(x), "rows\n")'
  epa_out crossscore Cartesian_Downsampling/Cross_Score
  epa_complete
}
