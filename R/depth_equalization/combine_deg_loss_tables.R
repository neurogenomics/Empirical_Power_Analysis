#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  combine_deg_loss_tables.R
# PURPOSE: Stack the per-cell-type outputs of deg_loss_analysis.R
#          (Datalink/Depth_Eq_DEG_Loss/<cell type>/) into one table per kind,
#          with a leading cell-type column, for Supplementary Table 3.
# OUTPUT:  Datalink/Depth_Eq_DEG_Loss/{deg_loss_table,deg_loss_counts,deg_loss_per_run}.csv
# USAGE:   Rscript R/depth_equalization/combine_deg_loss_tables.R
# ==============================================================================

root <- Sys.getenv("POWER_ANALYSIS_ROOT", ".")
out_dir <- file.path(root, "Datalink", "Depth_Eq_DEG_Loss")
cell_types <- c(Astro = "Astrocytes", Oligo = "Oligodendrocytes")   # reporting order

present <- cell_types[file.exists(file.path(out_dir, names(cell_types), "deg_loss_table.csv"))]
if (length(present) == 0) stop("No per-cell-type DEG-loss tables under ", out_dir)

for (kind in c("deg_loss_table", "deg_loss_counts", "deg_loss_per_run")) {
  parts <- lapply(names(present), function(ct) {
    d <- read.csv(file.path(out_dir, ct, paste0(kind, ".csv")),
                  stringsAsFactors = FALSE, check.names = FALSE)
    cbind(`Cell type` = unname(present[ct]), d, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (kind != "deg_loss_table") names(out)[1] <- "cell_type"
  write.csv(out, file.path(out_dir, paste0(kind, ".csv")), row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows; %s)\n", kind, nrow(out), paste(names(present), collapse = ", ")))
}
