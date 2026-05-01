#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Genes-per-cell histogram for Roussos microglia shards
# PURPOSE: For every cell across every Roussos micro shard, count the number of
#          genes with non-zero UMI counts. Plot a histogram.
# OUTPUT:  Datalink/Exploratory_Analysis/genes_per_cell_micro_roussos.png
#          Datalink/Exploratory_Analysis/genes_per_cell_micro_roussos_summary.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(Matrix)
  library(ggplot2)
  library(future.apply)
  library(here)
})

PROJECT_ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT", unset = here::here())
SHARD_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards", "Micro", "Roussos")
OUT_DIR <- file.path(PROJECT_ROOT, "Datalink", "Exploratory_Analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

shard_files <- list.files(SHARD_DIR, pattern = "^donor_.*\\.qs$", full.names = TRUE)
message("Found ", length(shard_files), " micro Roussos shards")

n_workers <- min(16L, parallel::detectCores() - 2L)
plan(multisession, workers = n_workers)
message("Using ", n_workers, " workers")

genes_per_cell_one_shard <- function(f) {
  sce <- qs::qread(f)
  cnt <- counts(sce)
  if (!inherits(cnt, "dgCMatrix")) cnt <- as(cnt, "CsparseMatrix")
  # Number of non-zero entries per column == unique genes detected per cell
  diff(cnt@p)
}

t0 <- Sys.time()
gpc_list <- future_lapply(shard_files, genes_per_cell_one_shard, future.seed = TRUE)
gpc <- unlist(gpc_list, use.names = FALSE)
message("Computed genes-per-cell for ", length(gpc), " cells in ",
        format(round(difftime(Sys.time(), t0, units = "secs"), 1)))

summary_df <- data.frame(
  metric = c("n_cells", "n_donors", "mean", "median", "sd", "min", "p1", "p5",
             "p25", "p75", "p95", "p99", "max"),
  value = c(length(gpc), length(shard_files), mean(gpc), median(gpc), sd(gpc),
            min(gpc),
            quantile(gpc, 0.01), quantile(gpc, 0.05), quantile(gpc, 0.25),
            quantile(gpc, 0.75), quantile(gpc, 0.95), quantile(gpc, 0.99),
            max(gpc))
)
print(summary_df)
write.csv(summary_df,
          file.path(OUT_DIR, "genes_per_cell_micro_roussos_summary.csv"),
          row.names = FALSE)

df <- data.frame(genes_per_cell = gpc)

p <- ggplot(df, aes(x = genes_per_cell)) +
  geom_histogram(bins = 100, fill = "#3182bd", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = median(gpc), colour = "red", linetype = "dashed") +
  annotate("text",
           x = median(gpc), y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("median = %.0f", median(gpc)),
           colour = "red", size = 3.5) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Unique genes detected per cell -- Roussos microglia shards",
    subtitle = sprintf("%s cells across %d donors",
                       formatC(length(gpc), format = "d", big.mark = ","),
                       length(shard_files)),
    x = "Number of unique genes detected per cell (non-zero UMI count)",
    y = "Number of cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

out_png <- file.path(OUT_DIR, "genes_per_cell_micro_roussos.png")
ggsave(out_png, p, width = 9, height = 5.5, dpi = 200)
message("Wrote ", out_png)

# Also save the raw vector so Chris can re-bin without re-running
qs::qsave(gpc, file.path(OUT_DIR, "genes_per_cell_micro_roussos.qs"))
message("Done.")
