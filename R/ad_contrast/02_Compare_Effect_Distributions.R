#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:  02_Compare_Effect_Distributions.R
# PURPOSE: Do AD vs control DEGs live in the same effect-size / expression
#          regime as the sex-biased DEGs the power model was trained on?
#          Compares the full-cohort DEGs of the two contrasts (from
#          01_Full_Cohort_DE.R): AD vs control (all chromosomes) against
#          AUTOSOMAL sex-biased DEGs only (sex-chromosome genes such as XIST
#          have extreme, mechanistically different effects and would distort
#          the comparison).
#
# DEG definition: padj < 0.05 and |log2FC| > 0.2 (the grid's truth-set rule).
# Every (gene, DE method) call counts once in the pooled figure; the by-method
# figure splits them.
#
# INPUT:   Datalink/AD_vs_Control/Full_Cohort_DE/{sex,ad_vs_control}_de_results.csv
# OUTPUT:  Datalink/AD_vs_Control/Effect_Distributions/
#            effect_distributions_pooled.{png,pdf}     |LFC| and mean expression, rows x cell type
#            effect_distributions_by_method.{png,pdf}
#            effect_distribution_summary.csv           counts, quantiles, LFC-bin shares, KS tests
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ROOT    <- Sys.getenv("POWER_ANALYSIS_ROOT", ".")
IN_DIR  <- file.path(ROOT, "Datalink", "AD_vs_Control", "Full_Cohort_DE")
OUT_DIR <- file.path(ROOT, "Datalink", "AD_vs_Control", "Effect_Distributions")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FDR <- 0.05
LFC_MIN <- 0.2
LFC_BREAKS <- c(0.2, 1, 2, 3, Inf)
LFC_LABELS <- c("Tiny (0.2-1)", "Small (1-2)", "Medium (2-3)", "Large (3+)")
GROUP_LEVELS <- c("Sex-biased (autosomal)", "AD vs control")
GROUP_COLOURS <- setNames(c("#0072B2", "#D55E00"), GROUP_LEVELS)
CT_LABELS <- c(Astro = "Astrocytes", Micro = "Microglia", Oligo = "Oligodendrocytes")
METHOD_LABELS <- c(deseq2 = "DESeq2", edger = "edgeR", limma = "limma-voom")

read_de <- function(contrast) {
  f <- file.path(IN_DIR, paste0(contrast, "_de_results.csv"))
  if (!file.exists(f)) stop("Missing ", f, " (run 01_Full_Cohort_DE.R)")
  read.csv(f, stringsAsFactors = FALSE)
}
sex <- read_de("sex") %>% filter(chr_class == "Autosomal") %>%
  mutate(group = GROUP_LEVELS[1])
ad  <- read_de("ad_vs_control") %>% mutate(group = GROUP_LEVELS[2])

degs <- bind_rows(sex, ad) %>%
  filter(!is.na(padj), padj < FDR, abs(log2FoldChange) > LFC_MIN,
         is.finite(base_mean), base_mean > 0) %>%
  mutate(abs_lfc   = abs(log2FoldChange),
         group     = factor(group, levels = GROUP_LEVELS),
         cell_type = factor(CT_LABELS[cell_type], levels = CT_LABELS),
         de_method = factor(METHOD_LABELS[de_method], levels = METHOD_LABELS),
         lfc_bin   = cut(abs_lfc, LFC_BREAKS, labels = LFC_LABELS, right = TRUE))

# ---- summary table --------------------------------------------------------------
ks_p <- function(d, var) {
  a <- d[[var]][d$group == GROUP_LEVELS[1]]; b <- d[[var]][d$group == GROUP_LEVELS[2]]
  if (length(a) < 2 || length(b) < 2) return(c(D = NA_real_, p = NA_real_))
  t <- suppressWarnings(stats::ks.test(a, b))
  c(D = unname(t$statistic), p = t$p.value)
}
summ <- degs %>%
  group_by(cell_type, de_method, group) %>%
  summarise(n_degs = n(),
            median_abs_lfc = median(abs_lfc),
            q90_abs_lfc = quantile(abs_lfc, 0.9),
            median_base_mean = median(base_mean),
            pct_tiny   = 100 * mean(lfc_bin == LFC_LABELS[1]),
            pct_small  = 100 * mean(lfc_bin == LFC_LABELS[2]),
            pct_medium = 100 * mean(lfc_bin == LFC_LABELS[3]),
            pct_large  = 100 * mean(lfc_bin == LFC_LABELS[4]),
            .groups = "drop")
ks <- degs %>%
  group_by(cell_type, de_method) %>%
  group_modify(~ {
    l <- ks_p(.x, "abs_lfc"); b <- ks_p(.x, "base_mean")
    data.frame(ks_D_abs_lfc = l[["D"]], ks_p_abs_lfc = l[["p"]],
               ks_D_base_mean = b[["D"]], ks_p_base_mean = b[["p"]])
  }) %>% ungroup()
summ <- left_join(summ, ks, by = c("cell_type", "de_method"))
write.csv(summ, file.path(OUT_DIR, "effect_distribution_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), row.names = FALSE)

# ---- figures ----------------------------------------------------------------------
# Histograms on log10 axes, each contrast normalised to its own total (share of
# that contrast's DEG calls), so the shapes compare despite very different counts.
long <- degs %>%
  select(cell_type, de_method, group, abs_lfc, base_mean) %>%
  pivot_longer(c(abs_lfc, base_mean), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("abs_lfc", "base_mean"),
                           labels = c("|log2 fold change|", "Gene counts (median raw count)")))

hist_plot <- function(d, facet_formula, ncol) {
  ggplot(d, aes(x = value, fill = group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, position = "identity",
                   alpha = 0.55, colour = NA) +
    scale_x_log10() +
    scale_fill_manual(values = GROUP_COLOURS, name = NULL) +
    facet_wrap(facet_formula, scales = "free", ncol = ncol) +
    labs(x = NULL, y = "Density (within contrast)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95"))
}

p_pooled <- hist_plot(long, ~ cell_type + variable, ncol = 2) +
  geom_vline(data = data.frame(variable = factor("|log2 fold change|", levels = levels(long$variable)),
                               x = c(1, 2)),
             aes(xintercept = x), linetype = "dashed", colour = "grey40", inherit.aes = FALSE)
ggsave(file.path(OUT_DIR, "effect_distributions_pooled.png"), p_pooled, width = 9, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "effect_distributions_pooled.pdf"), p_pooled, width = 9, height = 8)

p_method <- hist_plot(long, ~ cell_type + de_method + variable, ncol = 6)
ggsave(file.path(OUT_DIR, "effect_distributions_by_method.png"), p_method, width = 16, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "effect_distributions_by_method.pdf"), p_method, width = 16, height = 8)

cat("Done:", OUT_DIR, "\n")
