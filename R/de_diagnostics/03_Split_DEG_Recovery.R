# ==============================================================================
# SCRIPT: Split DEG Recovery — Compare Split GT vs Full GT DE Results
# PURPOSE: Assess how representative each 50% split ground truth is by
#          computing bidirectional overlap (recall + precision) against the
#          full-cohort ground truth, stratified by LFC magnitude and
#          chromosome type. Generates bar charts, UpSet plots, and LFC
#          scatter plots.
# INPUT: Full GT DE results from Datalink/Full_Ground_Truth/
#        Split GT DE results from Datalink/Cartesian_Downsampling/.../Ground_Truth.qs
#        BioMart cache for chromosome annotation
# OUTPUT: Recovery bar charts, UpSet plots, LFC scatter plots, and summary CSVs
#         in Datalink/DE_Diagnostics/Split_Recovery/
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(UpSetR)
  library(scales)
})

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Load config - skip if already loaded (e.g., sourced from run_pipeline.R)
if (!exists("SCRIPT_DIR") || !grepl("de_diagnostics", SCRIPT_DIR)) {
  possible_config_paths <- c(
    "R/de_diagnostics/config.R",
    "config.R",
    "../config.R"
  )

  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }

  if (!config_loaded) {
    stop("Could not find config.R. Run from project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

log_section("STEP 3: SPLIT DEG RECOVERY")
log_info("Thresholds: padj < %g, |LFC| > %g",
         RECOVERY_FDR_THRESHOLD, RECOVERY_LFC_THRESHOLD)
log_info("LFC bins: %s", paste(names(RECOVERY_LFC_BINS), collapse = ", "))
log_info("Modes: %s | Runs: %d | Cell types: %s",
         paste(CARTESIAN_MODES, collapse = "+"), N_ROBUSTNESS_RUNS,
         paste(CELL_TYPES, collapse = ", "))

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORIES
# ==============================================================================

top_output_dir <- get_split_recovery_output_dir()
if (!dir.exists(top_output_dir)) {
  dir.create(top_output_dir, recursive = TRUE)
  log_info("Created output directory: %s", top_output_dir)
}

for (ct in CELL_TYPES) {
  ct_dir <- get_split_recovery_output_dir(ct)
  if (!dir.exists(ct_dir)) {
    dir.create(ct_dir, recursive = TRUE)
  }
}

# Check skip-if-exists (look for summary CSV as sentinel)
if (!FORCE_RECALC) {
  summary_file <- file.path(top_output_dir, "Recovery_Summary.csv")
  summary_strat_file <- file.path(top_output_dir, "Recovery_Summary_Stratified.csv")
  if (file.exists(summary_file) && file.exists(summary_strat_file)) {
    log_info("[SKIP] Recovery outputs already exist. Set FORCE_RECALC=TRUE to regenerate.")
    .skip_step3 <- TRUE
  } else {
    .skip_step3 <- FALSE
  }
} else {
  .skip_step3 <- FALSE
}

if (!exists(".skip_step3") || !.skip_step3) {

# ==============================================================================
# 3. LOAD FULL GT DE RESULTS
# ==============================================================================

log_section("Loading Full Ground Truth DE results")

full_gt_list <- list()

for (ct in CELL_TYPES) {
  # Use the recovery thresholds (|LFC| > 1), not the density plot thresholds (|LFC| > 2)
  df <- load_de_results(ct,
                        padj_cutoff = RECOVERY_FDR_THRESHOLD,
                        lfc_min = RECOVERY_LFC_THRESHOLD)
  if (!is.null(df)) {
    full_gt_list[[ct]] <- df
  }
}

if (length(full_gt_list) == 0) {
  stop("No Full GT DE results loaded for any cell type. Check paths and filters.")
}

# Also load unfiltered Full GT for LFC scatter (need LFC for all genes)
full_gt_unfiltered <- list()
for (ct in CELL_TYPES) {
  csv_path <- get_de_results_path(ct)
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    df$gene_id <- clean_gene_ids(df$gene)
    rownames(df) <- df$gene_id
    # Remove rows with NA padj or LFC
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    full_gt_unfiltered[[ct]] <- df
    log_info("  %s: %d genes (unfiltered, non-NA)", ct, nrow(df))
  }
}

# ==============================================================================
# 4. GET CHROMOSOME ANNOTATIONS
# ==============================================================================

log_section("Loading chromosome annotations")

cache_path <- file.path(OUTPUT_BASE, BIOMART_CACHE_FILE)

# Collect all gene IDs across all cell types (from unfiltered for full coverage)
all_gene_ids <- unique(unlist(lapply(full_gt_unfiltered, function(df) df$gene_id)))
log_info("Total unique genes across cell types: %d", length(all_gene_ids))

chr_annotations <- get_chromosome_annotations(all_gene_ids, cache_path)

if (!is.null(chr_annotations) && nrow(chr_annotations) > 0) {
  log_info("Chromosome annotations: %d genes annotated", nrow(chr_annotations))
  chr_table <- table(chr_annotations$Chr_Type)
  log_info("  Breakdown: %s",
           paste(sprintf("%s=%d", names(chr_table), as.integer(chr_table)),
                 collapse = ", "))
} else {
  log_warn("Chromosome annotation failed — stratification by chromosome will use NA")
}

# ==============================================================================
# 5. COMPUTE OVERLAP STATISTICS
# ==============================================================================

log_section("Computing overlap statistics across all cell types, modes, and runs")

# Collectors for summary tables
recovery_rows <- list()
recovery_strat_rows <- list()

# Collector for UpSet data (per cell type × mode)
upset_data <- list()

# Collector for LFC scatter data (per cell type × mode)
scatter_data <- list()

for (ct in CELL_TYPES) {

  log_info("--- Cell Type: %s ---", ct)

  if (!(ct %in% names(full_gt_list))) {
    log_warn("  No Full GT results for %s — skipping", ct)
    next
  }

  full_gt_df <- full_gt_list[[ct]]
  full_gt_genes <- full_gt_df$gene_id

  # Annotate Full GT genes with LFC bin and chromosome type
  full_gt_df$Abs_LFC <- abs(full_gt_df$log2FoldChange)
  full_gt_df$LFC_Bin <- assign_lfc_bin(full_gt_df$Abs_LFC)

  if (!is.null(chr_annotations)) {
    chr_match <- chr_annotations[match(full_gt_df$gene_id,
                                        chr_annotations$ensembl_gene_id), ]
    full_gt_df$Chr_Type <- chr_match$Chr_Type
  } else {
    full_gt_df$Chr_Type <- NA_character_
  }

  for (mode in CARTESIAN_MODES) {
    log_info("  Mode: %s", mode)

    # Key for UpSet: "CellType_mode"
    upset_key <- paste(ct, mode, sep = "_")
    # Start with Full GT gene set
    upset_sets <- list("Full_GT" = full_gt_genes)

    # Collector for per-run LFC values (for averaging)
    run_lfc_list <- list()

    for (run_id in seq_len(N_ROBUSTNESS_RUNS)) {

      # Load split GT results
      gt_path <- get_cartesian_gt_path(ct, mode, run_id)
      split_df <- load_split_gt_results(gt_path)

      if (is.null(split_df)) {
        log_warn("    Run %d: failed to load — skipping", run_id)
        next
      }

      # Filter to significant genes using recovery thresholds
      split_sig <- split_df[!is.na(split_df$padj) &
                              split_df$padj < RECOVERY_FDR_THRESHOLD &
                              abs(split_df$log2FoldChange) > RECOVERY_LFC_THRESHOLD, ]
      split_genes <- split_sig$gene_id

      log_info("    Run %d: %d significant genes in split (vs %d in Full GT)",
               run_id, length(split_genes), length(full_gt_genes))

      # --- Overall overlap ---
      stats <- compute_overlap_stats(full_gt_genes, split_genes)

      recovery_rows[[length(recovery_rows) + 1]] <- data.frame(
        Cell_Type = ct,
        Mode = mode,
        Run = run_id,
        N_FullGT = stats$N_FullGT,
        N_Split = stats$N_Split,
        N_Overlap = stats$N_Overlap,
        N_FullOnly = stats$N_FullOnly,
        N_SplitOnly = stats$N_SplitOnly,
        Recall = stats$Recall,
        Precision = stats$Precision,
        Stratum = "Overall",
        stringsAsFactors = FALSE
      )

      # --- Stratified by LFC bin ---
      for (bin_name in names(RECOVERY_LFC_BINS)) {
        bin_bounds <- RECOVERY_LFC_BINS[[bin_name]]

        # Full GT genes in this bin
        bin_full <- full_gt_df$gene_id[full_gt_df$LFC_Bin == bin_name &
                                         !is.na(full_gt_df$LFC_Bin)]

        # Split genes in this bin (need to compute abs LFC for split)
        split_sig$Abs_LFC <- abs(split_sig$log2FoldChange)
        split_sig$LFC_Bin <- assign_lfc_bin(split_sig$Abs_LFC)
        bin_split <- split_sig$gene_id[split_sig$LFC_Bin == bin_name &
                                         !is.na(split_sig$LFC_Bin)]

        bin_stats <- compute_overlap_stats(bin_full, bin_split)

        recovery_strat_rows[[length(recovery_strat_rows) + 1]] <- data.frame(
          Cell_Type = ct,
          Mode = mode,
          Run = run_id,
          N_FullGT = bin_stats$N_FullGT,
          N_Split = bin_stats$N_Split,
          N_Overlap = bin_stats$N_Overlap,
          N_FullOnly = bin_stats$N_FullOnly,
          N_SplitOnly = bin_stats$N_SplitOnly,
          Recall = bin_stats$Recall,
          Precision = bin_stats$Precision,
          Stratum = paste0("LFC_", bin_name),
          stringsAsFactors = FALSE
        )
      }

      # --- Stratified by chromosome type ---
      if (!is.null(chr_annotations) && !all(is.na(full_gt_df$Chr_Type))) {
        for (chr_type in c("Autosomal", "SexChr")) {
          chr_full <- full_gt_df$gene_id[full_gt_df$Chr_Type == chr_type &
                                           !is.na(full_gt_df$Chr_Type)]

          # Annotate split genes with chromosome
          split_chr <- chr_annotations$Chr_Type[
            match(split_sig$gene_id, chr_annotations$ensembl_gene_id)
          ]
          chr_split <- split_sig$gene_id[split_chr == chr_type &
                                           !is.na(split_chr)]

          chr_stats <- compute_overlap_stats(chr_full, chr_split)

          recovery_strat_rows[[length(recovery_strat_rows) + 1]] <- data.frame(
            Cell_Type = ct,
            Mode = mode,
            Run = run_id,
            N_FullGT = chr_stats$N_FullGT,
            N_Split = chr_stats$N_Split,
            N_Overlap = chr_stats$N_Overlap,
            N_FullOnly = chr_stats$N_FullOnly,
            N_SplitOnly = chr_stats$N_SplitOnly,
            Recall = chr_stats$Recall,
            Precision = chr_stats$Precision,
            Stratum = paste0("Chr_", chr_type),
            stringsAsFactors = FALSE
          )
        }
      }

      # Add to UpSet sets
      upset_sets[[paste0("Run_", run_id)]] <- split_genes

      # Store LFC values for scatter (using unfiltered data for all genes)
      run_lfc_list[[run_id]] <- split_df[, c("gene_id", "log2FoldChange")]

    } # end run loop

    # Store UpSet data for this cell type × mode
    upset_data[[upset_key]] <- upset_sets

    # Compute mean LFC across runs for scatter
    if (length(run_lfc_list) > 0) {
      # Merge all runs by gene_id
      merged_lfc <- run_lfc_list[[1]]
      colnames(merged_lfc)[2] <- "LFC_Run1"

      for (r in seq_along(run_lfc_list)[-1]) {
        if (!is.null(run_lfc_list[[r]])) {
          tmp <- run_lfc_list[[r]]
          colnames(tmp)[2] <- paste0("LFC_Run", r)
          merged_lfc <- merge(merged_lfc, tmp, by = "gene_id", all = TRUE)
        }
      }

      # Compute mean LFC across runs
      lfc_cols <- grep("^LFC_Run", colnames(merged_lfc), value = TRUE)
      merged_lfc$Mean_Split_LFC <- rowMeans(merged_lfc[, lfc_cols, drop = FALSE],
                                             na.rm = TRUE)

      # Merge with Full GT unfiltered LFC
      if (ct %in% names(full_gt_unfiltered)) {
        full_uf <- full_gt_unfiltered[[ct]][, c("gene_id", "log2FoldChange", "padj")]
        colnames(full_uf)[2:3] <- c("Full_GT_LFC", "Full_GT_padj")

        scatter_df <- merge(merged_lfc[, c("gene_id", "Mean_Split_LFC")],
                            full_uf, by = "gene_id", all = TRUE)

        # Compute mean split padj for significance classification
        # (use per-run significance status instead)
        scatter_df$Full_Sig <- !is.na(scatter_df$Full_GT_padj) &
          scatter_df$Full_GT_padj < RECOVERY_FDR_THRESHOLD &
          !is.na(scatter_df$Full_GT_LFC) &
          abs(scatter_df$Full_GT_LFC) > RECOVERY_LFC_THRESHOLD

        # For split significance: gene is "split significant" if significant
        # in at least half the runs
        sig_count <- rep(0, nrow(scatter_df))
        for (r in seq_along(run_lfc_list)) {
          if (!is.null(run_lfc_list[[r]])) {
            gt_path <- get_cartesian_gt_path(ct, mode, r)
            split_df_r <- load_split_gt_results(gt_path)
            if (!is.null(split_df_r)) {
              sig_genes_r <- split_df_r$gene_id[
                !is.na(split_df_r$padj) &
                  split_df_r$padj < RECOVERY_FDR_THRESHOLD &
                  abs(split_df_r$log2FoldChange) > RECOVERY_LFC_THRESHOLD
              ]
              sig_count <- sig_count +
                as.integer(scatter_df$gene_id %in% sig_genes_r)
            }
          }
        }
        n_valid_runs <- sum(!sapply(run_lfc_list, is.null))
        scatter_df$Split_Sig <- sig_count >= ceiling(n_valid_runs / 2)

        # Classify points
        scatter_df$Category <- "Neither"
        scatter_df$Category[scatter_df$Full_Sig & scatter_df$Split_Sig] <- "Both Significant"
        scatter_df$Category[scatter_df$Full_Sig & !scatter_df$Split_Sig] <- "Full GT Only"
        scatter_df$Category[!scatter_df$Full_Sig & scatter_df$Split_Sig] <- "Split Only"

        scatter_data[[upset_key]] <- scatter_df
      }
    }

  } # end mode loop
} # end cell type loop

# ==============================================================================
# 6. BUILD SUMMARY DATA FRAMES
# ==============================================================================

log_section("Building summary tables")

recovery_summary <- do.call(rbind, recovery_rows)
recovery_strat_summary <- do.call(rbind, recovery_strat_rows)

# Combine overall + stratified for the stratified CSV
recovery_all_strat <- rbind(recovery_summary, recovery_strat_summary)

log_info("Overall summary: %d rows", nrow(recovery_summary))
log_info("Stratified summary: %d rows", nrow(recovery_strat_summary))

# Save CSVs
summary_file <- file.path(top_output_dir, "Recovery_Summary.csv")
write.csv(recovery_summary, summary_file, row.names = FALSE)
log_info("Saved: %s", basename(summary_file))

summary_strat_file <- file.path(top_output_dir, "Recovery_Summary_Stratified.csv")
write.csv(recovery_all_strat, summary_strat_file, row.names = FALSE)
log_info("Saved: %s", basename(summary_strat_file))

# Print summary statistics
for (ct in CELL_TYPES) {
  ct_data <- recovery_summary[recovery_summary$Cell_Type == ct, ]
  if (nrow(ct_data) > 0) {
    log_info("  %s: mean recall=%.3f (SD=%.3f), mean precision=%.3f (SD=%.3f)",
             ct,
             mean(ct_data$Recall, na.rm = TRUE),
             sd(ct_data$Recall, na.rm = TRUE),
             mean(ct_data$Precision, na.rm = TRUE),
             sd(ct_data$Precision, na.rm = TRUE))
  }
}

# ==============================================================================
# 7. GENERATE RECOVERY BAR CHARTS
# ==============================================================================

log_section("Generating recovery bar charts")

for (ct in CELL_TYPES) {
  for (mode in CARTESIAN_MODES) {

    ct_mode_data <- recovery_summary[recovery_summary$Cell_Type == ct &
                                       recovery_summary$Mode == mode, ]

    if (nrow(ct_mode_data) == 0) {
      log_warn("  No data for %s/%s — skipping bar chart", ct, mode)
      next
    }

    # --- Overall bar chart ---
    bar_df <- ct_mode_data %>%
      dplyr::select(Run, N_FullGT, N_Split, N_Overlap) %>%
      tidyr::pivot_longer(cols = c(N_FullGT, N_Split, N_Overlap),
                          names_to = "Category",
                          values_to = "Count") %>%
      dplyr::mutate(
        Category = factor(Category,
                          levels = c("N_FullGT", "N_Split", "N_Overlap"),
                          labels = c("Full GT", "Split GT", "Overlap")),
        Run = factor(Run)
      )

    # Add recall/precision as text annotation
    annot_df <- ct_mode_data %>%
      dplyr::mutate(
        Label = sprintf("R=%.0f%% P=%.0f%%", Recall * 100, Precision * 100),
        Run = factor(Run)
      )

    p_bar <- ggplot(bar_df, aes(x = Run, y = Count, fill = Category)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(data = annot_df,
                aes(x = Run, y = max(ct_mode_data$N_FullGT,
                                      ct_mode_data$N_Split) * 1.05,
                    label = Label),
                inherit.aes = FALSE, size = 3, vjust = 0) +
      scale_fill_manual(values = c("Full GT" = RECOVERY_BAR_COLOURS["Full_GT"],
                                    "Split GT" = RECOVERY_BAR_COLOURS["Split_GT"],
                                    "Overlap" = RECOVERY_BAR_COLOURS["Overlap"])) +
      labs(
        title = sprintf("%s — Split DEG Recovery (%s mode)", ct, mode),
        subtitle = sprintf("Thresholds: padj < %g, |LFC| > %g",
                           RECOVERY_FDR_THRESHOLD, RECOVERY_LFC_THRESHOLD),
        x = "Robustness Run",
        y = "Number of DEGs",
        fill = ""
      ) +
      theme_bw() +
      theme(legend.position = "bottom") +
      coord_cartesian(ylim = c(0, max(ct_mode_data$N_FullGT,
                                       ct_mode_data$N_Split) * 1.15))

    bar_file <- file.path(get_split_recovery_output_dir(ct),
                          sprintf("Recovery_Bars_%s.png", mode))
    ggsave(bar_file, p_bar, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    log_info("  Saved: %s", basename(bar_file))

    # --- Stratified bar chart (LFC bins + chromosome) ---
    ct_mode_strat <- recovery_all_strat[
      recovery_all_strat$Cell_Type == ct &
        recovery_all_strat$Mode == mode &
        recovery_all_strat$Stratum != "Overall", ]

    if (nrow(ct_mode_strat) > 0) {
      # Compute mean recall and precision per stratum
      strat_agg <- ct_mode_strat %>%
        dplyr::group_by(Stratum) %>%
        dplyr::summarise(
          Mean_Recall = mean(Recall, na.rm = TRUE),
          Mean_Precision = mean(Precision, na.rm = TRUE),
          SD_Recall = sd(Recall, na.rm = TRUE),
          SD_Precision = sd(Precision, na.rm = TRUE),
          Mean_N_FullGT = mean(N_FullGT, na.rm = TRUE),
          Mean_N_Overlap = mean(N_Overlap, na.rm = TRUE),
          .groups = "drop"
        )

      # Pivot for grouped bar chart
      strat_bar <- strat_agg %>%
        dplyr::select(Stratum, Mean_Recall, Mean_Precision) %>%
        tidyr::pivot_longer(cols = c(Mean_Recall, Mean_Precision),
                            names_to = "Metric",
                            values_to = "Value") %>%
        dplyr::mutate(
          Metric = factor(Metric,
                          levels = c("Mean_Recall", "Mean_Precision"),
                          labels = c("Recall", "Precision"))
        )

      # Error bars
      strat_err <- strat_agg %>%
        dplyr::select(Stratum, SD_Recall, SD_Precision) %>%
        tidyr::pivot_longer(cols = c(SD_Recall, SD_Precision),
                            names_to = "Metric",
                            values_to = "SD") %>%
        dplyr::mutate(
          Metric = factor(Metric,
                          levels = c("SD_Recall", "SD_Precision"),
                          labels = c("Recall", "Precision"))
        )

      strat_bar <- merge(strat_bar, strat_err, by = c("Stratum", "Metric"))

      p_strat <- ggplot(strat_bar, aes(x = Stratum, y = Value, fill = Metric)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_errorbar(aes(ymin = pmax(Value - SD, 0), ymax = pmin(Value + SD, 1)),
                      position = position_dodge(width = 0.7), width = 0.2) +
        scale_fill_manual(values = c("Recall" = "#2166AC", "Precision" = "#B2182B")) +
        scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
        labs(
          title = sprintf("%s — Stratified Recovery (%s mode)", ct, mode),
          subtitle = sprintf("Mean ± SD across %d runs", N_ROBUSTNESS_RUNS),
          x = "Stratum",
          y = "Rate",
          fill = ""
        ) +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 30, hjust = 1))

      strat_file <- file.path(get_split_recovery_output_dir(ct),
                              sprintf("Recovery_Bars_%s_Stratified.png", mode))
      ggsave(strat_file, p_strat,
             width = PLOT_WIDTH_STRATIFIED, height = PLOT_HEIGHT_STRATIFIED,
             dpi = PLOT_DPI)
      log_info("  Saved: %s", basename(strat_file))
    }

  } # end mode loop
} # end cell type loop

# ==============================================================================
# 8. GENERATE UPSET PLOTS
# ==============================================================================

log_section("Generating UpSet plots")

for (ct in CELL_TYPES) {
  for (mode in CARTESIAN_MODES) {

    upset_key <- paste(ct, mode, sep = "_")
    sets <- upset_data[[upset_key]]

    if (is.null(sets) || length(sets) < 2) {
      log_warn("  No UpSet data for %s/%s — skipping", ct, mode)
      next
    }

    # Convert to binary membership matrix
    all_genes <- unique(unlist(sets))
    membership <- data.frame(
      row.names = all_genes,
      stringsAsFactors = FALSE
    )
    for (set_name in names(sets)) {
      membership[[set_name]] <- as.integer(all_genes %in% sets[[set_name]])
    }

    # Generate UpSet plot
    upset_file <- file.path(get_split_recovery_output_dir(ct),
                            sprintf("UpSet_%s.png", mode))

    tryCatch({
      png(upset_file, width = PLOT_WIDTH * PLOT_DPI,
          height = PLOT_HEIGHT * PLOT_DPI, res = PLOT_DPI)

      print(UpSetR::upset(
        membership,
        nsets = UPSET_NSETS,
        order.by = UPSET_ORDER_BY,
        sets = rev(names(sets)),  # Full_GT first
        keep.order = TRUE,
        text.scale = c(1.3, 1.2, 1.0, 1.0, 1.3, 1.0),
        mainbar.y.label = "Intersection Size",
        sets.x.label = "DEGs per Set",
        main.bar.color = "#4DAF4A",
        sets.bar.color = c(RECOVERY_BAR_COLOURS["Full_GT"],
                           rep(RECOVERY_BAR_COLOURS["Split_GT"],
                               length(sets) - 1))
      ))

      dev.off()
      log_info("  Saved: %s", basename(upset_file))
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      log_warn("  UpSet plot failed for %s/%s: %s", ct, mode, e$message)
    })

  } # end mode loop
} # end cell type loop

# ==============================================================================
# 9. GENERATE LFC SCATTER PLOTS
# ==============================================================================

log_section("Generating LFC scatter plots")

for (ct in CELL_TYPES) {
  for (mode in CARTESIAN_MODES) {

    upset_key <- paste(ct, mode, sep = "_")
    scat_df <- scatter_data[[upset_key]]

    if (is.null(scat_df) || nrow(scat_df) == 0) {
      log_warn("  No scatter data for %s/%s — skipping", ct, mode)
      next
    }

    # Remove genes with NA in either LFC
    scat_clean <- scat_df[!is.na(scat_df$Full_GT_LFC) &
                            !is.na(scat_df$Mean_Split_LFC), ]

    if (nrow(scat_clean) == 0) {
      log_warn("  No non-NA scatter data for %s/%s — skipping", ct, mode)
      next
    }

    # Compute correlation
    cor_val <- cor(scat_clean$Full_GT_LFC, scat_clean$Mean_Split_LFC,
                   use = "complete.obs")

    # Count per category
    cat_counts <- table(scat_clean$Category)

    # Order categories for consistent layering (plot "Neither" first)
    scat_clean$Category <- factor(scat_clean$Category,
                                   levels = c("Neither", "Full GT Only",
                                              "Split Only", "Both Significant"))

    # Sort so "Both Significant" is plotted on top
    scat_clean <- scat_clean[order(scat_clean$Category), ]

    p_scatter <- ggplot(scat_clean,
                        aes(x = Full_GT_LFC, y = Mean_Split_LFC,
                            color = Category)) +
      geom_point(size = SCATTER_POINT_SIZE, alpha = SCATTER_ALPHA) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey40", linewidth = 0.5) +
      scale_color_manual(values = SCATTER_COLOURS) +
      labs(
        title = sprintf("%s — Full GT vs Mean Split LFC (%s mode)", ct, mode),
        subtitle = sprintf("r = %.3f | Both=%s, Full Only=%s, Split Only=%s",
                           cor_val,
                           ifelse("Both Significant" %in% names(cat_counts),
                                  cat_counts["Both Significant"], 0),
                           ifelse("Full GT Only" %in% names(cat_counts),
                                  cat_counts["Full GT Only"], 0),
                           ifelse("Split Only" %in% names(cat_counts),
                                  cat_counts["Split Only"], 0)),
        x = "Full GT log2FC",
        y = "Mean Split GT log2FC (across 6 runs)",
        color = ""
      ) +
      theme_bw() +
      theme(legend.position = "bottom") +
      coord_fixed(ratio = 1)

    scatter_file <- file.path(get_split_recovery_output_dir(ct),
                              sprintf("LFC_Scatter_%s.png", mode))
    ggsave(scatter_file, p_scatter,
           width = PLOT_WIDTH, height = PLOT_WIDTH,  # square
           dpi = PLOT_DPI)
    log_info("  Saved: %s", basename(scatter_file))

  } # end mode loop
} # end cell type loop

# ==============================================================================
# 10. COMPLETION REPORT
# ==============================================================================

log_section("STEP 3 COMPLETION REPORT")

# Count outputs
total_plots <- 0
for (ct in CELL_TYPES) {
  ct_dir <- get_split_recovery_output_dir(ct)
  if (dir.exists(ct_dir)) {
    n_plots <- length(list.files(ct_dir, pattern = "\\.png$"))
    total_plots <- total_plots + n_plots
    log_info("  %s: %d plots", ct, n_plots)
  }
}

log_info("Total plots generated: %d", total_plots)
log_info("Summary CSVs: Recovery_Summary.csv, Recovery_Summary_Stratified.csv")

# Print overall recovery statistics
log_info("Overall recovery statistics:")
for (ct in CELL_TYPES) {
  for (mode in CARTESIAN_MODES) {
    ct_mode <- recovery_summary[recovery_summary$Cell_Type == ct &
                                  recovery_summary$Mode == mode, ]
    if (nrow(ct_mode) > 0) {
      log_info("  %s/%s: recall=%.1f%% ± %.1f%%, precision=%.1f%% ± %.1f%%",
               ct, mode,
               mean(ct_mode$Recall, na.rm = TRUE) * 100,
               sd(ct_mode$Recall, na.rm = TRUE) * 100,
               mean(ct_mode$Precision, na.rm = TRUE) * 100,
               sd(ct_mode$Precision, na.rm = TRUE) * 100)
    }
  }
}

log_info("Step 3 complete")

} # end skip check
