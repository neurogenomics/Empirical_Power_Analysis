# ==============================================================================
# STEP 06: Recovery by Expression Level
# ==============================================================================
# PURPOSE: Analyse how bulk→SC gene recovery varies with expression level.
#   For each cell type, bulk-significant genes are binned into quintiles of
#   SC pseudobulk mean expression (baseMean from DESeq2) and recovery rates
#   are visualised as stacked bar charts matching the existing recovery palette.
#
# INPUT:
#   - Datalink/Bulk_Analysis/Filtered_Consistent_Genes.csv (from Step 02)
#   - Datalink/Full_Ground_Truth/{CT}/{Dataset}_Full_DESeq2_Results.qs
#
# OUTPUT:
#   - Per cell type × per scale (raw baseMean, log2 baseMean):
#     - Datalink/Bulk_Analysis/recovery_by_expression/{CT}_{Dataset}_{scale}.png
#   - CSV summaries:
#     - Datalink/Bulk_Analysis/recovery_by_expression/{CT}_{Dataset}_summary.csv
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/bulk_deg_processing/config.R",
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
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(glue)
  library(qs)
})

# ==============================================================================
# 1. LOAD BULK DATA
# ==============================================================================

message("\n", strrep("=", 60))
message("STEP 06: RECOVERY BY EXPRESSION LEVEL")
message(strrep("=", 60))

message("\nLoading bulk data...")
bulk_data <- read.csv(FILTERED_FILE, stringsAsFactors = FALSE)

if (nrow(bulk_data) == 0) {
  stop("Bulk data is empty. Run Step 02 first.")
}

message(sprintf("  Loaded %d bulk genes", nrow(bulk_data)))

# --- CHECK BULK DATA ORIENTATION VIA XIST ---
XIST_ID <- "ENSG00000229807"
xist_bulk <- bulk_data %>% filter(gene == XIST_ID)

if (nrow(xist_bulk) > 0 && !is.na(xist_bulk$log2FoldChange[1])) {
  xist_lfc_bulk <- xist_bulk$log2FoldChange[1]
  if (xist_lfc_bulk > 0) {
    message(sprintf("  Bulk XIST LFC = %.2f (positive) -> Flipping all bulk LFCs to match SC convention", xist_lfc_bulk))
    bulk_data$log2FoldChange <- -bulk_data$log2FoldChange
  } else {
    message(sprintf("  Bulk XIST LFC = %.2f (negative) -> No flip needed", xist_lfc_bulk))
  }
} else {
  message("  [WARNING] XIST not found in bulk data. Assuming correct orientation.")
}

# Filter for significant bulk genes
bulk_sig <- bulk_data %>%
  dplyr::filter(padj < P_VALUE_THRESHOLD & abs(log2FoldChange) > LFC_THRESHOLD)

message(sprintf("  %d significant bulk genes (padj < %.2f, |LFC| > %.2f)",
                nrow(bulk_sig), P_VALUE_THRESHOLD, LFC_THRESHOLD))

# --- ROBUST ID CLEANER ---
clean_ids <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*", "", x)
  x <- trimws(x)
  return(x)
}

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

RECOVERY_EXPR_DIR <- file.path(OUTPUT_DIR, "recovery_by_expression")
if (!dir.exists(RECOVERY_EXPR_DIR)) {
  dir.create(RECOVERY_EXPR_DIR, recursive = TRUE)
  message(sprintf("  Created output directory: %s", RECOVERY_EXPR_DIR))
}

# ==============================================================================
# 3. RECOVERY COLOURS (FROM CONFIG — 4-CATEGORY SCHEME)
# ==============================================================================
# Uses RECOVERY_COLORS and RECOVERY_LEVELS from config.R

# ==============================================================================
# 4. MAIN LOOP: PER CELL TYPE × PER SC DATASET
# ==============================================================================

CELL_TYPES <- c("Astro", "Micro", "Oligo")

for (Cell_Type in CELL_TYPES) {

  message(sprintf("\n=== PROCESSING CELL TYPE: %s ===", Cell_Type))

  for (Dataset in SC_DATASETS_FOR_CONCORDANCE) {

    tryCatch({
      message(sprintf("\n  Dataset: %s", Dataset))

      # --- A. LOAD SC DATA ---
      sc_file <- file.path(FULL_GT_BASE, Cell_Type,
                           paste0(Dataset, "_Full_DESeq2_Results.qs"))

      if (!file.exists(sc_file)) {
        message(sprintf("    [SKIP] SC file not found: %s", sc_file))
        next
      }

      message(sprintf("    Loading SC data from: %s", basename(sc_file)))
      sc_data_raw <- qread(sc_file)

      # Convert DESeqResults to data frame
      sc_data <- as.data.frame(sc_data_raw) %>%
        rownames_to_column("gene") %>%
        mutate(gene = clean_ids(gene))

      message(sprintf("    Loaded %d SC genes", nrow(sc_data)))

      # --- B. AUTO-DETECT LFC SIGN FLIP VIA XIST ---
      xist_row <- sc_data %>% filter(gene == XIST_ID)

      if (nrow(xist_row) > 0 && !is.na(xist_row$log2FoldChange[1])) {
        xist_lfc <- xist_row$log2FoldChange[1]
        if (xist_lfc > 0) {
          message(sprintf("    XIST LFC = %.2f (positive) -> Flipping all SC LFCs", xist_lfc))
          sc_data$log2FoldChange <- -sc_data$log2FoldChange
        } else {
          message(sprintf("    XIST LFC = %.2f (negative) -> No flip needed", xist_lfc))
        }
      } else {
        message("    [WARNING] XIST not found in SC data. Assuming correct orientation.")
      }

      # --- C. COMPUTE GENE UNIVERSE (INTERSECTION) ---
      # Genes expressed in both bulk and SC
      bulk_genes <- clean_ids(bulk_sig$gene)
      sc_genes   <- sc_data$gene
      universe   <- intersect(bulk_genes, sc_genes)

      message(sprintf("    Gene universe (bulk ∩ SC): %d genes", length(universe)))

      if (length(universe) == 0) {
        message("    [SKIP] No overlapping genes between bulk and SC")
        next
      }

      # --- D. BUILD RECOVERY TABLE (4-CATEGORY, DIRECTION-AWARE) ---
      # Start from ALL bulk-significant genes
      recovery_df <- bulk_sig %>%
        mutate(gene = clean_ids(gene)) %>%
        dplyr::select(gene, log2FoldChange) %>%
        dplyr::rename(lfc_bulk = log2FoldChange)

      # Join SC data (baseMean, LFC, padj) — left join keeps all bulk genes
      sc_lookup <- sc_data %>%
        dplyr::select(gene, baseMean, log2FoldChange, padj) %>%
        dplyr::rename(
          sc_baseMean = baseMean,
          lfc_sc      = log2FoldChange,
          padj_sc     = padj
        )

      recovery_df <- recovery_df %>%
        left_join(sc_lookup, by = "gene")

      # Classify recovery (4-category, direction-aware, matching Step 05)
      recovery_df <- recovery_df %>%
        mutate(
          Recovery = factor(case_when(
            !(gene %in% sc_genes)                                              ~ "Absent",
            is.na(padj_sc)                                                     ~ "Low-Count Filtered",
            padj_sc < P_VALUE_THRESHOLD & sign(lfc_sc) == sign(lfc_bulk)       ~ "Recovered",
            TRUE                                                               ~ "Insignificant"
          ), levels = RECOVERY_LEVELS)
        )

      message(sprintf("    Recovery breakdown: Recovered=%d, Insignificant=%d, NA-Filtered=%d, Absent=%d",
                      sum(recovery_df$Recovery == "Recovered"),
                      sum(recovery_df$Recovery == "Insignificant"),
                      sum(recovery_df$Recovery == "Low-Count Filtered"),
                      sum(recovery_df$Recovery == "Absent")))

      # --- E. QUINTILE BINNING ---
      # Only bin genes that have a valid SC baseMean (non-absent)
      has_expr <- !is.na(recovery_df$sc_baseMean)

      # Raw baseMean quintiles
      recovery_df$raw_quintile <- NA_character_
      if (sum(has_expr) >= 5) {
        breaks_raw <- quantile(recovery_df$sc_baseMean[has_expr],
                               probs = seq(0, 1, 0.2), na.rm = TRUE)
        # Make labels from break values
        raw_labels <- paste0("Q", 1:5, " [",
                             round(breaks_raw[1:5], 1), "-",
                             round(breaks_raw[2:6], 1), "]")
        recovery_df$raw_quintile[has_expr] <- as.character(cut(
          recovery_df$sc_baseMean[has_expr],
          breaks = breaks_raw,
          labels = raw_labels,
          include.lowest = TRUE
        ))
      }

      # Log2(baseMean + 1) quintiles
      recovery_df$log2_baseMean <- log2(recovery_df$sc_baseMean + 1)
      recovery_df$log2_quintile <- NA_character_
      if (sum(has_expr) >= 5) {
        breaks_log2 <- quantile(recovery_df$log2_baseMean[has_expr],
                                probs = seq(0, 1, 0.2), na.rm = TRUE)
        log2_labels <- paste0("Q", 1:5, " [",
                              round(breaks_log2[1:5], 2), "-",
                              round(breaks_log2[2:6], 2), "]")
        recovery_df$log2_quintile[has_expr] <- as.character(cut(
          recovery_df$log2_baseMean[has_expr],
          breaks = breaks_log2,
          labels = log2_labels,
          include.lowest = TRUE
        ))
      }

      # Assign absent genes to a special bin
      recovery_df$raw_quintile[is.na(recovery_df$raw_quintile)] <- "Not in SC"
      recovery_df$log2_quintile[is.na(recovery_df$log2_quintile)] <- "Not in SC"

      # --- F. SAVE CSV SUMMARY ---
      # Per-quintile summary for both scales
      make_summary <- function(df, quintile_col, scale_name) {
        df %>%
          group_by(quintile = .data[[quintile_col]], Recovery) %>%
          summarise(count = n(), .groups = "drop") %>%
          mutate(scale = scale_name)
      }

      summary_raw  <- make_summary(recovery_df, "raw_quintile", "raw_baseMean")
      summary_log2 <- make_summary(recovery_df, "log2_quintile", "log2_baseMean")
      summary_all  <- bind_rows(summary_raw, summary_log2)

      csv_file <- file.path(RECOVERY_EXPR_DIR,
                            paste0(Cell_Type, "_", Dataset, "_summary.csv"))
      write.csv(summary_all, csv_file, row.names = FALSE)
      message(sprintf("    Saved CSV summary: %s", basename(csv_file)))

      # --- G. PLOTTING FUNCTION ---
      plot_recovery_by_expression <- function(df, quintile_col, scale_label,
                                               cell_type, dataset) {
        # Prepare summary
        plot_summary <- df %>%
          group_by(quintile = .data[[quintile_col]], Recovery) %>%
          summarise(count = n(), .groups = "drop") %>%
          mutate(Recovery = factor(Recovery, levels = RECOVERY_LEVELS))

        # Get quintile order (Q1..Q5, then "Not in SC")
        q_levels <- sort(unique(as.character(plot_summary$quintile)))
        # Move "Not in SC" to end if present
        if ("Not in SC" %in% q_levels) {
          q_levels <- c(setdiff(q_levels, "Not in SC"), "Not in SC")
        }
        plot_summary$quintile <- factor(plot_summary$quintile, levels = q_levels)

        # Total per quintile for annotation
        totals <- plot_summary %>%
          group_by(quintile) %>%
          summarise(total = sum(count), .groups = "drop")

        p <- ggplot(plot_summary, aes(x = quintile, y = count, fill = Recovery)) +
          geom_bar(stat = "identity", position = "stack", width = 0.7) +
          geom_text(aes(label = count),
                    position = position_stack(vjust = 0.5),
                    color = "black", size = 4, fontface = "bold") +
          # Add total count above each bar
          geom_text(data = totals, aes(x = quintile, y = total, label = paste0("n=", total)),
                    inherit.aes = FALSE, vjust = -0.5, size = 3.5, color = "grey40") +
          scale_fill_manual(values = RECOVERY_COLORS) +
          theme_publication() +
          theme(
            axis.text.x = element_text(angle = 35, hjust = 1)
          ) +
          labs(
            title = glue("{prettify_label(cell_type)} — Bulk→SC Recovery by Expression Level"),
            subtitle = glue("Dataset: {prettify_label(dataset)} | Scale: {scale_label} | Recovered: FDR<0.05 + sign match"),
            x = glue("SC Mean Expression Quintile ({scale_label})"),
            y = "Number of Bulk-Significant Genes",
            fill = "Recovery Status"
          ) +
          # Give a bit of headroom for the n= labels
          scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

        return(p)
      }

      # --- H. GENERATE AND SAVE PLOTS ---

      # Raw baseMean plot
      p_raw <- plot_recovery_by_expression(
        recovery_df, "raw_quintile", "Raw baseMean",
        Cell_Type, Dataset
      )
      raw_file <- file.path(RECOVERY_EXPR_DIR,
                            paste0(Cell_Type, "_", Dataset, "_raw_baseMean.png"))
      ggsave(raw_file, p_raw, width = 8, height = 6, dpi = 150)
      message(sprintf("    Saved plot: %s", basename(raw_file)))

      # Log2 baseMean plot
      p_log2 <- plot_recovery_by_expression(
        recovery_df, "log2_quintile", "log2(baseMean + 1)",
        Cell_Type, Dataset
      )
      log2_file <- file.path(RECOVERY_EXPR_DIR,
                             paste0(Cell_Type, "_", Dataset, "_log2_baseMean.png"))
      ggsave(log2_file, p_log2, width = 8, height = 6, dpi = 150)
      message(sprintf("    Saved plot: %s", basename(log2_file)))

      message(sprintf("    Completed %s - %s", Cell_Type, Dataset))

    }, error = function(e) {
      message(sprintf("[ERROR] Processing %s - %s failed: %s", Dataset, Cell_Type, e$message))
    })
  }
}

# ==============================================================================
# 5. COMPLETION REPORT
# ==============================================================================

message("\n", strrep("=", 60))
message("STEP 06 COMPLETE")
message(strrep("=", 60))
message(sprintf("Generated recovery-by-expression plots for:"))
message(sprintf("  Cell types: %s", paste(CELL_TYPES, collapse = ", ")))
message(sprintf("  SC datasets: %s", paste(SC_DATASETS_FOR_CONCORDANCE, collapse = ", ")))
message(sprintf("  Output directory: %s", RECOVERY_EXPR_DIR))
message("\n")
