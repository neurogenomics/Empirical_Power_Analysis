# ==============================================================================
# STEP 05: Visualise Bulk vs Single-Cell Concordance
# ==============================================================================
# PURPOSE: Generate all bulk-vs-SC comparison plots including:
#   - Concordance scatter plots (bulk LFC vs SC LFC)
#   - Recovery statistics (Recovered/Insignificant/NA-Filtered/Absent)
#   - Volcano plots (SC DEGs with shared/unique translucency)
#   - LFC and SE distribution histograms (coding vs non-coding)
#   - Chromosome distribution plots
#   - Bulk volcano plot (global)
#
# INPUT:  
#   - Datalink/Bulk_Analysis/Filtered_Consistent_Genes.csv (from Step 02)
#   - Datalink/Full_Ground_Truth/{CT}/{Dataset}_Full_DESeq2_Results.qs
#
# OUTPUT: 
#   - Per cell type × per SC dataset:
#     - Datalink/Bulk_Analysis/{CT}/{Dataset}_*.png (15 plots per combination)
#   - Global:
#     - Datalink/Bulk_Analysis/Bulk_Results/Bulk_Volcano_Plot.png
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
  library(biomaRt)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(glue)
  library(httr)
  library(qs)
})

# Disable SSL verification for BioMart (as in archived script)
set_config(config(ssl_verifypeer = FALSE))

# ==============================================================================
# 1. LOAD BULK DATA
# ==============================================================================

message("\n", strrep("=", 60))
message("STEP 05: BULK VS SC VISUALIZATION")
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

# Filter for significant genes
bulk_sig <- bulk_data %>% 
  dplyr::filter(padj < P_VALUE_THRESHOLD & abs(log2FoldChange) > LFC_THRESHOLD)

message(sprintf("  %d significant bulk genes (padj < %.2f, |LFC| > %.2f)", 
                nrow(bulk_sig), P_VALUE_THRESHOLD, LFC_THRESHOLD))

# ==============================================================================
# 2. CONNECT TO BIOMART & PRE-FETCH GENE ANNOTATIONS (ONCE)
# ==============================================================================
# BioMart annotations (chromosome, gene name, biotype) are gene-level properties
# that don't change per cell type. We collect all gene IDs across all SC datasets
# and cell types, union with bulk genes, and query once.
# ==============================================================================

# --- ROBUST ID CLEANER (needed before gene collection) ---
clean_ids <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*", "", x)  # Remove everything after the first dot
  x <- trimws(x)            # Remove accidental whitespace
  return(x)
}

# --- COLLECT ALL GENE IDs ACROSS CELL TYPES & DATASETS ---
CELL_TYPES <- c("Astro", "Micro", "Oligo")
all_gene_ids <- unique(clean_ids(bulk_data$gene))

for (ct in CELL_TYPES) {
  for (ds in SC_DATASETS_FOR_CONCORDANCE) {
    sc_file_tmp <- file.path(FULL_GT_BASE, ct, paste0(ds, "_Full_DESeq2_Results.qs"))
    if (file.exists(sc_file_tmp)) {
      sc_tmp <- as.data.frame(qread(sc_file_tmp))
      all_gene_ids <- unique(c(all_gene_ids, clean_ids(rownames(sc_tmp))))
    }
  }
}

# --- BioMart with local cache fallback ---
BIOMART_CACHE <- file.path(OUTPUT_DIR, "biomart_gene_coords_cache.qs")

if (file.exists(BIOMART_CACHE)) {
  message("  Loading cached BioMart annotations...")
  raw_coords <- qread(BIOMART_CACHE)
  message(sprintf("  Loaded %d cached annotations", nrow(raw_coords)))
} else {
  message("\nConnecting to BioMart...")
  message(sprintf("  Ensembl version: %d", ENSEMBL_VERSION))

  ensembl <- tryCatch({
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = ENSEMBL_VERSION)
  }, error = function(e) {
    stop(sprintf("Failed to connect to BioMart: %s\nCheck internet connection.\nTo use cached annotations, run this script once with internet access.", e$message))
  })

  message("  Connected successfully")
  message(sprintf("  Querying BioMart for %d unique gene IDs (chr + biotype)...", length(all_gene_ids)))

  raw_coords <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name", "external_gene_name", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = all_gene_ids,
    mart = ensembl
  )

  # Cache for future runs
  qsave(raw_coords, BIOMART_CACHE)
  message(sprintf("  Cached %d annotations to %s", nrow(raw_coords), BIOMART_CACHE))
}

gene_coords <- raw_coords %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    external_gene_name = dplyr::first(external_gene_name),
    chr_name = dplyr::first(chromosome_name),
    biotype_raw = dplyr::first(gene_biotype)
  ) %>%
  mutate(
    # Fold Patch/Haplotype into Autosome (patches are alternate autosomal assemblies)
    robust_chr_type = case_when(
      chr_name %in% c("Y", "chrY") ~ "Y",
      chr_name %in% c("X", "chrX") ~ "X",
      TRUE ~ "Autosome"
    ),
    Biotype_Class = ifelse(biotype_raw == "protein_coding", "Coding", "Non-Coding")
  ) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, robust_chr_type, Biotype_Class)

message(sprintf("  Retrieved annotations for %d / %d genes (%.1f%%)",
                nrow(gene_coords), length(all_gene_ids),
                100 * nrow(gene_coords) / length(all_gene_ids)))

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================
# Note: clean_ids() is defined in Section 2 (needed before gene collection)

# --- VOLCANO PLOT FUNCTION ---
generate_volcano_plot <- function(df, title_prefix, fc_thresh, p_thresh, shared_list = NULL, unique_list = NULL) {
  
  # 1. Prepare Data
  df_plot <- df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      padj_safe = ifelse(padj == 0, .Machine$double.xmin, padj),
      logP = -log10(padj_safe),
      Sig_Cat = case_when(
        padj < p_thresh & log2FoldChange > fc_thresh ~ "Up",
        padj < p_thresh & log2FoldChange < -fc_thresh ~ "Down",
        TRUE ~ "NS"
      ),
      # Identify genes that are strictly unique for the "Opaque" setting
      Is_Shared = ifelse(!is.null(unique_list) & (gene %in% unique_list), "Unique", "Shared")
    )
  
  # 2. Select Labels (Conditional Logic)
  if (!is.null(unique_list)) {
    # CASE A: Unique List Provided (SC Plots) -> LABEL TOP 10 HIGHEST ABSOLUTE LFC
    top_labels <- df_plot %>%
      filter(Sig_Cat != "NS") %>%
      filter(gene %in% unique_list) %>%        # Filter for unique
      arrange(desc(abs(log2FoldChange))) %>%    # Sort by Magnitude (High LFC)
      slice(1:10)                               # Take Top 10
    
    subtitle_str <- glue("Labels = Top 10 Unique Genes (Highest LFC)\nShared genes are translucent")
    label_overlaps <- Inf # Force these 10 to show
    
  } else {
    # CASE B: No Unique List (Bulk Plot) -> Safety Limit (Top 10 by P-value)
    top_labels <- df_plot %>%
      filter(Sig_Cat != "NS") %>%
      group_by(Sig_Cat) %>%
      slice_min(order_by = padj, n = 10) %>%
      ungroup()
    
    subtitle_str <- glue("Labels = Top 10 Significant Genes (Low P-val)\nShared genes are translucent")
    label_overlaps <- 20 # Standard overlap limit
  }
  
  # 3. Plot
  p <- ggplot(df_plot, aes(x = log2FoldChange, y = logP, color = Sig_Cat, alpha = Is_Shared)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c("Up" = "#F8766D", "Down" = "#00BFC4", "NS" = "grey85")) +
    scale_alpha_manual(values = c("Shared" = 0.2, "Unique" = 1.0), guide = "none") +
    
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "grey50", alpha = 0.5) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "grey50", alpha = 0.5) +
    
    # Labels
    geom_text_repel(data = top_labels, aes(label = external_gene_name), 
                    color = "black", size = 3, 
                    max.overlaps = label_overlaps, 
                    box.padding = 0.5,
                    alpha = 1.0) + 
    
    theme_publication() +
    labs(title = paste0(prettify_label(title_prefix), ": Volcano Plot"),
         x = expression(Log[2]~Fold~Change),
         y = expression(-log[10](adjusted~italic(p)))) +
    guides(color = guide_legend(title = "Significance")) +
    theme(legend.position = "right")
  
  return(p)
}

# --- CORRELATION PLOT FUNCTION ---
generate_corr_plot <- function(df, title_prefix, cell_type, dataset_name) {
  df <- df %>% mutate(sc_sig_status = ifelse(!is.na(padj_SC) & padj_SC < P_VALUE_THRESHOLD, "Significant", "NS"))
  df_calc <- df %>% filter(sc_sig_status == "Significant")
  r_all <- round(cor(df_calc$log2FoldChange_SC, df_calc$log2FoldChange_Bulk, method = "pearson", use = "complete.obs"), 2)
  label_all <- glue("All Sig Genes (R = {r_all})")
  
  p <- ggplot(df, aes(x = log2FoldChange_Bulk, y = log2FoldChange_SC)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey80") +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_point(aes(color = chr_type, alpha = sc_sig_status), size = 2) + 
    geom_smooth(data = df_calc, method = "lm", se = FALSE, aes(linetype = "All Genes"), color = "blue", linewidth = 0.8)
  
  if (VIS_CALC_SUBGROUP_LINES) {
    df_auto <- df_calc %>% filter(chr_type == "Autosome")
    r_auto  <- round(cor(df_auto$log2FoldChange_SC, df_auto$log2FoldChange_Bulk, method = "pearson", use = "complete.obs"), 2)
    if(is.na(r_auto)) r_auto <- "NA"
    df_sex  <- df_calc %>% filter(chr_type %in% c("X", "Y"))
    r_sex   <- round(cor(df_sex$log2FoldChange_SC, df_sex$log2FoldChange_Bulk, method = "pearson", use = "complete.obs"), 2)
    if(is.na(r_sex)) r_sex <- "NA"
    label_auto <- glue("Sig Autosomes (R = {r_auto})")
    label_sex  <- glue("Sig Sex Chr (R = {r_sex})")
    p <- p + geom_smooth(data = df_auto, method = "lm", se = FALSE, aes(linetype = "Autosomes"), color = "grey50", linewidth = 0.8) +
      geom_smooth(data = df_sex, method = "lm", se = FALSE, aes(linetype = "Sex Chr"), color = "purple", linewidth = 0.8) +
      scale_linetype_manual(name = "Trend Lines (Sig Only)", values = c("All Genes" = "dotted", "Autosomes" = "solid", "Sex Chr" = "solid"), labels = c("All Genes" = label_all, "Autosomes" = label_auto, "Sex Chr" = label_sex))
  } else {
    p <- p + scale_linetype_manual(name = "Trend Lines (Sig Only)", values = c("All Genes" = "dotted"), labels = c("All Genes" = label_all))
  }
  
  p <- p + geom_text_repel(aes(label = label_gene), max.overlaps = 50, box.padding = 0.5, size=3.5) +
    scale_color_manual(name = "Chromosome", values = c("Autosome" = "grey60", "X" = "#f8766d", "Y" = "#00bfc4")) + 
    scale_alpha_manual(name = "SC Significance", values = c("Significant" = 1.0, "NS" = 0.4)) +
    theme_publication() +
    labs(title = paste0(prettify_label(cell_type), ": Bulk vs Single-Cell Concordance"),
         subtitle = NULL,
         x = expression(Bulk~Log[2]~Fold~Change),
         y = expression(Single-Cell~Log[2]~Fold~Change))
  return(p)
}

# --- LFC CATEGORIZATION FUNCTION ---
categorize_lfc <- function(df) {
  # Use bin_by_lfc_interval() from utils.R for flexible interval handling
  interval_names <- names(LFC_INTERVALS)
  below_label <- glue("< {LFC_INTERVALS[[1]][1]}")
  
  df %>% mutate(
    abs_lfc = abs(log2FoldChange),
    lfc_group = {
      labels <- bin_by_lfc_interval(abs_lfc, LFC_INTERVALS)
      labels[is.na(labels)] <- below_label
      factor(labels, levels = c(interval_names, below_label))
    }
  )
}

# ==============================================================================
# 4. PHASE 1: BULK VOLCANO PLOT (GLOBAL)
# ==============================================================================

message("\n", strrep("=", 60))
message("PHASE 1: Generating Bulk Volcano Plot")
message(strrep("=", 60))

tryCatch({
  # Use pre-fetched gene annotations (from Section 2)
  bulk_volcano_df <- bulk_data %>%
    left_join(gene_coords %>% dplyr::select(ensembl_gene_id, external_gene_name),
              by = c("gene" = "ensembl_gene_id")) %>%
    mutate(external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "", gene, external_gene_name))
  
  # For Bulk, unique_list is NULL, so it defaults to Top 10 safety limit (sorted by padj)
  p_bulk_volcano <- generate_volcano_plot(bulk_volcano_df, "Bulk Data", LFC_THRESHOLD, P_VALUE_THRESHOLD)
  
  BULK_OUT_DIR <- file.path(OUTPUT_DIR, "Bulk_Results")
  if(!dir.exists(BULK_OUT_DIR)) dir.create(BULK_OUT_DIR, recursive=TRUE)
  
  ggsave(file.path(BULK_OUT_DIR, "Bulk_Volcano_Plot.png"), p_bulk_volcano, width = 8, height = 6)
  message("  > Saved Bulk_Volcano_Plot.png")
  
}, error = function(e) {
  message(sprintf("[ERROR] Bulk Volcano generation failed: %s", e$message))
})

# ==============================================================================
# 5. PHASE 2: PRE-CALCULATE SHARED/UNIQUE GENE LISTS
# ==============================================================================

message("\n", strrep("=", 60))
message("PHASE 2: Pre-calculating Significance Lists")
message(strrep("=", 60))

# We'll calculate this per dataset
# For now, just initialize storage
sig_lists_by_dataset <- list()

# ==============================================================================
# 6. MAIN EXECUTION LOOP: PER CELL TYPE × PER SC DATASET
# ==============================================================================

message("\n", strrep("=", 60))
message("PHASE 3: Generating SC Plots")
message(strrep("=", 60))

CELL_TYPES <- c("Astro", "Micro", "Oligo")

for (Cell_Type in CELL_TYPES) {
  
  message(sprintf("\n=== PROCESSING CELL TYPE: %s ===", Cell_Type))
  
  for (Dataset in SC_DATASETS_FOR_CONCORDANCE) {
    
    tryCatch({
      message(sprintf("\n  Dataset: %s", Dataset))
      
      # --- A. LOAD SC DATA ---
      sc_file <- file.path(FULL_GT_BASE, Cell_Type, paste0(Dataset, "_Full_DESeq2_Results.qs"))
      
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
      XIST_ID <- "ENSG00000229807"
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
      
      # Filter for significant SC genes
      sc_sig <- sc_data %>% dplyr::filter(padj < P_VALUE_THRESHOLD)
      message(sprintf("    %d significant SC genes (padj < %.2f)", nrow(sc_sig), P_VALUE_THRESHOLD))
      
      # --- C. GENE ANNOTATIONS (PRE-FETCHED) ---
      # gene_coords was pre-fetched outside the loop in Section 2
      message(sprintf("    Using pre-fetched annotations (%d genes in lookup)", nrow(gene_coords)))
      
      # --- D. PREPARE DATA FOR PLOTTING ---
      lfc_palette <- scales::hue_pal()(length(LFC_INTERVALS))
      my_colors <- setNames(lfc_palette, names(LFC_INTERVALS))
      
      plot_data_sc <- sc_sig %>% 
        inner_join(gene_coords, by = c("gene"="ensembl_gene_id")) %>%
        mutate(chr_type = robust_chr_type) %>% 
        categorize_lfc() %>% 
        filter(lfc_group %in% names(my_colors))
      
      plot_data_bulk <- bulk_sig %>% 
        inner_join(gene_coords, by = c("gene"="ensembl_gene_id")) %>%
        mutate(chr_type = robust_chr_type) %>% 
        categorize_lfc() %>% 
        filter(lfc_group %in% names(my_colors))
      
      # --- E. PLOT COUNTS & PROPS ---
      chr_cols <- c("Autosome"="grey70", "X"="#f8766d", "Y"="#00bfc4")
      
      p1 <- ggplot(plot_data_sc, aes(x=lfc_group, fill=lfc_group)) + 
        geom_bar() + theme_publication() + scale_fill_manual(values=my_colors) +
        labs(title = paste0(prettify_label(Dataset), " — ", prettify_label(Cell_Type), " SC DEG Counts"), 
             x = "LFC Group", y = "Count")
      
      p2 <- ggplot(plot_data_sc, aes(x=lfc_group, fill=chr_type)) + 
        geom_bar(position="fill") + scale_y_continuous(labels=percent) + 
        theme_publication() + scale_fill_manual(values=chr_cols) +
        labs(title = paste0(prettify_label(Dataset), " — ", prettify_label(Cell_Type), " SC Chromosome Proportions"), 
             x = "LFC Group", y = "Proportion")
      
      p3 <- ggplot(plot_data_sc, aes(x=chr_type, fill=chr_type)) + 
        geom_bar() + facet_wrap(~lfc_group, scales="free_y") + 
        theme_publication() + scale_fill_manual(values=chr_cols) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste0(prettify_label(Dataset), " — ", prettify_label(Cell_Type), " SC Faceted by LFC"), 
             x = "Chromosome Type", y = "Count")
      
      p1_bulk <- ggplot(plot_data_bulk, aes(x=lfc_group, fill=lfc_group)) + 
        geom_bar() + theme_publication() + scale_fill_manual(values=my_colors) +
        labs(title = "Bulk DEG Counts", x = "LFC Group", y = "Count")
      
      p2_bulk <- ggplot(plot_data_bulk, aes(x=lfc_group, fill=chr_type)) + 
        geom_bar(position="fill") + scale_y_continuous(labels=percent) + 
        theme_publication() + scale_fill_manual(values=chr_cols) +
        labs(title = "Bulk Chromosome Proportions", x = "LFC Group", y = "Proportion")
      
      p3_bulk <- ggplot(plot_data_bulk, aes(x=chr_type, fill=chr_type)) + 
        geom_bar() + facet_wrap(~lfc_group, scales="free_y") + 
        theme_publication() + scale_fill_manual(values=chr_cols) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Bulk Faceted by LFC", x = "Chromosome Type", y = "Count")
      
      # --- F. CORRELATION PLOTS ---
      bulk_anchors <- bulk_data %>% 
        dplyr::filter(padj < P_VALUE_THRESHOLD & abs(log2FoldChange) > LFC_THRESHOLD) %>%
        dplyr::select(gene, log2FoldChange) %>% 
        dplyr::rename(log2FoldChange_Bulk = log2FoldChange)
      
      matched_df <- inner_join(bulk_anchors, sc_data, by = "gene") %>% 
        dplyr::rename(log2FoldChange_SC = log2FoldChange, padj_SC = padj)
      
      plot_data_comp <- matched_df %>% 
        inner_join(gene_coords, by = c("gene"="ensembl_gene_id")) %>%
        mutate(chr_type = robust_chr_type,
               label_gene = if (VIS_LABEL_GENES) {
                 case_when(
                   (!is.na(padj_SC) & padj_SC < P_VALUE_THRESHOLD & (sign(log2FoldChange_Bulk) != sign(log2FoldChange_SC))) ~ external_gene_name,
                   TRUE ~ NA_character_)
               } else { NA_character_ })
      
      p4 <- generate_corr_plot(plot_data_comp, "Bulk-Anchored Concordance", Cell_Type, Dataset)
      p4_zoom <- generate_corr_plot(filter(plot_data_comp, external_gene_name != "XIST"), "Bulk-Anchored (No XIST)", Cell_Type, Dataset)
      
      # --- G. RECOVERY STATS ---
      
      # Calculate shared/unique genes for this dataset
      if (!exists("sig_lists_by_dataset")) sig_lists_by_dataset <- list()
      if (is.null(sig_lists_by_dataset[[Dataset]])) sig_lists_by_dataset[[Dataset]] <- list()
      
      sig_lists_by_dataset[[Dataset]][[Cell_Type]] <- sc_sig$gene
      
      # Determine unique genes for this cell type in this dataset
      this_ct_sig <- sc_sig$gene
      other_ct_sig <- unique(unlist(sig_lists_by_dataset[[Dataset]][setdiff(CELL_TYPES, Cell_Type)]))
      unique_genes_for_ct <- setdiff(this_ct_sig, other_ct_sig)
      
      # Calculate shared genes across all cell types in this dataset
      if (length(sig_lists_by_dataset[[Dataset]]) > 1) {
        SHARED_SIG_GENES <- Reduce(intersect, sig_lists_by_dataset[[Dataset]])
      } else {
        SHARED_SIG_GENES <- c()
      }
      
      message(sprintf("    Found %d unique significant genes for %s", length(unique_genes_for_ct), Cell_Type))
      
      # SC Volcano (With Shared Translucency + Unique Labels)
      sc_volcano_df <- sc_data %>%
        inner_join(gene_coords, by = c("gene" = "ensembl_gene_id")) %>%
        mutate(external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "", gene, external_gene_name))
      
      # Pass unique list to trigger "Top 10 LFC Unique" logic
      p_sc_volcano <- generate_volcano_plot(sc_volcano_df, paste0(Dataset, " - ", Cell_Type), 
                                            LFC_THRESHOLD, P_VALUE_THRESHOLD, 
                                            shared_list = SHARED_SIG_GENES,
                                            unique_list = unique_genes_for_ct) 
      
      # Recovery Stats (4-category, direction-aware)
      sc_recovery_lookup <- sc_data %>%
        dplyr::select(gene, padj, log2FoldChange) %>%
        dplyr::rename(padj_SC = padj, lfc_SC = log2FoldChange)
      bulk_recovery_base <- bulk_sig %>% 
        filter(abs(log2FoldChange) > LFC_THRESHOLD) %>%
        left_join(gene_coords, by=c("gene"="ensembl_gene_id")) %>% 
        left_join(sc_recovery_lookup, by="gene") %>% 
        mutate(
          Recovery = factor(case_when(
            !(gene %in% sc_data$gene)                                          ~ "Absent",
            is.na(padj_SC)                                                     ~ "Low-Count Filtered",
            padj_SC < P_VALUE_THRESHOLD & sign(lfc_SC) == sign(log2FoldChange) ~ "Recovered",
            TRUE                                                               ~ "Insignificant"
          ), levels = RECOVERY_LEVELS),
          # Use pre-fetched chr_type; genes not in Ensembl 88 get Autosome as fallback
          robust_chr_type = ifelse(is.na(robust_chr_type), "Autosome", robust_chr_type),
          Chr_Cat = case_when(
            robust_chr_type == "Autosome" ~ "Autosome", 
            robust_chr_type %in% c("X", "Y") ~ "Sex Chr", 
            TRUE ~ "Autosome")
        ) %>% 
        categorize_lfc()
      
      recovery_summary <- bulk_recovery_base %>% 
        group_by(Chr_Cat, robust_chr_type, Recovery) %>% 
        summarise(Count=n(), .groups="drop")
      
      p5 <- ggplot(recovery_summary, aes(x=robust_chr_type, y=Count, fill=Recovery)) +
        geom_bar(stat="identity", position="stack") +
        geom_text(aes(label=Count), position=position_stack(vjust=0.5), color="black", size=4.5) +
        facet_grid(~Chr_Cat, scales="free_x", space="free") + 
        scale_fill_manual(values=RECOVERY_COLORS) + 
        theme_publication() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste0("Recovery of Bulk DEGs: ", prettify_label(Cell_Type)),
             x = "Chromosome", y = "Count")
      
      recovery_lfc_summary <- bulk_recovery_base %>% 
        group_by(robust_chr_type, lfc_group, Recovery) %>% 
        summarise(Count=n(), .groups="drop")
      
      p6 <- ggplot(recovery_lfc_summary, aes(x=lfc_group, y=Count, fill=Recovery)) + 
        geom_bar(stat="identity", position="stack") + 
        facet_wrap(~robust_chr_type, scales="fixed") + 
        geom_text(aes(label=Count), position=position_stack(vjust=0.5), color="black", size=3.5) + 
        scale_fill_manual(values=RECOVERY_COLORS) + 
        theme_publication() +
        labs(title = paste0("Recovery of Bulk DEGs: ", prettify_label(Cell_Type)),
             subtitle = "By LFC magnitude",
             x = "|LFC| Bin", y = "Count")
      
      # --- H. HIGH MAGNITUDE HISTOGRAM (COMBINED & STRATIFIED) ---
      
      # Combined
      hist_data <- sc_sig %>% 
        inner_join(gene_coords, by=c("gene"="ensembl_gene_id")) %>%
        filter(abs(log2FoldChange) > 2.0) %>%
        mutate(Magnitude = abs(log2FoldChange))
      
      if(nrow(hist_data) > 0) {
        mean_val <- mean(hist_data$Magnitude, na.rm=TRUE)
        med_val  <- median(hist_data$Magnitude, na.rm=TRUE)
        n_genes  <- nrow(hist_data)
        sub_str <- glue("N = {n_genes} | Mean |LFC|: {round(mean_val, 2)} | Median |LFC|: {round(med_val, 2)}")
        
        p7 <- ggplot(hist_data, aes(x=Magnitude)) +
          geom_histogram(binwidth=VIS_LFC_BIN_WIDTH, fill="#619CFF", color="black", alpha=0.7) +
          geom_vline(xintercept=mean_val, color="red", linetype="dashed", linewidth=1) +
          labs(title=glue("{prettify_label(Dataset)} — {prettify_label(Cell_Type)}: High Magnitude DEGs (|LFC| > 2)"),
               subtitle=sub_str, x="Absolute Log2 Fold Change", y="Count") +
          coord_cartesian(xlim=VIS_HIST_XLIM, ylim=VIS_HIST_YLIM) +
          theme_publication()
        
        # Stratified by Coding vs Non-Coding
        p7_strat <- ggplot(hist_data, aes(x=Magnitude, fill=Biotype_Class)) +
          geom_histogram(binwidth=VIS_LFC_BIN_WIDTH, color="black", alpha=0.7) +
          facet_wrap(~Biotype_Class, scales="fixed") +
          labs(title=glue("{prettify_label(Dataset)} — {prettify_label(Cell_Type)}: High Magnitude DEGs — Stratified"),
               x="Absolute Log2 Fold Change", y="Count") +
          coord_cartesian(xlim=VIS_HIST_XLIM) +
          theme_publication() + theme(legend.position="none")
        
      } else {
        p7 <- NULL
        p7_strat <- NULL
      }
      
      # --- I. SE HISTOGRAM (COMBINED & STRATIFIED) ---
      
      # Combined
      se_data <- sc_sig %>% 
        inner_join(gene_coords, by=c("gene"="ensembl_gene_id")) %>%
        filter(!is.na(lfcSE)) %>%
        filter(abs(log2FoldChange) > 2.0)
      
      if(nrow(se_data) > 0) {
        mean_se <- mean(se_data$lfcSE, na.rm=TRUE)
        med_se  <- median(se_data$lfcSE, na.rm=TRUE)
        n_se    <- nrow(se_data)
        sub_se <- glue("N = {n_se} | Mean SE: {round(mean_se, 3)} | Median SE: {round(med_se, 3)}")
        
        p8 <- ggplot(se_data, aes(x=lfcSE)) +
          geom_histogram(binwidth=VIS_SE_BIN_WIDTH, fill="#F8766D", color="black", alpha=0.7) +
          geom_vline(xintercept=mean_se, color="blue", linetype="dashed", linewidth=1) +
          labs(title=glue("{prettify_label(Dataset)} — {prettify_label(Cell_Type)}: DEG Standard Error Distribution"),
               subtitle=sub_se, x="lfcSE (Standard Error)", y="Count") +
          coord_cartesian(xlim=VIS_SE_XLIM, ylim=VIS_SE_YLIM) + 
          theme_publication()
        
        # Stratified by Coding vs Non-Coding
        p8_strat <- ggplot(se_data, aes(x=lfcSE, fill=Biotype_Class)) +
          geom_histogram(binwidth=VIS_SE_BIN_WIDTH, color="black", alpha=0.7) +
          facet_wrap(~Biotype_Class, scales="fixed") +
          labs(title=glue("{prettify_label(Dataset)} — {prettify_label(Cell_Type)}: DEG SE Distribution — Stratified"),
               x="lfcSE (Standard Error)", y="Count") +
          coord_cartesian(xlim=VIS_SE_XLIM) +
          theme_publication() + theme(legend.position="none")
        
      } else {
        p8 <- NULL
        p8_strat <- NULL
      }
      
      # --- J. SAVE ALL PLOTS ---
      OUT_DIR <- file.path(OUTPUT_DIR, Cell_Type)
      if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive=TRUE)
      
      message(sprintf("    Saving plots to: %s", OUT_DIR))
      
      # Prefix with dataset name if multiple datasets
      prefix <- if (length(SC_DATASETS_FOR_CONCORDANCE) > 1) paste0(Dataset, "_") else ""
      
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_counts_SC.png")), p1, width=6, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_props_SC.png")), p2, width=6, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_facets_SC.png")), p3, width=9, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_counts_Bulk.png")), p1_bulk, width=6, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_props_Bulk.png")), p2_bulk, width=6, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "plot_facets_Bulk.png")), p3_bulk, width=9, height=4)
      ggsave(file.path(OUT_DIR, paste0(prefix, "Compare_Correlation_BulkAnchored.png")), p4, width=7, height=6)
      ggsave(file.path(OUT_DIR, paste0(prefix, "Compare_Correlation_BulkAnchored_NoXIST.png")), p4_zoom, width=7, height=6)
      ggsave(file.path(OUT_DIR, paste0(prefix, "Recovery_Stats.png")), p5, width=9, height=5)
      ggsave(file.path(OUT_DIR, paste0(prefix, "Recovery_Stats_Stratified.png")), p6, width=12, height=6) 
      ggsave(file.path(OUT_DIR, paste0(prefix, "Volcano_Plot.png")), p_sc_volcano, width=8, height=6)
      
      if(!is.null(p7)) {
        ggsave(file.path(OUT_DIR, paste0(prefix, "High_Magnitude_LFC_Hist_Combined.png")), p7, width=7, height=5)
        ggsave(file.path(OUT_DIR, paste0(prefix, "High_Magnitude_LFC_Hist_Stratified.png")), p7_strat, width=9, height=5)
      }
      if(!is.null(p8)) {
        ggsave(file.path(OUT_DIR, paste0(prefix, "lfcSE_Distribution_Hist_Combined.png")), p8, width=7, height=5)
        ggsave(file.path(OUT_DIR, paste0(prefix, "lfcSE_Distribution_Hist_Stratified.png")), p8_strat, width=9, height=5)
      }
      
      message(sprintf("    Saved %d plots for %s - %s", 15, Dataset, Cell_Type))
      
    }, error = function(e) {
      message(sprintf("[ERROR] Processing %s - %s failed: %s", Dataset, Cell_Type, e$message))
    })
  }
}

# ==============================================================================
# 7. COMPLETION REPORT
# ==============================================================================

message("\n", strrep("=", 60))
message("STEP 05 COMPLETE")
message(strrep("=", 60))
message(sprintf("Generated bulk-vs-SC concordance plots for:"))
message(sprintf("  Cell types: %s", paste(CELL_TYPES, collapse = ", ")))
message(sprintf("  SC datasets: %s", paste(SC_DATASETS_FOR_CONCORDANCE, collapse = ", ")))
message(sprintf("  Output directory: %s", OUTPUT_DIR))
message("\n")
