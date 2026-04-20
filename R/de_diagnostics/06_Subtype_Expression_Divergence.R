# ==============================================================================
# SCRIPT: Subtype Expression Divergence — Per-Gene Split-Violin Figures
# PURPOSE: For each of the 23 shared truth genes, generate a single figure
#          showing all 3 cell types with split-violin geometry:
#            - Micro: 1 pooled violin (no subtypes, baseline)
#            - Oligo: 1 split violin (OPALIN | RBFOX1)
#            - Astro: 2 split violins (ADAMTSL3|GRIA1 + PLSCR1|WIF1)
#          When pseudobulking randomly samples cells regardless of subtype,
#          divergent subtype expression profiles introduce variance at low cell
#          counts — explaining power curve instability.
#
# INPUT:  - Shared truth genes from Power Curve Fitting pipeline
#           (Datalink/Power_Curve_Fitting/Fitting_Data/shared_truth_genes.qs)
#         - Donor shards for all 3 cell types (Datalink/Shards/{CT}/Roussos/)
# OUTPUT: - Per-gene PNG + PDF figures in Subtype_Expression_Divergence/
#         - Subtype_Expression_Summary.csv (per-gene × per-CT × per-subtype stats)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(scales)
  library(RColorBrewer)
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

log_section("STEP 6: SUBTYPE EXPRESSION DIVERGENCE")

# ==============================================================================
# 2. CREATE OUTPUT DIRECTORY
# ==============================================================================

output_dir <- get_subtype_div_output_dir()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  log_info("Created output directory: %s", output_dir)
}

# ==============================================================================
# 3. LOAD SHARED TRUTH GENES
# ==============================================================================

log_info("Loading shared truth genes from: %s", basename(SHARED_GENES_PATH))

if (!file.exists(SHARED_GENES_PATH)) {
  stop(sprintf("Shared truth genes file not found: %s\n  Run power_curve_fitting Step 01 first.",
               SHARED_GENES_PATH))
}

shared_genes <- qs::qread(SHARED_GENES_PATH)
shared_genes <- clean_gene_ids(shared_genes)
n_shared <- length(shared_genes)
log_info("Loaded %d shared truth genes", n_shared)

# ==============================================================================
# 4. EXTRACT DATA FOR ALL CELL TYPES
# ==============================================================================

# Helper function: extract per-cell UMI data from a single shard
# Returns a data.frame with columns: donor_id, gene_id, subtype, count
# For cell types with no subtypes (Micro), assigns subtype = "All"
extract_shard_data <- function(shard_path, expected_subtypes, subtype_col) {
  tryCatch({
    sce <- qs::qread(shard_path)
    cd <- as.data.frame(colData(sce))

    # Determine subtype assignments
    if (is.null(expected_subtypes)) {
      # No subtypes — assign all cells to "All"
      subtype_vec <- rep("All", nrow(cd))
      valid_cells <- rep(TRUE, nrow(cd))
    } else {
      subtype_vec <- as.character(cd[[subtype_col]])
      valid_cells <- subtype_vec %in% expected_subtypes
    }

    if (sum(valid_cells) == 0) return(NULL)

    # Get the count matrix (genes × cells)
    count_mat <- counts(sce)
    rownames(count_mat) <- sub("\\..*", "", rownames(count_mat))

    # Find shared genes present in this shard
    genes_present <- intersect(shared_genes, rownames(count_mat))
    if (length(genes_present) == 0) return(NULL)

    # Subset to shared genes and valid cells
    sub_mat <- count_mat[genes_present, valid_cells, drop = FALSE]
    sub_subtypes <- subtype_vec[valid_cells]
    donor_id <- as.character(cd[[DONOR_ID_COL]][1])

    # Convert to long-format data.frame
    rows <- vector("list", length(genes_present))
    for (i in seq_along(genes_present)) {
      gene <- genes_present[i]
      umi_counts <- as.numeric(sub_mat[gene, ])
      rows[[i]] <- data.frame(
        donor_id = donor_id,
        gene_id  = gene,
        subtype  = sub_subtypes,
        count    = umi_counts,
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, rows)

  }, error = function(e) {
    message(sprintf("[ERROR] Failed to process %s: %s",
                    basename(shard_path), e$message))
    return(NULL)
  })
}

# Loop over cell types and extract data
all_ct_data <- list()

for (ct in SUBTYPE_DIV_CT_ORDER) {
  log_info("")
  log_info("=== Processing %s ===", ct)

  ct_config <- SUBTYPE_DIV_CT_CONFIG[[ct]]
  shard_dir <- get_shard_dir(ct)

  if (!dir.exists(shard_dir)) {
    log_warn("Shard directory not found for %s: %s — skipping", ct, shard_dir)
    next
  }

  # List all donor shard files
  all_shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$",
                                 full.names = TRUE)
  n_total <- length(all_shard_files)
  log_info("Found %d donor shards in %s/%s", n_total, ct, "Roussos")

  if (n_total == 0) {
    log_warn("No shards found for %s — skipping", ct)
    next
  }

  # Sample donors reproducibly (same seed per CT for independence,
  # but offset by CT index for different samples)
  ct_idx <- which(SUBTYPE_DIV_CT_ORDER == ct)
  set.seed(SUBTYPE_DIV_SEED + ct_idx - 1)
  n_sample <- min(SUBTYPE_DIV_N_DONORS, n_total)
  sampled_shards <- sort(sample(all_shard_files, n_sample))
  log_info("Sampled %d donors (seed=%d)", n_sample, SUBTYPE_DIV_SEED + ct_idx - 1)

  # Parallel extraction
  extract_start <- Sys.time()
  log_info("Extracting per-cell UMI counts (%d cores)...", SUBTYPE_DIV_N_CORES)

  shard_results <- mclapply(
    sampled_shards,
    extract_shard_data,
    expected_subtypes = ct_config$subtypes,
    subtype_col = SUBTYPE_DIV_SUBTYPE_COL,
    mc.cores = SUBTYPE_DIV_N_CORES
  )

  shard_results <- Filter(Negate(is.null), shard_results)
  n_donors_ok <- length(shard_results)

  if (n_donors_ok == 0) {
    log_warn("No data extracted for %s — skipping", ct)
    next
  }

  ct_cell_data <- do.call(rbind, shard_results)
  ct_cell_data$cell_type <- ct

  extract_duration <- as.numeric(difftime(Sys.time(), extract_start, units = "secs"))
  log_info("Extraction complete: %d donors, %s cell-gene observations (%.1f sec)",
           n_donors_ok, format(nrow(ct_cell_data), big.mark = ","), extract_duration)

  # Report per-subtype counts
  subtype_counts <- table(ct_cell_data$subtype)
  for (st in names(subtype_counts)) {
    n_cells_approx <- subtype_counts[st] / n_shared
    log_info("  %s: ~%s cells", st, format(round(n_cells_approx), big.mark = ","))
  }

  all_ct_data[[ct]] <- ct_cell_data
}

# Combine all cell types
cell_data <- do.call(rbind, all_ct_data)
log_info("")
log_info("Combined data: %s cell-gene observations across %d cell types",
         format(nrow(cell_data), big.mark = ","), length(all_ct_data))

# Check gene coverage
genes_found <- unique(cell_data$gene_id)
n_genes_found <- length(genes_found)
log_info("Genes with data: %d/%d", n_genes_found, n_shared)

# ==============================================================================
# 5. COMPUTE SUMMARY STATISTICS
# ==============================================================================

log_info("Computing per-gene per-CT per-subtype summary statistics...")

summary_df <- cell_data %>%
  group_by(cell_type, gene_id, subtype) %>%
  summarise(
    n_cells      = n(),
    mean_count   = mean(count),
    median_count = median(count),
    sd_count     = sd(count),
    pct_zero     = 100 * mean(count == 0),
    q25          = quantile(count, 0.25),
    q75          = quantile(count, 0.75),
    max_count    = max(count),
    .groups = "drop"
  ) %>%
  arrange(cell_type, gene_id, subtype)

# Compute statistical tests per gene per cell type
test_results <- list()

for (ct in names(all_ct_data)) {
  ct_config <- SUBTYPE_DIV_CT_CONFIG[[ct]]
  ct_subset <- cell_data[cell_data$cell_type == ct, ]

  for (gene in genes_found) {
    gene_data <- ct_subset[ct_subset$gene_id == gene, ]

    p_val <- NA_real_
    test_type <- ct_config$test

    if (test_type == "wilcoxon" && length(ct_config$subtypes) == 2) {
      vals_a <- gene_data$count[gene_data$subtype == ct_config$subtypes[1]]
      vals_b <- gene_data$count[gene_data$subtype == ct_config$subtypes[2]]
      if (length(vals_a) >= 3 && length(vals_b) >= 3) {
        p_val <- tryCatch(
          wilcox.test(vals_a, vals_b, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
    } else if (test_type == "kruskal" && length(ct_config$subtypes) >= 3) {
      # Kruskal-Wallis omnibus test
      if (length(unique(gene_data$subtype)) >= 2) {
        p_val <- tryCatch(
          kruskal.test(count ~ subtype, data = gene_data)$p.value,
          error = function(e) NA_real_
        )
      }
    }
    # test_type == "none" (Micro) → p_val stays NA

    test_results[[length(test_results) + 1]] <- data.frame(
      cell_type = ct,
      gene_id   = gene,
      test_type = test_type,
      p_value   = p_val,
      stringsAsFactors = FALSE
    )
  }
}

test_df <- do.call(rbind, test_results)

# Merge test results into summary
summary_df <- summary_df %>%
  left_join(test_df, by = c("cell_type", "gene_id"))

# Save summary CSV
summary_path <- file.path(output_dir, "Subtype_Expression_Summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
log_info("Saved: %s (%d rows)", basename(summary_path), nrow(summary_df))

# ==============================================================================
# 6. COMPUTE DONOR-LEVEL MEANS (for overlay points)
# ==============================================================================

log_info("Computing donor-level means for overlay...")

donor_means <- cell_data %>%
  group_by(cell_type, gene_id, subtype, donor_id) %>%
  summarise(
    donor_mean = mean(count),
    .groups = "drop"
  )

log_info("Donor-level means: %s rows", format(nrow(donor_means), big.mark = ","))

# ==============================================================================
# 7. LOOK UP GENE SYMBOLS
# ==============================================================================

log_info("Looking up gene symbols for human-readable filenames...")

gene_symbols <- get_gene_symbols(genes_found, SUBTYPE_DIV_SYMBOL_CACHE)
log_info("Gene symbols resolved: %d/%d", sum(gene_symbols != genes_found), n_genes_found)

# ==============================================================================
# 8. BUILD SUBTYPE LABEL MAPPING (BUG FIX: use unname to avoid name mangling)
# ==============================================================================

# FIXED: unlist() on a named list prepends list-element names (e.g.,
# "Oligo.Oligo_OPALIN" instead of "Oligo_OPALIN"). Using unname() on the
# outer list prevents this, preserving the original inner names.
all_labels <- do.call(c, unname(lapply(SUBTYPE_DIV_CT_CONFIG, function(x) x$labels)))

log_info("Subtype label mapping: %s",
         paste(sprintf("%s->%s", names(all_labels), all_labels), collapse = ", "))

# ==============================================================================
# 9. BUILD COLOUR PALETTE
# ==============================================================================

# Collect all unique subtype keys across CTs
all_subtype_keys <- names(all_labels)
n_colours_needed <- length(all_subtype_keys)

# Use a qualitative palette with enough colours
if (n_colours_needed <= 8) {
  palette_colours <- brewer.pal(max(3, n_colours_needed), "Set2")[seq_len(n_colours_needed)]
} else {
  palette_colours <- rainbow(n_colours_needed, s = 0.7, v = 0.8)
}
names(palette_colours) <- all_labels[all_subtype_keys]

log_info("Colour palette: %d subtypes (%s)",
         n_colours_needed, paste(all_labels[all_subtype_keys], collapse = ", "))

# ==============================================================================
# 10. PREPARE P-VALUE ANNOTATIONS
# ==============================================================================

p_annot <- test_df %>%
  mutate(
    p_label = case_when(
      is.na(p_value)    ~ "",
      test_type == "none" ~ "",
      p_value < 1e-300  ~ "p < 1e-300",
      p_value < 0.001   ~ sprintf("p = %.1e", p_value),
      TRUE              ~ sprintf("p = %.3f", p_value)
    )
  )

# ==============================================================================
# 11. ORDER GENES BY MAX SUBTYPE DIVERGENCE
# ==============================================================================

gene_divergence <- summary_df %>%
  filter(cell_type != "Micro") %>%
  group_by(cell_type, gene_id) %>%
  summarise(
    mean_range = max(mean_count) - min(mean_count),
    .groups = "drop"
  ) %>%
  group_by(gene_id) %>%
  summarise(
    max_divergence = max(mean_range),
    .groups = "drop"
  ) %>%
  arrange(desc(max_divergence))

gene_order <- gene_divergence$gene_id

log_info("Gene order (by max subtype divergence): %s",
         paste(gene_symbols[gene_order[1:min(5, length(gene_order))]], collapse = ", "))

# ==============================================================================
# 12. GENERATE PER-GENE SPLIT-VIOLIN FIGURES
# ==============================================================================

log_info("")
log_info("Generating %d per-gene split-violin figures...", n_genes_found)

# Filter to non-zero counts for violin plotting
cell_data_nz <- cell_data[cell_data$count > 0, ]

# X-position labels (ordered as in config)
xpos_labels <- sapply(SUBTYPE_DIV_XPOS, function(x) x$label)

# Helper: build plot data for one gene across all x-positions
build_gene_figure <- function(gene_id) {

  gene_symbol <- gene_symbols[gene_id]
  title_text <- sprintf("%s (%s)", gene_symbol, gene_id)

  # Collect data for each x-position
  plot_rows <- list()
  donor_rows <- list()
  annot_rows <- list()

  for (pos_idx in seq_along(SUBTYPE_DIV_XPOS)) {
    pos <- SUBTYPE_DIV_XPOS[[pos_idx]]
    pos_label <- pos$label
    ct <- pos$ct

    # Get cell data for this gene + cell type
    ct_gene_nz <- cell_data_nz[cell_data_nz$gene_id == gene_id &
                                 cell_data_nz$cell_type == ct, ]
    ct_donor <- donor_means[donor_means$gene_id == gene_id &
                              donor_means$cell_type == ct, ]

    if (pos$type == "pooled") {
      # Micro: single pooled violin — use all cells, no split
      if (nrow(ct_gene_nz) > 0) {
        df_pos <- ct_gene_nz
        df_pos$xpos <- pos_label
        df_pos$split_group <- "All Micro"
        plot_rows[[length(plot_rows) + 1]] <- df_pos
      }
      if (nrow(ct_donor) > 0) {
        dm_pos <- ct_donor
        dm_pos$xpos <- pos_label
        dm_pos$split_group <- "All Micro"
        donor_rows[[length(donor_rows) + 1]] <- dm_pos
      }

      # Annotation: mean + dropout% for pooled
      ct_summary <- summary_df[summary_df$gene_id == gene_id &
                                 summary_df$cell_type == ct, ]
      if (nrow(ct_summary) > 0) {
        annot_rows[[length(annot_rows) + 1]] <- data.frame(
          xpos = pos_label,
          label = sprintf("\u03bc=%.2f\ndrop=%.0f%%",
                          ct_summary$mean_count[1], ct_summary$pct_zero[1]),
          stringsAsFactors = FALSE
        )
      }

    } else if (pos$type == "split") {
      # Split violin: filter to the 2 subtypes for this position
      st_a <- pos$subtypes[1]
      st_b <- pos$subtypes[2]

      for (st in c(st_a, st_b)) {
        st_nz <- ct_gene_nz[ct_gene_nz$subtype == st, ]
        if (nrow(st_nz) > 0) {
          df_st <- st_nz
          df_st$xpos <- pos_label
          df_st$split_group <- all_labels[st]
          plot_rows[[length(plot_rows) + 1]] <- df_st
        }

        st_dm <- ct_donor[ct_donor$subtype == st, ]
        if (nrow(st_dm) > 0) {
          dm_st <- st_dm
          dm_st$xpos <- pos_label
          dm_st$split_group <- all_labels[st]
          donor_rows[[length(donor_rows) + 1]] <- dm_st
        }
      }

      # Annotation: mean + dropout% per subtype
      st_a_summary <- summary_df[summary_df$gene_id == gene_id &
                                   summary_df$cell_type == ct &
                                   summary_df$subtype == st_a, ]
      st_b_summary <- summary_df[summary_df$gene_id == gene_id &
                                   summary_df$cell_type == ct &
                                   summary_df$subtype == st_b, ]

      label_parts <- c()
      if (nrow(st_a_summary) > 0) {
        label_parts <- c(label_parts, sprintf(
          "%s: \u03bc=%.2f, drop=%.0f%%",
          all_labels[st_a], st_a_summary$mean_count[1], st_a_summary$pct_zero[1]
        ))
      }
      if (nrow(st_b_summary) > 0) {
        label_parts <- c(label_parts, sprintf(
          "%s: \u03bc=%.2f, drop=%.0f%%",
          all_labels[st_b], st_b_summary$mean_count[1], st_b_summary$pct_zero[1]
        ))
      }

      if (length(label_parts) > 0) {
        annot_rows[[length(annot_rows) + 1]] <- data.frame(
          xpos = pos_label,
          label = paste(label_parts, collapse = "\n"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(plot_rows) == 0) {
    log_warn("No non-zero data for gene %s — skipping figure", gene_id)
    return(NULL)
  }

  plot_df <- do.call(rbind, plot_rows)
  donor_df <- if (length(donor_rows) > 0) do.call(rbind, donor_rows) else NULL
  annot_df <- if (length(annot_rows) > 0) do.call(rbind, annot_rows) else NULL

  # Set x-position as ordered factor
  plot_df$xpos <- factor(plot_df$xpos, levels = xpos_labels)
  if (!is.null(donor_df)) {
    donor_df$xpos <- factor(donor_df$xpos, levels = xpos_labels)
  }
  if (!is.null(annot_df)) {
    annot_df$xpos <- factor(annot_df$xpos, levels = xpos_labels)
  }

  # Force split_group factor levels to ensure deterministic left/right ordering

  # in geom_split_violin. The first level per x-position gets the odd group
  # (drawn on the left), the second gets the even group (drawn on the right).
  # We collect all subtype labels in the order they appear in SUBTYPE_DIV_XPOS
  # (left subtype first, right subtype second per position).
  split_level_order <- c()
  for (pos in SUBTYPE_DIV_XPOS) {
    if (pos$type == "pooled") {
      pooled_label <- SUBTYPE_DIV_CT_CONFIG[[pos$ct]]$labels
      split_level_order <- c(split_level_order, unname(pooled_label))
    } else {
      split_level_order <- c(split_level_order,
                              all_labels[pos$subtypes[1]],
                              all_labels[pos$subtypes[2]])
    }
  }
  split_level_order <- unique(split_level_order)
  plot_df$split_group <- factor(plot_df$split_group, levels = split_level_order)
  if (!is.null(donor_df)) {
    donor_df$split_group <- factor(donor_df$split_group, levels = split_level_order)
  }

  # --- Build p-value annotation for the CT-level tests ---
  # Oligo gets one p-value, Astro gets one p-value (omnibus), Micro gets none
  ct_pvals <- p_annot[p_annot$gene_id == gene_id & nchar(p_annot$p_label) > 0, ]

  # Place Oligo p-value at the Oligo x-position, Astro at Astro1 position
  pval_annot_rows <- list()
  for (i in seq_len(nrow(ct_pvals))) {
    ct_name <- ct_pvals$cell_type[i]
    p_text <- ct_pvals$p_label[i]
    # Find the first x-position for this CT
    for (pos in SUBTYPE_DIV_XPOS) {
      if (pos$ct == ct_name && pos$type == "split") {
        pval_annot_rows[[length(pval_annot_rows) + 1]] <- data.frame(
          xpos = pos$label,
          p_label = p_text,
          stringsAsFactors = FALSE
        )
        break
      }
    }
  }
  pval_df <- if (length(pval_annot_rows) > 0) do.call(rbind, pval_annot_rows) else NULL
  if (!is.null(pval_df)) {
    pval_df$xpos <- factor(pval_df$xpos, levels = xpos_labels)
  }

  # --- Determine which positions are pooled vs split ---
  pooled_positions <- sapply(SUBTYPE_DIV_XPOS, function(x) x$type == "pooled")
  pooled_labels <- xpos_labels[pooled_positions]
  split_labels <- xpos_labels[!pooled_positions]

  # Separate data for pooled vs split positions
  df_pooled <- plot_df[plot_df$xpos %in% pooled_labels, ]
  df_split <- plot_df[plot_df$xpos %in% split_labels, ]

  # --- Build the plot ---
  p <- ggplot()

  # Layer 1: Regular violins for pooled positions (Micro)
  if (nrow(df_pooled) > 0) {
    p <- p + geom_violin(
      data = df_pooled,
      aes(x = xpos, y = count, fill = split_group),
      alpha = SUBTYPE_DIV_VIOLIN_ALPHA,
      scale = SUBTYPE_DIV_VIOLIN_SCALE,
      trim = TRUE,
      linewidth = 0.3
    )
  }

  # Layer 2: Split violins for split positions (Oligo, Astro1, Astro2)
  if (nrow(df_split) > 0) {
    p <- p + geom_split_violin(
      data = df_split,
      aes(x = xpos, y = count, fill = split_group),
      alpha = SUBTYPE_DIV_VIOLIN_ALPHA,
      scale = SUBTYPE_DIV_VIOLIN_SCALE,
      trim = TRUE
    )
  }

  # Layer 3: Donor-level mean overlay
  # Use pmax to avoid log(0) — donors with mean=0 are placed at 0.5
  if (!is.null(donor_df) && nrow(donor_df) > 0) {
    donor_df$donor_mean_safe <- pmax(donor_df$donor_mean, 0.5)
    p <- p + geom_jitter(
      data = donor_df,
      aes(x = xpos, y = donor_mean_safe),
      size = SUBTYPE_DIV_POINT_SIZE,
      alpha = SUBTYPE_DIV_POINT_ALPHA,
      width = 0.15,
      colour = "grey20",
      shape = 16
    )
  }

  # Layer 4: Annotation text (mean + dropout% per subtype)
  if (!is.null(annot_df) && nrow(annot_df) > 0) {
    p <- p + geom_text(
      data = annot_df,
      aes(x = xpos, y = Inf, label = label),
      hjust = 0.5,
      vjust = 1.1,
      size = SUBTYPE_DIV_ANNOT_SIZE,
      colour = "grey30",
      lineheight = 0.85
    )
  }

  # Layer 5: P-value annotation (placed below annotation text, above violins)
  if (!is.null(pval_df) && nrow(pval_df) > 0) {
    p <- p + geom_label(
      data = pval_df,
      aes(x = xpos, y = Inf, label = p_label),
      hjust = 0.5,
      vjust = 2.8,
      size = SUBTYPE_DIV_PVAL_SIZE,
      colour = "firebrick",
      fontface = "bold",
      fill = "white",
      alpha = 0.8,
      label.size = 0.2,
      label.padding = unit(0.15, "lines")
    )
  }

  # Scales and theme
  p <- p +
    scale_y_log10(labels = comma) +
    scale_fill_manual(values = palette_colours, drop = FALSE) +
    labs(
      title = title_text,
      subtitle = sprintf(
        "%d donors per CT | Non-zero UMI counts (log10) | Oligo: Wilcoxon | Astro: Kruskal-Wallis",
        SUBTYPE_DIV_N_DONORS
      ),
      x = NULL,
      y = "Raw UMI Count (log10, non-zero only)",
      fill = "Subtype"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0.5),
      plot.margin = margin(t = 25, r = 10, b = 5, l = 10, unit = "pt"),
      axis.text.x = element_text(size = 9, lineheight = 0.9),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# --- Loop over genes and save figures ---
n_saved <- 0

for (gene_id in gene_order) {

  gene_symbol <- gene_symbols[gene_id]
  safe_symbol <- gsub("[^[:alnum:]_-]", "", gene_symbol)
  file_base <- sprintf("%s_%s", safe_symbol, gene_id)

  log_info("  [%d/%d] %s (%s)...", which(gene_order == gene_id), n_genes_found,
           gene_symbol, gene_id)

  p <- tryCatch({
    build_gene_figure(gene_id)
  }, error = function(e) {
    log_error("Failed to build figure for %s: %s", gene_id, e$message)
    NULL
  })

  if (is.null(p)) next

  # Save PNG
  png_path <- file.path(output_dir, paste0(file_base, ".png"))
  tryCatch({
    ggsave(png_path, p,
           width = SUBTYPE_DIV_FIG_WIDTH, height = SUBTYPE_DIV_FIG_HEIGHT,
           dpi = SUBTYPE_DIV_DPI, limitsize = FALSE)
  }, error = function(e) {
    log_error("Failed to save PNG for %s: %s", gene_id, e$message)
  })

  # Save PDF
  pdf_path <- file.path(output_dir, paste0(file_base, ".pdf"))
  tryCatch({
    ggsave(pdf_path, p,
           width = SUBTYPE_DIV_FIG_WIDTH, height = SUBTYPE_DIV_FIG_HEIGHT,
           limitsize = FALSE)
  }, error = function(e) {
    log_error("Failed to save PDF for %s: %s", gene_id, e$message)
  })

  n_saved <- n_saved + 1
}

# ==============================================================================
# 13. SUMMARY
# ==============================================================================

log_section("STEP 6 COMPLETE")
log_info("Cell types analysed: %s", paste(names(all_ct_data), collapse = ", "))
log_info("Donors sampled per CT: %d", SUBTYPE_DIV_N_DONORS)
log_info("Genes analysed: %d/%d shared truth genes", n_genes_found, n_shared)
log_info("Figures saved: %d (PNG + PDF each)", n_saved)
log_info("Total cell-gene observations: %s",
         format(nrow(cell_data), big.mark = ","))

# Print key findings per CT
log_info("")
log_info("=== KEY FINDINGS ===")

for (ct in names(all_ct_data)) {
  ct_config <- SUBTYPE_DIV_CT_CONFIG[[ct]]
  ct_tests <- test_df[test_df$cell_type == ct, ]

  if (ct_config$test == "none") {
    log_info("")
    log_info("%s: No subtypes — baseline (single pooled violin)", ct)
    # Report total cells
    ct_cells <- nrow(all_ct_data[[ct]]) / n_genes_found
    log_info("  ~%s cells total", format(round(ct_cells), big.mark = ","))
    next
  }

  sig_genes <- ct_tests[!is.na(ct_tests$p_value) & ct_tests$p_value < 0.05, ]
  log_info("")
  log_info("%s (%s test): %d/%d genes with significant subtype difference (p < 0.05)",
           ct, ct_config$test, nrow(sig_genes), n_genes_found)

  if (nrow(sig_genes) > 0) {
    top_n <- min(5, nrow(sig_genes))
    sig_genes <- sig_genes[order(sig_genes$p_value), ]
    log_info("  Top %d most significant:", top_n)
    for (i in seq_len(top_n)) {
      gene <- sig_genes$gene_id[i]
      p_val <- sig_genes$p_value[i]
      log_info("    %s (%s): p=%.1e", gene_symbols[gene], gene, p_val)
    }
  }

  # Report per-subtype cell counts
  ct_subset <- cell_data[cell_data$cell_type == ct, ]
  st_counts <- table(ct_subset$subtype)
  for (st in names(st_counts)) {
    n_cells_approx <- st_counts[st] / n_genes_found
    label <- if (st %in% names(all_labels)) all_labels[st] else st
    log_info("  %s: ~%s cells", label, format(round(n_cells_approx), big.mark = ","))
  }
}

log_info("")
log_info("Files saved:")
log_info("  %d PNG files", n_saved)
log_info("  %d PDF files", n_saved)
log_info("  %s", basename(summary_path))
log_info("Output directory: %s", output_dir)
