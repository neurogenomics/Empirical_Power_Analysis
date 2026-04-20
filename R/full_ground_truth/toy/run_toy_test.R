#!/usr/bin/env Rscript
# ==============================================================================
# STAGE 2: TOY INTEGRATION TEST
# Full Ground Truth Pipeline - Minimal End-to-End Run
# ==============================================================================
# Purpose: Verify entire DE + diagnostics + summary pipeline works on minimal
#          real data (2 datasets, 5 donors each, 1 cell type) before full run.
#
# Configuration:
#   Cell types:  1 (Micro only)
#   Datasets:    2 (Roussos + Allen)
#   Donors:      5 per dataset (first 5 shard files found)
#   Output:      {OUTPUT_BASE}/TOY_Micro/ (separate from production)
#   Max cores:   4 (for both AGG and DE)
#
# Run with: Rscript R/full_ground_truth/toy/run_toy_test.R
# ==============================================================================

toy_start_time <- Sys.time()

cat("\n========================================\n")
cat("STAGE 2: TOY INTEGRATION TEST\n")
cat("Full Ground Truth Pipeline\n")
cat("========================================\n\n")

# ==============================================================================
# 0. SOURCE CONFIG & OVERRIDE FOR TOY MODE
# ==============================================================================

# Source config first (defines PROJECT_ROOT, OUTPUT_BASE, etc.)
source(file.path(
  PROJECT_ROOT,
  "R", "full_ground_truth", "config.R"
))

# --- Toy overrides ---
original_OUTPUT_BASE <- OUTPUT_BASE
CELL_TYPES    <- "Micro"
RUN_DATASETS  <- c("Roussos", "Allen")
N_CORES_AGG   <- 4L
N_CORES_DE    <- 4L
OUTPUT_BASE   <- file.path(original_OUTPUT_BASE, "TOY_Micro")
FORCE_RECALC  <- TRUE
TOY_N_DONORS  <- 5L

log_info("=== TOY MODE OVERRIDES ACTIVE ===")
log_info("Cell types:  %s", paste(CELL_TYPES, collapse = ", "))
log_info("Datasets:    %s", paste(RUN_DATASETS, collapse = ", "))
log_info("Donors/dataset: %d (symlinked)", TOY_N_DONORS)
log_info("Cores (AGG/DE): %d / %d", N_CORES_AGG, N_CORES_DE)
log_info("Output base: %s", OUTPUT_BASE)
log_info("Force recalc: %s", FORCE_RECALC)

# Source utils AFTER overrides (so log functions exist)
source(file.path(SCRIPT_DIR, "utils.R"))

# ==============================================================================
# 0b. LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(qs)
  library(DESeq2)
  library(apeglm)
  library(Matrix)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(viridis)
  library(UpSetR)
  library(RColorBrewer)
})

# ==============================================================================
# 0c. VALIDATION TRACKER
# ==============================================================================

validation_results <- list()

validate <- function(description, condition) {
  status <- if (isTRUE(condition)) "PASS" else "FAIL"
  validation_results[[length(validation_results) + 1]] <<- list(
    description = description,
    status      = status
  )
  if (status == "PASS") {
    log_info("[PASS] %s", description)
  } else {
    log_warn("[FAIL] %s", description)
  }
  invisible(status == "PASS")
}

# ==============================================================================
# 0d. CREATE TOY OUTPUT DIRECTORY
# ==============================================================================

paths <- get_paths(CELL_TYPES[1], output_base = OUTPUT_BASE)
output_dir <- paths$output_dir
diag_dir   <- paths$diagnostic_dir

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(diag_dir))   dir.create(diag_dir, recursive = TRUE)

log_info("Output dir:     %s", output_dir)
log_info("Diagnostic dir: %s", diag_dir)

# ==============================================================================
# STEP 1: DE ANALYSIS (INLINE)
# ==============================================================================

log_section("STEP 1: DE Analysis (Toy - 2 datasets x 5 donors)")

de_param <- setup_parallel(N_CORES_DE)

# Temp dir for symlinked shards (cleaned up at end)
toy_temp_base <- tempfile(pattern = "toy_shards_")
dir.create(toy_temp_base, recursive = TRUE)
log_info("Temp shard base: %s", toy_temp_base)

for (d_name in RUN_DATASETS) {

  log_section(sprintf("Processing: %s / %s (TOY)", d_name, CELL_TYPES[1]))

  cfg <- DATASET_CONFIG[[d_name]]
  real_shard_dir <- get_shard_dir(CELL_TYPES[1], d_name)

  # ------------------------------------------------------------------
  # 1a. Check real shard directory exists
  # ------------------------------------------------------------------
  if (!dir.exists(real_shard_dir)) {
    log_error("Shard directory not found: %s", real_shard_dir)
    log_error("Cannot run toy test without real shard data.")
    next
  }

  # ------------------------------------------------------------------
  # 1b. Create temp dir with symlinks to first 5 shards
  # ------------------------------------------------------------------
  all_shard_files <- list.files(real_shard_dir, pattern = "\\.qs$", full.names = TRUE)
  all_shard_files <- all_shard_files[!grepl("MANIFEST", all_shard_files, ignore.case = TRUE)]
  all_shard_files <- all_shard_files[!grepl("donor_NA\\.qs$", all_shard_files)]
  all_shard_files <- sort(all_shard_files)  # deterministic order

  if (length(all_shard_files) < TOY_N_DONORS) {
    log_warn("Only %d shards found for %s (need %d). Using all available.",
             length(all_shard_files), d_name, TOY_N_DONORS)
    selected_files <- all_shard_files
  } else {
    selected_files <- all_shard_files[seq_len(TOY_N_DONORS)]
  }

  toy_shard_dir <- file.path(toy_temp_base, d_name)
  dir.create(toy_shard_dir, recursive = TRUE)

  for (sf in selected_files) {
    link_target <- file.path(toy_shard_dir, basename(sf))
    file.symlink(sf, link_target)
  }

  log_info("Symlinked %d shard files to: %s", length(selected_files), toy_shard_dir)

  # ------------------------------------------------------------------
  # 1c. Aggregate pseudobulk
  # ------------------------------------------------------------------
  tryCatch({
    vars_requested <- unique(c(
      cfg$cols$donor, cfg$cols$sex, cfg$cols$age,
      cfg$covariates, cfg$numeric_covariates
    ))

    pb_sce <- aggregate_shards_parallel(
      shard_dir     = toy_shard_dir,
      design_vars   = vars_requested,
      n_cores       = N_CORES_AGG
    )

    log_info("Pseudobulk: %d donors, %d genes", ncol(pb_sce), nrow(pb_sce))

    # ------------------------------------------------------------------
    # 1d. Run DESeq2
    # ------------------------------------------------------------------
    de_result <- run_deseq_for_dataset(
      pb_sce         = pb_sce,
      cfg            = cfg,
      de_param       = de_param,
      min_gene_count = MIN_GENE_COUNT
    )

    log_info("DESeq2 complete: %d genes tested, %d significant (padj < %g)",
             nrow(de_result$res_df),
             sum(de_result$res_df$padj < FDR_THRESHOLD, na.rm = TRUE),
             FDR_THRESHOLD)

    # ------------------------------------------------------------------
    # 1e. Save outputs
    # ------------------------------------------------------------------
    save_de_outputs(
      res_df       = de_result$res_df,
      dds          = de_result$dds,
      res          = de_result$results,
      pb_sce       = pb_sce,
      output_dir   = output_dir,
      dataset_name = d_name,
      cell_type    = CELL_TYPES[1]
    )

    log_info("[SUCCESS] %s / %s complete", d_name, CELL_TYPES[1])

  }, error = function(e) {
    log_error("[FAILED] %s / %s: %s", d_name, CELL_TYPES[1], e$message)
  })

  # Clean up memory
  if (exists("pb_sce", inherits = FALSE))   rm(pb_sce)
  if (exists("de_result", inherits = FALSE)) rm(de_result)
  gc(verbose = FALSE)
}

# ==============================================================================
# STEP 2: DIAGNOSTICS (INLINE)
# ==============================================================================

log_section("STEP 2: Cross-Dataset Diagnostics (Toy)")

tryCatch({

  all_results <- load_all_de_results(output_dir)
  log_info("Loaded %d dataset results: %s",
           length(all_results), paste(names(all_results), collapse = ", "))

  if (length(all_results) < 2) {
    log_warn("Need >= 2 datasets for diagnostics. Only %d found.", length(all_results))
  } else {

    # --- Extract significant genes ---
    sig_list <- extract_significant_genes(all_results,
                                          padj_cutoff = FDR_THRESHOLD,
                                          lfc_cutoff  = LFC_THRESHOLD)
    for (ds in names(sig_list)) {
      log_info("  %s: %d significant genes", ds, length(sig_list[[ds]]))
    }

    union_sig <- unique(unlist(sig_list))
    n_sets    <- length(sig_list)

    # --- Diagnostic 1: UpSet Plot ---
    log_info("[1/4] UpSet Plot...")
    tryCatch({
      if (length(union_sig) > 0 && n_sets > 1) {
        binary_mat <- matrix(0L, nrow = length(union_sig), ncol = n_sets)
        rownames(binary_mat) <- union_sig
        colnames(binary_mat) <- names(sig_list)
        for (nm in names(sig_list)) {
          binary_mat[sig_list[[nm]], nm] <- 1L
        }

        png(file.path(diag_dir, "Diagnostic_UpSet_Plot.png"),
            width = 1200, height = 800, res = 150)
        print(upset(as.data.frame(binary_mat),
                    sets = colnames(binary_mat), order.by = "freq",
                    mainbar.y.label = "Intersection Size",
                    sets.x.label = "Set Size"))
        dev.off()
        log_info("  Saved Diagnostic_UpSet_Plot.png")
      } else {
        log_warn("  Not enough significant genes for UpSet plot")
      }
    }, error = function(e) {
      log_error("  UpSet plot failed: %s", e$message)
      try(dev.off(), silent = TRUE)
    })

    # --- Diagnostic 2: LFC Correlation ---
    log_info("[2/4] LFC Correlation...")
    ref_name <- grep(paste0("^", REFERENCE_DATASET),
                     names(all_results), value = TRUE, ignore.case = TRUE)[1]

    tryCatch({
      if (!is.na(ref_name)) {
        ref_df <- all_results[[ref_name]]
        ref_sig <- ref_df[!is.na(ref_df$padj) & ref_df$padj < FDR_THRESHOLD &
                            abs(ref_df$log2FoldChange) > LFC_THRESHOLD, ]

        if ("baseMean" %in% colnames(ref_sig)) {
          ref_sig <- ref_sig[order(-ref_sig$baseMean), ]
        } else {
          ref_sig <- ref_sig[order(ref_sig$padj), ]
        }

        top_genes <- head(ref_sig$gene_id, N_TOP_GENES)

        if (length(top_genes) > 0) {
          lfc_mat <- build_lfc_matrix(top_genes, all_results)
          cor_mat <- cor(lfc_mat, use = "pairwise.complete.obs")
          cor_mat[is.na(cor_mat)] <- 0

          png(file.path(diag_dir, "LFC_Correlation_Top100.png"),
              width = 800, height = 700, res = 120)
          pheatmap(cor_mat,
                   display_numbers = TRUE, number_format = "%.2f",
                   main = sprintf("LFC Correlation (Top %d %s genes) - TOY",
                                  N_TOP_GENES, ref_name),
                   color = colorRampPalette(brewer.pal(11, "RdBu"))(100))
          dev.off()
          log_info("  Saved LFC_Correlation_Top100.png")
        } else {
          log_warn("  No significant genes in reference for LFC correlation")
        }
      } else {
        log_warn("  Reference dataset '%s' not found", REFERENCE_DATASET)
      }
    }, error = function(e) {
      log_error("  LFC correlation failed: %s", e$message)
      try(dev.off(), silent = TRUE)
    })

    # --- Diagnostic 3: Concordance Analysis ---
    log_info("[3/4] Concordance Analysis vs %s...", REFERENCE_DATASET)
    tryCatch({
      concordance_df <- compute_concordance(sig_list, reference = REFERENCE_DATASET)
      if (nrow(concordance_df) > 0) {
        conc_file <- file.path(diag_dir,
                               paste0("Concordance_Analysis_vs_", REFERENCE_DATASET, ".csv"))
        write.csv(concordance_df, conc_file, row.names = FALSE, quote = TRUE)
        log_info("  Saved Concordance_Analysis_vs_%s.csv", REFERENCE_DATASET)
      }
    }, error = function(e) {
      log_error("  Concordance analysis failed: %s", e$message)
    })

    # --- Diagnostic 4: XIST Detection ---
    log_info("[4/4] XIST Detection Summary...")
    tryCatch({
      xist_rows <- data.frame()
      for (ds in names(all_results)) {
        ds_df <- all_results[[ds]]
        idx   <- which(ds_df$gene_id == TARGET_GENE_ID)

        if (length(idx) > 0) {
          row <- ds_df[idx[1], ]
          xist_rows <- rbind(xist_rows, data.frame(
            Dataset        = ds,
            Found          = TRUE,
            log2FoldChange = row$log2FoldChange,
            padj           = row$padj,
            baseMean       = if ("baseMean" %in% colnames(row)) row$baseMean else NA_real_,
            stringsAsFactors = FALSE
          ))
        } else {
          xist_rows <- rbind(xist_rows, data.frame(
            Dataset        = ds,
            Found          = FALSE,
            log2FoldChange = NA_real_,
            padj           = NA_real_,
            baseMean       = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
      }
      xist_file <- file.path(diag_dir, "XIST_Detection_Summary.csv")
      write.csv(xist_rows, xist_file, row.names = FALSE, quote = TRUE)
      log_info("  Saved XIST_Detection_Summary.csv")

      for (i in seq_len(nrow(xist_rows))) {
        r <- xist_rows[i, ]
        if (r$Found) {
          log_info("    %s: %s found (LFC=%.3f, padj=%.2e)",
                   r$Dataset, TARGET_GENE_NAME, r$log2FoldChange, r$padj)
        } else {
          log_warn("    %s: %s NOT found", r$Dataset, TARGET_GENE_NAME)
        }
      }
    }, error = function(e) {
      log_error("  XIST detection failed: %s", e$message)
    })
  }

}, error = function(e) {
  log_error("Step 2 (Diagnostics) failed: %s", e$message)
})

# ==============================================================================
# STEP 3: SUMMARY (INLINE)
# ==============================================================================

log_section("STEP 3: Summary (Toy)")

tryCatch({

  summary_rows <- list()
  all_res <- load_all_de_results(output_dir)

  for (ds in names(all_res)) {
    df <- all_res[[ds]]
    sig <- df[!is.na(df$padj) & df$padj < FDR_THRESHOLD &
                abs(df$log2FoldChange) > LFC_THRESHOLD, ]
    n_sig    <- nrow(sig)
    n_tested <- nrow(df)

    df$gene_clean <- clean_gene_ids(df$gene)
    xist_row <- df[df$gene_clean == TARGET_GENE_ID, ]

    if (nrow(xist_row) > 0) {
      xist_row <- xist_row[1, ]
      xist_detected <- "yes"
      xist_padj <- xist_row$padj
      xist_lfc  <- xist_row$log2FoldChange
    } else {
      xist_detected <- "no"
      xist_padj <- NA_real_
      xist_lfc  <- NA_real_
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      CellType      = CELL_TYPES[1],
      Dataset       = ds,
      N_Significant = n_sig,
      N_Tested      = n_tested,
      XIST_Detected = xist_detected,
      XIST_padj     = xist_padj,
      XIST_LFC      = xist_lfc,
      stringsAsFactors = FALSE
    )
  }

  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    summary_csv <- file.path(OUTPUT_BASE, "TOY_Summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    log_info("Saved summary: %s", summary_csv)

    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      log_info("  %s: %d sig genes / %d tested, XIST=%s",
               r$Dataset, r$N_Significant, r$N_Tested, r$XIST_Detected)
    }
  }

}, error = function(e) {
  log_error("Step 3 (Summary) failed: %s", e$message)
})

# ==============================================================================
# VALIDATION CHECKS
# ==============================================================================

log_section("VALIDATION CHECKS")

# Check DE results CSVs exist for both datasets
for (d_name in RUN_DATASETS) {
  csv_path <- file.path(output_dir, paste0(d_name, "_DE_results.csv"))
  validate(
    sprintf("DE results CSV exists: %s", d_name),
    file.exists(csv_path)
  )
}

# Check pseudobulk SCE objects saved
for (d_name in RUN_DATASETS) {
  pb_path <- file.path(output_dir, paste0(d_name, "_pb_checkpoint.qs"))
  validate(
    sprintf("Pseudobulk SCE saved: %s", d_name),
    file.exists(pb_path)
  )
}

# Check DESeq2 dds objects saved
for (d_name in RUN_DATASETS) {
  dds_path <- file.path(output_dir, paste0(d_name, "_dds_checkpoint.qs"))
  validate(
    sprintf("DESeq2 dds object saved: %s", d_name),
    file.exists(dds_path)
  )
}

# Check UpSet plot generated
validate(
  "UpSet plot generated",
  file.exists(file.path(diag_dir, "Diagnostic_UpSet_Plot.png"))
)

# Check LFC correlation plot generated
validate(
  "LFC correlation plot generated",
  file.exists(file.path(diag_dir, "LFC_Correlation_Top100.png"))
)

# Check no NA values in critical columns (at least some non-NA)
for (d_name in RUN_DATASETS) {
  csv_path <- file.path(output_dir, paste0(d_name, "_DE_results.csv"))
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)

    n_non_na_padj <- sum(!is.na(df$padj))
    validate(
      sprintf("Non-NA padj values present: %s (%d non-NA)", d_name, n_non_na_padj),
      n_non_na_padj > 0
    )

    n_non_na_lfc <- sum(!is.na(df$log2FoldChange))
    validate(
      sprintf("Non-NA log2FoldChange values present: %s (%d non-NA)", d_name, n_non_na_lfc),
      n_non_na_lfc > 0
    )
  }
}

# Check XIST shows expected sex bias (if present)
for (d_name in RUN_DATASETS) {
  csv_path <- file.path(output_dir, paste0(d_name, "_DE_results.csv"))
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    df$gene_clean <- clean_gene_ids(df$gene)
    xist_row <- df[df$gene_clean == TARGET_GENE_ID, ]

    if (nrow(xist_row) > 0) {
      # XIST should be female-biased (negative LFC if male is reference)
      # or at least significant — just check it's present and non-NA
      validate(
        sprintf("XIST detected with non-NA LFC: %s (LFC=%.3f)",
                d_name, xist_row$log2FoldChange[1]),
        !is.na(xist_row$log2FoldChange[1])
      )
    } else {
      log_warn("XIST (%s) not found in %s results — cannot validate bias",
               TARGET_GENE_ID, d_name)
    }
  }
}

# Check significant gene counts are reasonable (not 0 and not all genes)
for (d_name in RUN_DATASETS) {
  csv_path <- file.path(output_dir, paste0(d_name, "_DE_results.csv"))
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    n_sig <- sum(!is.na(df$padj) & df$padj < FDR_THRESHOLD &
                   abs(df$log2FoldChange) > LFC_THRESHOLD, na.rm = TRUE)
    n_total <- nrow(df)

    validate(
      sprintf("Sig gene count reasonable: %s (%d sig / %d total)",
              d_name, n_sig, n_total),
      n_sig > 0 && n_sig < n_total
    )
  }
}

# ==============================================================================
# CLEANUP
# ==============================================================================

# Remove temp shard symlinks
if (dir.exists(toy_temp_base)) {
  unlink(toy_temp_base, recursive = TRUE)
  log_info("Cleaned up temp shard directory: %s", toy_temp_base)
}

# ==============================================================================
# FINAL REPORT
# ==============================================================================

toy_end_time <- Sys.time()
toy_duration <- difftime(toy_end_time, toy_start_time, units = "mins")

log_section("TOY TEST FINAL REPORT")

n_pass <- sum(sapply(validation_results, function(v) v$status == "PASS"))
n_fail <- sum(sapply(validation_results, function(v) v$status == "FAIL"))
n_total <- length(validation_results)

log_info("Duration: %.1f minutes", as.numeric(toy_duration))
log_info("Validations: %d PASS / %d FAIL / %d TOTAL", n_pass, n_fail, n_total)
log_info("")
log_info("Output paths:")
log_info("  DE results:  %s/", output_dir)
log_info("  Diagnostics: %s/", diag_dir)
log_info("  Summary CSV: %s", file.path(OUTPUT_BASE, "TOY_Summary.csv"))

if (n_fail > 0) {
  log_info("")
  log_warn("FAILED VALIDATIONS:")
  for (v in validation_results) {
    if (v$status == "FAIL") {
      log_warn("  [FAIL] %s", v$description)
    }
  }
}

log_info("")
if (n_fail == 0) {
  log_info("ALL %d VALIDATIONS PASSED", n_total)
  log_info("Toy test suggests full pipeline will work correctly.")
} else {
  log_warn("%d VALIDATION(S) FAILED — investigate before full run.", n_fail)
}

log_info("")
log_info("========================================")
log_info("AWAITING USER APPROVAL")
log_info("========================================")
log_info("Please review:")
log_info("  1. Validation results above")
log_info("  2. Diagnostic plots in: %s", diag_dir)
log_info("  3. DE result CSVs in: %s", output_dir)
log_info("If all looks correct, proceed to Stage 3 (full production run).")
