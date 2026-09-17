#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Sex Label Validation (v2 — pseudobulk CPM + logistic regression)
# PURPOSE:
#   Validates donor sex metadata labels by comparing pseudobulk CPM expression of:
#     - XIST       : X-inactivation gene, expected HIGH in females, ZERO in males
#     - DDX3Y      : Y-chromosome gene, expected HIGH in males, ZERO in females
#     - UTY        : Y-chromosome gene (note: has X paralog UTX/KDM6A — lower
#                    weight in regression; reported but not used as primary predictor)
#     - KDM5D      : Y-chromosome gene, definitively Y-only
#
#   APPROACH: Pseudobulk CPM per donor (sum gene counts / sum all UMIs * 1e6).
#   This avoids the noise inherent in per-cell percentages:
#     - A single stray XIST read across 500 cells = ~0.07 CPM (negligible)
#     - A genuine female across 500 cells = ~10,000 CPM (unambiguous)
#   The separation is ~4 orders of magnitude, making thresholds unnecessary.
#
#   A logistic regression is fitted per cell type on the log2-CPM values.
#   Donors whose predicted probability strongly disagrees with their label
#   are flagged as potential mislabels.
#
# USAGE:
#   Rscript R/eda_roussos/Validate_Sex_Labels.R [CellType]
#   CellType: one of Astro, Micro, Oligo, or "All" (default: All)
#
# OUTPUT:
#   Datalink/Sex_Validation/sex_validation_{CellType}.csv   — per-donor results
#   Datalink/Sex_Validation/sex_exclusion_{CellType}.csv    — flagged donors only
#   Datalink/Sex_Validation/sex_exclusion_combined.csv      — all cell types merged
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(dplyr)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Markers — Ensembl IDs (rownames in shards are Ensembl format)
X_GENE  <- "ENSG00000229807"   # XIST
Y_GENES <- c(
  DDX3Y = "ENSG00000067048",   # Y-linked RNA helicase (most robust, no X paralog)
  UTY   = "ENSG00000183878",   # Y-linked histone demethylase (has X paralog UTX — note)
  KDM5D = "ENSG00000012817"    # Y-linked chromatin regulator (no X paralog)
)
X_GENE_LABEL  <- "XIST"
Y_GENE_LABELS <- names(Y_GENES)

# Y genes used as the primary predictor in the logistic regression.
# UTY excluded because it shares sequence similarity with KDM6A/UTX on the X
# chromosome (the cross-mapping we observed in the prior analysis). DDX3Y and
# KDM5D have no X paralog and are clean predictors.
Y_PRIMARY_LABELS <- c("DDX3Y", "KDM5D")

# Logistic regression probability thresholds for flagging.
# A donor is flagged if:
#   - labelled male  but P(female | XIST_CPM) > XIST_PROB_THRESHOLD
#   - labelled female but P(female | XIST_CPM) < 1 - XIST_PROB_THRESHOLD
#   - labelled female but P(male   | Y_CPM)   > Y_PROB_THRESHOLD
#   - labelled male  but P(male   | Y_CPM)    < 1 - Y_PROB_THRESHOLD
XIST_PROB_THRESHOLD <- 0.80   # >80% probability of being female despite male label
Y_PROB_THRESHOLD    <- 0.80   # >80% probability of being male despite female label

# ==============================================================================
# RESOLVE PROJECT ROOT
# ==============================================================================

if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  PROJECT_ROOT <- Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` not installed. Run from repo root or set POWER_ANALYSIS_ROOT.")
  }
  PROJECT_ROOT <- here::here()
}

SHARDS_BASE <- file.path(PROJECT_ROOT, "Datalink", "Shards")
OUT_DIR     <- file.path(PROJECT_ROOT, "Datalink", "Sex_Validation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================

VALID_CELL_TYPES <- c("Astro", "Micro", "Oligo")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1]) && args[1] != "All") {
  if (!args[1] %in% VALID_CELL_TYPES) {
    stop(sprintf("Unknown cell type '%s'. Use one of: %s, All",
                 args[1], paste(VALID_CELL_TYPES, collapse = ", ")))
  }
  CELL_TYPES <- args[1]
} else {
  CELL_TYPES <- VALID_CELL_TYPES
}

# ==============================================================================
# HELPER: extract pseudobulk CPM per donor from one shard
# ==============================================================================
#
# Pseudobulk CPM for gene g in donor d:
#   CPM_g = sum(counts[g, ]) / sum(counts) * 1e6
#
# where sum(counts) is the total UMI count across ALL genes and ALL cells
# for that donor. This is the correct denominator: it is proportional to
# sequencing depth and makes CPM comparable across donors of different sizes.

extract_donor_pseudobulk <- function(shard_path, target_genes) {
  sce <- tryCatch(
    qs::qread(shard_path),
    error = function(e) {
      warning(sprintf("Cannot read shard %s: %s", basename(shard_path), e$message))
      return(NULL)
    }
  )
  if (is.null(sce)) return(NULL)

  cd         <- as.data.frame(colData(sce))
  donor_id   <- as.character(unique(cd[["donor_id"]]))[1]
  sex_label  <- as.character(unique(cd[["sex"]]))[1]
  n_cells    <- ncol(sce)

  # Total UMI across ALL genes and ALL cells for this donor.
  # Using colSums first (library size per cell) then summing is equivalent to
  # summing the whole count matrix, but slightly more memory-efficient.
  total_umi  <- sum(Matrix::colSums(counts(sce)))

  gene_names <- rownames(sce)
  # Map Ensembl ID -> display symbol (Y_GENES is defined symbol -> EnsemblID, so invert)
  all_labels <- c(setNames(X_GENE_LABEL, X_GENE),
                  setNames(names(Y_GENES), unname(Y_GENES)))

  out <- data.frame(
    donor_id  = donor_id,
    sex_label = sex_label,
    n_cells   = n_cells,
    total_umi = total_umi,
    stringsAsFactors = FALSE
  )

  for (g in target_genes) {
    lbl <- all_labels[[g]]
    col <- paste0(lbl, "_cpm")
    if (g %in% gene_names) {
      gene_total   <- sum(as.numeric(counts(sce)[g, ]))
      out[[col]]   <- if (total_umi > 0) round(gene_total / total_umi * 1e6, 4) else NA_real_
    } else {
      out[[col]]   <- NA_real_
    }
  }

  out
}

# ==============================================================================
# HELPER: fit logistic regression and return per-donor probabilities
# ==============================================================================
#
# fit_sex_glm fits:
#   P(is_positive ~ log2(cpm + 1), family = binomial)
#
# where is_positive = 1 if the label matches the "positive" class.
# Returns a vector of fitted probabilities (same length as nrow(df)).
# If the fit fails (e.g. singular — all CPM zero or all one class) it
# returns NA for all donors.

fit_sex_glm <- function(df, cpm_col, positive_label) {
  y   <- as.integer(df$sex_std == positive_label)
  x   <- log2(df[[cpm_col]] + 1)

  if (length(unique(y)) < 2 || all(is.na(x))) {
    warning(sprintf("Cannot fit GLM for %s ~ %s: not enough variation.",
                    positive_label, cpm_col))
    return(rep(NA_real_, nrow(df)))
  }

  fit <- tryCatch(
    suppressWarnings(
      glm(y ~ x, family = binomial(link = "logit"))
    ),
    error = function(e) {
      warning(sprintf("GLM fit failed for %s: %s", cpm_col, e$message))
      return(NULL)
    }
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(df)))

  predict(fit, type = "response")
}

# ==============================================================================
# HELPER: flag discordant donors based on logistic regression probabilities
# ==============================================================================

flag_donor_glm <- function(sex_std, p_female_xist, p_male_y,
                            xist_thresh = XIST_PROB_THRESHOLD,
                            y_thresh    = Y_PROB_THRESHOLD) {
  flags <- character(0)

  if (!is.na(p_female_xist)) {
    if (sex_std == "Female" && p_female_xist < (1 - xist_thresh)) {
      flags <- c(flags, "XIST_absent_in_female")
    }
    if (sex_std == "Male" && p_female_xist > xist_thresh) {
      flags <- c(flags, "XIST_present_in_male")
    }
  }

  if (!is.na(p_male_y)) {
    if (sex_std == "Male" && p_male_y < (1 - y_thresh)) {
      flags <- c(flags, "Y_absent_in_male")
    }
    if (sex_std == "Female" && p_male_y > y_thresh) {
      flags <- c(flags, "Y_present_in_female")
    }
  }

  if (length(flags) == 0) "OK" else paste(flags, collapse = " | ")
}

# ==============================================================================
# MAIN LOOP
# ==============================================================================

all_genes <- c(X_GENE, unname(Y_GENES))

for (ct in CELL_TYPES) {
  shard_dir <- file.path(SHARDS_BASE, ct, "Roussos")

  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  Cell type : %s\n", ct))
  cat(sprintf("  Shard dir : %s\n", shard_dir))
  cat(sprintf("  Approach  : pseudobulk CPM + logistic regression\n"))
  cat(sprintf("%s\n\n", strrep("=", 70)))

  if (!dir.exists(shard_dir)) {
    cat(sprintf("  [SKIP] Directory not found: %s\n", shard_dir))
    next
  }

  shard_files <- list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
  shard_files <- shard_files[!grepl("MANIFEST", shard_files, ignore.case = TRUE)]

  if (length(shard_files) == 0) {
    cat(sprintf("  [SKIP] No donor .qs shards found in %s\n", shard_dir))
    next
  }

  cat(sprintf("  Found %d donor shards. Computing pseudobulk CPM...\n", length(shard_files)))

  results_list <- lapply(seq_along(shard_files), function(i) {
    if (i %% 100 == 0) {
      cat(sprintf("    [%d / %d] processed...\n", i, length(shard_files)))
    }
    extract_donor_pseudobulk(shard_files[[i]], all_genes)
  })

  donor_df <- bind_rows(Filter(Negate(is.null), results_list))

  if (nrow(donor_df) == 0) {
    cat("  [ERROR] No donor data extracted.\n")
    next
  }

  cat(sprintf("  Successfully processed %d donors.\n\n", nrow(donor_df)))

  # --- Warn about completely missing genes ---
  for (lbl in c(X_GENE_LABEL, Y_GENE_LABELS)) {
    col <- paste0(lbl, "_cpm")
    if (col %in% names(donor_df) && all(is.na(donor_df[[col]]))) {
      cat(sprintf("  [WARNING] Gene '%s' not found in the %s dataset.\n", lbl, ct))
    }
  }

  # --- Standardise sex labels ---
  donor_df <- donor_df %>%
    mutate(sex_std = case_when(
      tolower(trimws(sex_label)) %in% c("male",   "m", "1") ~ "Male",
      tolower(trimws(sex_label)) %in% c("female", "f", "0") ~ "Female",
      TRUE ~ NA_character_
    ))

  n_unknown_sex <- sum(is.na(donor_df$sex_std))
  if (n_unknown_sex > 0) {
    cat(sprintf("  [WARNING] %d donors have unrecognised sex labels:\n", n_unknown_sex))
    print(table(donor_df$sex_label[is.na(donor_df$sex_std)]))
    cat("\n")
  }

  # --- Summary of pseudobulk CPM by sex ---
  cat("  --- Pseudobulk CPM summary by sex label ---\n")
  cpm_cols <- c(paste0(X_GENE_LABEL, "_cpm"), paste0(Y_GENE_LABELS, "_cpm"))
  cpm_cols_present <- intersect(cpm_cols, names(donor_df))

  for (col in cpm_cols_present) {
    gene_name <- sub("_cpm$", "", col)
    s <- donor_df %>%
      group_by(sex_std) %>%
      summarise(
        n      = n(),
        median_cpm = round(median(.data[[col]], na.rm = TRUE), 2),
        p95_cpm    = round(quantile(.data[[col]], 0.95, na.rm = TRUE), 2),
        pct_zero   = round(mean(.data[[col]] == 0, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      )
    cat(sprintf("    %s:\n", gene_name))
    print(as.data.frame(s))
    cat("\n")
  }

  # --- Fit logistic regressions ---
  #
  # XIST: P(Female) ~ log2(XIST_cpm + 1)
  # Y:    P(Male)   ~ log2(Y_combined_cpm + 1)
  #   where Y_combined = DDX3Y_cpm + KDM5D_cpm
  #   (UTY excluded from predictor due to X paralog cross-mapping noise)

  # XIST regression
  cat("  --- Fitting logistic regression (XIST) ---\n")
  donor_df$P_female_XIST <- fit_sex_glm(donor_df, "XIST_cpm", "Female")

  # Y combined regression
  y_primary_cpm_cols <- paste0(Y_PRIMARY_LABELS, "_cpm")
  y_present_cols     <- intersect(y_primary_cpm_cols, names(donor_df))

  if (length(y_present_cols) == 0) {
    cat("  [WARNING] No primary Y genes (DDX3Y, KDM5D) found — skipping Y regression.\n")
    donor_df$Y_combined_cpm <- NA_real_
    donor_df$P_male_Y       <- NA_real_
  } else {
    donor_df$Y_combined_cpm <- rowSums(
      donor_df[, y_present_cols, drop = FALSE],
      na.rm = TRUE
    )
    cat("  --- Fitting logistic regression (Y genes combined) ---\n")
    donor_df$P_male_Y <- fit_sex_glm(donor_df, "Y_combined_cpm", "Male")
  }

  # --- Flag donors ---
  donor_df$flag <- mapply(
    flag_donor_glm,
    sex_std      = donor_df$sex_std,
    p_female_xist = donor_df$P_female_XIST,
    p_male_y     = donor_df$P_male_Y,
    MoreArgs = list(
      xist_thresh = XIST_PROB_THRESHOLD,
      y_thresh    = Y_PROB_THRESHOLD
    )
  )

  # --- Concordance summary ---
  n_ok         <- sum(donor_df$flag == "OK",  na.rm = TRUE)
  n_discordant <- sum(donor_df$flag != "OK",  na.rm = TRUE)
  n_total      <- nrow(donor_df)

  cat(sprintf("\n  --- Concordance: %d / %d donors OK (%.1f%%)\n\n",
              n_ok, n_total, n_ok / n_total * 100))

  if (n_discordant > 0) {
    cat(sprintf("  --- DISCORDANT donors (%d) ---\n", n_discordant))

    display_cols <- c("donor_id", "sex_std", "n_cells",
                      "XIST_cpm", "P_female_XIST",
                      "Y_combined_cpm", "P_male_Y",
                      "flag")
    display_cols <- intersect(display_cols, names(donor_df))

    discordant_df <- donor_df %>%
      filter(flag != "OK") %>%
      select(all_of(display_cols)) %>%
      arrange(flag, sex_std, donor_id) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))

    print(as.data.frame(discordant_df))
    cat("\n")
  } else {
    cat("  All sex labels are concordant with XIST and Y-gene expression.\n\n")
  }

  # --- Flag breakdown ---
  cat("  --- Flag breakdown ---\n")
  print(table(donor_df$flag))
  cat("\n")

  # --- Probability distribution summary (diagnostic) ---
  cat("  --- XIST P(Female) distribution by sex label ---\n")
  prob_summary <- donor_df %>%
    group_by(sex_std) %>%
    summarise(
      n              = n(),
      p_median       = round(median(P_female_XIST, na.rm = TRUE), 3),
      p_min          = round(min(P_female_XIST,    na.rm = TRUE), 3),
      p_max          = round(max(P_female_XIST,    na.rm = TRUE), 3),
      n_prob_gt_0.8  = sum(P_female_XIST > 0.8,   na.rm = TRUE),
      n_prob_lt_0.2  = sum(P_female_XIST < 0.2,   na.rm = TRUE),
      .groups = "drop"
    )
  print(as.data.frame(prob_summary))
  cat("\n")

  # --- Save per-donor results ---
  out_file <- file.path(OUT_DIR, sprintf("sex_validation_%s.csv", ct))
  write.csv(donor_df, out_file, row.names = FALSE)
  cat(sprintf("  Per-donor results saved to:\n  %s\n\n", out_file))

  # --- Save exclusion list ---
  excl_df <- donor_df %>%
    filter(flag != "OK") %>%
    mutate(cell_type = ct) %>%
    select(donor_id, cell_type, sex_label = sex_std, flag,
           XIST_cpm, P_female_XIST,
           Y_combined_cpm, P_male_Y,
           n_cells, total_umi)

  excl_file <- file.path(OUT_DIR, sprintf("sex_exclusion_%s.csv", ct))
  write.csv(excl_df, excl_file, row.names = FALSE)

  if (nrow(excl_df) > 0) {
    cat(sprintf("  Exclusion list (%d donors) saved to:\n  %s\n\n",
                nrow(excl_df), excl_file))
  } else {
    cat(sprintf("  No exclusions for %s — empty exclusion list written.\n\n", ct))
  }

  if (!exists("all_exclusions")) {
    all_exclusions <- excl_df
  } else {
    all_exclusions <- bind_rows(all_exclusions, excl_df)
  }
}

# ==============================================================================
# COMBINED EXCLUSION LIST
# ==============================================================================

if (exists("all_exclusions")) {
  combined_file <- file.path(OUT_DIR, "sex_exclusion_combined.csv")
  all_exclusions <- all_exclusions %>%
    distinct(donor_id, cell_type, .keep_all = TRUE) %>%
    arrange(cell_type, donor_id)
  write.csv(all_exclusions, combined_file, row.names = FALSE)
  n_unique_donors <- n_distinct(all_exclusions$donor_id)
  cat(sprintf("  Combined exclusion list: %d unique donors (%d rows across cell types)\n",
              n_unique_donors, nrow(all_exclusions)))
  cat(sprintf("  Saved to: %s\n\n", combined_file))
} else {
  cat("  No cell types processed — no combined exclusion list written.\n\n")
}

cat(sprintf("\n%s\n", strrep("=", 70)))
cat("  Sex label validation complete.\n")
cat(sprintf("  Outputs saved in: %s\n", OUT_DIR))
cat(sprintf("%s\n\n", strrep("=", 70)))
