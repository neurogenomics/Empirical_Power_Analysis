# ==============================================================================
# SANITY CHECKS: Fundamental Data Validation
# PURPOSE: Catch data issues that indicate deeper problems (not just test failures)
# ==============================================================================

#' Sanity check helper function
#' 
#' @param data Data to check
#' @param description Description of what's being checked
#' @param checks Named list of check functions
#' @return TRUE if all checks pass, FALSE otherwise
sanity_check <- function(data, description, checks) {
  cat(sprintf("\n[SANITY CHECK] %s\n", description))
  issues <- c()
  
  for (check_name in names(checks)) {
    check_fn <- checks[[check_name]]
    result <- tryCatch({
      check_fn(data)
    }, error = function(e) {
      sprintf("ERROR: %s", e$message)
    })
    
    if (!isTRUE(result)) {
      issues <- c(issues, sprintf("⚠️  %s: %s", check_name, result))
    }
  }
  
  if (length(issues) == 0) {
    cat("  ✅ PASS - All sanity checks passed\n")
    return(TRUE)
  } else {
    cat("  🔴 SANITY CHECK FAILURES:\n")
    for (issue in issues) {
      cat(sprintf("     %s\n", issue))
    }
    cat("\n  ⚠️  These are NOT test failures - they indicate fundamental data issues.\n")
    cat("     Investigate before proceeding. Data may be normalized, corrupted, or wrong assay.\n")
    return(FALSE)
  }
}

#' Count matrix sanity checks (scRNA-seq)
#' 
#' @param counts Count matrix (genes x samples)
#' @return Named list of check functions
count_matrix_checks <- function() {
  list(
    "All integers" = function(counts) {
      if (all(counts == floor(counts), na.rm = TRUE)) TRUE
      else sprintf("%.2f%% non-integer values", mean(counts != floor(counts), na.rm = TRUE) * 100)
    },
    "Max value > 1000" = function(counts) {
      max_val <- max(counts, na.rm = TRUE)
      if (max_val > 1000) TRUE
      else sprintf("Max = %.0f (expect >1000 for real data)", max_val)
    },
    "Sparsity 70-95%" = function(counts) {
      zero_pct <- mean(counts == 0, na.rm = TRUE) * 100
      if (zero_pct >= 70 && zero_pct <= 95) TRUE
      else sprintf("%.1f%% zeros (expect 70-95%%)", zero_pct)
    },
    "No negatives" = function(counts) {
      if (!any(counts < 0, na.rm = TRUE)) TRUE
      else "CRITICAL: Negative counts detected"
    },
    "Non-empty" = function(counts) {
      if (nrow(counts) > 0 && ncol(counts) > 0) TRUE
      else sprintf("Empty matrix: %d x %d", nrow(counts), ncol(counts))
    }
  )
}

#' Pseudobulk sanity checks
#' 
#' @param pb_sce SingleCellExperiment pseudobulk object
#' @return Named list of check functions
pseudobulk_checks <- function() {
  list(
    "Counts are integers" = function(pb_sce) {
      counts <- assay(pb_sce, "counts")
      if (all(counts == floor(counts), na.rm = TRUE)) TRUE
      else sprintf("%.2f%% non-integer", mean(counts != floor(counts), na.rm = TRUE) * 100)
    },
    "Donor count reasonable" = function(pb_sce) {
      n_donors <- ncol(pb_sce)
      if (n_donors >= 5 && n_donors <= 10000) TRUE
      else sprintf("%d donors (expect 5-1000)", n_donors)
    },
    "Gene count reasonable" = function(pb_sce) {
      n_genes <- nrow(pb_sce)
      if (n_genes >= 10000 && n_genes <= 40000) TRUE
      else sprintf("%d genes (expect 10k-40k for human)", n_genes)
    },
    "No all-zero samples" = function(pb_sce) {
      counts <- assay(pb_sce, "counts")
      col_sums <- colSums(counts)
      n_zero <- sum(col_sums == 0)
      if (n_zero == 0) TRUE
      else sprintf("%d all-zero samples (aggregation failure)", n_zero)
    }
  )
}

#' Metadata sanity checks
#' 
#' @param metadata Data frame of metadata
#' @param required_cols Character vector of required column names
#' @return Named list of check functions
metadata_checks <- function(required_cols = c("donor_id", "sex", "age")) {
  list(
    "Required columns present" = function(metadata) {
      missing <- setdiff(required_cols, colnames(metadata))
      if (length(missing) == 0) TRUE
      else sprintf("Missing: %s", paste(missing, collapse = ", "))
    },
    "No completely empty columns" = function(metadata) {
      empty_cols <- colnames(metadata)[colSums(is.na(metadata)) == nrow(metadata)]
      if (length(empty_cols) == 0) TRUE
      else sprintf("Empty columns: %s", paste(empty_cols, collapse = ", "))
    },
    "Donor IDs unique" = function(metadata) {
      if ("donor_id" %in% colnames(metadata)) {
        n_dup <- sum(duplicated(metadata$donor_id))
        if (n_dup == 0) TRUE
        else sprintf("%d duplicate donor IDs", n_dup)
      } else TRUE
    },
    "Age values reasonable" = function(metadata) {
      if ("age" %in% colnames(metadata)) {
        ages <- as.numeric(metadata$age)
        ages <- ages[!is.na(ages)]
        if (length(ages) == 0) return("All ages NA")
        if (all(ages >= 0 & ages <= 120)) TRUE
        else sprintf("Ages out of range: [%.1f, %.1f]", min(ages), max(ages))
      } else TRUE
    }
  )
}

#' DESeq2 results sanity checks
#' 
#' @param res DESeqResults object
#' @return Named list of check functions
deseq2_results_checks <- function() {
  list(
    "padj in [0,1] range" = function(res) {
      padj <- res$padj[!is.na(res$padj)]
      if (length(padj) == 0) return("All padj are NA")
      if (all(padj >= 0 & padj <= 1)) TRUE
      else sprintf("padj out of range: [%.3f, %.3f]", min(padj), max(padj))
    },
    "log2FoldChange not all NA" = function(res) {
      lfc <- res$log2FoldChange
      n_na <- sum(is.na(lfc))
      pct_na <- (n_na / length(lfc)) * 100
      if (pct_na < 100) TRUE
      else "All log2FoldChange are NA"
    },
    "baseMean > 0" = function(res) {
      bm <- res$baseMean[!is.na(res$baseMean)]
      if (length(bm) == 0) return("All baseMean are NA")
      n_zero <- sum(bm == 0)
      if (n_zero == 0) TRUE
      else sprintf("%d genes with baseMean = 0", n_zero)
    },
    "Expected columns present" = function(res) {
      required <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
      missing <- setdiff(required, colnames(as.data.frame(res)))
      if (length(missing) == 0) TRUE
      else sprintf("Missing columns: %s", paste(missing, collapse = ", "))
    }
  )
}
