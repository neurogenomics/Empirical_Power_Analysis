# ==============================================================================
# TEST: DE_FIT_TYPE equivalence + speedup benchmark
#
# PURPOSE
#   Demonstrate that DESeq2::DESeq(fitType = "glmGamPoi") produces results
#   equivalent to the default "parametric" fit on sex-biased pseudobulk DE,
#   and quantify the speedup. This lets us set DE_FIT_TYPE <- "glmGamPoi" in
#   config.R with confidence that published findings are preserved.
#
#   Equivalence criteria (must ALL pass):
#     1. Pearson correlation of LFCs on shared tested genes      >= 0.99
#     2. Spearman correlation of -log10(padj) on shared genes    >= 0.99
#     3. Jaccard overlap of significant gene sets (padj < 0.05)  >= 0.95
#     4. Sign agreement on significant genes                     >= 0.99
#     5. Kendall's tau on top-100 gene ranking                   >= 0.95
#     6. No gene flips direction in the intersection of          == TRUE
#        significant-in-both sets with |LFC| > 0.5
#
#   Speedup criterion (informational, not gating):
#     wall_time_parametric / wall_time_glmgampoi > 3
#
# RUN ON:
#   A real aggregated pseudobulk DataSet (or a large enough synthetic one) so
#   that both fitTypes have meaningful dispersion estimation work to do. Small
#   toy data (6 donors x 100 cells) will NOT exercise glmGamPoi's strengths.
#
# USAGE:
#   module load R/4.5.1-gfbf-2025a
#   Rscript R/cartesian_downsampling/tests/test_fittype_equivalence.R /mnt/data/chrisc/Power_Analysis_Paper/Datalink/Cartesian_Downsampling/Oligo/Roussos/cells/Run_1/Ground_Truth.qs
#     /path/to/agg_result.qs       # optional; defaults to synthetic
#
#Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/cartesian_downsampling/tests/test_fittype_equivalence.R /mnt/data/chrisc/Power_Analysis_Paper/Datalink/Full_Ground_Truth/Oligo/Roussos_pb_checkpoint.qs /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/fittype_benchmark_report.md
#   Or from within testthat:
#     testthat::test_file("R/cartesian_downsampling/tests/test_fittype_equivalence.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(DESeq2)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(qs)
})

# Null-coalescing helper (avoids rlang dep)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# 0. Helper: build a realistic pseudobulk test input
# ------------------------------------------------------------------------------
#' @description Create a synthetic pseudobulk dataset large enough to exercise
#'   glmGamPoi's dispersion estimator (not a 6-donor toy). Injects 200
#'   sex-biased genes so both methods have real signal to recover.
#' @param n_donors  integer; 50-200 recommended for meaningful benchmarks
#' @param n_genes   integer; 5000-30000 recommended
#' @param n_biased  integer; number of "true" sex-biased genes to inject
#' @param seed      integer; deterministic
#' @return list(sce, truth_genes)
build_realistic_pseudobulk <- function(n_donors = 100, n_genes = 10000,
                                        n_biased = 200, seed = 42) {
  set.seed(seed)
  sex <- factor(rep(c("male", "female"), length.out = n_donors))
  Source <- factor(sample(c("MSSM", "RUSH", "HBCC"), n_donors, replace = TRUE))
  self_reported_ethnicity <- factor(sample(c("EUR", "AFR", "AMR"),
                                            n_donors, replace = TRUE,
                                            prob = c(0.7, 0.2, 0.1)))
  disease <- factor(sample(c("control", "AD"), n_donors, replace = TRUE))

  # Baseline negative-binomial counts — wide dispersion range so both fitTypes
  # have non-trivial work to do
  mu <- rgamma(n_genes, shape = 0.5, rate = 0.01)  # wide mean-expression range
  phi <- 0.1 + rgamma(n_genes, shape = 0.3, rate = 3)  # dispersion 0.1-3ish

  counts <- matrix(0L, nrow = n_genes, ncol = n_donors)
  for (g in seq_len(n_genes)) {
    counts[g, ] <- rnbinom(n_donors, mu = mu[g], size = 1 / phi[g])
  }

  # Inject sex-biased signal on n_biased genes — realistic LFC range for
  # sex-biased genes (mixture: some large like sex-chr, many small autosomal)
  biased_idx <- sample.int(n_genes, n_biased)
  lfc_biased <- c(
    rnorm(n_biased %/% 4, mean = 3, sd = 0.5),    # large LFC (sex-chr-like)
    rnorm(n_biased - n_biased %/% 4, mean = 1, sd = 0.5)  # small-moderate
  )
  male_mask <- sex == "male"
  for (k in seq_along(biased_idx)) {
    g <- biased_idx[k]
    fc_male <- 2^lfc_biased[k]
    counts[g, male_mask] <- rnbinom(sum(male_mask),
                                     mu = mu[g] * fc_male,
                                     size = 1 / phi[g])
  }

  rownames(counts) <- sprintf("gene%05d", seq_len(n_genes))
  colnames(counts) <- sprintf("donor%03d", seq_len(n_donors))

  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(sex = sex,
                        Source = Source,
                        self_reported_ethnicity = self_reported_ethnicity,
                        disease = disease,
                        row.names = colnames(counts))
  )
  list(sce = sce, truth_genes = rownames(counts)[biased_idx])
}

# ------------------------------------------------------------------------------
# 1. Run DESeq2 with one fitType, time it
# ------------------------------------------------------------------------------
run_de_timed <- function(sce, fit_type) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = assay(sce, "counts"),
    colData   = colData(sce),
    design = ~ Source + self_reported_ethnicity + disease + sex
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]

  t0 <- Sys.time()
  dds <- DESeq2::DESeq(dds, fitType = fit_type, parallel = FALSE, quiet = TRUE)
  res <- DESeq2::lfcShrink(dds, coef = "sex_male_vs_female", type = "apeglm",
                           quiet = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(
    res     = as.data.frame(res),
    gene_id = rownames(res),
    elapsed = elapsed
  )
}

# ------------------------------------------------------------------------------
# 2. Compute equivalence metrics on two result sets
# ------------------------------------------------------------------------------
compute_equivalence_metrics <- function(ref, alt) {
  shared <- intersect(ref$gene_id, alt$gene_id)
  stopifnot(length(shared) > 100)

  ref_df <- ref$res[shared, ]
  alt_df <- alt$res[shared, ]

  # Drop genes with NA in either
  keep <- complete.cases(ref_df[, c("log2FoldChange", "padj")]) &
          complete.cases(alt_df[, c("log2FoldChange", "padj")])
  ref_df <- ref_df[keep, ]
  alt_df <- alt_df[keep, ]

  sig_ref <- rownames(ref_df)[ref_df$padj < 0.05]
  sig_alt <- rownames(alt_df)[alt_df$padj < 0.05]
  sig_both <- intersect(sig_ref, sig_alt)
  sig_union <- union(sig_ref, sig_alt)

  top100_ref <- rownames(ref_df)[order(ref_df$padj)][1:min(100, nrow(ref_df))]
  top100_alt <- rownames(alt_df)[order(alt_df$padj)][1:min(100, nrow(alt_df))]
  top100_shared <- intersect(top100_ref, top100_alt)

  # Sign flips in the non-trivial range
  lfc_big <- abs(ref_df$log2FoldChange) > 0.5 & abs(alt_df$log2FoldChange) > 0.5
  sig_in_both <- rownames(ref_df) %in% sig_both
  check_flip <- lfc_big & sig_in_both
  any_flip <- any(
    sign(ref_df$log2FoldChange[check_flip]) !=
      sign(alt_df$log2FoldChange[check_flip])
  )

  kendall_tau <- if (length(top100_shared) > 5) {
    r1 <- match(top100_shared, top100_ref)
    r2 <- match(top100_shared, top100_alt)
    cor(r1, r2, method = "kendall")
  } else NA_real_

  list(
    n_genes_tested     = nrow(ref_df),
    n_sig_ref          = length(sig_ref),
    n_sig_alt          = length(sig_alt),
    n_sig_shared       = length(sig_both),
    jaccard_sig        = length(sig_both) / max(length(sig_union), 1),
    pearson_lfc        = cor(ref_df$log2FoldChange, alt_df$log2FoldChange,
                              method = "pearson", use = "complete"),
    spearman_log_padj  = cor(-log10(ref_df$padj), -log10(alt_df$padj),
                              method = "spearman", use = "complete"),
    sign_agreement_sig = if (length(sig_both) > 0) {
      mean(sign(ref_df[sig_both, "log2FoldChange"]) ==
           sign(alt_df[sig_both, "log2FoldChange"]))
    } else NA_real_,
    kendall_top100     = kendall_tau,
    any_sig_flip       = any_flip
  )
}

# ------------------------------------------------------------------------------
# 3. Main benchmark driver
# ------------------------------------------------------------------------------
#' @param sce optional SingleCellExperiment pseudobulk; synthetic if NULL
#' @param report_file optional path to write a markdown report to
  benchmark_fittype_equivalence <- function(sce = NULL, report_file = NULL,
                                            input_source = NULL) {
    if (is.null(sce)) {
      cat("[info] INPUT SOURCE: synthetic pseudobulk (100 donors, 10k genes)\n")
      if (!is.null(input_source)) cat(sprintf("[info] REASON: %s\n", input_source))
    } else {
      cat(sprintf("[info] INPUT SOURCE: %s\n",
                  input_source %||% "caller-provided sce (input_source not specified)"))
    }
    

  cat(sprintf("[info] pseudobulk shape: %d genes x %d donors\n",
              nrow(sce), ncol(sce)))
  cat(sprintf("[info] colData columns: %s\n",
              paste(colnames(SummarizedExperiment::colData(sce)), collapse = ", ")))

  # Warm up R so first-call overhead doesn't pollute baseline timing
  invisible(run_de_timed(sce[1:100, ], "parametric"))

  cat("[info] running DESeq2 with fitType='parametric' (baseline)...\n")
  ref <- run_de_timed(sce, "parametric")

  cat("[info] running DESeq2 with fitType='glmGamPoi' (candidate)...\n")
  alt <- run_de_timed(sce, "glmGamPoi")

  metrics <- compute_equivalence_metrics(ref, alt)
  speedup <- ref$elapsed / alt$elapsed

  report <- c(
    "# fitType Equivalence + Speedup Benchmark",
    sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("R: %s | DESeq2: %s | glmGamPoi: %s | apeglm: %s",
            R.Version()$version.string,
            packageVersion("DESeq2"),
            packageVersion("glmGamPoi"),
            packageVersion("apeglm")),
    "",
    "## Timing",
    sprintf("| Config | Wall (s) |"),
    sprintf("|---|---|"),
    sprintf("| fitType = parametric  | %.2f |", ref$elapsed),
    sprintf("| fitType = glmGamPoi   | %.2f |", alt$elapsed),
    sprintf("| **Speedup**           | **%.2fx** |", speedup),
    "",
    "## Equivalence",
    sprintf("| Metric | Value | Threshold | Pass |"),
    sprintf("|---|---|---|---|"),
    sprintf("| Pearson LFC        | %.4f | >= 0.99  | %s |",
            metrics$pearson_lfc,   tick(metrics$pearson_lfc >= 0.99)),
    sprintf("| Spearman -log10(padj) | %.4f | >= 0.99 | %s |",
            metrics$spearman_log_padj, tick(metrics$spearman_log_padj >= 0.99)),
    sprintf("| Jaccard sig (padj<0.05) | %.4f | >= 0.95 | %s |",
            metrics$jaccard_sig,   tick(metrics$jaccard_sig >= 0.95)),
    sprintf("| Sign agreement on sig | %.4f | >= 0.99 | %s |",
            metrics$sign_agreement_sig,
            tick(is.na(metrics$sign_agreement_sig) ||
                 metrics$sign_agreement_sig >= 0.99)),
    sprintf("| Kendall top-100 | %.4f | >= 0.95 | %s |",
            metrics$kendall_top100,
            tick(is.na(metrics$kendall_top100) ||
                 metrics$kendall_top100 >= 0.95)),
    sprintf("| No sig-set direction flips | %s | TRUE | %s |",
            !metrics$any_sig_flip, tick(!metrics$any_sig_flip)),
    "",
    "## Sig gene counts",
    sprintf("- Parametric: %d significant genes", metrics$n_sig_ref),
    sprintf("- glmGamPoi: %d significant genes", metrics$n_sig_alt),
    sprintf("- Shared significant: %d", metrics$n_sig_shared),
    ""
  )
  cat(paste(report, collapse = "\n"), "\n")

  if (!is.null(report_file)) {
    writeLines(report, report_file)
    cat(sprintf("[info] report written to %s\n", report_file))
  }

  invisible(list(metrics = metrics, ref = ref, alt = alt,
                 speedup = speedup, report = report))
}

tick <- function(x) if (isTRUE(x)) "PASS" else "FAIL"

# ------------------------------------------------------------------------------
# 4. testthat — enforces equivalence as CI gate
#
# Only runs under interactive sessions or when explicitly invoked via
# testthat::test_file() / RUN_TESTTHAT_INLINE=1 env var. This prevents the
# testthat block firing at file-scope during Rscript invocations, which
# previously made Rscript runs appear to use synthetic data (because
# testthat ran a synthetic benchmark BEFORE the Rscript argv block).
# ------------------------------------------------------------------------------
if (interactive() ||
    identical(Sys.getenv("RUN_TESTTHAT_INLINE"), "1") ||
    any(grepl("testthat", sys.calls(), fixed = TRUE))) {
  context("DESeq2 fitType equivalence")

  test_that("glmGamPoi fitType is equivalent to parametric on realistic pseudobulk", {
    skip_if_not_installed("glmGamPoi")

    built <- build_realistic_pseudobulk(n_donors = 50, n_genes = 3000,
                                         n_biased = 80, seed = 42)
    out <- benchmark_fittype_equivalence(
      built$sce,
      input_source = "synthetic fixture from build_realistic_pseudobulk (testthat)"
    )
    m <- out$metrics

    expect_gte(m$pearson_lfc,       0.99)
    expect_gte(m$spearman_log_padj, 0.99)
    expect_gte(m$jaccard_sig,       0.95)
    expect_gte(m$sign_agreement_sig, 0.99)
    expect_gte(m$kendall_top100,    0.95)
    expect_false(m$any_sig_flip)

    # Informational — not a gating assertion
    message(sprintf("[info] glmGamPoi speedup: %.2fx", out$speedup))
  })
}

# ------------------------------------------------------------------------------
# 5. Standalone invocation (Rscript)
# ------------------------------------------------------------------------------
if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  report <- if (length(args) >= 2) args[2] else "fittype_benchmark_report.md"

  # Load real-data input explicitly, and LOUDLY report any reason we fall back
  # to synthetic. Previously file-not-found, qread errors, and wrong object
  # shapes all silently produced synthetic data — this made debugging
  # a missing/renamed .qs file effectively impossible.
  sce_in <- NULL
  fallback_reason <- NULL

  if (length(args) < 1 || !nzchar(args[1])) {
    fallback_reason <- "no input path supplied as argv[1]"
    cat("[warn] ===========================================================\n")
    cat("[warn] FALLING BACK TO SYNTHETIC DATA\n")
    cat(sprintf("[warn] reason: %s\n", fallback_reason))
    cat("[warn] usage: Rscript test_fittype_equivalence.R <input.qs> [<report.md>]\n")
    cat("[warn] ===========================================================\n")
    input_source <- sprintf("synthetic (%s)", fallback_reason)
  } else {
    input_path <- args[1]
    cat(sprintf("[info] argv[1] = %s\n", input_path))
    cat(sprintf("[info] normalizePath: %s\n",
                tryCatch(normalizePath(input_path, mustWork = FALSE),
                         error = function(e) "<normalizePath failed>")))

    if (!file.exists(input_path)) {
      fallback_reason <- sprintf("file does not exist: %s", input_path)
      cat("[warn] ===========================================================\n")
      cat("[warn] FALLING BACK TO SYNTHETIC DATA\n")
      cat(sprintf("[warn] reason: %s\n", fallback_reason))
      cat("[warn] check that the path is correct and accessible from the\n")
      cat("[warn] current working directory / host.\n")
      cat("[warn] ===========================================================\n")
      input_source <- sprintf("synthetic (%s)", fallback_reason)
    } else {
      cat(sprintf("[info] file exists; size = %s bytes\n",
                  format(file.info(input_path)$size, big.mark = ",")))
      cat("[info] attempting qs::qread()...\n")

      obj <- tryCatch(
        qs::qread(input_path),
        error = function(e) {
          cat("[warn] ===========================================================\n")
          cat("[warn] qs::qread() FAILED — cannot use real data\n")
          cat(sprintf("[warn] error: %s\n", conditionMessage(e)))
          cat("[warn] ===========================================================\n")
          structure(NULL, qread_error = conditionMessage(e))
        }
      )

      if (is.null(obj) && !is.null(attr(obj, "qread_error"))) {
        fallback_reason <- sprintf("qs::qread failed: %s",
                                    attr(obj, "qread_error"))
        input_source <- sprintf("synthetic (%s)", fallback_reason)
      } else {
        cat(sprintf("[info] qread succeeded; object class = %s\n",
                    paste(class(obj), collapse = "/")))
        if (inherits(obj, "list")) {
          cat(sprintf("[info] list names: %s\n",
                      paste(names(obj), collapse = ", ")))
        }

        candidate <- if (inherits(obj, "list") && "sce" %in% names(obj)) {
          cat("[info] extracting obj$sce from list\n")
          obj$sce
        } else {
          obj
        }

        if (inherits(candidate, c("SingleCellExperiment",
                                   "SummarizedExperiment"))) {
          cat(sprintf(
            "[info] real-data object accepted (%s: %d genes x %d samples)\n",
            class(candidate)[1], nrow(candidate), ncol(candidate)))
          sce_in <- candidate
          input_source <- sprintf("real data from %s", input_path)
        } else {
          fallback_reason <- sprintf(
            "loaded object class '%s' is not a (Single)SummarizedExperiment; cannot use",
            paste(class(candidate), collapse = "/"))
          cat("[warn] ===========================================================\n")
          cat("[warn] FALLING BACK TO SYNTHETIC DATA\n")
          cat(sprintf("[warn] reason: %s\n", fallback_reason))
          cat("[warn] ===========================================================\n")
          input_source <- sprintf("synthetic (%s)", fallback_reason)
        }
      }
    }
  }

  benchmark_fittype_equivalence(sce = sce_in, report_file = report,
                                 input_source = input_source)
}
