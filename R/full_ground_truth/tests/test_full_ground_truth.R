# ==============================================================================
# STAGE 1 UNIT TESTS: Full Ground Truth Pipeline
# ==============================================================================
# Tests for config.R (age parsers) and utils.R (aggregation, DESeq2, diagnostics)
# Run with: Rscript -e "testthat::test_file('R/full_ground_truth/tests/test_full_ground_truth.R')"
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(Matrix)
  library(DESeq2)
  library(BiocParallel)
  library(qs)
})

# --- Source config.R (try multiple paths for flexibility) ---
config_paths <- c(
  "R/full_ground_truth/config.R",
  "config.R",
  "../config.R",
  file.path(Sys.getenv("PROJECT_ROOT", PROJECT_ROOT),
            "R", "full_ground_truth", "config.R")
)
config_found <- FALSE
for (cp in config_paths) {
  if (file.exists(cp)) {
    source(cp)
    config_found <- TRUE
    break
  }
}
if (!config_found) stop("Could not find config.R in any expected location")

# --- Source utils.R from SCRIPT_DIR ---
utils_path <- file.path(SCRIPT_DIR, "utils.R")
if (!file.exists(utils_path)) {
  stop(sprintf("utils.R not found at %s", utils_path))
}
source(utils_path)

# ==============================================================================
# HELPER: Create a minimal mock SCE with integer counts
# ==============================================================================
make_mock_sce <- function(n_genes, n_cells, donor_id, sex, age = "70",
                          extra_cols = list()) {
  set.seed(42 + nchar(donor_id))
  counts_mat <- matrix(
    as.integer(rpois(n_genes * n_cells, lambda = 5)),
    nrow = n_genes, ncol = n_cells
  )
  rownames(counts_mat) <- paste0("Gene", seq_len(n_genes))
  colnames(counts_mat) <- paste0(donor_id, "_cell", seq_len(n_cells))

  col_data <- DataFrame(
    donor_id   = rep(donor_id, n_cells),
    sex        = rep(sex, n_cells),
    development_stage = rep(age, n_cells)
  )

  # Add any extra metadata columns

for (nm in names(extra_cols)) {
    col_data[[nm]] <- rep(extra_cols[[nm]], n_cells)
  }

  SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = col_data
  )
}

# ==============================================================================
# 1. test_age_parsers
# ==============================================================================

test_that("parse_verbose_age handles year-old stage strings", {
  result <- parse_verbose_age("65 year-old stage")
  expect_true(is.numeric(result))
  expect_equal(result, 65)
})

test_that("parse_verbose_age handles month-old stage strings", {
  result <- parse_verbose_age("6 month-old stage")
  expect_true(is.numeric(result))
  expect_equal(result, 6 / 12, tolerance = 1e-4)
})

test_that("parse_verbose_age returns NA for 'unknown'", {
  result <- parse_verbose_age("unknown")
  expect_true(is.na(result))
})

test_that("parse_verbose_age returns NA for NA input", {
  result <- parse_verbose_age(NA)
  expect_true(is.na(result))
})

test_that("parse_verbose_age handles vector input correctly", {
  input <- c("65 year-old stage", "6 month-old stage", "unknown", "80 year-old stage")
  result <- parse_verbose_age(input)
  expect_length(result, 4)
  expect_equal(result[1], 65)
  expect_equal(result[2], 6 / 12, tolerance = 1e-4)
  expect_true(is.na(result[3]))
  expect_equal(result[4], 80)
})

test_that("parse_simple_numeric handles standard numeric strings", {
  expect_equal(parse_simple_numeric("80"), 80)
  expect_equal(parse_simple_numeric("80.5"), 80.5)
})

test_that("parse_simple_numeric returns NA for non-numeric strings", {
  result <- parse_simple_numeric("NA")
  expect_true(is.na(result))
})

test_that("parse_simple_numeric handles vector input", {
  input <- c("70", "80.5", "NA", "65")
  result <- parse_simple_numeric(input)
  expect_length(result, 4)
  expect_equal(result[1], 70)
  expect_equal(result[2], 80.5)
  expect_true(is.na(result[3]))
  expect_equal(result[4], 65)
})

# ==============================================================================
# 2. test_pseudobulk_aggregation (mock-based)
# ==============================================================================

test_that("aggregate_shards_parallel produces valid SCE from mock shards", {
  tmp_dir <- tempfile(pattern = "shards_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  n_genes <- 50
  n_cells_per_donor <- 10
  gene_names <- paste0("Gene", seq_len(n_genes))

  # Create 3 mock donor shards with known counts
  donor_info <- list(
    list(id = "donor_A", sex = "male", age = "65 year-old stage"),
    list(id = "donor_B", sex = "female", age = "70 year-old stage"),
    list(id = "donor_C", sex = "male", age = "75 year-old stage")
  )

  expected_sums <- list()
  for (d in donor_info) {
    sce <- make_mock_sce(n_genes, n_cells_per_donor, d$id, d$sex, d$age)
    rownames(sce) <- gene_names
    # Record expected sum per donor for verification
    expected_sums[[d$id]] <- as.integer(rowSums(assay(sce, "counts")))
    qs::qsave(sce, file.path(tmp_dir, paste0("donor_", d$id, ".qs")))
  }

  # Run aggregation
  design_vars <- c("donor_id", "sex", "development_stage")
  result <- aggregate_shards_parallel(
    shard_dir   = tmp_dir,
    design_vars = design_vars,
    n_cores     = 1
  )

  # Test 1: Output is a SingleCellExperiment
  expect_s4_class(result, "SingleCellExperiment")

  # Test 2: Correct dimensions (genes x donors)
  expect_equal(nrow(result), n_genes)
  expect_equal(ncol(result), 3)

  # Test 3: Counts are integers
  counts_out <- assay(result, "counts")
  expect_true(is.integer(counts_out))

  # Test 4: Sums match expected per donor
  # Each column should be the rowSums of that donor's cells
  for (i in seq_len(ncol(result))) {
    donor_col <- as.character(colData(result)$donor_id[i])
    if (donor_col %in% names(expected_sums)) {
      expect_equal(
        as.integer(counts_out[, i]),
        expected_sums[[donor_col]],
        info = paste("Sum mismatch for", donor_col)
      )
    }
  }
})

test_that("aggregate_shards_parallel errors on empty directory", {
  empty_dir <- tempfile(pattern = "empty_shards_")
  dir.create(empty_dir, recursive = TRUE)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  expect_error(
    aggregate_shards_parallel(empty_dir, c("donor_id", "sex"), n_cores = 1),
    "No valid shard files"
  )
})

test_that("aggregate_shards_parallel preserves metadata columns", {
  tmp_dir <- tempfile(pattern = "meta_shards_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  sce1 <- make_mock_sce(20, 5, "d1", "male", "65 year-old stage",
                         extra_cols = list(Source = "Study1"))
  sce2 <- make_mock_sce(20, 5, "d2", "female", "70 year-old stage",
                         extra_cols = list(Source = "Study1"))
  qs::qsave(sce1, file.path(tmp_dir, "donor_d1.qs"))
  qs::qsave(sce2, file.path(tmp_dir, "donor_d2.qs"))

  result <- aggregate_shards_parallel(
    shard_dir   = tmp_dir,
    design_vars = c("donor_id", "sex", "development_stage", "Source"),
    n_cores     = 1
  )

  expect_true("Source" %in% colnames(colData(result)))
  expect_true(all(colData(result)$Source == "Study1"))
})

test_that("aggregate_shards_parallel skips MANIFEST and donor_NA files", {
  tmp_dir <- tempfile(pattern = "skip_shards_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  sce <- make_mock_sce(20, 5, "d1", "male")
  qs::qsave(sce, file.path(tmp_dir, "donor_d1.qs"))

  # Create files that should be skipped
  qs::qsave(sce, file.path(tmp_dir, "MANIFEST.qs"))
  qs::qsave(sce, file.path(tmp_dir, "donor_NA.qs"))

  result <- aggregate_shards_parallel(
    shard_dir   = tmp_dir,
    design_vars = c("donor_id", "sex"),
    n_cores     = 1
  )

  # Only one real donor should be loaded
  expect_equal(ncol(result), 1)
})

# ==============================================================================
# 3. test_deseq2_setup (mock-based)
# ==============================================================================

test_that("run_deseq_for_dataset returns correct structure", {
  set.seed(123)
  n_genes  <- 100
  n_donors <- 6

  # Create pseudobulk with 3M, 3F and enough signal for DESeq2
  counts_mat <- matrix(
    as.integer(rpois(n_genes * n_donors, lambda = 100)),
    nrow = n_genes, ncol = n_donors
  )
  rownames(counts_mat) <- paste0("Gene", seq_len(n_genes))
  colnames(counts_mat) <- paste0("Donor", seq_len(n_donors))

  col_data <- DataFrame(
    donor_id          = paste0("Donor", seq_len(n_donors)),
    sex               = c("male", "male", "male", "female", "female", "female"),
    development_stage = paste0(c(65, 70, 75, 68, 72, 80), " year-old stage")
  )

  pb_sce <- SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = col_data
  )

  # Create a config matching Roussos format (simple, no extra covariates)
  mock_cfg <- list(
    cols = list(donor = "donor_id", sex = "sex", age = "development_stage"),
    covariates         = c(),
    numeric_covariates = c(),
    sex_values         = list(male = "male", female = "female"),
    age_parser         = parse_verbose_age
  )

  de_param <- BiocParallel::MulticoreParam(workers = 1, stop.on.error = TRUE)

  result <- run_deseq_for_dataset(pb_sce, mock_cfg, de_param, min_gene_count = 10)

  # Test 1: Returns a list with expected components
  expect_type(result, "list")
  expect_true(all(c("dds", "results", "res_df") %in% names(result)))

  # Test 2: res_df has required columns
  required_cols <- c("gene", "baseMean", "log2FoldChange", "padj", "pvalue")
  expect_true(all(required_cols %in% colnames(result$res_df)),
              info = paste("Missing cols:", paste(setdiff(required_cols, colnames(result$res_df)),
                                                  collapse = ", ")))

  # Test 3: padj values are in [0,1] or NA
  padj_vals <- result$res_df$padj
  valid_padj <- padj_vals[!is.na(padj_vals)]
  expect_true(all(valid_padj >= 0 & valid_padj <= 1),
              info = "All non-NA padj values should be in [0, 1]")
})

# ==============================================================================
# 4. test_output_structure
# ==============================================================================

test_that("mock res_df has required columns", {
  res_df <- data.frame(
    gene            = paste0("ENSG", 1:5),
    baseMean        = runif(5, 10, 1000),
    log2FoldChange  = rnorm(5),
    padj            = runif(5, 0, 0.1),
    pvalue          = runif(5, 0, 0.05),
    lfcSE           = abs(rnorm(5, 0, 0.1)),
    stringsAsFactors = FALSE
  )
  required_cols <- c("gene", "baseMean", "log2FoldChange", "padj", "pvalue", "lfcSE")
  expect_true(all(required_cols %in% colnames(res_df)))
})

test_that("clean_gene_ids strips version numbers", {
  input <- c("ENSG00000229807.5", "ENSG00000141510.12", "ENSG00000002330.1")
  result <- clean_gene_ids(input)
  expect_equal(result, c("ENSG00000229807", "ENSG00000141510", "ENSG00000002330"))
})

test_that("clean_gene_ids handles IDs without versions", {
  input <- c("ENSG00000229807", "GeneX")
  result <- clean_gene_ids(input)
  expect_equal(result, c("ENSG00000229807", "GeneX"))
})

test_that("extract_significant_genes filters correctly", {
  mock_results <- list(
    "DatasetA" = data.frame(
      gene            = paste0("Gene", 1:10),
      padj            = c(0.001, 0.01, 0.04, 0.06, 0.1, NA, 0.001, 0.03, 0.5, 0.02),
      log2FoldChange  = c(1.0, -0.5, 0.3, 0.8, -0.1, 0.5, -0.8, 0.26, 1.2, -0.3),
      row.names       = paste0("Gene", 1:10),
      stringsAsFactors = FALSE
    ),
    "DatasetB" = data.frame(
      gene            = paste0("Gene", 1:10),
      padj            = c(0.02, 0.03, NA, 0.001, 0.2, 0.3, 0.04, 0.01, 0.001, 0.5),
      log2FoldChange  = c(0.5, -1.0, 0.3, 0.6, -0.1, 0.2, 0.3, -0.4, 1.5, 0.0),
      row.names       = paste0("Gene", 1:10),
      stringsAsFactors = FALSE
    )
  )

  sig <- extract_significant_genes(mock_results, padj_cutoff = 0.05, lfc_cutoff = 0.25)

  expect_type(sig, "list")
  expect_length(sig, 2)
  expect_true("DatasetA" %in% names(sig))
  expect_true("DatasetB" %in% names(sig))

  # DatasetA: padj < 0.05 AND |lfc| > 0.25:
  #   Gene1:  0.001 & |1.0|  -> YES
  #   Gene2:  0.01  & |-0.5| -> YES
  #   Gene3:  0.04  & |0.3|  -> YES
  #   Gene7:  0.001 & |-0.8| -> YES
  #   Gene8:  0.03  & |0.26| -> YES
  #   Gene10: 0.02  & |-0.3| -> YES
  expect_true("Gene1" %in% sig$DatasetA)
  expect_true("Gene2" %in% sig$DatasetA)
  expect_true("Gene7" %in% sig$DatasetA)
  expect_false("Gene5" %in% sig$DatasetA)  # padj 0.1 > 0.05
  expect_false("Gene6" %in% sig$DatasetA)  # padj NA
})

test_that("build_overlap_matrix computes correct pairwise overlaps", {
  sig_list <- list(
    "A" = c("G1", "G2", "G3"),
    "B" = c("G2", "G3", "G4", "G5"),
    "C" = c("G1", "G5")
  )

  mat <- build_overlap_matrix(sig_list)

  expect_equal(nrow(mat), 3)
  expect_equal(ncol(mat), 3)
  expect_equal(rownames(mat), c("A", "B", "C"))

  # Diagonal = self-overlap = set size
  expect_equal(mat["A", "A"], 3)
  expect_equal(mat["B", "B"], 4)
  expect_equal(mat["C", "C"], 2)

  # A & B share G2, G3
  expect_equal(mat["A", "B"], 2)
  expect_equal(mat["B", "A"], 2)

  # A & C share G1
  expect_equal(mat["A", "C"], 1)

  # B & C share G5
  expect_equal(mat["B", "C"], 1)
})

test_that("compute_concordance returns correct metrics", {
  sig_list <- list(
    "Roussos" = c("G1", "G2", "G3", "G4", "G5"),
    "Allen"   = c("G1", "G2", "G6"),
    "DeJager" = c("G3", "G4", "G5", "G7")
  )

  conc <- compute_concordance(sig_list, reference = "Roussos")

  expect_s3_class(conc, "data.frame")
  expect_true(nrow(conc) == 2)  # Allen and DeJager (not Roussos itself)
  expect_true(all(c("Dataset", "Reference_Unique", "Dataset_Unique",
                     "Shared", "Concordance_Rate") %in% colnames(conc)))

  # Allen: shared = {G1, G2} = 2, ref_unique = {G3, G4, G5} = 3
  allen_row <- conc[conc$Dataset == "Allen", ]
  expect_equal(allen_row$Shared, 2)
  expect_equal(allen_row$Reference_Unique, 3)
  expect_equal(allen_row$Dataset_Unique, 1)  # G6
  expect_equal(allen_row$Concordance_Rate, 2 / 5, tolerance = 1e-10)

  # DeJager: shared = {G3, G4, G5} = 3
  dejager_row <- conc[conc$Dataset == "DeJager", ]
  expect_equal(dejager_row$Shared, 3)
  expect_equal(dejager_row$Concordance_Rate, 3 / 5, tolerance = 1e-10)
})

test_that("build_sign_agreement computes correct proportions", {
  lfc_matrix <- matrix(
    c(1.0, -0.5, 0.3, -0.2,   # Dataset A
      0.5, -1.0, -0.2, -0.8), # Dataset B
    nrow = 4, ncol = 2,
    dimnames = list(paste0("G", 1:4), c("A", "B"))
  )

  sa <- build_sign_agreement(lfc_matrix)

  expect_equal(nrow(sa), 2)
  expect_equal(ncol(sa), 2)

  # Diagonal should be 1 (perfect self-agreement)
  expect_equal(sa["A", "A"], 1.0)
  expect_equal(sa["B", "B"], 1.0)

  # A vs B: sign(A) = +, -, +, -; sign(B) = +, -, -, -
  # Agree on G1(+,+), G2(-,-), G4(-,-) = 3/4
  expect_equal(sa["A", "B"], 3 / 4, tolerance = 1e-10)
  expect_equal(sa["B", "A"], 3 / 4, tolerance = 1e-10)
})

test_that("build_lfc_matrix returns correct values", {
  mock_results <- list(
    "A" = data.frame(
      gene = c("G1", "G2", "G3"),
      log2FoldChange = c(1.5, -0.8, 0.3),
      padj = c(0.01, 0.05, 0.1),
      row.names = c("G1", "G2", "G3"),
      stringsAsFactors = FALSE
    ),
    "B" = data.frame(
      gene = c("G1", "G2", "G4"),
      log2FoldChange = c(0.9, -1.2, 0.5),
      padj = c(0.02, 0.03, 0.2),
      row.names = c("G1", "G2", "G4"),
      stringsAsFactors = FALSE
    )
  )

  gene_ids <- c("G1", "G2", "G3", "G4")
  mat <- build_lfc_matrix(gene_ids, mock_results)

  expect_equal(nrow(mat), 4)
  expect_equal(ncol(mat), 2)
  expect_equal(mat["G1", "A"], 1.5)
  expect_equal(mat["G2", "B"], -1.2)
  expect_true(is.na(mat["G3", "B"]))   # G3 not in B
  expect_true(is.na(mat["G4", "A"]))   # G4 not in A
  expect_equal(mat["G4", "B"], 0.5)
})

# ==============================================================================
# 5. test_no_silent_failures
# ==============================================================================

test_that("run_deseq_for_dataset errors on missing sex column", {
  set.seed(999)
  n_genes  <- 50
  n_donors <- 6

  counts_mat <- matrix(
    as.integer(rpois(n_genes * n_donors, lambda = 50)),
    nrow = n_genes, ncol = n_donors
  )
  rownames(counts_mat) <- paste0("Gene", seq_len(n_genes))
  colnames(counts_mat) <- paste0("D", seq_len(n_donors))

  # Missing sex column — only donor_id and age
  col_data <- DataFrame(
    donor_id          = paste0("D", seq_len(n_donors)),
    development_stage = rep("70 year-old stage", n_donors)
  )

  pb_sce <- SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = col_data
  )

  mock_cfg <- list(
    cols = list(donor = "donor_id", sex = "sex", age = "development_stage"),
    covariates         = c(),
    numeric_covariates = c(),
    sex_values         = list(male = "male", female = "female"),
    age_parser         = parse_verbose_age
  )

  de_param <- BiocParallel::MulticoreParam(workers = 1, stop.on.error = TRUE)

  expect_error(
    run_deseq_for_dataset(pb_sce, mock_cfg, de_param),
    "Sex column.*not found"
  )
})

test_that("run_deseq_for_dataset errors on insufficient replicates", {
  set.seed(888)
  n_genes  <- 50
  n_donors <- 2  # only 1M, 1F

  counts_mat <- matrix(
    as.integer(rpois(n_genes * n_donors, lambda = 50)),
    nrow = n_genes, ncol = n_donors
  )
  rownames(counts_mat) <- paste0("Gene", seq_len(n_genes))
  colnames(counts_mat) <- paste0("D", seq_len(n_donors))

  col_data <- DataFrame(
    donor_id          = paste0("D", seq_len(n_donors)),
    sex               = c("male", "female"),
    development_stage = c("65 year-old stage", "70 year-old stage")
  )

  pb_sce <- SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = col_data
  )

  mock_cfg <- list(
    cols = list(donor = "donor_id", sex = "sex", age = "development_stage"),
    covariates         = c(),
    numeric_covariates = c(),
    sex_values         = list(male = "male", female = "female"),
    age_parser         = parse_verbose_age
  )

  de_param <- BiocParallel::MulticoreParam(workers = 1, stop.on.error = TRUE)

  expect_error(
    run_deseq_for_dataset(pb_sce, mock_cfg, de_param),
    "Insufficient replicates"
  )
})

test_that("run_deseq_for_dataset errors on empty SCE", {
  # Zero genes, zero donors
  counts_mat <- matrix(integer(0), nrow = 0, ncol = 0)

  col_data <- DataFrame(
    donor_id          = character(0),
    sex               = character(0),
    development_stage = character(0)
  )

  pb_sce <- SingleCellExperiment(
    assays  = list(counts = counts_mat),
    colData = col_data
  )

  mock_cfg <- list(
    cols = list(donor = "donor_id", sex = "sex", age = "development_stage"),
    covariates         = c(),
    numeric_covariates = c(),
    sex_values         = list(male = "male", female = "female"),
    age_parser         = parse_verbose_age
  )

  de_param <- BiocParallel::MulticoreParam(workers = 1, stop.on.error = TRUE)

  # Should error — not return NA or NULL silently
  expect_error(
    run_deseq_for_dataset(pb_sce, mock_cfg, de_param)
  )
})

# ==============================================================================
# 6. test_mi_optimal_flip
# ==============================================================================

test_that("compute_mi_optimal_flip returns correct values for balanced case", {
  result <- compute_mi_optimal_flip(50, 50)

  # a* = round(50^2 / 100) = round(25) = 25

  expect_equal(result$a_star, 25)
  expect_equal(result$n_flip_males, 25)
  expect_equal(result$n_flip_females, 25)
  # P(AM|TM) = 25/50 = 0.5 = P(AM) = 50/100 = 0.5 -> deviation = 0
  expect_equal(result$mi_deviation, 0, tolerance = 1e-10)
})

test_that("compute_mi_optimal_flip returns correct values for imbalanced case", {
  result <- compute_mi_optimal_flip(40, 60)

  # a* = round(40^2 / 100) = round(16) = 16
  expect_equal(result$a_star, 16)
  expect_equal(result$n_flip_males, 24)
  expect_equal(result$n_flip_females, 24)
  # P(AM|TM) = 16/40 = 0.4 = P(AM) = 40/100 = 0.4 -> deviation = 0
  expect_equal(result$mi_deviation, 0, tolerance = 1e-10)
})

test_that("compute_mi_optimal_flip handles extreme imbalance", {
  result <- compute_mi_optimal_flip(10, 90)

  # a* = round(10^2 / 100) = round(1) = 1
  expect_equal(result$a_star, 1)
  expect_equal(result$n_flip_males, 9)
  expect_equal(result$n_flip_females, 9)
  # P(AM|TM) = 1/10 = 0.1 = P(AM) = 10/100 = 0.1 -> deviation = 0
  expect_equal(result$mi_deviation, 0, tolerance = 1e-10)
})

test_that("compute_mi_optimal_flip errors on zero donors", {
  expect_error(compute_mi_optimal_flip(0, 0), "No donors provided")
})

test_that("compute_mi_optimal_flip handles single-sex gracefully", {
  result <- compute_mi_optimal_flip(50, 0)
  expect_equal(result$n_flip_males, 0)
  expect_equal(result$n_flip_females, 0)
  expect_true(is.na(result$mi_deviation))
})

test_that("flip_sex_labels_mi_minimised preserves marginal ratio", {
  set.seed(42)
  meta_df <- data.frame(
    donor_id = paste0("D", 1:100),
    sex = c(rep("male", 40), rep("female", 60)),
    stringsAsFactors = FALSE
  )

  flipped <- flip_sex_labels_mi_minimised(
    meta_df    = meta_df,
    sex_col    = "sex",
    male_val   = "male",
    female_val = "female",
    seed       = 42
  )

  # Marginals preserved
  expect_equal(sum(flipped$sex == "male"), 40)
  expect_equal(sum(flipped$sex == "female"), 60)

  # Original sex saved
  expect_true("true_sex" %in% colnames(flipped))
  expect_equal(sum(flipped$true_sex == "male"), 40)
  expect_equal(sum(flipped$true_sex == "female"), 60)

  # Some labels actually changed
  n_changed <- sum(flipped$sex != flipped$true_sex)
  expect_true(n_changed > 0, info = "Some labels should have been flipped")

  # Check approximate independence: P(assigned_M | true_M) ≈ P(assigned_M)
  true_males <- flipped$true_sex == "male"
  p_am_given_tm <- mean(flipped$sex[true_males] == "male")
  p_am <- mean(flipped$sex == "male")
  expect_equal(p_am_given_tm, p_am, tolerance = 0.05)
})

test_that("flip_sex_labels_mi_minimised gives different results with different seeds", {
  meta_df <- data.frame(
    donor_id = paste0("D", 1:100),
    sex = c(rep("male", 40), rep("female", 60)),
    stringsAsFactors = FALSE
  )

  flipped1 <- flip_sex_labels_mi_minimised(meta_df, "sex", "male", "female", seed = 1)
  flipped2 <- flip_sex_labels_mi_minimised(meta_df, "sex", "male", "female", seed = 2)

  # Different seeds should produce different flips (with high probability)
  expect_false(identical(flipped1$sex, flipped2$sex),
               info = "Different seeds should produce different label assignments")
})

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n== Full Ground Truth Pipeline: Stage 1 Unit Tests Complete ==\n")
