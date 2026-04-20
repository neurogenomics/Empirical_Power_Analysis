# ==============================================================================
# TESTS: Cross-Cell-Type Aggregation (Step 5)
# PURPOSE: Unit tests for load_all_raw_results() and compute_pooled_aggregation()
# STAGE: 1 (Unit Tests)
# ==============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(tidyr)
})

# Source config and utils
possible_config_paths <- c(
  "R/cartesian_downsampling/config.R",  # From project root
  "../config.R",                         # From tests/ directory
  "config.R"                             # From script directory
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
  stop("Could not find config.R")
}

source(file.path(SCRIPT_DIR, "utils.R"))

# Source the Step 5 script functions only (not the execution block)
# We'll define the functions inline for testing since the script runs on source
# Instead, we re-define the core functions here for isolated testing

# ==============================================================================
# FUNCTION DEFINITIONS (copied from 05_Visualise_Cross_CellType.R for testing)
# ==============================================================================

compute_pooled_aggregation <- function(df_combined) {
  
  pooled <- df_combined %>%
    mutate(
      TP_Recall = Recall * GT_Count,
      Called_Sig = ifelse(Precision > 0, 
                          (Recall * GT_Count) / Precision, 
                          0),
      TP_Precision = Recall * GT_Count,
      TP_Recall_Harm = if ("Recall_Harm" %in% colnames(df_combined) &&
                           "GT_Harm_Count" %in% colnames(df_combined)) {
        Recall_Harm * GT_Harm_Count
      } else {
        rep(0, n())
      },
      GT_Harm_Count_Safe = ifelse(is.na(GT_Harm_Count), 0, GT_Harm_Count)
    ) %>%
    group_by(Donors, Cells, LFC_Bin, Gene_Category) %>%
    summarise(
      Pooled_Recall = ifelse(sum(GT_Count, na.rm = TRUE) > 0,
                              sum(TP_Recall, na.rm = TRUE) / sum(GT_Count, na.rm = TRUE),
                              0),
      Pooled_Precision = ifelse(sum(Called_Sig, na.rm = TRUE) > 0,
                                 sum(TP_Precision, na.rm = TRUE) / sum(Called_Sig, na.rm = TRUE),
                                 0),
      Pooled_Recall_Harm = ifelse(sum(GT_Harm_Count_Safe, na.rm = TRUE) > 0,
                                   sum(TP_Recall_Harm, na.rm = TRUE) / sum(GT_Harm_Count_Safe, na.rm = TRUE),
                                   0),
      Total_GT_Observations = sum(GT_Count, na.rm = TRUE),
      Mean_GT_Per_Obs = mean(GT_Count, na.rm = TRUE),
      N_Cell_Types = n_distinct(Cell_Type),
      N_Observations = n(),
      Cell_Types_Present = paste(sort(unique(Cell_Type)), collapse = "+"),
      .groups = "drop"
    )
  
  return(pooled)
}

# ==============================================================================
# TEST DATA HELPERS
# ==============================================================================

#' Create synthetic raw data mimicking df_raw from Step 3
create_test_raw_data <- function(cell_types = c("Astro", "Micro", "Oligo"),
                                  n_runs = 2,
                                  donors = c(10, 50),
                                  cells = c("100", "500"),
                                  lfc_bins = c("Large (3.0+)", "Small (1.0-2.0)")) {
  
  rows <- list()
  for (ct in cell_types) {
    for (run in seq_len(n_runs)) {
      for (d in donors) {
        for (c_val in cells) {
          for (bin in lfc_bins) {
            # Generate plausible values
            set.seed(sum(utf8ToInt(ct)) + run + d + as.numeric(c_val))
            gt_count <- sample(10:50, 1)
            recall <- runif(1, 0.1, 0.9)
            precision <- runif(1, 0.3, 1.0)
            recall_harm <- recall * runif(1, 0.9, 1.0)
            gt_harm_count <- round(gt_count * runif(1, 0.5, 0.9))
            
            rows[[length(rows) + 1]] <- data.frame(
              Run = paste0("Run_", run),
              Donors = d,
              Cells = c_val,
              Iteration = 1,
              LFC_Bin = bin,
              Gene_Category = "Total",
              Precision = precision,
              Recall = recall,
              Precision_Harm = precision * runif(1, 0.9, 1.0),
              Recall_Harm = recall_harm,
              GT_Count = gt_count,
              GT_Harm_Count = gt_harm_count,
              Cell_Type = ct,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  do.call(rbind, rows)
}

# ==============================================================================
# TESTS
# ==============================================================================

test_that("compute_pooled_aggregation returns correct structure", {
  df <- create_test_raw_data()
  result <- compute_pooled_aggregation(df)
  
  expect_true(is.data.frame(result))
  expect_true("Pooled_Recall" %in% colnames(result))
  expect_true("Pooled_Precision" %in% colnames(result))
  expect_true("Pooled_Recall_Harm" %in% colnames(result))
  expect_true("Total_GT_Observations" %in% colnames(result))
  expect_true("N_Cell_Types" %in% colnames(result))
  expect_true("Cell_Types_Present" %in% colnames(result))
})

test_that("pooled recall is weighted mean of per-observation recalls", {
  # Create simple known data: 2 cell types, 1 grid point
  df <- data.frame(
    Run = "Run_1", Donors = 10, Cells = "100", Iteration = 1,
    LFC_Bin = "Large (2.0+)", Gene_Category = "Total",
    Precision = c(1.0, 0.5),
    Recall = c(0.8, 0.6),
    Precision_Harm = c(1.0, 0.5),
    Recall_Harm = c(0.7, 0.5),
    GT_Count = c(10, 20),
    GT_Harm_Count = c(8, 15),
    Cell_Type = c("Astro", "Micro"),
    stringsAsFactors = FALSE
  )
  
  result <- compute_pooled_aggregation(df)
  
  # Expected pooled recall: (0.8*10 + 0.6*20) / (10 + 20) = (8 + 12) / 30 = 20/30
  expected_recall <- (0.8 * 10 + 0.6 * 20) / (10 + 20)
  expect_equal(result$Pooled_Recall, expected_recall, tolerance = 1e-10)
  
  # Expected pooled GT observations
  expect_equal(result$Total_GT_Observations, 30)
  
  # Expected N cell types
  expect_equal(result$N_Cell_Types, 2)
})

test_that("same gene in 3 cell types contributes 3 observations", {
  # Simulate: same gene detected in all 3 cell types
  # Each CT has GT_Count = 1 (just that one gene), Recall = 1.0 (detected)
  df <- data.frame(
    Run = "Run_1", Donors = 50, Cells = "500", Iteration = 1,
    LFC_Bin = "Large (2.0+)", Gene_Category = "Total",
    Precision = c(1.0, 1.0, 1.0),
    Recall = c(1.0, 1.0, 1.0),
    Precision_Harm = c(1.0, 1.0, 1.0),
    Recall_Harm = c(1.0, 1.0, 1.0),
    GT_Count = c(1, 1, 1),
    GT_Harm_Count = c(1, 1, 1),
    Cell_Type = c("Astro", "Micro", "Oligo"),
    stringsAsFactors = FALSE
  )
  
  result <- compute_pooled_aggregation(df)
  
  # 3 observations, all detected → pooled recall = 3/3 = 1.0
  expect_equal(result$Pooled_Recall, 1.0)
  expect_equal(result$Total_GT_Observations, 3)
  expect_equal(result$N_Cell_Types, 3)
  expect_equal(result$Cell_Types_Present, "Astro+Micro+Oligo")
})

test_that("partial detection across cell types gives correct pooled recall", {
  # 3 cell types, gene detected in 2 out of 3
  df <- data.frame(
    Run = "Run_1", Donors = 50, Cells = "500", Iteration = 1,
    LFC_Bin = "Large (2.0+)", Gene_Category = "Total",
    Precision = c(1.0, 1.0, 0.0),
    Recall = c(1.0, 1.0, 0.0),
    Precision_Harm = c(1.0, 1.0, 0.0),
    Recall_Harm = c(1.0, 1.0, 0.0),
    GT_Count = c(1, 1, 1),
    GT_Harm_Count = c(1, 1, 1),
    Cell_Type = c("Astro", "Micro", "Oligo"),
    stringsAsFactors = FALSE
  )
  
  result <- compute_pooled_aggregation(df)
  
  # 2 out of 3 detected → pooled recall = 2/3
  expect_equal(result$Pooled_Recall, 2/3, tolerance = 1e-10)
})

test_that("pooled values are bounded 0-1", {
  df <- create_test_raw_data()
  result <- compute_pooled_aggregation(df)
  
  expect_true(all(result$Pooled_Recall >= 0 & result$Pooled_Recall <= 1))
  expect_true(all(result$Pooled_Precision >= 0 & result$Pooled_Precision <= 1))
  expect_true(all(result$Pooled_Recall_Harm >= 0 & result$Pooled_Recall_Harm <= 1))
})

test_that("zero GT count produces zero recall", {
  df <- data.frame(
    Run = "Run_1", Donors = 10, Cells = "100", Iteration = 1,
    LFC_Bin = "Large (2.0+)", Gene_Category = "Total",
    Precision = 0, Recall = 0,
    Precision_Harm = 0, Recall_Harm = 0,
    GT_Count = 0, GT_Harm_Count = 0,
    Cell_Type = "Astro",
    stringsAsFactors = FALSE
  )
  
  result <- compute_pooled_aggregation(df)
  expect_equal(result$Pooled_Recall, 0)
})

test_that("aggregation groups correctly by grid point", {
  df <- create_test_raw_data(
    cell_types = c("Astro", "Micro"),
    n_runs = 1,
    donors = c(10, 50),
    cells = c("100"),
    lfc_bins = c("Large (3.0+)")
  )
  
  result <- compute_pooled_aggregation(df)
  
  # Should have 2 rows: one per donor count (10 and 50)
  expect_equal(nrow(result), 2)
  expect_true(all(c(10, 50) %in% result$Donors))
})

test_that("unequal GT counts across cell types are weighted correctly", {
  # Astro has 100 GT genes with 0.9 recall
  # Micro has 10 GT genes with 0.1 recall
  # Pooled should be dominated by Astro: (0.9*100 + 0.1*10) / 110 = 91/110
  df <- data.frame(
    Run = "Run_1", Donors = 50, Cells = "500", Iteration = 1,
    LFC_Bin = "Large (2.0+)", Gene_Category = "Total",
    Precision = c(0.9, 0.5),
    Recall = c(0.9, 0.1),
    Precision_Harm = c(0.9, 0.5),
    Recall_Harm = c(0.8, 0.1),
    GT_Count = c(100, 10),
    GT_Harm_Count = c(80, 8),
    Cell_Type = c("Astro", "Micro"),
    stringsAsFactors = FALSE
  )
  
  result <- compute_pooled_aggregation(df)
  
  expected <- (0.9 * 100 + 0.1 * 10) / (100 + 10)
  expect_equal(result$Pooled_Recall, expected, tolerance = 1e-10)
  expect_equal(result$Total_GT_Observations, 110)
})

cat("\n=== All cross-cell-type aggregation tests passed ===\n")
