#!/usr/bin/env Rscript
# ==============================================================================
# TEST: Job Planning Functions
# PURPOSE: Validate generate_job_plan() and sample_donors_for_job()
# STAGE: 1 (Unit Test)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 1 UNIT TEST: Job Planning Functions\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
})

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
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
    stop("Could not find config.R. Ensure working directory is project root or tests/ directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))

# Test counter
tests_passed <- 0
tests_failed <- 0

run_test <- function(name, expr) {
  cat(sprintf("[TEST] %s... ", name))
  result <- tryCatch({
    if (expr) {
      cat("PASSED\n")
      tests_passed <<- tests_passed + 1
      TRUE
    } else {
      cat("FAILED\n")
      tests_failed <<- tests_failed + 1
      FALSE
    }
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    tests_failed <<- tests_failed + 1
    FALSE
  })
  return(result)
}

# ==============================================================================
# 2. CREATE MOCK DATA
# ==============================================================================

cat("--- Creating mock validation data ---\n\n")

# Create mock validation donor data frame
# Simulates 50 males and 50 females with varying cell counts
set.seed(42)
mock_val_df <- data.frame(
  donor_id = paste0("DONOR_", sprintf("%03d", 1:100)),
  sex = rep(c(MALE_VAL, FEMALE_VAL), each = 50),
  cell_count = c(
    sample(100:5000, 50, replace = TRUE),  # Males
    sample(100:5000, 50, replace = TRUE)   # Females
  ),
  stringsAsFactors = FALSE
)

cat(sprintf("Mock data: %d donors (%d M, %d F)\n", 
            nrow(mock_val_df),
            sum(mock_val_df$sex == MALE_VAL),
            sum(mock_val_df$sex == FEMALE_VAL)))
cat(sprintf("Cell count range: %d - %d\n\n", 
            min(mock_val_df$cell_count), 
            max(mock_val_df$cell_count)))

# Mock paths
mock_run_dir <- tempdir()
mock_shard_dir <- "/fake/shard/dir"
mock_prop_male <- 0.5

# ==============================================================================
# 3. TEST generate_job_plan() - BASIC FUNCTIONALITY
# ==============================================================================

cat("--- Testing generate_job_plan() ---\n\n")

# Generate plan with small grid for testing
test_donor_steps <- c(5, 10, 20)
test_cell_steps <- c(100, 500, NA)
test_n_subsamples <- 2

plan <- generate_job_plan(
  val_df = mock_val_df,
  donor_steps = test_donor_steps,
  cell_steps = test_cell_steps,
  n_subsamples = test_n_subsamples,
  prop_male = mock_prop_male,
  run_dir = mock_run_dir,
  shard_dir = mock_shard_dir
)

run_test("Returns a list", {
  is.list(plan)
})

run_test("Plan is non-empty", {
  length(plan) > 0
})

# Expected: 3 donor steps x 3 cell steps x 2 iterations = 18 jobs max
# Some may be filtered out due to cell count requirements
run_test("Plan has reasonable number of jobs", {
  length(plan) <= 18 && length(plan) > 0
})

# ==============================================================================
# 4. TEST JOB STRUCTURE
# ==============================================================================

cat("\n--- Testing job structure ---\n\n")

first_job <- plan[[1]]

run_test("Job has n_donors field", {
  "n_donors" %in% names(first_job)
})

run_test("Job has n_cells field", {
  "n_cells" %in% names(first_job)
})

run_test("Job has iter field", {
  "iter" %in% names(first_job)
})

run_test("Job has n_male field", {
  "n_male" %in% names(first_job)
})

run_test("Job has n_female field", {
  "n_female" %in% names(first_job)
})

run_test("Job has male_pool field", {
  "male_pool" %in% names(first_job)
})

run_test("Job has female_pool field", {
  "female_pool" %in% names(first_job)
})

run_test("Job has out field (output path)", {
  "out" %in% names(first_job)
})

run_test("Job has log_out field (log path)", {
  "log_out" %in% names(first_job)
})

run_test("Job has shard_dir field", {
  "shard_dir" %in% names(first_job)
})

# ==============================================================================
# 5. TEST SEX RATIO PRESERVATION
# ==============================================================================

cat("\n--- Testing sex ratio preservation ---\n\n")

run_test("n_male + n_female equals n_donors", {
  all(sapply(plan, function(j) j$n_male + j$n_female == j$n_donors))
})

run_test("Sex ratio approximately matches population (within 1 donor)", {
  all(sapply(plan, function(j) {
    expected_male <- round(j$n_donors * mock_prop_male)
    abs(j$n_male - expected_male) <= 1
  }))
})

run_test("At least 1 male when n_donors > 0 and prop_male > 0", {
  all(sapply(plan, function(j) {
    if (j$n_donors > 0 && mock_prop_male > 0) j$n_male >= 1 else TRUE
  }))
})

# ==============================================================================
# 6. TEST CELL COUNT FILTERING
# ==============================================================================

cat("\n--- Testing cell count filtering ---\n\n")

# Jobs with high cell requirements should have fewer eligible donors
high_cell_jobs <- Filter(function(j) !is.na(j$n_cells) && j$n_cells >= 1000, plan)
low_cell_jobs <- Filter(function(j) !is.na(j$n_cells) && j$n_cells <= 200, plan)

if (length(high_cell_jobs) > 0 && length(low_cell_jobs) > 0) {
  run_test("High cell jobs have smaller donor pools than low cell jobs", {
    avg_high_pool <- mean(sapply(high_cell_jobs, function(j) length(j$male_pool) + length(j$female_pool)))
    avg_low_pool <- mean(sapply(low_cell_jobs, function(j) length(j$male_pool) + length(j$female_pool)))
    avg_high_pool <= avg_low_pool
  })
}

# NA cell jobs should have full pool
na_cell_jobs <- Filter(function(j) is.na(j$n_cells), plan)
if (length(na_cell_jobs) > 0) {
  run_test("NA cell jobs have full donor pool", {
    all(sapply(na_cell_jobs, function(j) {
      length(j$male_pool) == sum(mock_val_df$sex == MALE_VAL) &&
      length(j$female_pool) == sum(mock_val_df$sex == FEMALE_VAL)
    }))
  })
}

# ==============================================================================
# 7. TEST FILE NAMING
# ==============================================================================

cat("\n--- Testing file naming ---\n\n")

run_test("Output filename follows N{d}_C{c}_I{i}.qs pattern", {
  all(sapply(plan, function(j) {
    fname <- basename(j$out)
    grepl("^N\\d+_C(\\d+|NA)_I\\d+\\.qs$", fname)
  }))
})

run_test("Log filename follows N{d}_C{c}_I{i}_sampling_metadata.qs pattern", {
  all(sapply(plan, function(j) {
    fname <- basename(j$log_out)
    grepl("^N\\d+_C(\\d+|NA)_I\\d+_sampling_metadata\\.qs$", fname)
  }))
})

run_test("Filenames encode correct parameters", {
  # Check first job
  j <- plan[[1]]
  fname <- basename(j$out)
  cell_str <- ifelse(is.na(j$n_cells), "NA", as.character(j$n_cells))
  expected <- sprintf("N%d_C%s_I%d.qs", j$n_donors, cell_str, j$iter)
  fname == expected
})

# ==============================================================================
# 8. TEST sample_donors_for_job()
# ==============================================================================

cat("\n--- Testing sample_donors_for_job() ---\n\n")

run_test("Returns correct number of donors", {
  sampled <- sample_donors_for_job(first_job, seed = 42)
  length(sampled) == first_job$n_donors
})

run_test("Sampled donors are from the pools", {
  sampled <- sample_donors_for_job(first_job, seed = 42)
  all(sampled %in% c(first_job$male_pool, first_job$female_pool))
})

run_test("Same seed produces same sample", {
  sample1 <- sample_donors_for_job(first_job, seed = 42)
  sample2 <- sample_donors_for_job(first_job, seed = 42)
  identical(sample1, sample2)
})

run_test("Different seeds produce different samples", {
  sample1 <- sample_donors_for_job(first_job, seed = 42)
  sample2 <- sample_donors_for_job(first_job, seed = 123)
  !identical(sample1, sample2)
})

run_test("No duplicate donors in sample", {
  sampled <- sample_donors_for_job(first_job, seed = 42)
  length(sampled) == length(unique(sampled))
})

# ==============================================================================
# 9. TEST EDGE CASES
# ==============================================================================

cat("\n--- Testing edge cases ---\n\n")

# Test with insufficient donors
small_val_df <- mock_val_df[1:5, ]  # Only 5 donors total

small_plan <- generate_job_plan(
  val_df = small_val_df,
  donor_steps = c(10, 20, 50),  # All require more than available
  cell_steps = c(NA),
  n_subsamples = 1,
  prop_male = mock_prop_male,
  run_dir = mock_run_dir,
  shard_dir = mock_shard_dir
)

run_test("Returns empty plan when insufficient donors", {
  length(small_plan) == 0
})

# Test with very high cell requirement
high_cell_plan <- generate_job_plan(
  val_df = mock_val_df,
  donor_steps = c(5),
  cell_steps = c(10000),  # Higher than any donor has
  n_subsamples = 1,
  prop_male = mock_prop_male,
  run_dir = mock_run_dir,
  shard_dir = mock_shard_dir
)

run_test("Returns empty plan when cell requirement too high", {
  length(high_cell_plan) == 0
})

# ==============================================================================
# 10. TEST READS MODE
# ==============================================================================

cat("\n--- Testing reads mode ---\n\n")

# Add total_reads column to mock data for reads mode
set.seed(42)
mock_val_df$total_reads <- mock_val_df$cell_count * runif(nrow(mock_val_df), 3000, 7000)

# Test reads mode job planning
test_read_steps <- c(100000, 500000, NA)

plan_reads <- generate_job_plan(
  val_df = mock_val_df,
  donor_steps = c(5, 10),
  read_steps = test_read_steps,
  mode = "reads",
  n_subsamples = 2,
  prop_male = mock_prop_male,
  run_dir = mock_run_dir,
  shard_dir = mock_shard_dir
)

run_test("Reads mode returns a plan", {
  is.list(plan_reads) && length(plan_reads) > 0
})

run_test("Reads mode plan has correct structure", {
  first_job <- plan_reads[[1]]
  "target_reads" %in% names(first_job) && 
  "mode" %in% names(first_job) &&
  first_job$mode == "reads"
})

run_test("Reads mode uses R in filename", {
  first_job <- plan_reads[[1]]
  grepl("_R", basename(first_job$out))
})

run_test("Reads mode filters by total_reads", {
  # High read requirement should reduce eligible donors
  high_read_plan <- generate_job_plan(
    val_df = mock_val_df,
    donor_steps = c(5),
    read_steps = c(max(mock_val_df$total_reads) * 0.9),  # Very high requirement
    mode = "reads",
    n_subsamples = 1,
    prop_male = mock_prop_male,
    run_dir = mock_run_dir,
    shard_dir = mock_shard_dir
  )
  
  low_read_plan <- generate_job_plan(
    val_df = mock_val_df,
    donor_steps = c(5),
    read_steps = c(min(mock_val_df$total_reads) * 0.5),  # Low requirement
    mode = "reads",
    n_subsamples = 1,
    prop_male = mock_prop_male,
    run_dir = mock_run_dir,
    shard_dir = mock_shard_dir
  )
  
  # Low requirement should have more eligible donors (or at least same)
  length(low_read_plan) >= length(high_read_plan)
})

run_test("Reads mode NA uses all donors", {
  na_read_plan <- generate_job_plan(
    val_df = mock_val_df,
    donor_steps = c(5),
    read_steps = c(NA),
    mode = "reads",
    n_subsamples = 1,
    prop_male = mock_prop_male,
    run_dir = mock_run_dir,
    shard_dir = mock_shard_dir
  )
  
  length(na_read_plan) > 0 && 
  length(na_read_plan[[1]]$male_pool) + length(na_read_plan[[1]]$female_pool) == nrow(mock_val_df)
})

run_test("Reads mode filename encodes target_reads", {
  job <- plan_reads[[1]]
  fname <- basename(job$out)
  # Should match N{d}_R{reads}_I{iter}.qs pattern
  grepl("^N\\d+_R[0-9NA]+_I\\d+\\.qs$", fname)
})

# ==============================================================================
# 11. SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed\n", tests_passed, tests_failed))
cat("================================================================================\n")
cat("\n")

if (tests_failed > 0) {
  cat("  STATUS: FAILED\n")
  quit(status = 1)
} else {
  cat("  STATUS: PASSED\n")
  quit(status = 0)
}
