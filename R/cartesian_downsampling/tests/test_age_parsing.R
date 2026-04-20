#!/usr/bin/env Rscript
# ==============================================================================
# TEST: Age Parsing Functions
# PURPOSE: Validate parse_age_stage_to_years() and is_unknown_string()
# STAGE: 1 (Unit Test)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  STAGE 1 UNIT TEST: Age Parsing Functions\n")
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# 1. SETUP
# ==============================================================================

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

# Source dependencies
source(file.path(SCRIPT_DIR, "utils.R"))

# Test counter
tests_passed <- 0
tests_failed <- 0

# Test helper
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
# 2. TEST is_unknown_string()
# ==============================================================================

cat("\n--- Testing is_unknown_string() ---\n\n")

run_test("Detects 'unknown'", 
         is_unknown_string("unknown") == TRUE)

run_test("Detects 'Unknown' (case insensitive)", 
         is_unknown_string("Unknown") == TRUE)

run_test("Detects 'not available'", 
         is_unknown_string("not available") == TRUE)

run_test("Detects 'NA'", 
         is_unknown_string("NA") == TRUE)

run_test("Detects 'n/a'", 
         is_unknown_string("n/a") == TRUE)

run_test("Detects empty string", 
         is_unknown_string("") == TRUE)

run_test("Does not flag '65 year-old stage'", 
         is_unknown_string("65 year-old stage") == FALSE)

run_test("Does not flag numeric string '42'", 
         is_unknown_string("42") == FALSE)

run_test("Handles vector input", 
         all(is_unknown_string(c("unknown", "65 year-old stage")) == c(TRUE, FALSE)))

# ==============================================================================
# 3. TEST parse_age_stage_to_years()
# ==============================================================================

cat("\n--- Testing parse_age_stage_to_years() ---\n\n")

# Year parsing
run_test("Parses '65 year-old stage' to 65", 
         abs(parse_age_stage_to_years("65 year-old stage") - 65) < 0.01)

run_test("Parses '25 year-old stage' to 25", 
         abs(parse_age_stage_to_years("25 year-old stage") - 25) < 0.01)

run_test("Parses '100 year-old stage' to 100", 
         abs(parse_age_stage_to_years("100 year-old stage") - 100) < 0.01)

# Month parsing
run_test("Parses '6 month-old stage' to 0.5", 
         abs(parse_age_stage_to_years("6 month-old stage") - 0.5) < 0.01)

run_test("Parses '18 month-old stage' to 1.5", 
         abs(parse_age_stage_to_years("18 month-old stage") - 1.5) < 0.01)

# Week parsing
run_test("Parses '52 week-old stage' to ~1 year", 
         abs(parse_age_stage_to_years("52 week-old stage") - 1) < 0.05)

# Day parsing
run_test("Parses '365 day-old stage' to ~1 year", 
         abs(parse_age_stage_to_years("365 day-old stage") - 1) < 0.01)

# Unknown values
run_test("Returns NA for 'unknown'", 
         is.na(parse_age_stage_to_years("unknown")))

run_test("Returns NA for 'not available'", 
         is.na(parse_age_stage_to_years("not available")))

# Numeric fallback
run_test("Parses plain number '42' as years", 
         abs(parse_age_stage_to_years("42") - 42) < 0.01)

# Vector input
run_test("Handles vector input correctly", {
  input <- c("65 year-old stage", "unknown", "6 month-old stage")
  result <- parse_age_stage_to_years(input)
  abs(result[1] - 65) < 0.01 && is.na(result[2]) && abs(result[3] - 0.5) < 0.01
})

# Edge cases
run_test("Handles decimal ages '65.5 year-old stage'", 
         abs(parse_age_stage_to_years("65.5 year-old stage") - 65.5) < 0.01)

run_test("Handles whitespace '  65 year-old stage  '", 
         abs(parse_age_stage_to_years("  65 year-old stage  ") - 65) < 0.01)

# ==============================================================================
# 4. TEST standardize_age()
# ==============================================================================

cat("\n--- Testing standardize_age() ---\n\n")

run_test("Standardizes to mean 0", {
  ages <- c(20, 40, 60, 80)
  std <- standardize_age(ages, impute_na = FALSE)
  abs(mean(std)) < 0.01
})

run_test("Standardizes to variance 1", {
  ages <- c(20, 40, 60, 80)
  std <- standardize_age(ages, impute_na = FALSE)
  abs(sd(std) - 1) < 0.01
})

run_test("Imputes NA with mean when impute_na=TRUE", {
  ages <- c(20, 40, NA, 80)
  std <- standardize_age(ages, impute_na = TRUE)
  !any(is.na(std))
})

run_test("Preserves NA when impute_na=FALSE", {
  ages <- c(20, 40, NA, 80)
  std <- standardize_age(ages, impute_na = FALSE)
  is.na(std[3])
})

# ==============================================================================
# 5. SUMMARY
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
