# ==============================================================================
# UNIT TESTS: Shard Validation Framework
# DESC: Tests to verify shard integrity and identify issues
# USAGE: source("R/sharding/tests/test_sharding.R")
# ==============================================================================

library(testthat)
library(qs)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

SHARD_BASE_DIR <- file.path(PROJECT_ROOT, "Datalink", "Shards")
DATASETS <- c("Roussos", "Roussos_Full", "Allen", "DeJager", "Johansen", "Rexach")
CELL_TYPES <- c("Micro", "Oligo", "Astro")

# Expected standardized metadata columns
REQUIRED_COLS <- c("donor_id")  # Must exist in all shards
RECOMMENDED_COLS <- c("age", "sex", "cell_type")  # Should exist if available in source

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Find all shard files for a dataset
find_shards <- function(cell_type, dataset) {
  shard_dir <- file.path(SHARD_BASE_DIR, cell_type, dataset)
  if (!dir.exists(shard_dir)) return(character(0))
  list.files(shard_dir, pattern = "^donor_.*\\.qs$", full.names = TRUE)
}

#' Find manifest for a dataset
find_manifest <- function(cell_type, dataset) {
  manifest_path <- file.path(SHARD_BASE_DIR, cell_type, dataset, "00_MASTER_MANIFEST.qs")
  if (file.exists(manifest_path)) return(manifest_path)
  return(NULL)
}

#' Load a single shard safely
load_shard <- function(shard_path) {
  tryCatch({
    qs::qread(shard_path)
  }, error = function(e) {
    stop(sprintf("Failed to load shard %s: %s", shard_path, e$message))
  })
}

#' Load manifest safely
load_manifest <- function(manifest_path) {
  tryCatch({
    qs::qread(manifest_path)
  }, error = function(e) {
    stop(sprintf("Failed to load manifest %s: %s", manifest_path, e$message))
  })
}

# ==============================================================================
# TEST SUITE 1: SHARD STRUCTURE & INTEGRITY
# ==============================================================================

context("Shard Structure & Integrity")

test_that("All shards are loadable SingleCellExperiment objects", {
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      # Test first shard
      test_shard <- load_shard(shards[1])
      expect_is(test_shard, "SingleCellExperiment",
                info = sprintf("%s/%s: Shard is not SingleCellExperiment", cell_type, dataset))
      
      # Clean up
      rm(test_shard); gc()
    }
  }
})

test_that("Shards contain required metadata column: donor_id", {
  issues <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      for (shard_file in shards) {
        test_shard <- load_shard(shard_file)
        col_names <- names(colData(test_shard))
        
        if (!"donor_id" %in% col_names) {
          donor_cols <- col_names[grepl("donor|individual|subject", col_names, ignore.case = TRUE)]
          issues[[sprintf("%s/%s/%s", cell_type, dataset, basename(shard_file))]] <- list(
            missing = "donor_id",
            available = donor_cols
          )
        }
        
        rm(test_shard); gc()
      }
    }
  }
  
  if (length(issues) > 0) {
    # Print detailed info for debugging
    message("\n=== DONOR_ID COLUMN ISSUES ===")
    for (name in names(issues)) {
      message(sprintf("%s:", name))
      message(sprintf("  Missing: 'donor_id'"))
      message(sprintf("  Similar columns found: %s", 
                      ifelse(length(issues[[name]]$available) > 0,
                             paste(issues[[name]]$available, collapse = ", "),
                             "NONE")))
    }
    message("============================\n")
  }
  
  expect_equal(length(issues), 0, 
               info = sprintf("Found %d shards without 'donor_id' column. See details above.", 
                             length(issues)))
})

test_that("Count matrices contain integer values (raw counts, not normalized)", {
  non_integer_shards <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      # Test sample of shards (first 3 per dataset)
      test_shards <- head(shards, 3)
      
      for (shard_file in test_shards) {
        test_shard <- load_shard(shard_file)
        
        if (!"counts" %in% assayNames(test_shard)) {
          non_integer_shards[[basename(shard_file)]] <- "Missing 'counts' assay"
          rm(test_shard); gc()
          next
        }
        
        counts_mat <- assay(test_shard, "counts")
        
        # Check for non-integer values
        if (inherits(counts_mat, "dgCMatrix")) {
          # For sparse matrix, check the @x slot
          non_int_rate <- mean(counts_mat@x != floor(counts_mat@x), na.rm = TRUE)
        } else if (is.matrix(counts_mat) || inherits(counts_mat, "Matrix")) {
          non_int_rate <- mean(counts_mat != floor(counts_mat), na.rm = TRUE)
        } else {
          non_integer_shards[[basename(shard_file)]] <- sprintf("Unknown matrix type: %s", 
                                                                  class(counts_mat))
          rm(test_shard); gc()
          next
        }
        
        if (non_int_rate > 0.01) {  # More than 1% non-integer
          non_integer_shards[[sprintf("%s/%s/%s", cell_type, dataset, basename(shard_file))]] <-
            sprintf("%.2f%% non-integer values", non_int_rate * 100)
        }
        
        rm(test_shard, counts_mat); gc()
      }
    }
  }
  
  if (length(non_integer_shards) > 0) {
    message("\n=== NON-INTEGER COUNT MATRIX ISSUES ===")
    for (name in names(non_integer_shards)) {
      message(sprintf("%s: %s", name, non_integer_shards[[name]]))
    }
    message("======================================\n")
  }
  
  expect_equal(length(non_integer_shards), 0,
               info = sprintf("Found %d shards with non-integer counts (likely normalized, not raw). See details above.",
                             length(non_integer_shards)))
})

test_that("Shard dimensions are reasonable for scRNA-seq", {
  dimension_issues <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      # Test sample of shards
      test_shards <- head(shards, 3)
      
      for (shard_file in test_shards) {
        test_shard <- load_shard(shard_file)
        
        n_genes <- nrow(test_shard)
        n_cells <- ncol(test_shard)
        
        # Reasonable ranges for single-cell
        if (n_genes < 1000) {
          dimension_issues[[basename(shard_file)]] <- sprintf("Too few genes: %d", n_genes)
        } else if (n_genes > 80000) {
          dimension_issues[[basename(shard_file)]] <- sprintf("Too many genes: %d", n_genes)
        }
        
        if (n_cells < 1) {
          dimension_issues[[basename(shard_file)]] <- sprintf("No cells: %d", n_cells)
        } else if (n_cells > 100000) {
          dimension_issues[[basename(shard_file)]] <- sprintf("Suspiciously many cells: %d", n_cells)
        }
        
        rm(test_shard); gc()
      }
    }
  }
  
  if (length(dimension_issues) > 0) {
    message("\n=== DIMENSION ISSUES ===")
    for (name in names(dimension_issues)) {
      message(sprintf("%s: %s", name, dimension_issues[[name]]))
    }
    message("========================\n")
  }
  
  expect_equal(length(dimension_issues), 0,
               info = sprintf("Found %d shards with suspicious dimensions. See details above.",
                             length(dimension_issues)))
})

# ==============================================================================
# TEST SUITE 2: MANIFEST INTEGRITY
# ==============================================================================

context("Manifest Integrity")

test_that("Manifest files exist for all datasets with shards", {
  missing_manifests <- character(0)
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      manifest <- find_manifest(cell_type, dataset)
      if (is.null(manifest)) {
        missing_manifests <- c(missing_manifests, sprintf("%s/%s", cell_type, dataset))
      }
    }
  }
  
  expect_equal(length(missing_manifests), 0,
               info = sprintf("Missing manifests for: %s", paste(missing_manifests, collapse = ", ")))
})

test_that("Manifest donor_ids match actual shard files", {
  mismatches <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      manifest_path <- find_manifest(cell_type, dataset)
      if (is.null(manifest_path)) next
      
      manifest <- load_manifest(manifest_path)
      if (!"manifest" %in% names(manifest)) {
        mismatches[[sprintf("%s/%s", cell_type, dataset)]] <- "Manifest missing 'manifest' element"
        next
      }
      
      manifest_donors <- manifest$manifest$donor_id
      shards <- find_shards(cell_type, dataset)
      shard_files <- basename(shards)
      
      # Extract donor IDs from shard filenames (donor_<safe_id>.qs)
      expected_files <- paste0("donor_", gsub("[^[:alnum:]]", "_", manifest_donors), ".qs")
      
      # Check for files in manifest but missing on disk
      missing_on_disk <- setdiff(expected_files, shard_files)
      # Check for files on disk but not in manifest
      missing_in_manifest <- setdiff(shard_files, expected_files)
      
      if (length(missing_on_disk) > 0 || length(missing_in_manifest) > 0) {
        mismatches[[sprintf("%s/%s", cell_type, dataset)]] <- list(
          manifest_donors = length(manifest_donors),
          shard_files = length(shard_files),
          missing_on_disk = missing_on_disk,
          missing_in_manifest = missing_in_manifest
        )
      }
    }
  }
  
  if (length(mismatches) > 0) {
    message("\n=== MANIFEST/SHARD MISMATCHES ===")
    for (name in names(mismatches)) {
      info <- mismatches[[name]]
      if (is.character(info)) {
        message(sprintf("%s: %s", name, info))
      } else {
        message(sprintf("%s:", name))
        message(sprintf("  Manifest donors: %d, Shard files: %d", 
                        info$manifest_donors, info$shard_files))
        if (length(info$missing_on_disk) > 0) {
          message(sprintf("  Missing on disk: %s", paste(info$missing_on_disk, collapse = ", ")))
        }
        if (length(info$missing_in_manifest) > 0) {
          message(sprintf("  Missing in manifest: %s", paste(info$missing_in_manifest, collapse = ", ")))
        }
      }
    }
    message("==================================\n")
  }
  
  expect_equal(length(mismatches), 0,
               info = sprintf("Found %d manifest/shard mismatches. See details above.",
                             length(mismatches)))
})

test_that("Manifest cell counts match actual shard dimensions", {
  count_mismatches <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      manifest_path <- find_manifest(cell_type, dataset)
      if (is.null(manifest_path)) next
      
      manifest <- load_manifest(manifest_path)
      if (!"manifest" %in% names(manifest)) next
      
      manifest_df <- manifest$manifest
      
      # Test subset of donors
      test_donors <- head(manifest_df$donor_id, 5)
      
      for (donor_id in test_donors) {
        safe_id <- gsub("[^[:alnum:]]", "_", donor_id)
        shard_file <- file.path(SHARD_BASE_DIR, cell_type, dataset, 
                                 paste0("donor_", safe_id, ".qs"))
        
        if (!file.exists(shard_file)) {
          count_mismatches[[sprintf("%s/%s/%s", cell_type, dataset, donor_id)]] <- 
            "Shard file missing"
          next
        }
        
        test_shard <- load_shard(shard_file)
        actual_cells <- ncol(test_shard)
        manifest_cells <- manifest_df$cell_count[manifest_df$donor_id == donor_id]
        
        if (length(manifest_cells) == 0) {
          count_mismatches[[sprintf("%s/%s/%s", cell_type, dataset, donor_id)]] <-
            "Donor not in manifest"
        } else if (actual_cells != manifest_cells[1]) {
          count_mismatches[[sprintf("%s/%s/%s", cell_type, dataset, donor_id)]] <-
            sprintf("Manifest: %d, Actual: %d", manifest_cells[1], actual_cells)
        }
        
        rm(test_shard); gc()
      }
    }
  }
  
  if (length(count_mismatches) > 0) {
    message("\n=== CELL COUNT MISMATCHES ===")
    for (name in names(count_mismatches)) {
      message(sprintf("%s: %s", name, count_mismatches[[name]]))
    }
    message("=============================\n")
  }
  
  expect_equal(length(count_mismatches), 0,
               info = sprintf("Found %d cell count mismatches. See details above.",
                             length(count_mismatches)))
})

# ==============================================================================
# TEST SUITE 3: DATASET CONSISTENCY
# ==============================================================================

context("Dataset Consistency")

test_that("Donor ID naming conventions are consistent within datasets", {
  naming_issues <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      # Extract donor IDs from first few shards
      test_shards <- head(shards, 5)
      donor_ids <- character(0)
      
      for (shard_file in test_shards) {
        test_shard <- load_shard(shard_file)
        if ("donor_id" %in% names(colData(test_shard))) {
          donor_ids <- c(donor_ids, unique(colData(test_shard)$donor_id))
        }
        rm(test_shard); gc()
      }
      
      # Check for inconsistent naming patterns
      has_prefix <- grepl("^(Donor|donor|Individual|Subject)", donor_ids)
      if (any(has_prefix) && !all(has_prefix)) {
        naming_issues[[sprintf("%s/%s", cell_type, dataset)]] <- 
          sprintf("Mixed naming: some with prefix, some without. Examples: %s", 
                  paste(head(donor_ids, 3), collapse = ", "))
      }
    }
  }
  
  if (length(naming_issues) > 0) {
    message("\n=== DONOR ID NAMING ISSUES ===")
    for (name in names(naming_issues)) {
      message(sprintf("%s: %s", name, naming_issues[[name]]))
    }
    message("==============================\n")
  }
  
  expect_equal(length(naming_issues), 0,
               info = sprintf("Found %d datasets with inconsistent donor ID naming. See details above.",
                             length(naming_issues)))
})

test_that("All datasets have minimum required donors per cell type", {
  low_donor_datasets <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      manifest_path <- find_manifest(cell_type, dataset)
      if (is.null(manifest_path)) {
        low_donor_datasets[[sprintf("%s/%s", cell_type, dataset)]] <- "No manifest"
        next
      }
      
      manifest <- load_manifest(manifest_path)
      n_donors <- nrow(manifest$manifest)
      
      if (n_donors < 10) {
        low_donor_datasets[[sprintf("%s/%s", cell_type, dataset)]] <- 
          sprintf("Only %d donors", n_donors)
      }
    }
  }
  
  if (length(low_donor_datasets) > 0) {
    message("\n=== LOW DONOR COUNT ===")
    for (name in names(low_donor_datasets)) {
      message(sprintf("%s: %s", name, low_donor_datasets[[name]]))
    }
    message("=======================\n")
  }
  
  # Warning only - some datasets may legitimately have few donors
  if (length(low_donor_datasets) > 0) {
    message("WARNING: Found datasets with < 10 donors. This may be expected for some cell types.")
  }
})

# ==============================================================================
# TEST SUITE 4: METADATA PRESERVATION
# ==============================================================================

context("Metadata Preservation")

test_that("Critical metadata columns are preserved", {
  metadata_issues <- list()
  
  for (cell_type in CELL_TYPES) {
    for (dataset in DATASETS) {
      shards <- find_shards(cell_type, dataset)
      if (length(shards) == 0) next
      
      # Test one shard per dataset
      test_shard <- load_shard(shards[1])
      col_names <- names(colData(test_shard))
      
      # Check for recommended columns
      missing_recommended <- setdiff(RECOMMENDED_COLS, col_names)
      
      if (length(missing_recommended) > 0) {
        metadata_issues[[sprintf("%s/%s", cell_type, dataset)]] <- 
          sprintf("Missing recommended: %s", paste(missing_recommended, collapse = ", "))
      }
      
      rm(test_shard); gc()
    }
  }
  
  if (length(metadata_issues) > 0) {
    message("\n=== METADATA PRESERVATION ISSUES ===")
    for (name in names(metadata_issues)) {
      message(sprintf("%s: %s", name, metadata_issues[[name]]))
    }
    message("=====================================\n")
  }
  
  # Warning only - not all metadata may be available in all sources
  if (length(metadata_issues) > 0) {
    message("INFO: Some recommended metadata columns missing. This may be expected depending on source data.")
  }
})

# ==============================================================================
# RUN TESTS
# ==============================================================================

message("\n========================================")
message("SHARD VALIDATION TEST SUITE")
message("========================================\n")

message("Testing datasets: ", paste(DATASETS, collapse = ", "))
message("Testing cell types: ", paste(CELL_TYPES, collapse = ", "))
message("")

# Run all tests
test_results <- test_dir("R/sharding/tests", reporter = "progress")

message("\n========================================")
message("TEST COMPLETE")
message("========================================")
