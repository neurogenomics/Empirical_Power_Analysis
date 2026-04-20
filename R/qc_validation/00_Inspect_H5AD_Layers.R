# ==============================================================================
# SCRIPT: Inspect PsychAD H5AD File Layers
# PURPOSE: Determine whether adata.X contains raw counts, CellBender-corrected
#          counts, or normalized data. Check for adata.raw.X and named layers.
# INPUT: PsychAD h5ad files (PsychAD-HBCC.h5ad, PsychAD-MSSM.h5ad, PsychAD-RADC.h5ad)
# OUTPUT: Datalink/QC_Validation/h5ad_layer_report.txt
# ==============================================================================

if (!exists("PROJECT_ROOT")) source("config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

suppressPackageStartupMessages({
  library(reticulate)
})

message("=" |> rep(70) |> paste0(collapse = ""))
message("STEP 00: Inspect PsychAD H5AD File Layers")
message("=" |> rep(70) |> paste0(collapse = ""))

# Import Python modules
ad <- import("anndata")
np <- import("numpy")

# Output report
report_file <- file.path(OUTPUT_BASE, "h5ad_layer_report.txt")
dir.create(OUTPUT_BASE, recursive = TRUE, showWarnings = FALSE)
report_lines <- c(
  "=" |> rep(70) |> paste0(collapse = ""),
  "H5AD LAYER INSPECTION REPORT",
  sprintf("Generated: %s", Sys.time()),
  "=" |> rep(70) |> paste0(collapse = ""),
  ""
)

for (h5ad_name in names(H5AD_FILES)) {
  h5ad_path <- H5AD_FILES[[h5ad_name]]
  message(sprintf("\n>>> Inspecting: %s (%s)", h5ad_name, basename(h5ad_path)))

  if (!file.exists(h5ad_path)) {
    msg <- sprintf("[SKIP] File not found: %s", h5ad_path)
    message(msg)
    report_lines <- c(report_lines, msg, "")
    next
  }

  report_lines <- c(report_lines,
    "-" |> rep(60) |> paste0(collapse = ""),
    sprintf("FILE: %s (%s)", h5ad_name, basename(h5ad_path)),
    "-" |> rep(60) |> paste0(collapse = "")
  )

  tryCatch({
    # Open in backed mode (memory-efficient)
    adata <- ad$read_h5ad(h5ad_path, backed = "r")

    # Basic shape
    shape <- adata$shape
    report_lines <- c(report_lines,
      sprintf("  Shape: %d cells x %d genes", shape[[1]], shape[[2]])
    )
    message(sprintf("  Shape: %d cells x %d genes", shape[[1]], shape[[2]]))

    # Check adata.X properties
    report_lines <- c(report_lines, "", "  --- adata.X ---")

    # Sample a small chunk from X
    n_sample <- as.integer(min(1000L, shape[[1]]))
    x_sample <- adata$X[0:(n_sample - 1L), ]

    # Convert to dense if sparse
    scipy_sparse <- import("scipy.sparse")
    if (scipy_sparse$issparse(x_sample)) {
      x_dense <- x_sample$toarray()
      report_lines <- c(report_lines, "  Storage: sparse")
    } else {
      x_dense <- np$array(x_sample)
      report_lines <- c(report_lines, "  Storage: dense")
    }

    x_vals <- as.numeric(x_dense)
    x_vals <- x_vals[!is.na(x_vals)]

    x_dtype <- as.character(adata$X$dtype)
    x_min <- min(x_vals)
    x_max <- max(x_vals)
    x_mean <- mean(x_vals)
    x_nonzero <- x_vals[x_vals != 0]
    x_is_integer <- all(x_nonzero == floor(x_nonzero))
    x_sparsity <- mean(x_vals == 0) * 100

    report_lines <- c(report_lines,
      sprintf("  dtype: %s", x_dtype),
      sprintf("  Min: %.4f", x_min),
      sprintf("  Max: %.4f", x_max),
      sprintf("  Mean: %.4f", x_mean),
      sprintf("  Sparsity: %.1f%% zeros", x_sparsity),
      sprintf("  All non-zero values are integers: %s", x_is_integer),
      sprintf("  Non-zero value range: [%.4f, %.4f]", min(x_nonzero), max(x_nonzero)),
      sprintf("  Non-zero mean: %.4f", mean(x_nonzero))
    )

    # Interpretation
    if (x_is_integer && x_max > 100) {
      x_interpretation <- "LIKELY RAW/CELLBENDER-CORRECTED COUNTS (integer, high max)"
    } else if (!x_is_integer && x_max < 20) {
      x_interpretation <- "LIKELY LOG-NORMALIZED (non-integer, low max)"
    } else if (!x_is_integer && x_max > 100) {
      x_interpretation <- "LIKELY NORMALIZED BUT NOT LOG-TRANSFORMED (non-integer, high max)"
    } else {
      x_interpretation <- "AMBIGUOUS — needs further investigation"
    }
    report_lines <- c(report_lines, sprintf("  INTERPRETATION: %s", x_interpretation))
    message(sprintf("  adata.X interpretation: %s", x_interpretation))

    # Check for adata.raw
    report_lines <- c(report_lines, "", "  --- adata.raw ---")
    has_raw <- tryCatch({
      raw_obj <- adata$raw
      !is.null(raw_obj)
    }, error = function(e) FALSE)

    if (has_raw) {
      report_lines <- c(report_lines, "  adata.raw EXISTS")
      message("  adata.raw: EXISTS")

      tryCatch({
        raw_x_sample <- adata$raw$X[0:(n_sample - 1L), ]
        if (scipy_sparse$issparse(raw_x_sample)) {
          raw_dense <- raw_x_sample$toarray()
        } else {
          raw_dense <- np$array(raw_x_sample)
        }
        raw_vals <- as.numeric(raw_dense)
        raw_vals <- raw_vals[!is.na(raw_vals)]
        raw_nonzero <- raw_vals[raw_vals != 0]
        raw_is_integer <- all(raw_nonzero == floor(raw_nonzero))

        report_lines <- c(report_lines,
          sprintf("  raw.X shape: %d x %d", nrow(raw_dense), ncol(raw_dense)),
          sprintf("  raw.X Min: %.4f", min(raw_vals)),
          sprintf("  raw.X Max: %.4f", max(raw_vals)),
          sprintf("  raw.X Mean: %.4f", mean(raw_vals)),
          sprintf("  raw.X All integers: %s", raw_is_integer),
          sprintf("  raw.X Sparsity: %.1f%% zeros", mean(raw_vals == 0) * 100)
        )

        # Compare X vs raw.X
        if (x_is_integer && raw_is_integer) {
          are_identical <- all(x_vals == raw_vals)
          report_lines <- c(report_lines,
            sprintf("  X vs raw.X identical: %s", are_identical)
          )
          if (are_identical) {
            report_lines <- c(report_lines,
              "  CONCLUSION: X and raw.X are identical — no separate normalized layer"
            )
          } else {
            report_lines <- c(report_lines,
              "  CONCLUSION: X and raw.X differ — X may be CellBender-corrected, raw.X may be original"
            )
          }
        } else if (!x_is_integer && raw_is_integer) {
          report_lines <- c(report_lines,
            "  CONCLUSION: X is normalized (non-integer), raw.X is raw counts (integer)",
            "  ⚠️  SHARDING EXTRACTED NORMALIZED DATA, NOT RAW COUNTS"
          )
        } else {
          report_lines <- c(report_lines,
            sprintf("  CONCLUSION: Unusual — X integers=%s, raw.X integers=%s", x_is_integer, raw_is_integer)
          )
        }
      }, error = function(e) {
        report_lines <<- c(report_lines,
          sprintf("  [ERROR] Could not read raw.X: %s", e$message)
        )
      })
    } else {
      report_lines <- c(report_lines, "  adata.raw DOES NOT EXIST")
      message("  adata.raw: DOES NOT EXIST")
      report_lines <- c(report_lines,
        "  CONCLUSION: No separate raw layer — X is the only matrix"
      )
    }

    # Check for named layers
    report_lines <- c(report_lines, "", "  --- adata.layers ---")
    tryCatch({
      layer_keys <- adata$layers$keys()
      layer_names <- as.character(reticulate::iterate(layer_keys))
      if (length(layer_names) > 0) {
        report_lines <- c(report_lines,
          sprintf("  Named layers: %s", paste(layer_names, collapse = ", "))
        )
        message(sprintf("  Named layers: %s", paste(layer_names, collapse = ", ")))
      } else {
        report_lines <- c(report_lines, "  No named layers")
        message("  No named layers")
      }
    }, error = function(e) {
      report_lines <- c(report_lines,
        sprintf("  [ERROR] Could not read layers: %s", e$message)
      )
    })

    # Check obs columns for QC-related metadata
    report_lines <- c(report_lines, "", "  --- QC-related obs columns ---")
    tryCatch({
      obs_cols <- as.character(adata$obs$columns$tolist())
      qc_patterns <- c("n_genes", "n_counts", "n_UMIs", "mito", "doublet",
                        "scrublet", "cellbender", "qc", "QC", "total_counts",
                        "pct_counts", "is_cell")
      qc_cols <- obs_cols[grepl(paste(qc_patterns, collapse = "|"), obs_cols, ignore.case = TRUE)]
      if (length(qc_cols) > 0) {
        report_lines <- c(report_lines,
          sprintf("  QC-related columns found: %s", paste(qc_cols, collapse = ", "))
        )
      } else {
        report_lines <- c(report_lines, "  No QC-related columns found in obs")
      }
    }, error = function(e) {
      report_lines <- c(report_lines,
        sprintf("  [ERROR] Could not read obs columns: %s", e$message)
      )
    })

    # Close
    tryCatch(adata$file$close(), error = function(e) NULL)

    report_lines <- c(report_lines, "")

  }, error = function(e) {
    msg <- sprintf("[ERROR] Failed to inspect %s: %s", h5ad_name, e$message)
    message(msg)
    report_lines <<- c(report_lines, msg, "")
  })
}

# Write report
report_lines <- c(report_lines,
  "=" |> rep(70) |> paste0(collapse = ""),
  "END OF REPORT",
  "=" |> rep(70) |> paste0(collapse = "")
)

writeLines(report_lines, report_file)
message(sprintf("\n[DONE] Report written to: %s", report_file))
