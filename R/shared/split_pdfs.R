# ==============================================================================
# SCRIPT: Split Dataset PDFs into Per-Page Files
# PURPOSE: Split large PDFs in DatasetPapers/ into individual per-page PDFs
#          so they can be read page-by-page by AI tools
# INPUT: DatasetPapers/*.pdf (7 dataset papers)
# OUTPUT: DatasetPapers/extracted/{Author}_p{NN}.pdf (one file per page)
# ==============================================================================

suppressPackageStartupMessages({
  library(qpdf)
})

# --- Configuration -----------------------------------------------------------
INPUT_DIR  <- "DatasetPapers"
OUTPUT_DIR <- file.path(INPUT_DIR, "extracted")

# --- Main --------------------------------------------------------------------
main <- function() {
  # Validate input directory
  if (!dir.exists(INPUT_DIR)) {
    stop(sprintf("[ERROR] Input directory '%s' does not exist. Run from project root.", INPUT_DIR))
  }

  # Find all PDFs
  pdf_files <- list.files(INPUT_DIR, pattern = "\\.pdf$", full.names = TRUE)
  if (length(pdf_files) == 0) {
    stop(sprintf("[ERROR] No PDF files found in '%s'", INPUT_DIR))
  }

  # Create output directory
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
    message(sprintf("[INFO] Created output directory: %s", OUTPUT_DIR))
  }

  message(sprintf("[INFO] Found %d PDF files to split", length(pdf_files)))
  message(sprintf("[INFO] Output directory: %s", OUTPUT_DIR))
  message("")

  total_pages <- 0

  for (pdf_path in pdf_files) {
    author <- tools::file_path_sans_ext(basename(pdf_path))

    # Get page count
    n_pages <- pdf_length(pdf_path)
    message(sprintf("[INFO] %s: %d pages (%.1f MB)",
                    author, n_pages, file.size(pdf_path) / 1e6))

    # Split each page
    for (page in seq_len(n_pages)) {
      out_file <- file.path(OUTPUT_DIR, sprintf("%s_p%02d.pdf", author, page))
      pdf_subset(pdf_path, pages = page, output = out_file)
    }

    total_pages <- total_pages + n_pages
    message(sprintf("  -> Split into %d page files", n_pages))
  }

  message("")
  message(sprintf("[DONE] Split %d PDFs into %d individual page files in '%s'",
                  length(pdf_files), total_pages, OUTPUT_DIR))
}

main()
