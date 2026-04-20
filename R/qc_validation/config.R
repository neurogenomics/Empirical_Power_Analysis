# ==============================================================================
# SCRIPT: QC Validation Pipeline - Configuration
# PURPOSE: Shared configuration for validating Roussos shard QC status
# ==============================================================================

# ==============================================================================
# PROJECT ROOT (HARDCODED FOR RELIABILITY)
# ==============================================================================
PROJECT_ROOT <- if (nzchar(Sys.getenv("POWER_ANALYSIS_ROOT"))) {
  Sys.getenv("POWER_ANALYSIS_ROOT")
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("The `here` package is required. Install with: install.packages(\"here\")")
  }
  here::here()
}
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "qc_validation")

project_path <- function(...) {
  file.path(PROJECT_ROOT, ...)
}

# ==============================================================================
# INPUT PATHS
# ==============================================================================

# Shard directories (Roussos_Full = unfiltered, all donors)
SHARD_BASE <- project_path("Datalink", "Shards")
CELL_TYPES <- c("Astro", "Micro", "Oligo")
SHARD_DATASET <- "Roussos_Full"

# PsychAD h5ad source files
H5AD_DIR <- project_path("Datalink", "raw_data")
H5AD_FILES <- c(
  HBCC = file.path(H5AD_DIR, "PsychAD-HBCC.h5ad"),
  MSSM = file.path(H5AD_DIR, "PsychAD-MSSM.h5ad"),
  RADC = file.path(H5AD_DIR, "PsychAD-RADC.h5ad")
)

# ==============================================================================
# OUTPUT PATHS
# ==============================================================================

OUTPUT_BASE <- project_path("Datalink", "QC_Validation")

# ==============================================================================
# QC THRESHOLDS (from Fullard et al. 2025, Nature Scientific Data)
# ==============================================================================

# Cell-level thresholds (stringent QC, post-demultiplexing)
QC_MIN_UMIS   <- 1500
QC_MAX_UMIS   <- 110000
QC_MIN_GENES  <- 1100
QC_MAX_GENES  <- 12500
QC_MAX_MITO_PCT <- 5.0

# Donor-level threshold
QC_MIN_CELLS_PER_DONOR <- 50

# ==============================================================================
# FAIL-EARLY THRESHOLDS
# ==============================================================================
# If a high proportion of cells violate thresholds, stop early rather than
# processing all ~1,490 donors per cell type.

# After processing 10 donors: if >50% of cells violate ANY threshold, STOP
FAIL_EARLY_10_DONORS_PCT <- 0.50

# After processing 50 donors: if >20% of cells violate ANY threshold, STOP
FAIL_EARLY_50_DONORS_PCT <- 0.20

# ==============================================================================
# BIOMART CONFIGURATION (for mitochondrial gene identification)
# ==============================================================================

# Ensembl dataset for human genes
BIOMART_DATASET <- "hsapiens_gene_ensembl"

# Cache file for mito gene list (avoid repeated biomaRt queries)
MITO_GENE_CACHE <- file.path(OUTPUT_BASE, "mito_genes_ensembl.qs")

# ==============================================================================
# PLOT CONFIGURATION
# ==============================================================================

PLOT_WIDTH  <- 10
PLOT_HEIGHT <- 7
PLOT_DPI    <- 300

# Threshold line colours for plots
THRESHOLD_COLOUR <- "red"
THRESHOLD_LINETYPE <- "dashed"

# ==============================================================================
# PARALLELISATION
# ==============================================================================

N_CORES <- 24
