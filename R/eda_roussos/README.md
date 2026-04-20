# EDA Roussos Pipeline

**Purpose:** Exploratory data analysis for Roussos dataset across multiple cell types.

**Version:** 1.0.0  
**Last Updated:** 2026-02-10  
**Status:** Active

---

## Overview

This pipeline performs comprehensive exploratory data analysis on the Roussos single-cell RNA-seq dataset, generating:

1. **PCA Diagnostics** - Principal component analysis colored by metadata variables
2. **Cell Count Analysis** - Distribution of cells per donor
3. **Read Depth Analysis** - Distribution of reads (UMIs) per cell
4. **Gene Expression Summaries** - Mean expression, detection rate, zero-inflated genes
5. **QC Metrics** - Mitochondrial %, library size, genes detected per cell
6. **DE Filtering Summary** - DESeq2 gene filtering breakdown by category and LFC magnitude

---

## Pipeline Structure

```
R/eda_roussos/
├── README.md                          # This file
├── config.R                           # Configuration (paths, parameters)
├── utils.R                            # Shared utility functions
├── 01_PCA_Diagnostics.R              # Step 1: PCA analysis
├── 02_Cell_Count_Analysis.R          # Step 2: Cells per donor
├── 03_Read_Depth_Analysis.R          # Step 3: Reads per cell
├── 04_Gene_Expression_Summary.R      # Step 4: Gene-level statistics
├── 05_QC_Metrics.R                   # Step 5: Cell-level QC
├── 06_DE_Filtering_Summary.R         # Step 6: DE filtering breakdown
├── run_pipeline.R                    # Master runner
├── tests/                            # Unit tests (Stage 1)
│   ├── test_age_parsing.R
│   ├── test_pca_data_prep.R
│   ├── test_gene_stats.R
│   └── test_qc_metrics.R
└── toy/                              # Toy integration test (Stage 2)
    └── run_toy_test.R
```

---

## Input Data

**Location:** `Datalink/Shards/{CellType}/Roussos/`

**Format:** Donor-level shards (`.qs` files) containing SingleCellExperiment objects

**Cell Types:** Astro, Micro, Oligo

**Required Metadata Columns:**
- `donor_id` (required)
- `sex` (required)
- `development_stage` (optional, for age parsing)
- `Source` (optional, for batch effects)
- `disease` (optional)
- `self_reported_ethnicity` (optional)
- `subtype` (optional)

---

## Output Structure

```
Datalink/Exploratory_Analysis/
├── Astro/
│   └── Roussos/
│       ├── PCA/
│       │   ├── PCA_FullCohort_sex.png
│       │   ├── PCA_FullCohort_age.png
│       │   ├── PCA_FullCohort_source.png
│       │   ├── PCA_FullCohort_disease.png
│       │   ├── PCA_FullCohort_ethnicity.png
│       │   ├── PCA_FullCohort_subtype.png
│       │   └── PCA_FullCohort_cell_count.png
│       ├── Histograms/
│       │   ├── Histogram_Cells_Per_Donor.png
│       │   └── Histogram_Reads_Per_Cell.png
│       ├── QC/
│       │   ├── QC_Mitochondrial_Percent.png
│       │   ├── QC_Genes_Detected.png
│       │   └── QC_Library_Size.png
│       └── Tables/
│           ├── Gene_Expression_Summary.csv
│           ├── Top_Expressed_Genes.csv
│           ├── QC_Summary_Stats.csv
│           └── DE_Filtering_Summary.csv
├── Micro/
│   └── Roussos/
│       └── [same structure]
├── Oligo/
│   └── Roussos/
│       └── [same structure]
├── DE_Filtering_Summary_Combined.csv      # All cell types combined
├── DE_Filtering_Summary_BarChart.png      # Stacked bar (LFC bins, faceted)
└── DE_Filtering_Summary_Totals.png        # Totals bar (all genes)
```

---

## Running the Pipeline

### Stage 1: Unit Tests

```bash
# Test individual functions
Rscript R/eda_roussos/tests/test_age_parsing.R
Rscript R/eda_roussos/tests/test_pca_data_prep.R
Rscript R/eda_roussos/tests/test_gene_stats.R
Rscript R/eda_roussos/tests/test_qc_metrics.R
```

### Stage 2: Toy Integration Test

```bash
# Run on small subset (Micro, 10 donors)
Rscript R/eda_roussos/toy/run_toy_test.R
```

### Stage 3: Production Run

```bash
# Run full pipeline on all cell types
Rscript R/eda_roussos/run_pipeline.R

# Or run individual steps
Rscript R/eda_roussos/01_PCA_Diagnostics.R
Rscript R/eda_roussos/02_Cell_Count_Analysis.R
# ... etc
```

### RStudio Background Job

```r
# Source as background job (works due to config-first design)
source("R/eda_roussos/run_pipeline.R")
```

---

## Configuration

Edit `config.R` to customize:

- **Cell types:** `CELL_TYPES` (default: Astro, Micro, Oligo)
- **PCA genes:** `N_PCA_GENES` (default: 10000)
- **Histogram bins:** `HISTOGRAM_BINS` (default: 100)
- **Parallelization:** `N_CORES` (default: 12)
- **Zero-inflated threshold:** `ZERO_INFLATED_THRESHOLD` (default: 0.90)
- **Top genes:** `TOP_GENES_N` (default: 50)

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| **N_PCA_GENES** | 10000 | Top variable genes for PCA |
| **MIN_GENE_COUNTS** | 10 | Minimum rowSums for gene filtering |
| **HISTOGRAM_BINS** | 100 | Number of bins for histograms |
| **ZERO_INFLATED_THRESHOLD** | 0.90 | Proportion of zeros to flag gene |
| **TOP_GENES_N** | 50 | Number of top expressed genes to report |
| **MT_PATTERN** | "^MT-" | Regex for mitochondrial genes |
| **N_CORES** | 12 | Cores for parallel processing |
| **DE_DATASET** | "Roussos" | Dataset for DE filtering summary |
| **DE_FILTER_LFC_BINS** | c(1,2,3,Inf) | \|LFC\| bin boundaries |

---

## Dependencies

**R Packages:**
- SingleCellExperiment
- DESeq2
- ggplot2
- dplyr
- tidyr
- Matrix
- qs
- viridis
- BiocParallel
- parallel

**Install:**
```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "Matrix", "qs", "viridis", "parallel"))

# Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "DESeq2", "BiocParallel"))
```

---

## Outputs for Paper

The pipeline automatically generates figure blocks in `docs/paper_skeleton.Rmd` under:

**Section:** Supplementary Materials → Dataset Characterization

**Figure blocks per cell type:**
- PCA plots (7 figures)
- Histograms (2 figures)
- QC plots (3 figures)

---

## Troubleshooting

### No mitochondrial genes found
- Check `MT_PATTERN` in config.R
- Verify gene names in your data (e.g., "MT-" vs "mt-" vs "Mt-")

### Age parsing fails
- Check `AGE_COL` column format
- Expected format: "80 year-old stage", "6 month-old stage", etc.
- Use `looks_like_age_stage()` to validate

### PCA fails with "too few samples"
- Ensure at least 3 donors per cell type
- Check manifest loading with `load_donor_ids_from_manifest()`

### Memory issues
- Reduce `N_CORES` in config.R
- Process cell types sequentially instead of in parallel
- Reduce `N_PCA_GENES` to 5000

---

## Version History

- **1.4.0** (2026-02-20) - DE Filtering Summary
  - Step 6: DESeq2 gene filtering breakdown by category and LFC magnitude
  - Classifies genes into: Low-count pre-filter, Cook's distance, Independent filtering, Fully tested
  - Stacked bar chart faceted by cell type + totals bar chart
  - Combined and per-cell-type summary CSVs

- **1.3.0** (2026-02-19) - Filtered cohort support
  - Steps 1-3 now process both Roussos_Full and Roussos (filtered) datasets
  - Skip-if-exists logic, cohort labels, dataset-suffixed CSVs

- **1.0.0** (2026-02-10) - Initial release
  - PCA diagnostics with 7 metadata variables
  - Cell count and read depth histograms
  - Gene expression summaries
  - QC metrics (MT%, library size, genes detected)
  - Automatic paper integration

---

## Contact

