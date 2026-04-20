# EDA Roussos Pipeline

**Purpose:** Exploratory data analysis for the Roussos dataset across Astro/Micro/Oligo.

---

## Overview

Retained steps used by the paper:

1. **PCA Diagnostics** — Principal component analysis coloured by metadata variables (feeds the Figure 1 PCA supplement).
3. **Read Depth Analysis** — Distribution of reads (UMIs) per cell; the resulting summary (`Read_Depth_Summary_All_CellTypes_Roussos.csv`) is referenced by the `power_curve_fitting` and `depth_equalization` configs as the documented source of hardcoded UMI-per-cell values.
6. **DE Filtering Summary** — DESeq2 gene filtering breakdown by category and LFC magnitude (source of the 37% / 4% / 27% figures in the paper).

Step numbering is preserved from the original pipeline; gaps are intentional.

---

## Pipeline Structure

```
R/eda_roussos/
├── README.md                       # This file
├── config.R                        # Configuration (paths, parameters)
├── utils.R                         # Shared utility functions
├── 01_PCA_Diagnostics.R            # Step 1: PCA analysis
├── 03_Read_Depth_Analysis.R        # Step 3: Reads per cell
├── 06_DE_Filtering_Summary.R       # Step 6: DE filtering breakdown
├── run_pipeline.R                  # Master runner
├── tests/                          # Unit tests
└── toy/                            # Toy integration test
```

---

## Input Data

**Location:** `Datalink/Shards/{CellType}/Roussos/`
**Format:** Donor-level shards (`.qs` files) containing SingleCellExperiment objects
**Cell Types:** Astro, Micro, Oligo

**Required metadata columns:** `donor_id`, `sex`; optional: `development_stage`, `Source`, `disease`, `self_reported_ethnicity`, `subtype`.

---

## Outputs

```
Datalink/Exploratory_Analysis/
├── {CellType}/Roussos/
│   ├── PCA/PCA_FullCohort_*.png          # Step 1
│   └── Histograms/Histogram_Reads_Per_Cell.png   # Step 3
├── Read_Depth_Summary_All_CellTypes_Roussos.csv  # Step 3 (consumed by power_curve_fitting/depth_equalization configs)
├── DE_Filtering_Summary_Combined.csv     # Step 6
├── DE_Filtering_Summary_BarChart.png     # Step 6
└── DE_Filtering_Summary_Totals.png       # Step 6
```

---

## Running the Pipeline

```bash
# Run all retained steps (1, 3, 6)
Rscript R/eda_roussos/run_pipeline.R

# Run a specific step
Rscript R/eda_roussos/run_pipeline.R 1
Rscript R/eda_roussos/run_pipeline.R 3
Rscript R/eda_roussos/run_pipeline.R 6
```

---

## Configuration

Edit `config.R` for cell types, PCA gene counts, histogram bins, cores, and DE filtering bins. Default `STEPS_TO_RUN` is `c(1, 3, 6)`.

---

## Dependencies

`SingleCellExperiment`, `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `Matrix`, `qs`, `viridis`, `BiocParallel`, `parallel`.
