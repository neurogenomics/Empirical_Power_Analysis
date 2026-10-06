# EDA Roussos Pipeline

**Purpose:** Exploratory data analysis for the PsychAD dataset (`Roussos` in code and paths) across Astro/Micro/Oligo, on both the unfiltered (`Roussos_Full`) and filtered (`Roussos`) cohorts.

---

## Overview

Steps run by `run_pipeline.R` (default `STEPS_TO_RUN <- c(1, 3, 6)` in `config.R`):

1. **PCA Diagnostics**: principal component analysis coloured by metadata variables, for the unfiltered and filtered cohorts (Supplementary PCA galleries).
3. **Read Depth Analysis**: distribution of UMIs per cell. The summary (`Read_Depth_Summary_All_CellTypes_Roussos.csv`) is cited in the `power_curve_fitting` and `depth_equalization` configs as the source of their hardcoded UMI-per-cell values.
6. **DE Filtering Summary**: DESeq2 gene filtering breakdown by category and |LFC| magnitude (bins at 1, 2, 3), from the full ground truth results.

Step numbering is preserved from the original pipeline; gaps are intentional.

Standalone scripts:

- `Validate_Sex_Labels.R`: checks donor sex labels against pseudobulk XIST and Y-gene (e.g. DDX3Y, UTY) expression (Supplementary Table 4).
- `Analyse_Sex_Validation_CrossCellType.R`: follow-up on flagged donors across cell types. Not run by the reproduction pipeline.
- `Run_Genes_Per_Cell_Histogram_Micro.R`: genes-per-cell histogram for microglia. Not run by the reproduction pipeline.

---

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Stage `eda` runs:

```bash
Rscript R/eda_roussos/run_pipeline.R 1
Rscript R/eda_roussos/Validate_Sex_Labels.R All
Rscript R/figure_data_generation/run_pipeline.R 1
```

Steps 3 and 6 are not run by the reproduction pipeline; their outputs are not read by the manuscript.

---

## Pipeline Structure

```
R/eda_roussos/
├── README.md                                 # This file
├── config.R                                  # Configuration (paths, parameters)
├── utils.R                                   # Shared utility functions
├── 01_PCA_Diagnostics.R                      # Step 1: PCA analysis
├── 03_Read_Depth_Analysis.R                  # Step 3: UMIs per cell
├── 06_DE_Filtering_Summary.R                 # Step 6: DE filtering breakdown
├── Validate_Sex_Labels.R                     # Sex label validation
├── Analyse_Sex_Validation_CrossCellType.R    # Cross-cell-type follow-up
├── Run_Genes_Per_Cell_Histogram_Micro.R      # Microglia genes-per-cell histogram
├── run_pipeline.R                            # Master runner
├── tests/                                    # Unit tests
└── toy/                                      # Toy integration test
```

---

## Input Data

**Location:** `Datalink/Shards/{CellType}/Roussos_Full/` and `Datalink/Shards/{CellType}/Roussos/`
**Format:** Donor-level shards (`.qs` files) containing SingleCellExperiment objects
**Cell Types:** Astro, Micro, Oligo

**Required metadata columns:** `donor_id`, `sex`; optional: `development_stage`, `Source`, `disease`, `self_reported_ethnicity`, `subtype`.

Step 6 additionally reads `Datalink/Full_Ground_Truth/{CellType}/` DE results and pseudobulk checkpoints.

---

## Outputs

```
Datalink/Exploratory_Analysis/
├── {CellType}/Roussos_Full/PCA/PCA_FullCohort_*.png       # Step 1 (unfiltered)
├── {CellType}/Roussos/PCA/PCA_FilteredCohort_*.png        # Step 1 (filtered)
├── {CellType}/{Dataset}/Histograms/Histogram_Reads_Per_Cell.png   # Step 3
├── {CellType}/{Dataset}/Read_Depth_Summary.csv             # Step 3
├── Read_Depth_Summary_All_CellTypes_{Roussos,Roussos_Full}.csv   # Step 3
├── DE_Filtering_Summary_Combined.csv                       # Step 6
├── DE_Filtering_Summary_BarChart.png                       # Step 6
└── DE_Filtering_Summary_Totals.png                         # Step 6

Datalink/Sex_Validation/
├── sex_validation_{CellType}.csv      # Per-donor results
├── sex_exclusion_{CellType}.csv       # Flagged donors
└── sex_exclusion_combined.csv         # All cell types merged
```

---

## Running the Pipeline

```bash
# Run the default steps (1, 3, 6)
Rscript R/eda_roussos/run_pipeline.R

# Run a specific step
Rscript R/eda_roussos/run_pipeline.R 1

# Sex label validation, all cell types
Rscript R/eda_roussos/Validate_Sex_Labels.R All
```

---

## Configuration

Edit `config.R` for cell types, cohorts (`SHARD_DATASETS`), PCA gene counts, histogram bins, cores, and DE filtering bins.

---

## Dependencies

`SingleCellExperiment`, `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `Matrix`, `qs`, `viridis`, `BiocParallel`, `parallel`.
