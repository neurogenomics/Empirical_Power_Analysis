# EDA Roussos Pipeline

**Purpose:** Exploratory data analysis for the PsychAD dataset (`Roussos` in code and paths) across Astro/Micro/Oligo, on both the unfiltered (`Roussos_Full`) and filtered (`Roussos`) cohorts.

---

## Overview

`run_pipeline.R` has one step:

1. **PCA Diagnostics**: principal component analysis coloured by metadata variables, for the unfiltered and filtered cohorts (Supplementary PCA galleries).

Standalone script:

- `Validate_Sex_Labels.R`: checks donor sex labels against pseudobulk XIST and Y-gene (e.g. DDX3Y, UTY) expression (Supplementary Table 4).

The former steps 3 (read depth) and 6 (DE filtering summary) were not used by the paper and are no longer in the repository. The per-cell-type mean UMIs per cell hardcoded in the `power_curve_fitting` and `depth_equalization` configs (Astro 9,456, Micro 4,426, Oligo 5,742) come from step 3's `Read_Depth_Summary_All_CellTypes_Roussos.csv`.

---

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Stage `eda` runs:

```bash
Rscript R/eda_roussos/run_pipeline.R 1
Rscript R/eda_roussos/Validate_Sex_Labels.R All
Rscript R/figure_data_generation/run_pipeline.R 1
```

---

## Pipeline Structure

```
R/eda_roussos/
├── README.md                                 # This file
├── config.R                                  # Configuration (paths, parameters)
├── utils.R                                   # Shared utility functions
├── 01_PCA_Diagnostics.R                      # Step 1: PCA analysis
├── Validate_Sex_Labels.R                     # Sex label validation
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

---

## Outputs

```
Datalink/Exploratory_Analysis/
├── {CellType}/Roussos_Full/PCA/PCA_FullCohort_*.png       # Step 1 (unfiltered)
└── {CellType}/Roussos/PCA/PCA_FilteredCohort_*.png        # Step 1 (filtered)

Datalink/Sex_Validation/
├── sex_validation_{CellType}.csv      # Per-donor results
├── sex_exclusion_{CellType}.csv       # Flagged donors
└── sex_exclusion_combined.csv         # All cell types merged
```

---

## Running the Pipeline

The supported way to run this is `pipeline/run.sh` (`--only <stage>` reruns one stage). The commands below are what the stage runs, for reading and debugging inside the container.

```bash
# PCA diagnostics (step 1)
Rscript R/eda_roussos/run_pipeline.R 1

# Sex label validation, all cell types
Rscript R/eda_roussos/Validate_Sex_Labels.R All
```

---

## Configuration

Edit `config.R` for cell types, cohorts (`SHARD_DATASETS`), PCA gene counts and cores. `config.R` still sets `STEPS_TO_RUN <- c(1, 3, 6)` and the histogram and DE-filtering settings of the removed steps; with no arguments, `run_pipeline.R` runs step 1 only.

---

## Dependencies

`SingleCellExperiment`, `DESeq2`, `ggplot2`, `dplyr`, `tidyr`, `Matrix`, `qs`, `viridis`, `BiocParallel`, `parallel`.
