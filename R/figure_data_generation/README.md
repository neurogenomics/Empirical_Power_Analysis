# Figure Data Generation Pipeline

Generates the CSV source tables read by `nature_comms_draft.rmd`,
`nature_comms_supplementary.rmd` and the numeric-reference layer
`_paper_values.R` (`val()`), so that numbers quoted in the paper are computed
from pipeline outputs rather than typed by hand.

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Only
steps 1, 6 and 8 are run:

| Stage | Command | Used for |
|-------|---------|----------|
| `eda` | `Rscript R/figure_data_generation/run_pipeline.R 1` | Cohort statistics (donor counts quoted in the text, Supplementary Data 1) |
| `tables` | `Rscript R/figure_data_generation/run_pipeline.R 6,8` | Supplementary Table 2 and truth-set sizes; tidy recall/precision grid for numbers quoted in the Results |

## Steps

| Step | Script | Status | Output |
|------|--------|--------|--------|
| 1 | `01_Dataset_Statistics.R` | Paper | `Paper_Source_Data/dataset_statistics.csv` (PsychAD x 3 cell types) |
| 2 | `02_Truth_Gene_Characteristics.R` | Obsolete (not used for the paper) | `Paper_Source_Data/truth_gene_characteristics.csv` (uses the legacy \|LFC\| > 1 cut-off) |
| 3 | `03_Collect_Model_Data.R` | Obsolete (not used for the paper) | `model_coefficients.csv`, `model_summary.csv` from the removed logistic-regression power model |
| 4 | `04_Collect_FDP_Data.R` | Obsolete (not used for the paper) | `fdp_summary.csv` from the removed cross-dataset FDP analysis |
| 5 | `05_Validate_Exports.R` | Obsolete (not used for the paper) | Checks the step 1-4 files; expects the old five-dataset layout |
| 6 | `06_Autosomal_vs_SexChr_Summary.R` | Paper | `Datalink/_scratch_autosomal_vs_sexchr_summary.csv` and `_pivot.csv` |
| 7 | `07_Depth_Eq_Truth_Summary.R` | Not run by the pipeline | `depth_eq_truth_summary.csv`, `depth_eq_truth_summary_sd.csv` |
| 8 | `08_Recall_Precision_Grid.R` | Paper | `Paper_Source_Data/recall_precision_grid.csv` |

Notes:

- Step 1 loads every filtered donor shard (`Datalink/Shards/{CellType}/Roussos/`) and reports N donors, N cells, and per-donor and per-cell summaries (cells, UMIs, genes). It sources the shard helpers in `R/glial_characteristics/utils.R`. `DATASETS` in `config.R` contains only `Roussos` (label `PsychAD`).
- Step 6 reads the truth-half DESeq2 checkpoints (`Power_Curve_Fitting/GT_DDS_Checkpoints/`) and uses the paper truth set, padj < 0.05 and \|LFC\| >= 0.2.
- Step 8 flattens the cartesian `Aggregated_Summary_Data.qs` checkpoints into one row per (cell type, mode, method, donors, cells/UMI level, \|LFC\| bin, gene category).

## Usage

```bash
# Steps used for the paper
Rscript R/figure_data_generation/run_pipeline.R 1
Rscript R/figure_data_generation/run_pipeline.R 6,8

# With more cores
Rscript R/figure_data_generation/run_pipeline.R 1 24
```

`run_pipeline.R` with no arguments runs all steps, including the obsolete ones,
which fail if their legacy inputs are absent.

## Step 1 runtime

Step 1 is the heavyweight step: it reads every filtered donor shard for three
cell types. The other steps read existing pipeline outputs and are fast.

## Why This Exists

The paper previously hardcoded dataset statistics and other tables in the
RMarkdown, which let quoted values drift from the pipeline outputs. Generating
them here means a change upstream propagates to the paper on re-knitting.
