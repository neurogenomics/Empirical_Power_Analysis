# Output

This document describes the outputs produced by the Empirical Power Analysis
pipeline.  All outputs are written under `params.outdir` (default: `results/`).
When running against the published paper data, `outdir` is set to `Datalink/`
so that the in-repository data tree is updated in place.

---

## Pipeline information

`results/pipeline_info/`

Nextflow automatically generates execution metadata for every run:

| File | Description |
|---|---|
| `execution_timeline.html` | Interactive Gantt chart of process runtimes |
| `execution_report.html` | Summary report (CPU/memory/I-O per process) |
| `execution_trace.txt` | TSV with per-task resource usage |
| `dag.svg` | Pipeline DAG |

---

## Shards

`Datalink/Shards/{CellType}/{Dataset}/`

Donor-level `SingleCellExperiment` objects serialised as `.qs` files.  One file
per donor, plus a `00_MASTER_MANIFEST.qs` index.  Inputs to all downstream
steps.

---

## Exploratory analysis

`Datalink/Exploratory_Analysis/`

PCA gallery and filtering diagnostics for the Roussos (PsychAD) dataset.
Informational only; no downstream dependencies.

---

## Full ground truth

`Datalink/Full_Ground_Truth/`

| File / directory | Description |
|---|---|
| `{CellType}/Roussos/` | Per-dataset DESeq2 DE results (`.qs`) |
| `FDP_Summary.csv` | False-discovery-proportion summary across 5 datasets |
| `FDP_Heatmap.png` | Cross-dataset FDP heatmap |
| `Gene_Audit_Full_Universe_Combined.csv` | Full gene universe with biomaRt annotation |

---

## Cartesian downsampling

`Datalink/Cartesian_Downsampling/{CellType}/{Dataset}/{mode}/Run_{1-6}/`

Empirical power surfaces: DESeq2 results for every grid point (donors ×
cells-per-donor × sequencing depth × robustness run × subsample).  ~22 million
per-gene binary detection outcomes.

---

## Depth equalisation

`Datalink/Depth_Equalization/`

Sensitivity analysis: astrocyte Cartesian grid at microglia-equivalent sequencing
depth.  Produces one supplementary figure in the paper.

---

## Power curve fitting

`Datalink/Power_Curve_Fitting/`

| Path | Description |
|---|---|
| `Exports/power_equation_coefficients.csv` | Logistic regression coefficients |
| `Exports/power_lookup_table.csv` | Pre-computed 11 250-row prediction grid |
| `Exports/predict_power.R` | Standalone callable R function |
| `Models/` | Serialised GLM objects (`.qs`) |
| `Plots/` | Power curves, depth decomposition, cross-validation |

---

## DE diagnostics

`Datalink/DE_Diagnostics/`

Density plots, chromosome-stratified recall surfaces, split-recovery curves and
shared-gene diagnostics.  Informational; no downstream dependencies.

---

## Bulk DEG processing

`Datalink/Bulk_Analysis/`

GTEx v8 DLPFC sex-biased DE concordance with the single-cell ground truth.
Produces the bulk-validation figure in the paper.

---

## Paper source data

`Datalink/Paper_Source_Data/`

Canonical CSVs imported directly by `paper/paper.Rmd`:

| File | Description |
|---|---|
| `dataset_statistics.csv` | Cohort size / cell counts per dataset × cell type |
| `fdp_summary.csv` | Cross-dataset FDP summary |
| `model_coefficients.csv` | Power equation coefficients (paper Table 2) |
| `model_summary.csv` | Model performance metrics |
| `truth_gene_characteristics.csv` | Ground-truth gene-level statistics |

---

## Paper

`Datalink/Paper/` (if `params.render_paper = true`)

| File | Description |
|---|---|
| `paper.html` | Rendered manuscript (HTML) |
| `paper.pdf` | Rendered manuscript (PDF, if TeX is available) |
