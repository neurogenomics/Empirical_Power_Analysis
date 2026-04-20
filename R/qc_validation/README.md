# QC Validation Pipeline

## Pipeline Purpose

Validates whether the Roussos (PsychAD) donor shard files contain raw (pre-QC) or post-QC data by computing per-cell QC metrics and comparing against the thresholds published in Fullard et al. 2025 (Nature Scientific Data). Also inspects the original PsychAD h5ad files to determine which matrix layer (adata.X vs adata.raw.X) was extracted during sharding.

## Inputs Expected

| Input | Path | Format |
|-------|------|--------|
| Donor shards | `Datalink/Shards/{CellType}/Roussos_Full/donor_*.qs` | SingleCellExperiment (.qs) |
| PsychAD h5ad files | `Datalink/raw_data/PsychAD-{HBCC,MSSM,RADC}.h5ad` | AnnData h5ad |

**Prerequisites:**
- Sharding pipeline must have been run (`R/sharding/Run_Shard_Roussos_Full.R`)
- PsychAD h5ad files must be accessible for Step 00

## Outputs Produced

| Output | Path | Description |
|--------|------|-------------|
| H5AD layer report | `Datalink/QC_Validation/h5ad_layer_report.txt` | Definitive answer on adata.X contents |
| Per-donor summary | `Datalink/QC_Validation/{CellType}/qc_summary_per_donor.csv` | QC stats per donor |
| Per-cell metrics | `Datalink/QC_Validation/{CellType}/qc_all_cells.qs` | Full per-cell QC data |
| UMI histogram | `Datalink/QC_Validation/{CellType}/Histogram_UMIs.png` | With threshold lines |
| Gene histogram | `Datalink/QC_Validation/{CellType}/Histogram_Genes.png` | With threshold lines |
| Mito histogram | `Datalink/QC_Validation/{CellType}/Histogram_Mito.png` | With 5% threshold |
| Scatter plot | `Datalink/QC_Validation/{CellType}/Scatter_Genes_vs_UMIs.png` | Genes vs UMIs, coloured by mito% |
| Verdict | `Datalink/QC_Validation/{CellType}/QC_Verdict.txt` | Pass/fail assessment |
| Global summary | `Datalink/QC_Validation/Global_QC_Summary.csv` | Cross-cell-type comparison |

## Tests Implemented

### Stage 1: Unit Tests (`tests/test_qc_utils.R`)
- `identify_mito_genes()` — Ensembl list matching, ^MT- fallback, no-match handling
- `compute_cell_qc_metrics()` — correct structure, value ranges, no-mito handling
- `check_thresholds()` — all 5 threshold types correctly flagged
- `summarise_violations()` — correct counts and percentages
- `summarise_donor()` — single-row summary with correct values
- `check_fail_early()` — triggers at 10/50 donor checkpoints, does not trigger with low violations
- `extract_donor_id()` — filename parsing
- `sanity_check_counts()` — passes for valid counts, fails for normalized data

### Stage 2: Toy Integration (`toy/run_toy_test.R`)
- 10 donors per cell type, full pipeline end-to-end
- Sanity checks on first donor per cell type
- Mito gene detection verification
- Histogram generation for visual inspection
- Per-donor summary CSV output

### Stage 3: Production (`run_pipeline.R`)
- All ~1,490 donors × 3 cell types
- Fail-early logic (stops if >50% violations after 10 donors, >20% after 50)
- Full report with plots and verdicts

## What Is NOT Tested

- **CellBender correction verification**: We can determine if counts are integer (raw-like) vs non-integer (normalized), but cannot distinguish pre-CellBender raw counts from post-CellBender corrected counts (both are integer-valued). The h5ad layer inspection helps but is not definitive.
- **Scrublet doublet detection**: We do not re-run doublet detection. We only check if the cell-level QC thresholds have been applied.
- **Gene-level QC**: We do not verify the 0.05% expression threshold from the paper.
- **Donor-level QC**: We check cell counts per donor but do not verify the 50-cell minimum was applied (this is a sharding-level filter, not a cell-level one).
- **Batch correction**: We do not verify Harmony was applied (this is an embedding-level operation, not visible in count matrices).

## Running

```bash
# Stage 1: Unit tests
Rscript R/qc_validation/tests/test_qc_utils.R

# Stage 2: Toy integration test
Rscript R/qc_validation/toy/run_toy_test.R

# Stage 3: Full production run
Rscript R/qc_validation/run_pipeline.R
```

## Reference

Fullard et al. 2025. "Population-scale cross-disorder atlas of the human prefrontal cortex at single-cell resolution." Nature Scientific Data 12:954.
https://www.nature.com/articles/s41597-025-04687-5
