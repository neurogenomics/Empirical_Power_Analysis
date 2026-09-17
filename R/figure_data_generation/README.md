# Figure Data Generation Pipeline

Generates all CSV source data tables needed by `docs/paper_final_draft.rmd`. Replaces hardcoded `data.frame()` calls in the paper with `read.csv()` from canonical, pipeline-generated CSV files.

## Output

All files written to `Datalink/Paper_Source_Data/`:

| File | Rows | Source | Replaces in Paper |
|------|------|--------|-------------------|
| `dataset_statistics.csv` | 15 (5 datasets × 3 CTs) | Computed from shard files | `source_data_1` (hardcoded, had placeholder values) |
| `truth_gene_characteristics.csv` | 3 | Computed from Roussos DE results | `truth_chars_median` and `truth_chars_mean` |
| `model_coefficients.csv` | 6 | Copy from PCF exports | `model_coeffs` |
| `model_summary.csv` | 6 | Computed from PCF coefficients | `model_summary` |
| `fdp_summary.csv` | 39-40 | Copy from Full GT diagnostics | `fdp_summary` |

## Usage

```bash
# Run all steps
Rscript run_pipeline.R

# Run only dataset statistics (Step 1 — requires shard access, slow)
Rscript run_pipeline.R 1

# Run steps 2-5 (fast, reads existing CSVs)
Rscript run_pipeline.R 2,3,4,5

# Run with more cores
Rscript run_pipeline.R all 24
```

## Step 1: Dataset Statistics

This is the heavyweight step — it loads every donor shard for all 5 datasets × 3 cell types to compute real statistics. Takes ~10-30 minutes depending on I/O and cores.

Steps 2-5 are fast (seconds) as they read from existing pipeline outputs.

## Why This Exists

The paper previously hardcoded dataset statistics, model coefficients, and other tables directly in the RMarkdown. This caused:
- Placeholder values (e.g., `98765`, `52341`) for datasets never processed
- Inconsistencies between hardcoded values and actual pipeline outputs
- No reproducibility — values couldn't be regenerated from data

This pipeline makes the paper fully reproducible: change a parameter upstream, rerun this pipeline, and the paper picks up the new values automatically.
