# Full Ground Truth Pipeline

## Pipeline Purpose

Performs sex-biased differential expression on the full PsychAD cohort (`Roussos` in code and paths) for 3 cell types (Astro, Micro, Oligo), using pseudobulk aggregation and DESeq2 with all filtered donors. It also builds a gene annotation audit (chromosome, biotype, PAR status) of the tested gene universe, which feeds Supplementary Data 1.

`DATASET_CONFIG` in `config.R` now contains only `Roussos`. The pipeline was originally written for five datasets (PsychAD, SEA-AD/Allen, ROSMAP/DeJager, Johansen, Rexach) with cross-dataset IoU, LFC and sign-agreement diagnostics. That cross-dataset analysis has been removed from the paper; the scripts that implement it are kept in the directory but are obsolete (not used for the paper).

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Stage `fgt` runs:

```bash
Rscript R/full_ground_truth/run_pipeline.R 1,4 <n_cores_agg> <n_cores_de>
Rscript R/de_diagnostics/run_pipeline.R 4
```

with `BIOMART_OFFLINE=true`, using the shipped biomaRt caches (`gene_audit_biomart_cache.csv`). Steps 2, 3 and 5 are not run.

## Quick Start

```bash
# Run unit tests (Stage 1)
Rscript R/full_ground_truth/tests/test_full_ground_truth.R

# Run the steps used for the paper (DE and gene audit)
Rscript R/full_ground_truth/run_pipeline.R 1,4

# Run with explicit core counts
Rscript R/full_ground_truth/run_pipeline.R 1,4 12 12
```

`run_pipeline.R` with no arguments runs steps 1-5, including the obsolete cross-dataset steps.

## Pipeline Steps

| Step | Script | Status | Description |
|------|--------|--------|-------------|
| 1 | `01_Run_DE_All_Datasets.R` | Paper | Pseudobulk aggregation + DESeq2 per cell type (PsychAD only) |
| 2 | `02_Cross_Dataset_Diagnostics.R` | Obsolete (not used for the paper) | UpSet plots, IoU, LFC/sign agreement, concordance, forest plots, FDP |
| 3 | `03_Visualise_Summary.R` | Obsolete (not used for the paper) | Cross-dataset summary tables and overview figures |
| 4 | `04_Gene_Audit_Matrix.R` | Paper | Gene universe audit by chromosome, biotype and PAR status |
| 5 | `02b_Regenerate_Heatmaps.R` | Obsolete (not used for the paper) | Re-renders IoU and FDP heatmaps from cached CSVs |

## DE Model

- Pseudobulk counts summed per donor; genes with fewer than `MIN_GENE_COUNT` (10) total counts removed.
- Design: `~ Source + self_reported_ethnicity + disease + pmi_scaled + internal_age_scaled + sex` (sex last). PMI is annotated onto the shards by `R/sharding/Annotate_Roussos_PMI.R`; age is parsed from `development_stage` and standardised.
- Covariates with a single level or zero variance are dropped with a warning.
- LFCs are shrunk with `DESeq2::lfcShrink(type = "ashr")`.

## File Structure

```
R/full_ground_truth/
  config.R                          # Configuration (paths, DATASET_CONFIG, thresholds)
  utils.R                           # Shared utility functions
  run_pipeline.R                    # Master runner (steps 1-5)
  01_Run_DE_All_Datasets.R          # Step 1: DE analysis
  02_Cross_Dataset_Diagnostics.R    # Step 2: cross-dataset diagnostics (obsolete)
  02b_Regenerate_Heatmaps.R         # Step 5: heatmap re-render (obsolete)
  03_Visualise_Summary.R            # Step 3: summary visualisations (obsolete)
  04_Gene_Audit_Matrix.R            # Step 4: gene audit
  _run_step2.R, _run_step2_shared_datalink.R   # Helper wrappers for step 2 (obsolete)
  tests/
    test_full_ground_truth.R        # Stage 1: unit tests
  toy/
    run_toy_test.R                  # Stage 2: integration test (Roussos + Allen; see below)
  README.md                         # This file
```

## Inputs Expected

### Source Data
- **Location:** `Datalink/Shards/{CellType}/Roussos/`
- **Format:** Individual donor shard files (`donor_*.qs`) containing SingleCellExperiment objects with a `counts` assay
- **Manifest:** `Datalink/Shards/{CellType}/Roussos/00_MASTER_MANIFEST.qs`

### Prerequisites
- Sharding pipeline (stage `shard`) completed, including PMI annotation
- For step 4 offline: `Datalink/Full_Ground_Truth/gene_audit_biomart_cache.csv`

## Outputs Produced

### Per Cell Type (Step 1)
- **Location:** `Datalink/Full_Ground_Truth/{CellType}/`
- `Roussos_DE_results.csv`: full DE results table
- `Roussos_Full_DESeq2_Results.qs`: DESeqResults object
- `Roussos_dds_checkpoint.qs`: DESeqDataSet object
- `Roussos_pb_checkpoint.qs`: pseudobulk SingleCellExperiment

### Gene Audit (Step 4)
- `Datalink/Full_Ground_Truth/{CellType}/Gene_Audit_Full_Universe.csv` and `Gene_Audit_Significant.csv`
- `Datalink/Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv` (used for Supplementary Data 1, gene universe composition)
- `Datalink/Full_Ground_Truth/Gene_Audit_Significant_Combined.csv`

### Obsolete outputs (steps 2, 3, 5)
`{CellType}/Diagnostic_Results/` (IoU, LFC/sign agreement, UpSet, forest plots), `FDP_Summary.csv`, `FDP_Heatmap.png` and the `Summary_*` files. These are not referenced by the paper.

## Configuration

All configuration is in `config.R`. Key settings:

| Setting | Default | Description |
|---------|---------|-------------|
| `RUN_DATASETS` | `Roussos` | Datasets to process (all of `DATASET_CONFIG`) |
| `CELL_TYPES` | Astro, Micro, Oligo | Cell types |
| `N_CORES_AGG` | 4 | Cores for pseudobulk aggregation (overridden by argument 2) |
| `N_CORES_DE` | 4 | Cores for DESeq2 (overridden by argument 3) |
| `MIN_GENE_COUNT` | 10 | Minimum total counts per gene |
| `FDR_THRESHOLD` | 0.05 | Significance threshold |
| `LFC_THRESHOLD` | 1 | Legacy \|LFC\| > 1 cut-off, used only by the obsolete diagnostics and `Gene_Audit_Significant*.csv`; the paper truth set uses \|LFC\| > 0.2 (`GT_LFC_THRESHOLD` in `R/power_curve_fitting/config.R`) |
| `FORCE_RECALC` | TRUE | Recompute even if outputs exist (`POWER_ANALYSIS_FORCE_RECALC`) |

## Tests Implemented

### Stage 1: Unit Tests (`tests/test_full_ground_truth.R`)

| Test Group | What's Validated |
|------------|-----------------|
| Age parsers | `parse_verbose_age()`, `parse_simple_numeric()` |
| Pseudobulk aggregation | `aggregate_shards_parallel()` with mock shards, metadata preservation, manifest skipping |
| DESeq2 setup | `run_deseq_for_dataset()` with mock data |
| Output structure | `clean_gene_ids()`, `extract_significant_genes()`, `build_overlap_matrix()`, `compute_concordance()` |
| No silent failures | Error propagation (not silent NA returns) |

### Stage 2: Toy Integration Test (`toy/run_toy_test.R`)

Configured for Micro, Roussos + Allen, 5 donors each. Allen is no longer in `DATASET_CONFIG` and its shards are not produced by the reproduction pipeline, so this test is obsolete in its current form.

### Stage 3: Full Production Run

Stage `fgt` of `pipeline/`: 3 cell types, PsychAD, full filtered donor pool, steps 1 and 4.

## What Is NOT Tested

1. **Shard Integrity**: assumes the sharding pipeline validated counts
2. **Covariate Selection**: domain knowledge decision, not testable programmatically
3. **Statistical Power**: handled by `cartesian_downsampling`
4. **Multiple Testing Correction**: assumes BH FDR is appropriate
5. **Outlier Detection**: would require PCA-based outlier detection
6. **Reproducibility Across R Versions**: the pinned `renv.lock` / Docker image is the control

## History

This pipeline merged two earlier analyses (full ground truth DE and a cross-set correlation analysis; the original scripts are not part of this repository). The cross-dataset component was dropped from the paper in 2026; the code remains for reference.

---

**Last Updated:** 2026-10-06
