# Full Ground Truth Pipeline

## Pipeline Purpose

Performs sex-biased differential expression analysis across 5 datasets (Roussos/PsychAD, Allen/SEA-AD, DeJager/ROSMAP, Johansen, Rexach) and 3 cell types (Astro, Micro, Oligo) using pseudobulk aggregation and DESeq2 with ALL donors per dataset. Then performs cross-dataset comparison diagnostics including UpSet plots, **IoU (Intersection over Union) matrices**, LFC correlation analyses, concordance metrics, forest plots, and gene overlap heatmaps.

**Note:** IoU has been elevated to the primary reproducibility metric (replacing Pearson correlation as the main figure) because it directly measures whether the same genes are called significant across datasets, providing a more interpretable assessment of gene set reproducibility.

This pipeline merges the former "Full Ground Truth" DE analysis and "Cross-Set Analysis" diagnostics into a single unified three-stage validated pipeline.

## Quick Start

```bash
# Run unit tests (Stage 1)
Rscript R/full_ground_truth/tests/test_full_ground_truth.R

# Run toy integration test (Stage 2)
Rscript R/full_ground_truth/toy/run_toy_test.R

# Run full pipeline (Stage 3)
Rscript R/full_ground_truth/run_pipeline.R

# Run specific steps
Rscript R/full_ground_truth/run_pipeline.R 1,2    # Steps 1 and 2 only
Rscript R/full_ground_truth/run_pipeline.R 2       # Step 2 only (diagnostics)
Rscript R/full_ground_truth/run_pipeline.R all 8 8 # All steps, 8 cores each
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_Run_DE_All_Datasets.R` | Pseudobulk aggregation + DESeq2 for all 7 datasets x 3 cell types |
| 2 | `02_Cross_Dataset_Diagnostics.R` | UpSet plots, correlation matrices, concordance, forest plots |
| 3 | `03_Visualise_Summary.R` | Cross-cell-type summary tables and overview figures |

## File Structure

```
R/full_ground_truth/
  config.R                          # All configuration (paths, datasets, thresholds)
  utils.R                           # Shared utility functions
  run_pipeline.R                    # Master runner (executes steps 1-3)
  01_Run_DE_All_Datasets.R          # Step 1: DE analysis
  02_Cross_Dataset_Diagnostics.R    # Step 2: Cross-dataset diagnostics
  03_Visualise_Summary.R            # Step 3: Summary visualisations
  tests/
    test_full_ground_truth.R        # Stage 1: Unit tests (24 tests)
  toy/
    run_toy_test.R                  # Stage 2: Integration test (Micro, 2 datasets, 5 donors)
  README.md                         # This file
```

## Inputs Expected

### Source Data
- **Location:** `Datalink/Shards/{CellType}/{Dataset}/`
- **Format:** Individual donor shard files (`donor_*.qs`) containing SingleCellExperiment objects
- **Structure Required:**
  - `counts` assay: Sparse integer matrix (genes x cells)
  - Metadata columns vary by dataset (see DATASET_CONFIG in config.R)

### Manifests
- **Location:** `Datalink/Shards/{CellType}/{Dataset}/00_MASTER_MANIFEST.qs`

### External Metadata (DeJager only)
- **Location:** `Datalink/raw_data/DeJagerMeta.csv`

### Prerequisites
- Sharding pipeline must have completed successfully
- All donor shards must be valid SingleCellExperiment objects with integer counts

## Outputs Produced

### Per Dataset x Cell Type (Step 1)
- **Location:** `Datalink/Full_Ground_Truth/{CellType}/`
- `{Dataset}_DE_results.csv` — Full DE results table
- `{Dataset}_Full_DESeq2_Results.qs` — Complete DESeqResults object
- `{Dataset}_dds_checkpoint.qs` — DESeqDataSet object
- `{Dataset}_pb_checkpoint.qs` — Pseudobulk SingleCellExperiment

### Per Cell Type Diagnostics (Step 2)
- **Location:** `Datalink/Full_Ground_Truth/{CellType}/Diagnostic_Results/`
- `Diagnostic_UpSet_Plot.png` — Gene overlap UpSet plot
- `UpSet_Male_Biased_Genes.png` — Male-biased gene overlaps
- `UpSet_Female_Biased_Genes.png` — Female-biased gene overlaps
- `Diagnostic_Overlap_Heatmap.png` — Pairwise overlap heatmap
- `Expanded_IoU_Matrix.png` — **IoU matrix (primary reproducibility metric)**
- `Expanded_IoU_Matrix.csv` — Raw IoU values
- `IoU_Summary_Statistics.csv` — Per-pair IoU statistics
- `LFC_Correlation_Top100.png` — LFC correlation for top 100 genes (supplementary)
- `LFC_Correlation_Matrix.csv` — Pairwise LFC correlation matrix
- `LFC_Agreement_Matrix_Roussos_Sig_Genes.png` — Pearson LFC agreement (supplementary)
- `Sign_Agreement_Matrix_Roussos_Sig_Genes.png` — Direction agreement
- `Forest_Plot_Top20_Genes.png` — Cross-dataset forest plots
- `Concordance_Analysis_vs_Roussos.csv` — Concordance metrics
- `XIST_Detection_Summary.csv` — XIST sanity check
- `Significant_Gene_Counts.csv` — Per-dataset sig gene counts

### Cross-Cell-Type Summary (Step 3)
- **Location:** `Datalink/Full_Ground_Truth/`
- `Summary_All_CellTypes.csv` — Master summary table
- `Summary_SigGene_Heatmap.png` — Sig gene count heatmap
- `Summary_XIST_Across_Datasets.png` — XIST LFC bar plot

## Configuration

All configuration is in `config.R`. Key settings:

| Setting | Default | Description |
|---------|---------|-------------|
| `RUN_DATASETS` | All 7 | Which datasets to process |
| `CELL_TYPES` | Astro, Micro, Oligo | Which cell types |
| `N_CORES_AGG` | 12 | Cores for pseudobulk aggregation |
| `N_CORES_DE` | 12 | Cores for DESeq2 |
| `FDR_THRESHOLD` | 0.05 | Significance threshold |
| `LFC_THRESHOLD` | 1 | Effect size threshold (|LFC| > 1.0) |
| `FORCE_RECALC` | FALSE | Skip existing results |
| `REFERENCE_DATASET` | Roussos | Reference for concordance |

## Tests Implemented

### Stage 1 — Unit Tests (`tests/test_full_ground_truth.R`)

| Test Group | Tests | What's Validated |
|------------|-------|-----------------|
| `test_age_parsers` | 8 | parse_verbose_age(), parse_simple_numeric() |
| `test_pseudobulk_aggregation` | 4 | aggregate_shards_parallel() with mock shards |
| `test_deseq2_setup` | 3 | run_deseq_for_dataset() with mock data |
| `test_output_structure` | 5+ | clean_gene_ids(), extract_significant_genes(), build_overlap_matrix(), compute_concordance() |
| `test_no_silent_failures` | 3 | Error propagation (not silent NA returns) |

### Stage 2 — Toy Integration Test (`toy/run_toy_test.R`)

**Configuration:** Micro only, Roussos + Allen, 5 donors each, 4 cores

**Validation Checklist:**
- [ ] DE results CSV files created for both datasets
- [ ] Pseudobulk SCE objects saved
- [ ] DESeq2 dds objects saved
- [ ] UpSet plot generated
- [ ] Correlation plot generated
- [ ] Non-NA padj and log2FoldChange values present
- [ ] XIST detected with expected sex bias
- [ ] Significant gene counts reasonable

### Stage 3 — Full Production Run

**Configuration:** 3 cell types, 7 datasets, full donor pools, 12 cores

## What Is NOT Tested

1. **Shard Integrity** — Assumes sharding pipeline validated counts
2. **Gene Annotation Consistency** — Assumes Ensembl IDs consistent across datasets
3. **Biological Replicability** — Correlation plots provide qualitative check only
4. **Metadata Completeness** — XIST analysis provides indirect sex label validation
5. **Covariate Selection** — Domain knowledge decision, not testable programmatically
6. **Statistical Power** — Separate pipeline (cartesian_downsampling) handles this
7. **Multiple Testing Correction** — Assumes BH FDR appropriate
8. **Cross-Platform Batch Effects** — Correlation analyses provide indirect assessment
9. **Outlier Detection** — Would require PCA-based outlier detection
10. **Reproducibility Across R Versions** — Would require CI/CD infrastructure

## History

This pipeline merges two previously separate analyses:
- **Full Ground Truth DE** (archived: `R/_archived/RoussosHalfnHalf/Run_All_Dataset_Full_Ground_Truths_Updated.R`)
- **Cross-Set Analysis** (archived: `R/_archived/Cross_Set_Analysis/Run_Correlation_Matrix.R`)

Previous output location: `Datalink/Cross-Dataset-Analysis/Full_Set_Outputs/`
Current output location: `Datalink/Full_Ground_Truth/`

---

**Last Updated:** 2026-02-10
**Framework Version:** Three-Stage Pipeline Validation v1.0
**Pipeline Version:** 1.0.0
**Status:** Active
