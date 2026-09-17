# DE Diagnostics Pipeline

**Location:** `R/de_diagnostics/`

## Purpose

Generates diagnostic visualisations and statistical tests for the DE analysis and cross-cell-type gene characterisation:

4. **DEG Chromosome Grid** — Counts of significant DEGs stratified by chromosome type (autosomal/sex) and LFC magnitude bin.

5. **Shared Gene Diagnostics** — Paired cross-cell-type comparisons of the 23 shared truth genes (from the power curve fitting pipeline). Computes paired Wilcoxon tests, Pearson correlations, and summary statistics for 5 metrics.

## Pipeline Steps

### Step 4: DEG Chromosome Grid (`04_DEG_Chromosome_Grid.R`)

**Inputs:**
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_DE_results.csv` (3 files)
- `Datalink/DE_Diagnostics/gene_biotype_cache.csv` (with chromosome_name)

**Outputs:**
- `Datalink/DE_Diagnostics/Chromosome_Grid/DEG_Chr_Grid_{CellType}.png` (3 per-cell-type grids)
- `Datalink/DE_Diagnostics/Chromosome_Grid/DEG_Chr_Grid_Combined.png` (1 faceted grid)
- `Datalink/DE_Diagnostics/Chromosome_Grid/DEG_Chr_Grid_Summary.csv`

**Total plots:** 3 per-cell-type + 1 combined = **4 PNG files** + 1 CSV

**Filters:** padj < 0.05 AND |LFC| > 1

**Grid layout:**
|  | Small (1-2) | Medium (2-3) | Large (3+) |
|--|-------------|--------------|------------|
| **SexChr** | count | count | count |
| **Autosomal** | count | count | count |

### Step 5: Shared Gene Diagnostics (`05_Shared_Gene_Diagnostics.R`)

**Inputs:**
- `Datalink/Power_Curve_Fitting/Fitting_Data/shared_truth_genes.qs` (23 gene IDs)
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_DE_results.csv` (3 files)
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_dds_checkpoint.qs` (3 files, for dispersion)

**Outputs:**
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Shared_Gene_Paired_Tests.csv` (all test statistics)
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Shared_Gene_Comparison_Matrix.csv` (long format)
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Shared_Gene_Comparison_Wide.csv` (wide format)
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Shared_Gene_Summary.csv` (shared gene medians)
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Full_Truth_Set_Summary.csv` (full truth set medians)
- `Datalink/DE_Diagnostics/Shared_Gene_Diagnostics/Scatter_Plots/Scatter_{metric}_{CT_A}_vs_{CT_B}.png`

**Total outputs:** 5 CSVs + 15 scatter plots (5 metrics × 3 CT pairs)

**Metrics tested (paired Wilcoxon signed-rank):**
| Metric | Description | Pearson correlation? |
|--------|-------------|---------------------|
| `dispersion` | Gene-level dispersion from DDS checkpoint | Yes |
| `abs_lfc` | Absolute log2 fold-change | No |
| `lfcSE` | Standard error of LFC estimate | No |
| `snr` | Signal-to-noise ratio (\|LFC\| / lfcSE) | No |
| `baseMean` | Mean normalised counts | Yes |

**Cell type pairs:** Astro–Micro, Astro–Oligo, Micro–Oligo

## Configuration

All settings in `config.R`. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DE_DATASET` | `"Roussos"` | Dataset for DE analyses |
| `PADJ_CUTOFF` | `0.05` | FDR threshold |
| `RECOVERY_LFC_THRESHOLD` | `1.0` | Minimum |LFC| for recovery comparison |
| `FORCE_RECALC` | `FALSE` | Skip existing outputs |

## Dependencies

- `ggplot2`, `dplyr`, `tidyr` — plotting and data wrangling
- `DESeq2`, `SummarizedExperiment` — dispersion
- `qs` — data loading
- `biomaRt` — chromosome annotation (with CSV caching)

## File Structure

```
R/de_diagnostics/
├── config.R                       # Shared configuration
├── utils.R                        # Shared utility functions
├── run_pipeline.R                 # Master runner (steps 4-5)
├── 04_DEG_Chromosome_Grid.R       # Step 4: DEG chromosome grid
├── 05_Shared_Gene_Diagnostics.R   # Step 5: shared gene cross-CT comparisons
├── tests/
│   └── test_de_diagnostics.R      # Unit tests
├── toy/
│   └── run_toy_test.R             # Toy integration test
└── README.md                      # This file
```
