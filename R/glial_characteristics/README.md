# Glial Characteristics Pipeline

**Purpose:** Compare cell-level and donor-level characteristics between Roussos and DeJager datasets (both DLPFC) across the three glial cell types to determine whether inter-cell-type ratios are biological invariants or dataset-specific artifacts.

## Scientific Question

Are patterns like:
- "Astrocytes always have ~2× the UMIs per cell of microglia"
- "Oligodendrocytes always have ~7× more cells per donor than microglia"

...biological invariants (driven by cell biology and brain composition), or dataset-specific artifacts (driven by sequencing protocol, QC thresholds, or annotation differences)?

## Pipeline Structure

```
R/glial_characteristics/
├── config.R                        # Shared configuration
├── utils.R                         # Shared utility functions  
├── run_pipeline.R                  # Master runner (3 steps)
├── 01_Extract_Metrics.R            # Step 1: Compute per-dataset metrics
├── 02_Cross_Dataset_Comparison.R   # Step 2: Comparison plots + ratio table
├── 03_Bootstrap_Ratio_Tests.R      # Step 3: Bootstrap CIs on ratios
├── tests/test_all_utils.R          # Unit tests
├── toy/run_toy_test.R              # Toy integration test
└── README.md                       # This file
```

## Usage

```bash
# Run full pipeline
Rscript run_pipeline.R

# Run specific steps
Rscript run_pipeline.R 1         # Only Step 1 (Extract Metrics)
Rscript run_pipeline.R 2 3       # Steps 2 and 3 only

# Run tests
Rscript tests/test_all_utils.R

# Run toy test with synthetic data
Rscript toy/run_toy_test.R
```

## Datasets Compared

| Dataset | Brain Region | Approx. Donors | Source |
|---------|-------------|-----------------|--------|
| **Roussos** (filtered) | DLPFC | ~1,200 | PsychAD |
| **DeJager** | DLPFC | ~450 | De Jager et al. |

Both datasets are filtered to age ≥40 and cells ≥30 per donor.

## Metrics Computed

1. **UMIs per cell** — Library size (sequencing depth)
2. **Genes detected per cell** — Transcriptomic complexity
3. **Cells per donor** — Cell abundance/capture efficiency

## Outputs

```
Datalink/Glial_Characteristics/
├── Metrics/
│   ├── Summary_Stats_{Dataset}_{CellType}.csv
│   ├── Donor_Summaries_{Dataset}_{CellType}.csv
│   └── Cell_Metrics_{Dataset}_{CellType}.qs
├── Plots/
│   ├── Histogram_UMIs_Per_Cell_Comparison.png
│   ├── Histogram_Genes_Per_Cell_Comparison.png
│   ├── Histogram_Cells_Per_Donor_Comparison.png
│   ├── Violin_UMIs_Per_Cell_Comparison.png
│   ├── Violin_Genes_Per_Cell_Comparison.png
│   ├── Violin_Cells_Per_Donor_Comparison.png
│   ├── Ratio_Comparison_BarChart.png
│   └── Combined_Violin_Summary.png
├── Tests/
│   ├── Bootstrap_Ratio_CIs.csv
│   ├── Ratio_Difference_Tests.csv
│   └── Forest_Plot_Ratio_CIs.png
└── Summary/
    ├── All_Summary_Stats.csv
    ├── All_Donor_Summaries.csv
    ├── Inter_CellType_Ratios.csv
    ├── Ratio_Comparison.csv
    ├── Median_Comparison_Wide.csv
    └── Ratio_Interpretation_Summary.csv
```

## Interpretation

For each inter-cell-type ratio (e.g., Astro:Micro UMIs/cell):
- **PRESERVED**: 95% bootstrap CIs overlap between datasets → biological invariant
- **DIVERGENT**: CIs do not overlap → dataset-specific pattern

## Configuration

Key settings in `config.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_BOOTSTRAP` | 10,000 | Bootstrap resamples for CIs |
| `CI_LEVEL` | 0.95 | Confidence level |
| `MAX_CELLS_PER_CT` | 100,000 | Max cells for plotting (subsampled) |
| `N_CORES` | 12 | Parallel cores for shard loading |
| `REFERENCE_CELL_TYPE` | "Micro" | Denominator for ratios |

## Dependencies

- SingleCellExperiment
- qs
- ggplot2
- dplyr
- tidyr
- scales
- patchwork
- BiocParallel
- Matrix

## Version

- **v1.0.0** (2026-02-27) — Initial release

## Related Pipelines

- `R/eda_roussos/` — Single-dataset EDA (Roussos only)
- `R/full_ground_truth/` — Cross-dataset DE comparison
