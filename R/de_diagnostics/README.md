# DE Diagnostics Pipeline

**Version:** 1.5.0  
**Created:** 2026-02-19  
**Location:** `R/de_diagnostics/`

## Purpose

Generates diagnostic visualisations and statistical tests for the DE analysis, cartesian downsampling validation, and cross-cell-type gene characterisation:

1. **DE Density Plots** — Density plots of differential expression metric distributions (|LFC|, SE, SNR, baseMean, dispersion) overlaid across all 3 cell types, with optional coding/non-coding biotype stratification.

2. **GT-vs-Validation PCA** — PCA plots of the cartesian downsampling 50:50 donor splits, coloured by split membership and key covariates (sex, age, source), to verify no systematic bias was introduced by the random split.

3. **Split DEG Recovery** — Compares the 50% split ground truth DE results against the full-cohort ground truth to assess how representative each split is. Computes bidirectional overlap (recall + precision), stratified by LFC magnitude and chromosome type.

4. **DEG Chromosome Grid** — Counts of significant DEGs stratified by chromosome type (autosomal/sex) and LFC magnitude bin.

5. **Shared Gene Diagnostics** — Paired cross-cell-type comparisons of the 23 shared truth genes (from the power curve fitting pipeline). Computes paired Wilcoxon tests, Pearson correlations, and summary statistics for 5 metrics. Provides the reproducible evidence chain for the Subtype Heterogeneity Hypothesis.

6. **Subtype Expression Divergence** — For each of the 23 shared truth genes, generates a per-gene figure showing all 3 cell types with split-violin geometry: Micro (pooled baseline), Oligo (OPALIN|RBFOX1 split), and Astro (2 split violins for 4 subtypes). Tests whether subtypes have divergent expression that would introduce pseudobulk variance at low cell counts.

## Pipeline Steps

### Step 1: DE Density Plots (`01_DE_Density_Plots.R`)

**Inputs:**
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_DE_results.csv` (3 files)
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_dds_checkpoint.qs` (3 files, for dispersion)
- BioMart (Ensembl, with CSV cache)

**Outputs:**
- `Datalink/DE_Diagnostics/Density_Plots/Cross_CellType/Density_{Metric}.png` (5 aggregated)
- `Datalink/DE_Diagnostics/Density_Plots/Cross_CellType/Density_{Metric}_Stratified.png` (5 stratified)
- `Datalink/DE_Diagnostics/Density_Plots/Cross_CellType/Density_Summary.csv`
- `Datalink/DE_Diagnostics/gene_biotype_cache.csv`

**Metrics plotted:**
| Metric | Column | Description |
|--------|--------|-------------|
| LFC | `Abs_LFC` | Absolute log2 fold-change |
| SE | `lfcSE` | Standard error of LFC |
| LogSNR | `Log10_SNR` | Log10 signal-to-noise ratio (|LFC|/SE) |
| BaseMean | `Log10_BaseMean` | Log10 mean normalised counts |
| Dispersion | `Log10_Dispersion` | Log10 gene-level dispersion |

**Filters:** padj < 0.05 AND |LFC| > 2

### Step 2: Split PCA Validation (`02_Split_PCA_Validation.R`)

**Inputs:**
- `Datalink/Shards/{CellType}/Roussos/` (donor shard files)
- `Datalink/Cartesian_Downsampling/{CellType}/Roussos/{mode}/Run_{N}/Ground_Truth_sampling_metadata.qs`

**Outputs:**
- `Datalink/DE_Diagnostics/PCA/{mode}/{CellType}/Run_{N}_PCA_{Covariate}.png`

**Total plots:** 3 cell types × 2 modes × 6 runs × 4 covariates = **144 PNG files**

**PCA covariates:**
| Covariate | Colour | Shape | Density Contours |
|-----------|--------|-------|-----------------|
| Split (GT/Val) | GT=orange, Val=blue | — | Yes |
| Sex | Categorical | GT/Val | No |
| Age | Viridis (magma) | GT/Val | No |
| Source | Categorical | GT/Val | No |

### Step 3: Split DEG Recovery (`03_Split_DEG_Recovery.R`)

**Inputs:**
- `Datalink/Full_Ground_Truth/{CellType}/Roussos_DE_results.csv` (3 files)
- `Datalink/Cartesian_Downsampling/{CellType}/Roussos/{mode}/Run_{N}/Ground_Truth.qs` (36 files)
- `Datalink/DE_Diagnostics/gene_biotype_cache.csv` (extended with chromosome_name)

**Outputs:**
- `Datalink/DE_Diagnostics/Split_Recovery/{CellType}/Recovery_Bars_{mode}.png` (6 bar charts)
- `Datalink/DE_Diagnostics/Split_Recovery/{CellType}/Recovery_Bars_{mode}_Stratified.png` (6 stratified bar charts)
- `Datalink/DE_Diagnostics/Split_Recovery/{CellType}/UpSet_{mode}.png` (6 UpSet plots)
- `Datalink/DE_Diagnostics/Split_Recovery/{CellType}/LFC_Scatter_{mode}.png` (6 scatter plots)
- `Datalink/DE_Diagnostics/Split_Recovery/Recovery_Summary.csv`
- `Datalink/DE_Diagnostics/Split_Recovery/Recovery_Summary_Stratified.csv`

**Total plots:** 3 cell types × 2 modes × 4 plot types = **24 PNG files** + 2 CSVs

**Filters:** padj < 0.05 AND |LFC| > 1

**Stratification:**
| Stratum | Description |
|---------|-------------|
| Overall | All significant DEGs |
| LFC Small (1-2) | Genes with 1 ≤ \|LFC\| < 2 |
| LFC Medium (2-3) | Genes with 2 ≤ \|LFC\| < 3 |
| LFC Large (3+) | Genes with \|LFC\| ≥ 3 |
| Chr Autosomal | Chromosomes 1-22 |
| Chr SexChr | Chromosomes X, Y |

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
| `dispersion` | Gene-level dispersion from DDS checkpoint | ✅ Yes |
| `abs_lfc` | Absolute log2 fold-change | No |
| `lfcSE` | Standard error of LFC estimate | No |
| `snr` | Signal-to-noise ratio (\|LFC\| / lfcSE) | No |
| `baseMean` | Mean normalised counts | ✅ Yes |

**Cell type pairs:** Astro–Micro, Astro–Oligo, Micro–Oligo

**Key question answered:** Do the 23 shared truth genes have similar intrinsic characteristics across cell types? If so, power differences must come from a cell-type-level mechanism (subtype heterogeneity), not gene-level properties.

### Step 6: Subtype Expression Divergence (`06_Subtype_Expression_Divergence.R`)

**Inputs:**
- `Datalink/Power_Curve_Fitting/Fitting_Data/shared_truth_genes.qs` (23 gene IDs)
- `Datalink/Shards/{Oligo,Astro,Micro}/Roussos/donor_*.qs` (500 randomly sampled per CT)
- BioMart (Ensembl, for gene symbol lookup with CSV cache)

**Outputs:**
- `Datalink/DE_Diagnostics/Subtype_Expression_Divergence/{Symbol}_{EnsemblID}.png` (23 per-gene figures)
- `Datalink/DE_Diagnostics/Subtype_Expression_Divergence/{Symbol}_{EnsemblID}.pdf` (23 per-gene figures)
- `Datalink/DE_Diagnostics/Subtype_Expression_Divergence/Subtype_Expression_Summary.csv`
- `Datalink/DE_Diagnostics/gene_symbol_cache.csv` (BioMart cache)

**Total outputs:** 23 PNG + 23 PDF (one per shared truth gene) + 1 CSV

**Method:**
1. Load 23 shared truth genes (Ensembl IDs, version-stripped).
2. For each cell type (Oligo, Astro, Micro), randomly sample 500 Roussos donor shards (seed=42+CT_offset). Parallel extraction with 40 cores.
3. For each shard, extract per-cell raw UMI counts for the 23 genes, split by subtype:
   - **Oligo:** 2 subtypes (OPALIN, RBFOX1) — Wilcoxon rank-sum test
   - **Astro:** 4 subtypes (ADAMTSL3, GRIA1, PLSCR1, WIF1) — Kruskal-Wallis omnibus test
   - **Micro:** no subtypes (single pooled violin as baseline)
4. Pool all cells across donors. Compute per-gene × per-CT × per-subtype summary statistics (n_cells, mean, median, SD, %zero, quartiles) and appropriate statistical tests.
5. Look up gene symbols via BioMart (cached) for human-readable filenames.
6. For each gene, generate a single figure with all 3 cell types using split-violin geometry:
   - **X-axis positions (left to right):** Micro (pooled) | Oligo (OPALIN|RBFOX1 split) | Astro₁ (ADAMTSL3|GRIA1 split) | Astro₂ (PLSCR1|WIF1 split)
   - **Micro:** Standard violin (no subtypes, baseline control)
   - **Oligo:** Split violin — left half OPALIN, right half RBFOX1
   - **Astro₁:** Split violin — left half ADAMTSL3, right half GRIA1
   - **Astro₂:** Split violin — left half PLSCR1, right half WIF1
   - Zero counts excluded from violins; dropout probability annotated per subtype
   - Mean expression annotated per subtype
   - Donor-level mean points overlaid (jittered, semi-transparent)
   - Black crossbar at median
   - Statistical test p-value annotation (Oligo: Wilcoxon, Astro: Kruskal-Wallis omnibus)
   - Genes ordered by max mean expression range across subtypes

**Figure layout (per gene):**
```
  Micro        Oligo           Astro₁            Astro₂
  (All)    (OPALIN|RBFOX1)  (ADAMTSL3|GRIA1)  (PLSCR1|WIF1)
   ┃          ┃┃              ┃┃                ┃┃
   ┃          ┃┃              ┃┃                ┃┃
  violin    split-violin    split-violin      split-violin
```

**Key question answered:** Do cell subtypes have divergent expression for the shared truth genes? If so, random subtype composition at low cell counts introduces pseudobulk variance, explaining power curve instability. Microglia (no subtypes) serves as a baseline control.

## Configuration

All settings in `config.R`. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DE_DATASET` | `"Roussos"` | Dataset for DE density plots |
| `PADJ_CUTOFF` | `0.05` | FDR threshold |
| `LFC_MIN_HIST` | `2.0` | Minimum |LFC| for density plots |
| `N_PCA_GENES` | `10000` | Top variable genes for PCA |
| `N_ROBUSTNESS_RUNS` | `6` | Number of cartesian runs |
| `CARTESIAN_MODES` | `c("cells", "reads")` | Downsampling modes |
| `RECOVERY_FDR_THRESHOLD` | `0.05` | FDR threshold for recovery comparison |
| `RECOVERY_LFC_THRESHOLD` | `1.0` | Minimum |LFC| for recovery comparison |
| `FORCE_RECALC` | `FALSE` | Skip existing outputs |
| `SUBTYPE_DIV_N_DONORS` | `500` | Random donors to sample per CT for Step 6 |
| `SUBTYPE_DIV_SEED` | `42` | Seed for donor sampling |
| `SUBTYPE_DIV_N_CORES` | `40` | Parallel cores for shard loading |

## Dependencies

- `ggplot2`, `dplyr`, `tidyr` — plotting and data wrangling
- `DESeq2`, `SingleCellExperiment`, `SummarizedExperiment` — VST, PCA, dispersion
- `qs`, `Matrix` — data loading
- `biomaRt` — gene biotype and chromosome annotation (with CSV caching)
- `BiocParallel`, `parallel` — parallel shard loading
- `UpSetR` — UpSet intersection plots (Step 3)
- `RColorBrewer`, `viridis`, `scales` — colour palettes
- Custom `geom_split_violin` — split-violin geometry for Step 6 (defined in `utils.R`)

## Estimated Runtime

- **Step 1** (density plots): ~2–5 minutes
- **Step 2** (PCA): ~30–60 minutes (36 pseudobulk aggregations + VST)
- **Step 3** (split recovery): ~5–15 minutes (loading 36 Ground_Truth.qs files + BioMart)
- **Step 4** (chromosome grid): ~1–3 minutes (loads DE CSVs + BioMart cache)
- **Step 5** (shared gene diagnostics): ~2–5 minutes (loads DE CSVs + DDS checkpoints)
- **Step 6** (subtype expression divergence): ~10–20 minutes (loads 500 donor shards × 3 CTs in parallel, 40 cores; generates 23 per-gene figures)

## File Structure

```
R/de_diagnostics/
├── config.R                    # Shared configuration
├── utils.R                     # Shared utility functions
├── run_pipeline.R              # Master runner (steps 1-6)
├── 01_DE_Density_Plots.R       # Step 1: density plots
├── 02_Split_PCA_Validation.R   # Step 2: PCA diagnostics
├── 03_Split_DEG_Recovery.R     # Step 3: split DEG recovery
├── 04_DEG_Chromosome_Grid.R    # Step 4: DEG chromosome grid
├── 05_Shared_Gene_Diagnostics.R # Step 5: shared gene cross-CT comparisons
├── 06_Subtype_Expression_Divergence.R # Step 6: per-gene split-violin subtype divergence
├── tests/
│   └── test_de_diagnostics.R   # Unit tests
├── toy/
│   └── run_toy_test.R          # Toy integration test
└── README.md                   # This file
```

## Archived Predecessors

This pipeline replaces functionality from:
- `R/_archived/RoussosHalfnHalf/Run_CellType_Predictors.R` (density plots)
- `R/_archived/RoussosHalfnHalf/Run_GTvsVal_PCA_Validation.R` (PCA split validation)
- `R/_archived/RoussosHalfnHalf/Run_PCA_All.R` (PCA with density contours)

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2026-02-19 | Initial pipeline creation (Steps 1-2) |
| 1.1.0 | 2026-02-19 | Added Step 3: Split DEG Recovery |
| 1.2.0 | 2026-02-19 | Added Step 4: DEG Chromosome Grid |
| 1.3.0 | 2026-02-22 | Added Step 5: Shared Gene Diagnostics (paired cross-CT comparisons) |
| 1.4.0 | 2026-02-23 | Added Step 6: Subtype Expression Divergence (oligo subtype UMI violin plots) |
| 1.4.1 | 2026-02-23 | Expanded Step 6 to all 3 cell types (Oligo/Astro/Micro), 500 donors, zeros excluded, combined mega-figure |
| 1.5.0 | 2026-02-23 | Redesigned Step 6: per-gene split-violin figures (Micro pooled + Oligo split + 2× Astro split). Fixed annotation NA bug from `unlist()` name mangling. Added `geom_split_violin` custom geom and `get_gene_symbols()` to utils.R. |
