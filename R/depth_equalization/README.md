# Depth Equalization Pipeline

## Purpose

Tests how much of the difference in DEG yield between microglia and the other two glial cell types is due to sequencing depth and cells per donor. Astrocytes have about 2.1x, and oligodendrocytes about 1.3x, the mean UMIs per cell of microglia (9,456 and 5,742 vs 4,426). For astrocytes and for oligodendrocytes, this pipeline builds shard variants matched to microglia on depth, cell count, or both, and counts how many ground-truth DEGs each variant loses relative to the native cell type. The result is Supplementary Table 3, which has a row block per cell type.

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Stage `deg_loss` runs, with `N_ROBUSTNESS_RUNS=10`, once for each of `DEPTH_EQ_CELL_TYPE=Astro` and `DEPTH_EQ_CELL_TYPE=Oligo`:

```bash
Rscript R/depth_equalization/run_pipeline.R
Rscript R/depth_equalization/deg_loss_analysis.R <n_cores>
```

then deletes that cell type's shard variants and `GT_Cache/`, and finally stacks the two cell types:

```bash
Rscript R/depth_equalization/combine_deg_loss_tables.R
```

## Method

### Depth equalisation (`run_pipeline.R`, `Roussos_DepthEq`)

For each cell in each `{CellType}/Roussos` donor shard (`CellType` = `DEPTH_EQ_CELL_TYPE`):

```
new_total = round(original_total × f)
new_counts = Multinomial(new_total, prob = original_counts / original_total)
```

The fraction `f` is the Micro mean over the source cell type's mean (`TARGET_MEAN_RPC` / `SOURCE_MEAN_RPC`): 4,426 / 9,456 = 0.4679 for Astro, 4,426 / 5,742 = 0.7708 for Oligo.

- Gene-level proportions preserved within each cell (multinomial sampling)
- All cells downsampled uniformly; cell count unchanged
- Deterministic per-shard seeds (`SEED_BASE` + donor index)

### Depth + cell-count matching (`run_pipeline.R`, `Roussos_ShallowFewer`)

Starting from `Roussos_DepthEq`, cells per donor are scaled by the ratio of the Micro mean cell count to the source cell type's (`FEWER_EXTENSION`, `scaling_mode = "mean_ratio"`); donors left with fewer than 30 cells are dropped.

### Cell-count matching only (`deg_loss_analysis.R`, `Roussos_CellsOnly`)

Built inside `deg_loss_analysis.R` if missing: the same cell removal as `ShallowFewer`, applied to full-depth shards.

### DEG loss (`deg_loss_analysis.R`)

For each of the four configurations (native, depth only, cells only, depth and cells) and each of 10 runs, the same 50% ground-truth donor split as `R/cartesian_downsampling/01_Generate_Ground_Truth.R` is aggregated and analysed with DESeq2, edgeR and limma-voom. DEGs are counted at padj < 0.05 with |LFC| > 0.2 and |LFC| > 1; the 1 threshold is reported alongside 0.2 because limma-voom is unstable at the permissive cut-off.

## Input

- `Datalink/Shards/{Astro,Oligo}/Roussos/` (filtered shards + manifest)
- `Datalink/Shards/Micro/Roussos/00_MASTER_MANIFEST.qs` (target cell-count distribution)

## Output

- `Datalink/Shards/{CellType}/Roussos_DepthEq/`: depth-equalised shards, manifest, `depth_equalization_summary.csv`
- `Datalink/Shards/{CellType}/Roussos_ShallowFewer/`: depth- and cell-matched shards
- `Datalink/Shards/{CellType}/Roussos_CellsOnly/`: cell-matched shards
- `Datalink/Depth_Eq_DEG_Loss/{CellType}/`: `deg_loss_table.csv`, `deg_loss_counts.csv`, `deg_loss_per_run.csv` and `GT_Cache/` (cached per-run DE tables) for one cell type
- `Datalink/Depth_Eq_DEG_Loss/deg_loss_table.csv`: Supplementary Table 3 (mean DEG count and % of native, per cell type, method and threshold), from `combine_deg_loss_tables.R`
- `Datalink/Depth_Eq_DEG_Loss/deg_loss_counts.csv`, `deg_loss_per_run.csv`: numeric backing, both cell types

The `deg_loss` stage keeps only the CSVs; the shard variants and `GT_Cache/` are deleted after each cell type.

The shard variants use the same file naming, SCE structure and manifest format as the source.

## Usage

The supported way to run this is `pipeline/run.sh` (`--only <stage>` reruns one stage). The commands below are what the stage runs, for reading and debugging inside the container.

```bash
DEPTH_EQ_CELL_TYPE=Oligo Rscript R/depth_equalization/run_pipeline.R
DEPTH_EQ_CELL_TYPE=Oligo Rscript R/depth_equalization/deg_loss_analysis.R 24
Rscript R/depth_equalization/combine_deg_loss_tables.R
```

## Tests

```bash
Rscript R/depth_equalization/tests/test_depth_equalization.R
Rscript R/depth_equalization/tests/test_astro_shallow_fewer.R
```

## Post-Downsampling Statistics (Expected, Astro `Roussos_DepthEq`)

| Metric | Source (Roussos) | Target (Roussos_DepthEq) |
|--------|-----------------|--------------------------|
| Mean UMIs/cell | 9,456 | ~4,426 |
| Median UMIs/cell | 8,197 | ~3,836 |
| Min UMIs/cell | 1,240 | ~580 |
| Max UMIs/cell | 60,867 | ~28,484 |
| Total cells | 675,930 | 675,930 (unchanged) |
| N donors | 1,293 | 1,293 (unchanged) |

## History

The earlier workflow ran the cartesian downsampling and power curve fitting on the depth-equalised shards and compared fitted models (summarised by a `figure_data_generation` step 7 that has since been removed). That route is not used for the paper; `deg_loss_analysis.R` replaces it.
