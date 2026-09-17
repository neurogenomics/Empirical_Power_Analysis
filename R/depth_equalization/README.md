# Depth Equalization Pipeline

## Purpose

Sensitivity analysis for sequencing depth confounding in the power curve fitting results. Astrocytes have ~2.1× the mean reads per cell of microglia (9,456 vs 4,426 UMIs). In cells mode of the cartesian downsampling pipeline, full cells are taken with all reads intact, meaning astrocytes get ~2× the total sequencing information per donor at the same cell count grid point. This pipeline multinomially downsamples every astrocyte cell proportionally to match microglia's mean depth, producing depth-equalized shards.

## Method

For each cell in each Astro/Roussos donor shard:

```
new_total = round(original_total × 0.4679)
new_counts = Multinomial(new_total, prob = original_counts / original_total)
```

The fraction 0.4679 = 4,426 / 9,456 (Micro mean / Astro mean).

**Key properties:**
- Gene-level proportions preserved within each cell (multinomial sampling)
- All cells downsampled uniformly (proportional, not threshold-based)
- Cell count unchanged; only per-cell depth reduced
- Distribution shape preserved, shifted down uniformly
- Deterministic per-shard seeds for reproducibility

## Input

- `Datalink/Shards/Astro/Roussos/` — 1,293 donor shards + manifest

## Output

- `Datalink/Shards/Astro/Roussos_DepthEq/` — 1,293 depth-equalized shards + manifest + summary CSV

The output directory is a drop-in replacement for the source: same file naming, same SCE structure, same manifest format. Downstream pipelines (cartesian downsampling, power curve fitting) can use it by changing `SHARD_DATASET <- "Roussos_DepthEq"`.

## Usage

```bash
Rscript R/depth_equalization/run_pipeline.R
```

## Tests

```bash
Rscript R/depth_equalization/tests/test_depth_equalization.R
```

## Post-Downsampling Statistics (Expected)

| Metric | Source (Roussos) | Target (Roussos_DepthEq) |
|--------|-----------------|--------------------------|
| Mean reads/cell | 9,456 | ~4,426 |
| Median reads/cell | 8,197 | ~3,836 |
| Min reads/cell | 1,240 | ~580 |
| Max reads/cell | 60,867 | ~28,484 |
| Total cells | 675,930 | 675,930 (unchanged) |
| N donors | 1,293 | 1,293 (unchanged) |

## Downstream Workflow

1. **Run this pipeline** → produces `Roussos_DepthEq` shards
2. **Run cartesian downsampling** with `SHARD_DATASET = "Roussos_DepthEq"`, reusing original GT splits
3. **Run power curve fitting** on depth-equalized results
4. **Compare coefficients** between original and depth-equalized models

## Version

- v1.0.0 (2026-02-25): Initial implementation
