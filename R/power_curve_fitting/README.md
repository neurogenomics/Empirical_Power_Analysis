# Power Curve Fitting Pipeline (v4.0.0)

## 1. Pipeline Purpose

Fits continuous power equations to the empirical cartesian downsampling results using logistic regression. For each gene in the ground truth (|LFC| > 1, padj < 0.05), computes per-iteration binary detection outcomes, then models detection probability as a function of experimental design parameters and gene-level characteristics.

### v4.0.0: Logistic Regression (Replacing Beta Regression)

v4.0.0 switches from beta regression (`betareg`) to logistic regression (`glm(family=binomial)`) because the response is binary per-iteration detection (0/1), not a continuous proportion. Logistic regression uses the same logit link, produces identical coefficient interpretation, and handles binary data natively without the Smithson-Verkuilen squeeze transformation. This is a simplification: removes the `betareg` dependency, squeeze functions, `BETA_SQUEEZE` config flag, and the precision parameter (phi).

### v3.1.0: Expanded Truth Set & Simplified Aggregation

v3.1.0 lowers the truth set threshold from |LFC| > 2 to |LFC| > 1, capturing smaller effect sizes. It also removes the "aggregated" aggregation level (mean detection across iterations), keeping only per-iteration binary detection data. This simplifies the pipeline from 16 to 8 primary models.

### v3.0.0: Unified Model with Actual Reads-Per-Cell

v3.0.0 merges cells-mode and reads-mode cartesian downsampling data into unified models by parsing all ~83K sampling metadata files to extract **actual** n_cells and mean reads_per_cell per iteration. This breaks the collinearity that made v2.0.0's decomposed formula uninformative (where reads_per_cell was a per-cell-type constant).

**Unified formula:**
```
detection ~ log_donors + log_n_cells + log_reads_per_cell +
            log_abs_gt_lfc + log_gt_basemean + log_gt_dispersion
```

**Key design decisions:**
- Combined model has **NO cell_type term** — tests whether 6 predictors are sufficient
- LRT compares: combined (no CT) vs +CT (additive) vs CT×predictors (interaction)
- 8 models total: (3 per-CT + 1 combined) × 2 gene sets × 1 aggregation level

## 2. Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_Generate_GT_DDS_Checkpoints.R` | Re-run DESeq2 on GT halves for dispersion extraction |
| 0.5 | `00b_Cache_Shard_Read_Counts.R` | Pre-compute per-cell read counts for all donor shards |
| 1 | `01_Build_Fitting_Dataset.R` | Parse metadata, merge cells+reads modes, build unified datasets |
| 2 | `02_Fit_Models.R` | Fit 8 model variants, LRT tests, coefficient comparison |
| 3 | `03_Visualise.R` | Power curves, depth decomposition, cell type comparisons |
| 4 | `04_Cross_Validate_Modes.R` | Validate unified model on mode subsets, combined vs per-CT |
| 5 | `05_Export_Equations.R` | Export coefficients, lookup tables, standalone predict_power.R |

## 3. Inputs Expected

| Input | Path | Format | Prerequisites |
|-------|------|--------|---------------|
| Cartesian downsampling DE results | `Datalink/Cartesian_Downsampling/{CT}/Roussos/{mode}/Run_{1-6}/N*_C*_I*.qs` | DESeqResults (.qs) | Cartesian downsampling complete |
| Sampling metadata | `Datalink/Cartesian_Downsampling/{CT}/Roussos/{mode}/Run_{1-6}/*_sampling_metadata.qs` | Named list (.qs) | Cartesian downsampling complete |
| Ground truth results | `Datalink/Cartesian_Downsampling/{CT}/Roussos/{mode}/Run_{1-6}/Ground_Truth.qs` | DESeqResults (.qs) | Cartesian downsampling complete |
| Donor shards | `Datalink/Shards/{CT}/Roussos/donor_*.qs` | SingleCellExperiment (.qs) | Sharding pipeline complete |
| BioMart cache | `Datalink/Full_Ground_Truth/gene_audit_biomart_cache.csv` | CSV | Full ground truth pipeline complete |

## 4. Outputs Produced

### Step 0: GT DDS Checkpoints
| Output | Path |
|--------|------|
| GT DDS checkpoints | `GT_DDS_Checkpoints/{CT}/GT_DDS_Run_{1-6}.qs` |

### Step 0.5: Shard Cache
| Output | Path |
|--------|------|
| Per-cell read counts | `Shard_Cache/{CT}/donor_{ID}_reads.qs` |

### Step 1: Fitting Datasets
| Output | Path |
|--------|------|
| Per-iteration fitting data (full) | `Fitting_Data/{CT}/unified_fitting_data_per_iteration.qs` |
| Per-iteration fitting data (shared) | `Fitting_Data/{CT}/unified_fitting_data_per_iteration_shared.qs` |
| Shared gene list | `Fitting_Data/shared_truth_genes.qs` |

### Step 2: Models
| Output | Path |
|--------|------|
| Per-CT models | `Models/{gene_set}_per_iteration/{CT}_model.qs` |
| Combined model (no CT) | `Models/{gene_set}_per_iteration/combined_model.qs` |
| Combined + CT (additive) | `Models/{gene_set}_per_iteration/combined_additive_model.qs` |
| Combined + CT×pred | `Models/{gene_set}_per_iteration/combined_interaction_model.qs` |
| Model checkpoint | `Models/{gene_set}_per_iteration/all_models_checkpoint.qs` |
| Coefficients | `Models/{gene_set}_per_iteration/all_coefficients.csv` |
| Model comparison | `Models/{gene_set}_per_iteration/model_comparison.csv` |
| **v3 comparison** | `Models/model_comparison_v3.csv` |

### Steps 3-5: Plots, Cross-Validation, Export
| Output | Path |
|--------|------|
| Per-variant plots | `Plots/{gene_set}_per_iteration/` |
| Comparison plots | `Plots/Comparison/` |
| Cross-validation | `Plots/Cross_Validation/` |
| Equation coefficients | `Exports/power_equation_coefficients.csv` |
| Lookup table | `Exports/power_lookup_table.csv` |
| Prediction function | `Exports/predict_power.R` |

## 5. Model Variants

### Gene sets
- **full**: All truth genes per cell type (different genes per CT)
- **shared**: Intersection of truth genes across all CTs

### Aggregation level
- **per_iteration**: No aggregation → 1 row per gene × grid point × run × iteration

### Models per variant
- 3 per-CT models (Astro, Micro, Oligo) — identical formula, different data
- 1 combined model (no cell_type term) — tests generalisability
- 1 combined + cell_type (additive) — for LRT comparison
- 1 combined + cell_type × predictors (interaction) — for LRT comparison

### LRT tests
- **LRT1**: combined (no CT) vs combined + CT → does cell type add explanatory power?
- **LRT2**: combined + CT vs combined + CT × predictors → do cell types differ in slopes?

## 6. Usage

```bash
# Run all steps
Rscript R/power_curve_fitting/run_pipeline.R all 48

# Run specific steps
Rscript R/power_curve_fitting/run_pipeline.R 0.5,1,2 48

# Step IDs: 0, 0.5, 1, 2, 3, 4, 5
```

## 7. Tests

### Stage 1: Unit Tests (`tests/test_power_curve_fitting.R`)

| Test | What it validates |
|------|-------------------|
| Config consistency | Grid dimensions, thresholds, v3.1.0 settings |
| Path helpers | get_paths(), get_model_dir(), shard_cache_dir |
| File pattern parsing | DE results + metadata files, both modes |
| Logistic regression fit | GLM convergence and coefficient signs |
| Model formula (v3.1.0) | Unified 6-predictor formula, cell_type variants |
| Shared gene identification | Intersection logic |
| Shard cache paths | Path generation for cache files |
| Prediction grid (v3.1.0) | n_cells and reads_per_cell axes |
| Coefficient formatting | NULL handling, column names |
| Model variant count | Confirms 8 total models |

### Stage 2: Toy Integration Test (`toy/run_toy_test.R`)

- Tests shard cache loading
- Tests metadata parsing (compute_iteration_depth)
- Builds mini unified dataset from both modes (per-iteration)
- Fits v3.1.0 unified formula
- Validates predictions in [0,1]
- Tests prediction grid generation

## 8. What Is NOT Tested

- **Dispersion accuracy**: DESeq2 dispersions from GT half assumed reasonable
- **Model family comparison**: Logistic regression chosen for binary data; not compared to mixed-effects logistic (glmer) or quasi-binomial
- **Extrapolation**: Model fitted on 5-600 donors, 5-3000 cells; outside range unreliable
- **Cross-dataset generalisability**: Fitted on Roussos data only
- **Gene random effects**: Deliberately excluded for generalisability
- **Covariate interactions**: Per-CT models don't include gene × design interactions
