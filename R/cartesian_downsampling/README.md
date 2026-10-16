# Cartesian Downsampling Pipeline

This pipeline implements empirical power analysis through systematic downsampling of PsychAD single-cell RNA-seq data (`Roussos` in code and paths) across a Cartesian grid of donor counts and cells (or UMIs) per donor.

**Pipeline Version:** 2.3.0 (`config.R`)
**Framework:** Three-Stage Validation

---

## 1. Pipeline Purpose

The Cartesian Downsampling Pipeline estimates statistical power for detecting sex-biased differential expression by:
1. Splitting the filtered donors 50:50 into a ground-truth half and a validation half, repeated over `N_ROBUSTNESS_RUNS` (10) independent splits
2. Systematically subsampling the validation half across a grid of (N donors) x (N cells per donor), or (N donors) x (UMI budget per donor) in reads mode
3. Running pseudobulk DE at each grid point with DESeq2, edgeR and limma-voom
4. Computing precision and recall against the ground-truth half, per |LFC| bin, expression bin and chromosome category
5. Pooling the three methods and visualising power as discrete and interpolated heatmaps

---

## 2. Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). The stages call `run_pipeline.R` with `SHARD_DATASET=Roussos`, `DE_METHODS=deseq2,edger,limma` and `N_ROBUSTNESS_RUNS=10`:

| Stage | Steps | Notes |
|-------|-------|-------|
| `gt_manifest` | 1 (cells), then 1.5 for cells and reads | Ground truth is mode-independent; it is written under `cells/` and copied to `reads/` |
| `grid` | 2 | PBS array: one task per (mode, chunk), `JOB_N_CHUNKS` chunks per mode |
| `cartesian_cells`, `cartesian_reads` | 3, 3.5, 4, 5, 6 | With `CARTESIAN_INCLUDE_ALL_METHODS=true` so the pooled `all_methods` outputs are plotted |
| `crossscore` | `10_CrossScore_Methods.R` | Run directly, both modes, all runs |
| `ad_gt_manifest`, `ad_grid`, `ad_cartesian` | 1 and 1.5, then 2 (PBS array), then 3, 3.5, 4, 5, 6 | AD vs control contrast: `EPA_CONTRAST=ad_vs_control`, `SHARD_DATASET=Roussos_ADvsCtrl`, cells mode only, same seeds. Outputs under `Cartesian_Downsampling/{CT}/Roussos_ADvsCtrl/` |

---

## 3. Inputs Expected

### Data Sources

| Input | Path | Format | Description |
|-------|------|--------|-------------|
| Donor Shards | `Datalink/Shards/{CellType}/Roussos/` | `.qs` files | One SCE per donor |
| Manifest | `Datalink/Shards/{CellType}/Roussos/00_MASTER_MANIFEST.qs` | `.qs` | Donor metadata and file paths |

### Data Structure Requirements

**Donor Shard (SingleCellExperiment):**
- `assay(sce, "counts")` - Raw integer count matrix (genes x cells)
- `colData(sce)$donor_id` - Donor identifier
- `colData(sce)$sex` - "male" or "female"
- `colData(sce)$development_stage` - Age string (e.g., "65 year-old stage")
- `colData(sce)$Source`, `self_reported_ethnicity`, `disease`, `pmi` - Design covariates

**Manifest Structure:**
```r
manifest <- qs::qread("00_MASTER_MANIFEST.qs")
# manifest$manifest is a data.frame with columns:
#   - donor_id: character
#   - cell_count: integer
#   - file_path: character (path to donor shard)
```

### Prerequisites

- Sharding pipeline (including PMI annotation) must have completed successfully
- Donor shards must contain raw counts (not normalised)

---

## 4. Method

- **Grid:** `DONOR_STEPS` = 5, 10, 15, 25, 50, 75, 100, 150, 200, 300, 500, 600 (12 levels). Cells mode: `CELL_STEPS` = 5, 10, 20, 50, 75, 100, 250, 500, 750, 1000, 1500, 2000, 3000, Max (14 levels). Reads mode: `READ_STEPS` = 25,000 to 15,000,000 total UMIs per donor plus Max (13 levels).
- **Iterations:** `N_SUBSAMPLES` = 20 per grid point per run.
- **Design:** `~ sex + age_scaled + pmi_scaled + Source + self_reported_ethnicity + disease`, falling back to `~ sex` if rank-deficient.
- **DE:** DESeq2 (`fitType = "glmGamPoi"`), edgeR quasi-likelihood (`glmQLFit`, robust), limma-voom (`eBayes`, robust). All LFCs are shrunk with ashr so effect sizes are comparable across methods.
- **Truth set:** ground-truth genes with padj < 0.05 (BH) and |LFC| >= 0.2. Precision and recall are reported in `LFC_INTERVALS`: Tiny (0.2-1.0), Small (1.0-2.0), Medium (2.0-3.0), Large (3.0+), and in `BASEMEAN_INTERVALS` for the Tiny and Small bins.
- **Pooling:** step 3.5 pools the three methods trial-weighted into `all_methods`; cell types stay separate.

---

## 5. Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_Generate_Ground_Truth.R` | 50:50 donor split and ground-truth DE per run and method |
| 1.5 | `01b_Build_HPC_Job_Manifest.R` | Cost-balanced job manifest for PBS array chunks |
| 2 | `02_Run_Downsampling_Grid.R` | Executes the donor x cells (or UMIs) grid |
| 3 | `03_Aggregate_Results.R` | Precision/recall per grid point, bin and chromosome category |
| 3.5 | `08_Aggregate_CrossMethod.R` | Trial-weighted pooling of methods into `all_methods` |
| 4 | `04_Visualise_Power.R` | Discrete heatmaps and summary plots |
| 5 | `07_Visualise_Continuous_Grid.R` | Bilinear-interpolated continuous surfaces |
| 6 | `09_Visualise_Figure_TS.R` | Tiny/Small bins x four expression categories (Figure T / Figure S) |
| (direct) | `10_CrossScore_Methods.R` | 3x3 cross-method scoring (subsample method x truth method) |

---

## 6. Outputs Produced

### Per Cell Type, Mode, Method, Run

Under `Datalink/Cartesian_Downsampling/{CellType}/Roussos/{cells|reads}/{deseq2|edger|limma}/Run_{i}/`:

| Output | File | Description |
|--------|------|-------------|
| Ground Truth | `Ground_Truth.qs` | DE results from the ground-truth half |
| Ground Truth Log | `Ground_Truth_sampling_metadata.qs` | Donors/cells used |
| Validation Donors | `Validation_Donors.qs` | Donors available for subsampling |
| Subsample Results | `N{d}_C{c}_I{i}.qs` (cells) or `N{d}_R{r}_I{i}.qs` (reads) | DE results per grid point and iteration |
| Subsample Logs | `..._sampling_metadata.qs` | Cells/UMIs sampled |

### Aggregated Outputs

Under `{CellType}/Roussos/{mode}/{method or all_methods}/`:

| Output | Path | Description |
|--------|------|-------------|
| Summary Data | `Plots/Aggregated_Summary_Data.qs` | Precision, recall, read stats (input to `R/figure_data_generation` step 8 and the composers) |
| Heatmaps | `Plots/Mean_{Recall,Precision}_{Total,Autosomal,SexChr}.png` | Discrete heatmaps by \|LFC\| bin |
| P-value quartiles | `Plots/Mean_{Recall,Precision}_PValQuartile_AllMethods.png` | Recall/precision by ground-truth padj quartile |
| Figure T / S | `Plots/Figure_TS/` | Tiny/Small bin x expression panels |
| Continuous surfaces | `Continuous_Viz/`, `Continuous_Viz_Log/` | Interpolated surfaces (linear and semi-log) |

### Cross-method scoring

`Datalink/Cartesian_Downsampling/Cross_Score/cross_score_matrix.csv`, rendered by `R/figure_composition/compose_crossscore.R`.

### HPC manifests

`Datalink/Cartesian_Downsampling/_hpc_job_manifests/manifest_Roussos_{cells,reads}_runsXX-YY_NNN_chunks.csv`.

---

## 7. Tests Implemented

### Stage 1: Unit Tests (`tests/`)

| Test File | What's Validated |
|-----------|------------------|
| `test_age_parsing.R` | Age string parsing, NA handling, edge cases |
| `test_aggregation.R` | SCE loading, cell sampling, pseudobulk aggregation, gene alignment |
| `test_deseq_worker.R` | DESeq2 execution, formula construction, LFC shrinkage, result saving |
| `test_job_planning.R` | Job grid generation, filtering by cell count, file naming |
| `test_reads_sampling.R` | Reads-mode sampling (multinomial overflow, exact targeting) |
| `test_ad_sex_balance.R` | AD vs control: AD status x sex ground-truth split, sex-balanced subsamples, sex contrast unchanged |
| `test_fittype_equivalence.R` | `glmGamPoi` vs `parametric` DESeq2 fits agree (LFC and -log10 padj correlation >= 0.99) |

**Run Stage 1:**
```bash
Rscript R/cartesian_downsampling/tests/test_age_parsing.R
Rscript R/cartesian_downsampling/tests/test_aggregation.R
Rscript R/cartesian_downsampling/tests/test_deseq_worker.R
Rscript R/cartesian_downsampling/tests/test_job_planning.R
Rscript R/cartesian_downsampling/tests/test_reads_sampling.R
```

### Stage 2: Toy Integration Test (`toy/run_toy_test.R`)

- Cell type: Micro (smallest)
- Robustness runs: 1
- Donor steps: 5, 10, 25
- Cell steps: 50, 100, Max
- Iterations: 2 per grid point

A reduced end-to-end run is also available with `MINI_RUN=true` (writes to `Datalink/MINI_Cartesian_Downsampling/`).

### Stage 3: Production Run

The `gt_manifest`, `grid`, `cartesian_cells`, `cartesian_reads` and `crossscore` stages of `pipeline/`: 3 cell types, 2 modes, 3 methods, 10 runs, full grid.

---

## 8. What Is NOT Tested

### Data Quality Assumptions
- [ ] **Donor shards contain raw counts**: if shards contain TPM/CPM, DE results are invalid.
- [ ] **Metadata is accurate**: sex labels, age, PMI and donor IDs are assumed correct (sex labels are checked separately by `R/eda_roussos/Validate_Sex_Labels.R`).

### Statistical Assumptions
- [ ] **Ground truth is "true"**: the ground-truth half is an estimate with its own error; false negatives in the truth lower apparent precision.
- [ ] **LFC bins map to effect sizes**: bin boundaries are a convention.

### Computational Assumptions
- [ ] **Parallel jobs are independent**: no shared state between workers is assumed.
- [ ] **Memory is sufficient**: large grid points (3000 cells x 600 donors) need substantial RAM.

### Known Limitations
- [ ] **No multiple testing correction across the grid**: each grid point is treated independently.
- [ ] **Chromosome mapping failures**: genes without a chromosome annotation are excluded from the Autosomal/SexChr categories.

---

## Usage

The supported way to run this is `pipeline/run.sh` (`--only <stage>` reruns one stage). The commands below are what the stage runs, for reading and debugging inside the container.

### Run Specific Steps
```bash
# Ground truth and job manifest (cells mode)
DOWNSAMPLE_MODE=cells Rscript R/cartesian_downsampling/run_pipeline.R 1,1.5 24

# One grid chunk
DOWNSAMPLE_MODE=reads JOB_CHUNK=3 Rscript R/cartesian_downsampling/run_pipeline.R 2 24

# Aggregation and figures, including pooled all_methods
CARTESIAN_INCLUDE_ALL_METHODS=true Rscript R/cartesian_downsampling/run_pipeline.R 3,3.5,4,5,6 24
```

### Configuration

`config.R` settings (most overridable by environment variable):
- `SHARD_DATASET` - Dataset (default `Roussos`)
- `EPA_CONTRAST` - `sex` (default) or `ad_vs_control` (see `R/shared/contrast.R`; must match `SHARD_DATASET`)
- `DOWNSAMPLE_MODE` - `cells` or `reads`
- `DONOR_STEPS` / `CELL_STEPS` / `READ_STEPS` - Grid dimensions
- `N_ROBUSTNESS_RUNS` (10), `RUN_START` / `RUN_END` - Runs to process
- `N_SUBSAMPLES` (20) - Iterations per grid point
- `DE_METHODS` - Subset of `deseq2`, `edger`, `limma`
- `JOB_N_CHUNKS` / `JOB_CHUNK` - Array chunking
- `FORCE_RECALC` - Whether to overwrite existing results

---

## References

- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
- Stephens M (2017). "False discovery rates: a new deal." Biostatistics.

---

*Last Updated: 2026-10-06*
