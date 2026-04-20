# Cartesian Downsampling Pipeline

This pipeline implements empirical power analysis through systematic downsampling of single-cell RNA-seq data across a Cartesian grid of donor counts and cells per donor.

**Last Updated:** 2026-03-31  
**Pipeline Version:** 2.2.1  
**Framework:** Three-Stage Validation

---

## 1. Pipeline Purpose

The Cartesian Downsampling Pipeline estimates statistical power for detecting sex-biased differential expression by:
1. Creating a "ground truth" from 50% of available donors
2. Systematically subsampling the remaining donors across a grid of (N donors) x (N cells per donor)
3. Running DESeq2 differential expression at each grid point
4. Computing precision and recall relative to ground truth
5. Visualizing power as heatmaps

This provides empirical power estimates that can be compared against simulation-based approaches (e.g., hierarchicell).

---

## 2. Inputs Expected

### Data Sources

| Input | Path | Format | Description |
|-------|------|--------|-------------|
| Donor Shards | `Datalink/Shards/{CellType}/{Dataset}/` | `.qs` files | One SCE per donor |
| Manifest | `Datalink/Shards/{CellType}/{Dataset}/00_MASTER_MANIFEST.qs` | `.qs` | Donor metadata and file paths |

### Data Structure Requirements

**Donor Shard (SingleCellExperiment):**
- `assay(sce, "counts")` - Raw integer count matrix (genes x cells)
- `colData(sce)$donor_id` - Donor identifier
- `colData(sce)$sex` - "male" or "female"
- `colData(sce)$development_stage` - Age string (e.g., "65 year-old stage")

**Manifest Structure:**
```r
manifest <- qs::qread("00_MASTER_MANIFEST.qs")
# manifest$manifest is a data.frame with columns:
#   - donor_id: character
#   - file_path: character (path to donor shard)
```

### Prerequisites

- Sharding pipeline must have completed successfully
- Donor shards must contain raw counts (not normalized)
- At least 10 donors per sex required for meaningful analysis

---

## 3. Outputs Produced

### Per Cell Type, Per Run

| Output | Path | Format | Description |
|--------|------|--------|-------------|
| Ground Truth | `{output}/Run_{i}/Ground_Truth.qs` | DESeqResults | DE results from 50% of donors |
| Ground Truth Log | `{output}/Run_{i}/Ground_Truth_sampling_metadata.qs` | list | Cell indices used |
| Subsample Results | `{output}/Run_{i}/N{d}_C{c}_I{i}.qs` | DESeqResults | DE results per grid point |
| Subsample Logs | `{output}/Run_{i}/N{d}_C{c}_I{i}_sampling_metadata.qs` | list | Cell indices used |

### Aggregated Outputs

| Output | Path | Format | Description |
|--------|------|--------|-------------|
| Summary Data | `{plots}/Aggregated_Summary_Data.qs` | list | Precision, recall, read stats |
| Precision Heatmap | `{plots}/Mean_Precision_{category}.png` | PNG | Heatmap by LFC bin |
| Recall Heatmap | `{plots}/Mean_Recall_{category}.png` | PNG | Heatmap by LFC bin |
| XIST Check | `{plots}/XIST_Sanity_Check.png` | PNG | XIST detection rate |
| Read Depth | `{plots}/Total_Reads_Per_Matrix.png` | PNG | Read depth heatmap |

### Astro Variant Outputs (Step 6)

When Step 6 is run, additional comparison plots are generated for Astrocyte variant datasets (if available):

| Output | Path | Format | Description |
|--------|------|--------|-------------|
| Individual Heatmaps | `Astro/Variants/{mode}/Plots/{Dataset}_Recall_{category}.png` | PNG | Per-variant recall heatmaps |
| Comparison Grid | `Astro/Variants/{mode}/Plots/Comparison_{metric}_Total.png` | PNG | Side-by-side variant comparison |
| Difference Heatmaps | `Astro/Variants/{mode}/Plots/Difference_{Variant}_vs_Standard_Recall.png` | PNG | Delta plots (Variant - Standard) |
| Summary Table | `Astro/Variants/{mode}/Plots/Variant_Summary_Table.csv` | CSV | Max recall per variant/bin |

**Astro Variants:**
- `Roussos` - Standard (deep sequencing, ~9,456 reads/cell)
- `Roussos_DepthEq` - Depth-equalized to match Microglia (~4,426 reads/cell)
- `Roussos_ShallowFewer` - Depth-equalized + cell count matched to Microglia

### Output Paths

```
Datalink/Roussos_Analysis/{CellType}/
├── Roussos_Empirical_Downsampling/
│   ├── Run_1/
│   │   ├── Ground_Truth.qs
│   │   ├── Ground_Truth_sampling_metadata.qs
│   │   ├── N5_C50_I1.qs
│   │   ├── N5_C50_I1_sampling_metadata.qs
│   │   └── ...
│   ├── Run_2/
│   └── ...
└── Roussos_Empirical_Downsampling_Plots/
    ├── Aggregated_Summary_Data.qs
    ├── Mean_Precision_Total.png
    ├── Mean_Recall_Total.png
    └── ...
```

---

## 4. Tests Implemented

### Stage 1: Unit Tests (`tests/`)

| Test File | Functions Tested | What's Validated |
|-----------|------------------|------------------|
| `test_age_parsing.R` | `parse_age_stage_to_years()`, `is_unknown_string()` | Correct parsing of age strings, NA handling, edge cases |
| `test_aggregation.R` | `aggregate_from_shards()` | SCE loading, cell sampling, pseudobulk aggregation, gene alignment |
| `test_deseq_worker.R` | `run_deseq_worker()` | DESeq2 execution, formula construction, LFC shrinkage, result saving |
| `test_job_planning.R` | Job grid generation | Correct filtering by cell count, sex ratio preservation, file naming |

**Run Stage 1:**
```bash
Rscript R/cartesian_downsampling/tests/test_age_parsing.R
Rscript R/cartesian_downsampling/tests/test_aggregation.R
Rscript R/cartesian_downsampling/tests/test_deseq_worker.R
Rscript R/cartesian_downsampling/tests/test_job_planning.R
```

### Stage 2: Toy Integration Test (`toy/`)

| Test | Parameters | What's Validated |
|------|------------|------------------|
| `run_toy_test.R` | 1 cell type, 1 run, 3x3 grid, 2 iterations | End-to-end pipeline execution, output file creation, heatmap generation |

**Toy Test Parameters:**
- Cell type: Micro (smallest dataset)
- Robustness runs: 1
- Donor steps: 5, 10, 25
- Cell steps: 50, 100, NA (Max)
- Iterations: 2 per grid point

**Expected Outputs:**
- 1 ground truth file
- 18 subsample results (3 donors x 3 cells x 2 iterations)
- Aggregated summary
- Heatmaps for visual inspection

**Run Stage 2:**
```bash
Rscript R/cartesian_downsampling/toy/run_toy_test.R
```

### Stage 3: Production Run

| Test | Parameters | What's Validated |
|------|------------|------------------|
| `run_pipeline.R` | All cell types, 6 runs, full grid | Complete analysis at scale |

**Run Stage 3:**
```bash
Rscript R/cartesian_downsampling/run_pipeline.R
```

---

## 5. What Is NOT Tested

**CRITICAL:** The following assumptions are made but NOT explicitly validated by tests:

### Data Quality Assumptions
- [ ] **Donor shards contain raw counts** - We assume counts are integers, not normalized. If shards contain TPM/CPM, DESeq2 will produce incorrect results.
- [ ] **No duplicate donors across datasets** - If the same donor appears in multiple datasets, they may be sampled multiple times.
- [ ] **Metadata is accurate** - Sex labels, age values, and donor IDs are assumed correct.

### Statistical Assumptions
- [ ] **Ground truth is "true"** - The 50% split ground truth is treated as the reference, but it's still an estimate with uncertainty.
- [ ] **DESeq2 model is appropriate** - The `~ sex + age_scaled + covariates` model assumes these are the relevant factors.
- [ ] **LFC bins map to effect sizes** - The mapping of fold changes to "Small/Medium/Large" bins is arbitrary.

### Computational Assumptions
- [ ] **Parallel jobs are independent** - We assume no race conditions or shared state between parallel workers.
- [ ] **Memory is sufficient** - Large cell counts (3000 cells x 500 donors) may exceed available RAM.

### Edge Cases NOT Tested
- [ ] **Very small donor counts** - Behavior with <3 donors per sex is undefined.
- [ ] **Missing metadata columns** - Pipeline may fail silently if expected columns are missing.
- [ ] **Manifest-shard mismatches** - If manifest lists files that don't exist, errors may not be clear.

### Known Limitations
- [ ] **No cross-validation** - Ground truth uses a single 50/50 split, not k-fold.
- [ ] **No multiple testing correction across grid** - Each grid point is treated independently.
- [ ] **Chromosome mapping failures** - Genes not in org.Hs.eg.db are silently excluded from autosomal/sex chromosome categories.

---

## Usage

### Run Complete Pipeline
```bash
cd R/cartesian_downsampling
Rscript run_pipeline.R
```

### Run Specific Steps
```bash
# Run only steps 1 and 2
Rscript run_pipeline.R 1,2

# Run only visualization
Rscript run_pipeline.R 4

# Run Astro variants comparison (requires Steps 1-3 for all variants first)
Rscript run_pipeline.R 6
```

### Configuration

Edit `config.R` to modify:
- `SHARD_DATASET` - Which dataset to use (default: "Roussos")
- `DONOR_STEPS` / `CELL_STEPS` - Grid dimensions
- `N_ROBUSTNESS_RUNS` - Number of independent runs
- `N_SUBSAMPLES` - Iterations per grid point
- `N_CORES` - Parallel workers
- `FORCE_RECALC` - Whether to overwrite existing results

---

## References

- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology.
- Zhu A, Ibrahim JG, Love MI (2019). "Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences." Bioinformatics.

---

*Pipeline Version: 2.2.1*  
*Framework: Three-Stage Validation*  
*Created: 2026-02-06*  
*Updated: 2026-03-31*  
