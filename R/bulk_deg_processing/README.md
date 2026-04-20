# Bulk DEG Processing Pipeline

## Purpose

This pipeline processes bulk RNA-seq differential expression data from multi-tissue studies. It filters genes for tissue consistency, annotates with gene names/descriptions, and analyzes sex chromosome representation across log fold-change intervals.

**Input:** Raw bulk DEG data (tissue-level statistics)  
**Output:** Filtered consistent gene lists, annotated tables, sex chromosome analysis

---

## Pipeline Steps

1. **Extract Bulk DEGs** (`01_Extract_Bulk_DEGs.R`)
   - Load raw bulk DEG file
   - Filter genes by tissue coverage threshold
   - Save to `.qs` format

2. **Collect Consistent DEGs** (`02_Collect_Consistent_DEGs.R`)
   - Filter for genes significant in ≥90% of tissues
   - Clean Ensembl gene IDs (remove version suffixes)
   - Calculate average log2FoldChange, lfcSE, padj across tissues
   - Output filtered gene list

3. **Annotate DEGs** (`03_Annotate_DEGs.R`)
   - Query BioMart for gene symbols and descriptions
   - Merge annotations with filtered gene list
   - Output annotated table

4. **Sex Chromosome Analysis** (`04_Sex_Chr_Analysis.R`)
   - Annotate genes with chromosome locations
   - Calculate proportion of sex chromosome genes by LFC interval
   - Output summary statistics

5. **Visualise Bulk vs SC Concordance** (`05_Visualise_Bulk_vs_SC.R`)
   - Pre-fetch gene annotations from BioMart (Ensembl 109) once for all cell types
   - Generate concordance scatter plots (bulk LFC vs SC LFC)
   - Generate 4-category recovery statistics (Recovered/Insignificant/NA-Filtered/Absent)
   - Direction-aware recovery: requires FDR < 0.05 AND matching LFC sign
   - Generate volcano plots (SC DEGs with shared/unique translucency)
   - Generate LFC and SE distribution histograms (coding vs non-coding)
   - Generate chromosome distribution plots (Autosome/X/Y; patches folded into Autosome)
   - Generate bulk volcano plot (global)
   - Output: 15 plots per cell type × per SC dataset + 1 global bulk volcano

---

## Inputs Expected

### Primary Input
- **File:** `Datalink/Bulk_Data/Bulk_Significant_DEGs.txt`
- **Format:** Tab/space-delimited text file
- **Required columns:**
  - `gene` - Ensembl gene ID (e.g., ENSG00000123456.12)
  - `tissue` - Tissue name
  - `effsize` - Effect size (log2 fold change)
  - `effsize_se` - Standard error of effect size
  - `lfsr` - Local false sign rate (p-value analog)

### Prerequisites
- Steps 1-4: None (standalone pipeline)
- Step 5: Requires Full Ground Truth pipeline outputs (`Datalink/Full_Ground_Truth/{CT}/{Dataset}_Full_DESeq2_Results.qs`)

---

## Outputs Produced

### Intermediate Files (cached)
- `Datalink/Bulk_Analysis/Bulk_DE_genes_All.qs` - All genes passing tissue coverage threshold
- `Datalink/Bulk_Analysis/Filtered_Consistent_Genes.csv` - Genes significant in ≥90% of tissues

### Final Outputs
- `Datalink/Bulk_Analysis/Filtered_Consistent_Genes_with_Descriptions.csv` - Annotated gene list
- `Datalink/Bulk_Analysis/Sex_Chr_Proportions_Intervals.csv` - Sex chromosome analysis summary
- `Datalink/Bulk_Analysis/{CT}/*.png` - Bulk-vs-SC concordance plots (15 per cell type × per SC dataset)
- `Datalink/Bulk_Analysis/Bulk_Results/Bulk_Volcano_Plot.png` - Bulk volcano plot

### Toy Test Outputs
- `Datalink/Bulk_Analysis/toy/` - All outputs from minimal test run (first 1000 genes)

---

## Configuration

Edit `config.R` to customize:

```r
# Filtering thresholds
TISSUE_CONSISTENCY_THRESHOLD <- 0.90  # Genes must be significant in ≥90% of tissues
P_VALUE_THRESHOLD <- 0.05             # lfsr threshold for significance
LFC_THRESHOLD <- 0                    # Minimum absolute log2 fold change (0 = no filter)

# Sex chromosome analysis intervals
LFC_INTERVALS <- list(
  "1.0 to 2.0" = c(1.0, 2.0),
  "2.0 to 3.0" = c(2.0, 3.0),
  "> 3.0"      = c(3.0, Inf)
)

# Visualization parameters (Step 05)
SC_DATASETS_FOR_CONCORDANCE <- c("Roussos")  # Which SC datasets to compare against bulk
VIS_LABEL_GENES <- TRUE                      # Label genes in plots
VIS_CALC_SUBGROUP_LINES <- FALSE             # Separate trend lines for autosomes/sex chr
VIS_HIST_XLIM <- c(2, 8)                     # Histogram X-axis limits
VIS_HIST_YLIM <- c(0, 5)                     # Histogram Y-axis limits

# BioMart version (Ensembl 109 — best available for GTEx v8 gene IDs)
ENSEMBL_VERSION <- 109

# Recovery classification (Step 05)
# 4 categories: Recovered, Insignificant, NA-Filtered, Absent
# Direction-aware: Recovered requires FDR < 0.05 AND matching LFC sign
RECOVERY_COLORS <- c(Recovered="#00BA38", Insignificant="#E7B800",
                     `NA-Filtered`="#F8766D", Absent="#999999")
```

---

## Running the Pipeline

### Full Pipeline (Stage 3 - Production)
```bash
# Run all steps
Rscript R/bulk_deg_processing/run_pipeline.R

# Run specific steps (e.g., steps 1-2 only)
Rscript R/bulk_deg_processing/run_pipeline.R 1,2

# Run with custom parameters (override N_CORES)
Rscript R/bulk_deg_processing/run_pipeline.R all 8
```

**RStudio:** Open `run_pipeline.R` and click "Source" or "Source as Background Job"

### Toy Test (Stage 2 - Integration Validation)
```bash
# Command line
Rscript R/bulk_deg_processing/toy/run_toy_test.R

# RStudio
# Open toy/run_toy_test.R and click "Source"
```

### Unit Tests (Stage 1 - Component Validation)
```bash
# Run all unit tests
Rscript R/bulk_deg_processing/tests/test_gene_id_cleaning.R
Rscript R/bulk_deg_processing/tests/test_tissue_counting.R
Rscript R/bulk_deg_processing/tests/test_consistency_filter.R
Rscript R/bulk_deg_processing/tests/test_chromosome_annotation.R
Rscript R/bulk_deg_processing/tests/test_interval_binning.R
Rscript R/bulk_deg_processing/tests/test_biomart_query.R
```

---

## Tests Implemented

### Stage 1: Unit Tests (Component Validation)

**Test: Gene ID Cleaning** (`tests/test_gene_id_cleaning.R`)
- Validates removal of Ensembl version suffixes
- Tests: `ENSG00000123456.12` → `ENSG00000123456`
- Edge cases: IDs without versions, empty input, NA values

**Test: Tissue Counting** (`tests/test_tissue_counting.R`)
- Validates calculation of tissues per gene
- Tests: Correct counts for genes appearing in multiple tissues
- Edge cases: Single tissue, all tissues, duplicates

**Test: Consistency Filtering** (`tests/test_consistency_filter.R`)
- Validates 90% threshold logic
- Tests: Genes meeting/not meeting threshold
- Edge cases: Exactly 90%, boundary conditions

**Test: Chromosome Annotation** (`tests/test_chromosome_annotation.R`)
- Validates gene-to-chromosome mapping via org.Hs.eg.db
- Tests: Known genes map to correct chromosomes
- Edge cases: Unmapped genes, sex chromosomes (X/Y)

**Test: Interval Binning** (`tests/test_interval_binning.R`)
- Validates LFC interval assignment
- Tests: Correct binning for 0.5-1.0, 1.0-2.0, >2.0 intervals
- Edge cases: Boundary values, negative LFC, zero LFC

**Test: BioMart Querying** (`tests/test_biomart_query.R`)
- Validates annotation fetching from BioMart
- Tests: Known genes return expected symbols/descriptions
- Edge cases: Connection failures, unmapped genes

### Stage 2: Toy Integration Test

**Test: End-to-End Pipeline** (`toy/run_toy_test.R`)
- Runs all 4 steps on first 1000 genes from raw data
- Validates: All outputs generated, no NA in critical columns
- Manual inspection: Check CSV files for reasonable values
- Execution contexts: Command line, RStudio Source, RStudio Background Job

### Stage 3: Production Run

**Validation:** End-of-run report with success/failure counts
- Logs to `Datalink/Run_Output_Logs/bulk_deg_processing_*.log`
- Reports: Total genes processed, filtering statistics, annotation success rate

---

## What Is NOT Tested

### Assumptions Made
1. **Input file format is correct** - Pipeline assumes `Bulk_Significant_DEGs.txt` has required columns with correct data types. No validation of raw input file structure.

2. **Ensembl gene IDs are valid** - Pipeline assumes gene IDs are real Ensembl IDs. No validation that IDs exist in current Ensembl release.

3. **Tissue names are consistent** - Pipeline assumes tissue names are standardized (no typos, no case variations). No normalization of tissue names.

4. **Effect sizes are on log2 scale** - Pipeline assumes `effsize` is log2 fold change. No validation of scale or transformation.

5. **BioMart is accessible** - Pipeline assumes internet connection and BioMart server availability. Fails immediately if connection fails (no retry logic, no cached fallback).

6. **org.Hs.eg.db is up-to-date** - Chromosome annotations depend on installed package version. No version checking or update prompts.

### Edge Cases Not Covered
1. **Empty input file** - Pipeline will fail if input file has no data rows
2. **All genes filtered out** - Pipeline will produce empty outputs if no genes meet thresholds
3. **Duplicate gene-tissue pairs** - Pipeline does not check for or handle duplicates in raw data
4. **Mixed species** - Pipeline assumes human data (uses org.Hs.eg.db)
5. **Very large files** - No memory optimization for files >10GB

### Known Limitations
1. **BioMart version drift** - Gene symbols/descriptions may change between Ensembl releases
2. **No provenance tracking** - Output files do not record which version of pipeline/parameters produced them
3. **No incremental updates** - Re-running pipeline regenerates all outputs (no smart caching of unchanged data)

---

## Troubleshooting

### BioMart Connection Failures
```
Error: biomaRt query failed
```
**Solution:** Check internet connection, try again later (BioMart servers occasionally down)

### Missing Input File
```
Error: Input file not found: Datalink/Bulk_Data/Bulk_Significant_DEGs.txt
```
**Solution:** Verify file path in `config.R`, ensure raw data file exists

### No Genes Pass Filter
```
Warning: 0 genes passed consistency threshold
```
**Solution:** Lower `TISSUE_CONSISTENCY_THRESHOLD` in `config.R` (e.g., 0.80 instead of 0.90)

### Checkpoint Cache Invalid
```
Warning: Cached data failed validation. Regenerating...
```
**Solution:** This is normal - pipeline automatically regenerates. If persists, delete cache files in `Datalink/Bulk_Analysis/` and re-run.

---

## Version History

- **v1.2.0** (2026-02-17) - Recovery category refinement + Ensembl version fix
  - 4-category recovery: Recovered, Insignificant, NA-Filtered, Absent (was 3-category)
  - Direction-aware recovery in both Step 05 and Step 06 (Step 05 was direction-agnostic)
  - BioMart switched to Ensembl 109 (resolves 77/79 bulk genes vs 74/79 with latest)
  - BioMart query hoisted out of cell-type loop (was querying ~35k genes 3× redundantly)
  - Patch/Haplotype chromosomes folded into Autosome; "ID Not Found" category eliminated
  - Consistent recovery classification across Steps 05 and 06

- **v1.1.0** (2026-02-17) - Added Step 06: Recovery by Expression Level
  - Expression-stratified recovery analysis (baseMean quintiles)
  - Two scales: raw baseMean and log2(baseMean + 1)
  - Stacked bar charts per quintile
  - CSV summaries per cell type × dataset

- **v1.0.0** (2026-02-10) - Initial pipeline conversion from standalone scripts
  - Three-stage validation framework
  - Checkpoint caching with validation
  - Auto-logging for background jobs
  - RStudio sourcing compatibility

---

## References

- Original scripts: `R/Bulk Analysis/` (preserved for reference)
