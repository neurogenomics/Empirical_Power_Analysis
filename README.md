# Single-Cell Study Designs Are Systematically Underpowered for Small-Effect Genes

![Study overview](Figure1_Diagram.png)

<sub>Study overview diagram designed with Claude.</sub>

Code and paper-ready figures accompanying:

> **Single-Cell Study Designs Are Systematically Underpowered for Small-Effect Genes**
> Chris Crampton, Salman Fawad, Toby Clark, Hiru Dash, Bence Kövér, Garry Cotton, Donghoon Lee,  Eugene Duff, Leonardo Bottolo, Paul Matthews, Nathan Skene
> bioRxiv preprint (2026), doi: [10.64898/2026.09.17.752063](https://doi.org/10.64898/2026.09.17.752063)

## Overview

This repository contains the complete analysis pipeline used to empirically measure detection power for pseudobulk differential expression in human brain single-nucleus RNA-seq. The PsychAD cohort (1,494 donors before filtering) is split ten times into sex-stratified halves: one half defines the ground truth, and the other is systematically downsampled across a grid of donor count × cells per donor (and, separately, donor count × UMI budget). Every pseudo-experiment is analysed with three DE packages (DESeq2, edgeR and limma-voom) on identical samples, and recall and precision are scored against the ground-truth half. Sex-biased expression is used as the case-control variable because cohorts split close to 50:50 and sex-biased genes give a reproducible ground truth.

Outputs include:

- The empirical recall/precision dataset across the design space, per cell type and DE method, for three glial cell types (astrocytes, microglia, oligodendrocytes), with continuous (interpolated) surfaces.
- A gradient-boosted (XGBoost) per-gene power model with a dependency-light `predict_power()` wrapper for study design.
- A gene-level detectability classification (transcription-factor families, PANTHER protein classes, loss-of-function intolerant genes).
- A one-command pipeline (`pipeline/`) that regenerates every figure, table, quoted number and supplementary data file from the raw data.

## Repository layout

```
.
├── pipeline/                  # One-command reproduction (run.sh, stage graph, stages, tests); see pipeline/README.md
├── R/                         # Analysis code called by the pipeline stages
│   ├── sharding/              #   Raw PsychAD h5ad -> donor-level shards (plus AD vs control shards)
│   ├── eda_roussos/           #   PCA gallery and sex-label validation
│   ├── full_ground_truth/     #   Full-cohort PsychAD DE and gene-universe audit
│   ├── de_diagnostics/        #   DEG chromosome grid and shared-gene diagnostics
│   ├── cartesian_downsampling/#   The empirical power grid (cells and UMI-budget modes), cross-method scoring
│   ├── depth_equalization/    #   Depth- and cell-count-equalised truth sets and DEG loss
│   ├── power_curve_fitting/   #   XGBoost power model: training table, search, selection, SHAP, export
│   ├── gene_classification/   #   Per-gene power by TF family / protein class (Figure 6)
│   ├── ad_contrast/           #   AD vs control contrast and the external test of the power model
│   ├── figure_data_generation/#   Source CSVs the manuscript reads
│   ├── figure_composition/    #   Multi-panel figures
│   ├── supplementary_data/    #   Supplementary Data files
│   ├── glial_characteristics/ #   utils.R only (shard readers used by figure_data_generation)
│   └── shared/                #   Helpers used across pipelines
├── Dockerfile, docker/        # The pinned container image; see docs/container.md
├── renv.lock                  # Exact R package versions
├── nature_comms_draft.rmd     # Manuscript source (knitted by the pipeline)
├── nature_comms_supplementary.rmd
├── _paper_values.R            # val("key"): numbers quoted in the text, read from output tables at knit time
├── _supp_fig_index.R          # SF("key"): supplementary figure numbering for both documents
├── references.bib
├── Datalink/                  # Output tree; git tracks only the files the manuscript renders
├── docs/                      # container.md, reproducing-the-paper.md, architecture.md, model card
├── AGENTS.md                  # Orientation for coding agents and new contributors
└── Figure1_Diagram.png        # Study overview (above)
```

Code that is not on the reproduction path (superseded runners, other cohorts,
the archived logistic-regression model, exploratory analyses) is listed in
`.gitignore` and is not part of this repository.

## Paper asset lookup

Every figure, table, and source-data file rendered in the manuscript is present in this repository. Paths below are relative to the repo root. The `{Astro,Micro,Oligo}` brace notation expands to one file per cell type; `{cells,reads}` to the cells-per-donor and UMI-budget downsampling modes.

### Main-text figures

| Figure | Chunk label | File |
| --- | --- | --- |
| Figure 1: Study overview | `FigOverview` | `Figure1_Diagram.png` |
| Figure 2: Oligodendrocyte recall (continuous + discrete) | `FigRecall` | `Datalink/Composed_Figures/figure_cartesian_oligo_recall.png` |
| Figure 3: Recall by adjusted-p-value quartile and expression | `FigRecallPanels` | `Datalink/Composed_Figures/figure_recall_panels_oligo_cells.png` |
| Figure 4: Oligodendrocyte precision (continuous + discrete) | `FigPrecision` | `Datalink/Composed_Figures/figure_cartesian_oligo_precision.png` |
| Figure 5: XGBoost power model (calibration, SHAP, F1) | `FigXGB` | `Datalink/Power_Curve_Fitting/Models/xgboost/figure_xgb_composed.png` |
| Figure 6: Per-gene power by TF family / protein class | `FigGeneGroups` | `Datalink/Gene_Classification/oligo_paper_expr_disp_panels.png` |

### Supplementary figures

Supplementary figures are numbered by `_supp_fig_index.R`, the single source of truth shared by the main text and the supplement (`nature_comms_supplementary.rmd`); re-knitting renumbers both documents automatically. The families below cover every supplementary figure the manuscript renders, in canonical order.

| Supplement family | Symbolic key(s) | Files / glob |
| --- | --- | --- |
| PCA gallery (unfiltered + filtered PsychAD cohort, 6 covariates × 3 cell types) | `pca_gallery` | `Datalink/Exploratory_Analysis/{Astro,Micro,Oligo}/{Roussos_Full,Roussos}/PCA/PCA_{FullCohort,FilteredCohort}_{Log2_Cells,sex,Age_Years,Source,disease,self_reported_ethnicity}_scaled.png` |
| Cross-method (truth-method × downsampling-method) recall/precision grids | `crossmethod_{recall,precision}_{cells,reads}` | `Datalink/Composed_Figures/Cross_Score/crossscore_{recall,precision}_mega_{Astro,Micro,Oligo}_{cells,reads}_discrete_raw.png` |
| Per-method downsampling grids (recall + precision; surface, tiny-effect, small-effect) | `bymethod_{astro,micro,oligo}_{recall,precision}_{surface,tiny,small}` | `Datalink/Composed_Figures/Supplementary/supp_{astro,micro,oligo}_{recall,precision}_{surface,tiny,small}_bymethod_{cells,reads}.png` |
| Method-pooled continuous recall surfaces | `recall_continuous_surfaces` | `Datalink/Cartesian_Downsampling/{Astro,Micro,Oligo}/Roussos/{cells,reads}/all_methods/Continuous_Viz_Log/Mean_Recall_Total_Continuous_SemiLog.png` |
| Method-pooled precision heatmaps | `precision_heatmaps` | `Datalink/Cartesian_Downsampling/{Astro,Micro,Oligo}/Roussos/{cells,reads}/all_methods/Plots/Mean_Precision_Total.png` |
| Method-pooled discrete recall heatmaps | `recall_discrete_heatmaps` | `Datalink/Cartesian_Downsampling/{Astro,Micro,Oligo}/Roussos/{cells,reads}/all_methods/Plots/Mean_Recall_Total.png` |
| Chromosome-stratified recall (autosomal vs sex-chromosome) | `recall_chrom_stratified` | `Datalink/Cartesian_Downsampling/{Astro,Micro,Oligo}/Roussos/cells/all_methods/Plots/Mean_Recall_{Autosomal,SexChr}.png` |
| DEG chromosome-distribution grid | `deg_chrom_distribution` | `Datalink/DE_Diagnostics/Chromosome_Grid/DEG_Chr_Grid_Combined.png` |
| XGBoost SHAP dependence / violin (winner + reweighting configs) | `xgb_shap_dependence`, `xgb_shap_violin`, `xgb_shap_dependence_by_config` | `Datalink/Power_Curve_Fitting/Models/xgboost/shap_figures/{median,winner_alpha_0,winner_alpha_0.5,winner_alpha_1}/shap_{dependence,violin}.png` |

("Roussos" in paths is the lab's legacy name for the PsychAD cohort.)

### Supplementary tables

| Table | File |
| --- | --- |
| Supplementary Table 1: brain sc/snRNA-seq DE studies surveyed | Typed into `nature_comms_supplementary.rmd` |
| Supplementary Table 2: autosomal vs sex-chromosome truth-gene characteristics | `Datalink/_scratch_autosomal_vs_sexchr_summary.csv` |
| Supplementary Table 3: DEG loss under depth and cell-count equalisation | `Datalink/Depth_Eq_DEG_Loss/deg_loss_table.csv` |
| Supplementary Table 4: sex-label validation by chromosomal marker expression | `Datalink/Sex_Validation/sex_exclusion_combined.csv` |

### Supplementary Data

The manuscript cites Supplementary Data 1 and 4. Both are built by `R/supplementary_data/build_supplementary_data.R` into `Datalink/Supplementary_Data/`:

| File | Contents | Built from |
| --- | --- | --- |
| `Supplementary_Data_1_dataset_summary.xlsx` | Dataset characteristics + gene-universe composition | `Datalink/Paper_Source_Data/dataset_statistics.csv`, `Datalink/Full_Ground_Truth/Gene_Audit_Full_Universe_Combined.csv` |
| `Supplementary_Data_4_xgboost_hyperparameters.csv` | Ranking of all 108 XGBoost configurations | `Datalink/Power_Curve_Fitting/Models/xgboost/model_selection.csv` |

Numbers quoted in the manuscript text are read at knit time by `_paper_values.R`, mainly from `Datalink/Paper_Source_Data/` (`dataset_statistics.csv`, `recall_precision_grid.csv`) and the XGBoost model directory.

## Scientific scope

- **Dataset:** PsychAD (Fullard et al., 2025), dorsolateral prefrontal cortex, 1,494 donors before filtering. After filtering (age ≥ 40 years, ≥ 30 cells per donor): 1,293 astrocyte donors, 1,214 microglia donors and 1,303 oligodendrocyte donors. See `Datalink/Paper_Source_Data/dataset_statistics.csv`.
- **Case-control variable:** Sex. Chosen because (i) sex is close to 50:50 in brain cohorts; (ii) sex-biased genes are well documented and reproducible; (iii) the small-effect sex-biased genes are predominantly autosomal, which supports extrapolation to other contrasts.
- **DE methods:** DESeq2, edgeR and limma-voom on pseudobulk aggregates, design `~ source + ethnicity + disease + pmi_scaled + age_scaled + sex`. All log fold changes are shrunk with ashr so effect sizes are comparable across methods. Significance is BH-adjusted p < 0.05; the ground-truth set additionally requires |LFC| > 0.2.
- **Grid:** 12 donor counts (5-600) × 14 cells-per-donor levels (5-3,000 plus all available cells), 20 iterations per grid point, 10 robustness runs. The paper's primary grids pool detection outcomes across all three methods (`all_methods/`); per-method grids are in the supplement.
- **Power model:** Gradient-boosted decision-tree ensemble (XGBoost), trained on the design rows with ≥ 100 cells per donor. It predicts each gene's detection probability from |LFC|, mean expression, dispersion, donor count, cells per donor, UMIs per cell, DE method and cell type. A dependency-light `predict_power()` wrapper, the fitted model (`power_model_bundle.rds`) and an example power curve are shipped in `Datalink/Power_Curve_Fitting/Exports/xgboost/`, so the wrapper runs without refitting. See [`docs/model-card-predict_power.md`](docs/model-card-predict_power.md).

## Key findings

Effect-size bins: tiny (|LFC| 0.2-1), small (1-2), medium (2-3), large (> 3).

1. **Genes with |LFC| < 2 did not reach 80% pooled recall anywhere in the design space**, even with 600 donors and all available cells. Tiny-effect genes needed 75 donors to reach 10% recall, reached 21% at 200 donors and never exceeded 42%.
2. **Large-effect genes reached 80% recall with as few as 15 donors** when at least 1,500 cells per donor were available.
3. **Only the most significant ground-truth genes were robustly recovered.** Recall falls sharply from the top to the bottom adjusted-p-value quartile, although every truth gene passes padj < 0.05, and mean expression strongly determines recall.
4. **Recall ranks microglia > astrocytes > oligodendrocytes**, anticorrelated with the amount of data per cell type. Reducing astrocytes to microglial cells per donor and UMIs per cell lost roughly 55% of ground-truth DEGs for DESeq2 and edgeR. Reducing cell numbers cost more DEGs than matching UMI depth, so maximising cells per donor matters more than deeper sequencing.
5. **The power model is reliable at the extremes but not in the contested middle.** 97.6% and 98.0% of predictions in the 0-0.1 and 0.9-1 bands fall within 0.10 of the observed detection rate, while the 0.2-0.8 band has a mean absolute error of 0.22-0.24. Only 7.8% of grid points fall between 0.1 and 0.9, so a skewed three-class (underpowered / contested / powered) output is recommended.

## Reproducing the analysis

`pipeline/` rebuilds every figure, table, quoted number and supplementary data
file from the raw PsychAD files, then knits the manuscript and the supplement.
One command submits the whole chain to a PBS Pro cluster (written for Imperial
CX3); `--local` runs the same stages one after another on a single large
machine. Everything runs inside the pinned container.

```bash
WORK=/path/to/empty/dir
mkdir -p $WORK/raw && cp PsychAD-{HBCC,MSSM,RADC}.h5ad media-1.csv $WORK/raw/
pipeline/reference/fetch_reference.sh $WORK/reference      # offline annotation bundle (needs internet once)
cp empirical-power-analysis.sif $WORK/epa.sif               # see docs/container.md
pipeline/run.sh $WORK                                       # check inputs, freeze code, submit every stage
pipeline/run.sh $WORK --status
pipeline/unpack_results.sh $WORK/results /path/to/checkout  # afterwards, wherever you want the Datalink/ tree
```

- **Full guide:** [`pipeline/README.md`](pipeline/README.md) (stages, resources,
  resuming, what differs from the preprint figures).
- **Container:** [`docs/container.md`](docs/container.md) (what is pinned, the
  smoke test, building or pulling the `.sif`).
- **Agents and new contributors:** [`AGENTS.md`](AGENTS.md).

### Inputs

Every input is public. The raw data are the three PsychAD h5ad files
(`PsychAD-HBCC.h5ad`, `PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad`) from
[CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f)
and the PsychAD donor table with post-mortem interval (`media-1.csv`, a
supplementary table of Fullard et al. 2025). The first job checks them against
`pipeline/inputs.sha256`. Annotation resources (biomaRt caches, HGNC, PANTHER,
the Lambert transcription factor catalogue, gnomAD pLI, the citation style) are
fetched once into a reference bundle and checked against
`pipeline/reference/MANIFEST.sha256` before every run, so compute nodes need no
internet access and annotation versions cannot drift. To use newer versions,
regenerate the bundle and its manifest deliberately.

We do not deposit copies of the raw data or of intermediates (donor shards, DE
checkpoints, the 10-run grid): the pipeline regenerates all of them from the
public inputs.

### Compute

The design grid dominates: 80 array tasks of roughly 6.5 hours on 24 cores
each, plus the AD vs control grid and a 108-task hyperparameter search. Stages
request up to 24 cores and 200 GB of memory, and peak node-local scratch is
about 330 GB. A full run takes 3-4 days on CX3. With `--local` the same stages
run sequentially on one large workstation, much more slowly.

### Reproducibility guarantees

- **Software:** base image by digest, apt packages by version, R packages by
  `renv.lock`, Python wheels by hash. The smoke test fails if any installed
  version differs from its pin.
- **Inputs:** raw files and the reference bundle are verified by SHA-256 before
  anything runs.
- **Code:** the first invocation freezes the commit with `git archive`; jobs
  run that snapshot, never the live checkout. Runs resume after interruption.
- **Randomness:** every stochastic step seeds from the run, design point and
  iteration, never from the job or chunk it lands in.
- **Outputs:** checksummed tars per stage, `MANIFEST.sha256` and
  `PROVENANCE.txt` (commit, image checksum, settings, input checksums, stage
  completion times).
- **Manuscript:** the `knit` stage renders both documents from the pipeline's
  outputs, with quoted numbers read from output tables via `_paper_values.R`.

Multithreaded XGBoost (8 threads) is deterministic on a given machine and
thread count; bitwise identity across different CPUs is not guaranteed.

## Citation

Please cite the accompanying preprint:

> Crampton C, Fawad S, Clark T, Dash H, Kövér B, Cotton G, Lee D, Duff E, Bottolo L, Matthews PM, Skene NG. *Single-Cell Study Designs Are Systematically Underpowered for Small-Effect Genes.* bioRxiv (2026). doi: [10.64898/2026.09.17.752063](https://doi.org/10.64898/2026.09.17.752063)

Machine-readable citation metadata is in `CITATION.cff`. The preprint DOI is the persistent identifier for this work; we do not mint a separate archive DOI for the code or data (see [Data prerequisites](#data-prerequisites)).

## License

Code is released under the MIT licence (see `LICENSE`).
Figures and source data derive from the publicly available PsychAD dataset (Fullard et al., 2025). Please refer to the original study for dataset-specific licensing terms.

## Contact

Questions and issues are welcome through the GitHub issue tracker.
