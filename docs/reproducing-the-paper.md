# Reproducing the paper

The supported way to reproduce the paper is `pipeline/run.sh`. It runs every
stage in the pinned Apptainer image, checks inputs against checksums, and ends
by knitting the manuscript and supplement. This reproduces every figure, table,
quoted number and Supplementary Data file. It submits the stages to a PBS Pro
cluster, or, with `--local`, runs them one after another on a single host (this
needs a large machine, and the full grid is impractical without a cluster). See
[`pipeline/README.md`](../pipeline/README.md).

The table below lists, for each stage, the commands it runs. Use it as a map
for reading or debugging the code. To rerun one stage, use
`pipeline/run.sh WORK --only <stage>` rather than calling the scripts directly.

## Inputs

- Raw data: `PsychAD-HBCC.h5ad`, `PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad`
  ([CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f))
  and the PsychAD donor table with post-mortem interval (`media-1.csv`, from
  Fullard et al. 2025). Checksums: `pipeline/inputs.sha256`. Place them under
  `WORK/raw/`; each job mounts that directory read-only at `Datalink/raw_data/`
  inside the container. Every input is publicly available, so no data deposit
  is needed to reproduce the paper.
- Reference bundle (biomaRt caches, HGNC, PANTHER, Lambert TFs, gnomAD pLI,
  citation style): `pipeline/reference/fetch_reference.sh WORK/reference`.
- Environment: the Apptainer image built from `Dockerfile` with
  `docker/docker-build.sh`, checked with `docker/docker-smoke.sh` and converted
  with `docker/build-sif.sh` (tag `empirical-power-analysis:v1.0.0`). Once
  published, it can be pulled from `ghcr.io/neurogenomics/empirical_power_analysis`
  instead. See [`container.md`](container.md).

## Stage-by-stage commands

Paths are relative to the repository root. Each `pipeline/stages/<stage>.sh` is
the authoritative version: it also sets the environment variables noted here
(defaults in `pipeline/epa.env.defaults`, e.g. 10 robustness runs and the three
DE methods).

| Stage | Commands | Main outputs |
| --- | --- | --- |
| `shard` | `Rscript R/sharding/Run_Shard_Roussos.R Astro,Micro,Oligo` | `Datalink/Shards/` |
| `eda` | `Rscript R/eda_roussos/run_pipeline.R 1`; `Rscript R/eda_roussos/Validate_Sex_Labels.R All`; `Rscript R/figure_data_generation/run_pipeline.R 1` | PCA gallery, sex validation, dataset statistics |
| `fgt` | `Rscript R/full_ground_truth/run_pipeline.R 1,4`; `Rscript R/de_diagnostics/run_pipeline.R 4` | Full-cohort DE, gene-universe audit, chromosome grid |
| `deg_loss` | `Rscript R/depth_equalization/run_pipeline.R` and `Rscript R/depth_equalization/deg_loss_analysis.R`, once with `DEPTH_EQ_CELL_TYPE=Astro` and once with `Oligo`; then `Rscript R/depth_equalization/combine_deg_loss_tables.R` | `Datalink/Depth_Eq_DEG_Loss/` (Supp. Table 3) |
| `gt_manifest` | `R/cartesian_downsampling/run_pipeline.R` steps `1` and `1.5` (`DOWNSAMPLE_MODE=cells`, then `1.5` with `reads`) | Ground-truth splits and grid manifests |
| `grid` | `R/cartesian_downsampling/run_pipeline.R 2` per mode and `JOB_CHUNK` (array job) | Per-iteration DE results for the whole grid |
| `cartesian_cells`, `cartesian_reads` | `R/cartesian_downsampling/run_pipeline.R 3,3.5,4,5,6` with `DOWNSAMPLE_MODE` and `CARTESIAN_INCLUDE_ALL_METHODS=true` | Recall/precision grids and surfaces |
| `crossscore` | `Rscript R/cartesian_downsampling/10_CrossScore_Methods.R` | `Datalink/Cartesian_Downsampling/Cross_Score/` |
| `pcf_aggregate` | `R/power_curve_fitting/` `00`, `00b`, `01` (once per DE method), `01b`, `01e_Filter_Design.R` | Power-model training table (`aggregated_fitting_data.qs`, rows with `n_cells >= 100`) and truth-half DESeq2 checkpoints |
| `xgb_search` (array of 108) | `R/power_curve_fitting/02_Fit_XGBoost.R` with `XGB_CONFIG_INDEX` = 1 to 108 | One hyperparameter configuration per task |
| `xgb_select` | `R/power_curve_fitting/03_Select_XGBoost.R` | The power model, config ranking (`model_selection.csv`, per-fold `model_selection_folds.csv`), five SHAP refits |
| `xgb_shap` | `R/power_curve_fitting/04_SHAP_Figures.R` | SHAP figures |
| `xgb_post` | `04_Export_XGBoost.R`, `07_LowPower_Diversity.R`, `08_Compare_Reweighted.R`, `09_Trinary_Classification.R` (`XGB_TRINARY_MODE=oof`), `_regenerate_xgb_panels.R`, `R/figure_composition/compose_xgb_panels.R` | Exported model bundle, Figure 5 |
| `gene_class` | `R/gene_classification/classify_oligo_genes.R`, `annotate_oligo_genes.R`, `classify_oligo_multi_lfc.R`, `plot_oligo_expr_disp_panels.py` | Figure 6 |
| `tables` | `Rscript R/figure_data_generation/run_pipeline.R 6,8` | Paper source CSVs |
| `compose` | `R/figure_composition/compose_figures.py --only figure_cartesian_oligo`, `compose_recall_panels.R --cell Oligo --mode cells`, `compose_crossscore.R`, `compose_supplementary.py` | Figures 2-4, composed supplementary figures |
| `supp_data` | `Rscript R/supplementary_data/build_supplementary_data.R` | `Datalink/Supplementary_Data/` |
| `knit` | `rmarkdown::render()` on `nature_comms_draft.rmd` and `nature_comms_supplementary.rmd` | Manuscript and supplement (HTML) |
| `ad_shard` | `Rscript R/sharding/Build_AD_Control_Shards.R` | AD vs control donor shards (`Shards/*/Roussos_ADvsCtrl/`) |
| `ad_gt_manifest`, `ad_grid`, `ad_cartesian` | `R/cartesian_downsampling/run_pipeline.R` steps `1, 1.5`, then `2` (array job), then `3,3.5,4,5,6`, with `EPA_CONTRAST=ad_vs_control`, `SHARD_DATASET=Roussos_ADvsCtrl`, cells mode only; `ad_cartesian` also runs `compose_recall_panels.R` and `compose_precision_panels.R` | AD vs control recall/precision grids and figures |
| `ad_de` | `R/ad_contrast/01_Full_Cohort_DE.R` (sex, then AD vs control), `02_Compare_Effect_Distributions.R` | Full-cohort DE for both contrasts, effect distributions |
| `finalise` | (shell only) | `results/MANIFEST.sha256`, `results/PROVENANCE.txt` |

The order above is a valid topological order; `pipeline/stages.tsv` gives the
exact dependencies. The AD vs control stages write to `Datalink/AD_vs_Control/`;
where that analysis goes in the manuscript is not yet decided.

## Checking a reproduction

- A `pipeline/` run writes `results/MANIFEST.sha256` (every result tar) and
  `results/PROVENANCE.txt` (commit, container digest, input checksums,
  settings).
- Quoted numbers in the manuscript are read at knit time from the regenerated
  CSVs (`_paper_values.R`), so the knitted text reflects the rerun directly.
- Every random step (splits, donor and cell subsampling, model folds) is
  seeded deterministically from the run, grid point and iteration.
- Intentional differences from the figures in the preprint are listed in
  [`pipeline/README.md`](../pipeline/README.md#what-differs-from-the-figures-in-the-preprint).

## Notes

- "Roussos" in paths and code is the lab's legacy name for the PsychAD cohort.
- Code that `pipeline/run.sh` does not execute, and outputs the manuscript does
  not read, are not in the repository (they were untracked in October 2026 and
  remain in the git history). This includes the earlier single-host drivers,
  the other cohorts and the GTEx comparison, and the retired logistic-regression
  power model.
