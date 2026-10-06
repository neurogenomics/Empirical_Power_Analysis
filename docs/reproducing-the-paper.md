# Reproducing the paper

There are two routes. Both use the same analysis code in `R/` and the same pinned
environment (`renv.lock`, or the container built from `Dockerfile`).

1. **One command on a PBS Pro cluster (recommended).** `pipeline/run.sh` submits
   every stage, runs it in the pinned Apptainer image, checks inputs against
   checksums, and ends by knitting the manuscript and supplement. This is the
   route that reproduces every figure, table, quoted number and Supplementary
   Data file. See [`pipeline/README.md`](../pipeline/README.md).
2. **By hand, one pipeline at a time.** Useful for rerunning or inspecting a
   single analysis. The table below lists, for each stage of the cluster
   pipeline, the commands it runs. The full grid (`grid`) is impractical without
   a cluster.

## Inputs

- Raw data: `PsychAD-HBCC.h5ad`, `PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad`
  ([CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f))
  and the PsychAD donor table with post-mortem interval (`media-1.csv`, from
  Fullard et al. 2025). Checksums: `pipeline/inputs.sha256`. For the manual
  route, place them under `Datalink/raw_data/`.
- Reference bundle (biomaRt caches, HGNC, PANTHER, Lambert TFs, gnomAD pLI,
  citation style): `pipeline/reference/fetch_reference.sh <dir>`.
- Environment: build the image with `./docker/docker-build.sh` (see
  [`docker.md`](docker.md)), or `renv::restore()` under R 4.5.1 / Bioconductor 3.22.

## Stage-by-stage commands

Run from the repository root. Each `pipeline/stages/<stage>.sh` is the
authoritative version: it also sets the environment variables noted here
(defaults in `pipeline/epa.env.defaults`, e.g. 10 robustness runs and the three
DE methods).

| Stage | Commands | Main outputs |
| --- | --- | --- |
| `shard` | `Rscript R/sharding/Run_Shard_Roussos.R Astro,Micro,Oligo` | `Datalink/Shards/` |
| `eda` | `Rscript R/eda_roussos/run_pipeline.R 1`; `Rscript R/eda_roussos/Validate_Sex_Labels.R All`; `Rscript R/figure_data_generation/run_pipeline.R 1` | PCA gallery, sex validation, dataset statistics |
| `fgt` | `Rscript R/full_ground_truth/run_pipeline.R 1,4`; `Rscript R/de_diagnostics/run_pipeline.R 4` | Full-cohort DE, gene-universe audit, chromosome grid |
| `deg_loss` | `Rscript R/depth_equalization/run_pipeline.R`; `Rscript R/depth_equalization/deg_loss_analysis.R` | `Datalink/Depth_Eq_DEG_Loss/` (Supp. Table 3) |
| `gt_manifest` | `R/cartesian_downsampling/run_pipeline.R` steps `1` and `1.5` (`DOWNSAMPLE_MODE=cells`, then `1.5` with `reads`) | Ground-truth splits and grid manifests |
| `grid` | `R/cartesian_downsampling/run_pipeline.R 2` per mode and `JOB_CHUNK` (array job) | Per-iteration DE results for the whole grid |
| `cartesian_cells`, `cartesian_reads` | `R/cartesian_downsampling/run_pipeline.R 3,3.5,4,5,6` with `DOWNSAMPLE_MODE` and `CARTESIAN_INCLUDE_ALL_METHODS=true` | Recall/precision grids and surfaces |
| `crossscore` | `Rscript R/cartesian_downsampling/10_CrossScore_Methods.R` | `Datalink/Cartesian_Downsampling/Cross_Score/` |
| `pcf_aggregate` | `R/power_curve_fitting/` `00`, `00b`, `01` (once per DE method), `01b`, `01e` | Aggregated fitting data (baseline and realistic) |
| `xgb_base_search`, `xgb_base_select` | `02_Fit_XGBoost.R` per `XGB_CONFIG_INDEX`, then `03_Select_XGBoost.R` | Baseline model (comparator) |
| `xgb_real_search`, `xgb_real_select`, `xgb_real_shap` | the same on the realistic (`n_cells >= 100`) data, then `04_SHAP_Figures.R` | Production model, config ranking, SHAP figures |
| `xgb_post` | `08b_Baseline_OOF.R`, `04_Export_XGBoost.R`, `07_LowPower_Diversity.R`, `08_Compare_Reweighted.R`, `09_Trinary_Classification.R` (`XGB_TRINARY_MODE=oof`), `_regenerate_xgb_panels.R`, `R/figure_composition/compose_xgb_panels.R` | Exported model bundle, Figure 5 |
| `gene_class` | `R/gene_classification/classify_oligo_genes.R`, `annotate_oligo_genes.R`, `classify_oligo_multi_lfc.R`, `plot_oligo_expr_disp_panels.py` | Figure 6 |
| `tables` | `Rscript R/figure_data_generation/run_pipeline.R 6,8`; `Rscript R/de_diagnostics/run_pipeline.R 5` | Paper source CSVs |
| `compose` | `R/figure_composition/compose_figures.py --only figure_cartesian_oligo`, `compose_recall_panels.R --cell Oligo --mode cells`, `compose_crossscore.R`, `compose_supplementary.py` | Figures 2-4, composed supplementary figures |
| `supp_data` | `Rscript R/supplementary_data/build_supplementary_data.R` | `Datalink/Supplementary_Data/` |
| `knit` | `rmarkdown::render()` on `nature_comms_draft.rmd` and `nature_comms_supplementary.rmd` | Manuscript and supplement (HTML) |

The order above is a valid topological order; `pipeline/stages.tsv` gives the
exact dependencies.

## Checking a reproduction

- A `pipeline/` run writes `results/MANIFEST.sha256` (every result tar) and
  `results/PROVENANCE.txt` (commit, container digest, input checksums,
  settings).
- Quoted numbers in the manuscript are read at knit time from the regenerated
  CSVs (`_paper_values.R`), so the knitted text reflects the rerun directly.
- Every random step (splits, donor and cell subsampling, model folds) is
  seeded deterministically from the run, grid point and iteration.

## Notes

- "Roussos" in paths and code is the lab's legacy name for the PsychAD cohort.
- `R/bulk_deg_processing/` and `R/glial_characteristics/` belong to an earlier
  multi-cohort version of the study and are not part of either route.
- `run_all.sh` and `run_power_curve_fitting.sh` are earlier single-host drivers
  that run pipelines serially. They are fine for rerunning one pipeline, but the
  full grid and the XGBoost search are impractical on one host; `pipeline/` is
  the tested route for the whole paper.
