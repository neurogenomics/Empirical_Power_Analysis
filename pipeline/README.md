# One-command reproduction on an HPC cluster

`pipeline/` regenerates every figure, table, number and supplementary data file
in the paper from the raw PsychAD data, and then knits the manuscript and the
supplement. One command submits the whole chain to a PBS Pro cluster (written
for Imperial CX3) and you can walk away. A few days later you get a set of
result tars that unpack into the `Datalink/` tree the manuscript reads.

Everything runs in one container image, and every random step is seeded. The
code is frozen at the commit you start from. Every intermediate is checksummed,
and the run records exactly what it used (`results/PROVENANCE.txt`).

## In short

```bash
# on a CX3 login node, inside a clone of this repository
WORK=$EPHEMERAL/epa-repro              # any empty directory on the cluster filesystem
mkdir -p $WORK/raw
cp /path/to/PsychAD-{HBCC,MSSM,RADC}.h5ad /path/to/media-1.csv $WORK/raw/
pipeline/reference/fetch_reference.sh $WORK/reference       # needs internet; login node is fine
cp /path/to/empirical-power-analysis.sif $WORK/epa.sif       # see docs/container.md

pipeline/run.sh $WORK          # check inputs, freeze code, submit all jobs, return
pipeline/run.sh $WORK --status # progress, any time
```

When the last stage (`finalise`) completes, `$WORK/results/` holds one directory
of tars per stage, plus `MANIFEST.sha256` and `PROVENANCE.txt`. Copy that
directory wherever you want (e.g. with Globus) and unpack it there:

```bash
pipeline/unpack_results.sh $WORK/results /path/to/checkout   # -> /path/to/checkout/Datalink/
```

Unpacking is deliberately not part of the pipeline. The cluster filesystem only
ever sees a few hundred large files, never the ~100k small ones inside them.

## Inputs

| Path in `WORK/` | What | Notes |
|---|---|---|
| `raw/PsychAD-HBCC.h5ad`, `raw/PsychAD-MSSM.h5ad`, `raw/PsychAD-RADC.h5ad` | PsychAD snRNA-seq ([CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f), Fullard et al. 2025), 57 GB | checked against `pipeline/inputs.sha256` by the first job |
| `raw/media-1.csv` | PsychAD donor table with post-mortem interval (supplementary table of Fullard et al.) | same |
| `reference/` | Offline annotation bundle (biomaRt caches, HGNC, PANTHER 19.0, Lambert TFs, gnomAD pLI, citation style) | built by `pipeline/reference/fetch_reference.sh`, checked against `pipeline/reference/MANIFEST.sha256` on every invocation |
| `epa.sif` | Apptainer image (or set `EPA_SIF=` in `WORK/epa.env`) | its sha256 is recorded; pin it with `EPA_SIF_SHA256=` |
| `epa.env` | Run settings | created from `pipeline/epa.env.defaults` on first use; the defaults reproduce the paper |

Compute nodes have no internet access, so nothing downloads anything at run
time. The raw directory is mounted read-only.

## Container image

The image holds software only, never data. Build it from this repository, or
pull the published image by digest once it exists; either way you end up with
a `.sif` and its checksum:

```bash
./docker/docker-build.sh       # Docker, ~1 h cold; ends with the smoke test
./docker/build-sif.sh --out empirical-power-analysis.sif   # converts, then smoke-tests the .sif
```

Both default to `empirical-power-analysis:v1.0.0`. The base image, apt
packages, R packages (`renv.lock`) and Python wheels are all pinned, and the
smoke test fails if any installed version differs from its pin. Routes without
Docker on the cluster, pinning a run to the image checksum, and how to change
the software are in [`docs/container.md`](../docs/container.md).

## How it works

- **Stages.** `pipeline/stages.tsv` defines 29 stages, their dependencies and
  their resources. Each stage is one PBS job (three of them, `grid`, `xgb_search`
  and `ad_grid`, are job arrays) and
  is implemented in `pipeline/stages/<stage>.sh`, which calls the analysis code
  in `R/`.
- **One command.** `run.sh` submits every stage up front, chained with
  `qsub -W depend=afterok:<id>`. If a stage fails, PBS cancels everything that
  depends on it.
- **No `qsub -v`.** All settings are baked into the generated job scripts
  (`WORK/jobs/<stage>.pbs`). Values with commas, such as
  `deseq2,edger,limma`, reach the container through an env file.
- **Frozen code.** The first invocation freezes the code with `git archive` of
  `HEAD` into `WORK/code/code.tar`. Every job runs that snapshot, so editing
  the checkout mid-run changes nothing. `run.sh` refuses to start from a
  repository with uncommitted changes unless you pass `--allow-dirty`.

**Anatomy of a stage.** Every job runs the same five steps, using the
functions in `pipeline/lib/common.sh`:
1. Unpack the code snapshot into node-local scratch (`$TMPDIR`).
2. Unpack the upstream result tars it needs, verifying each checksum.
3. Run R or Python inside the container, bound to that local copy.
4. Tar the declared outputs back to `WORK/results/<stage>/`, written
   atomically with a `.sha256`.
5. Write `results/<stage>/.complete`.

A declared output that is missing fails the job. This matters because several
R runners report success even when a step skipped work.

**Resuming.** Run the same command again. It does the following:
- **Complete stages** (`.complete` present, or every array task's tar present)
  are skipped and dropped from dependency lists.
- **Stages still queued or running** are reused as dependencies, not
  resubmitted.
- **Failed array stages** are resubmitted as an array of only the missing
  tasks; a single missing task goes as a plain job.

All of this is decided on the login node, so a large allocation is never
requested just to find there is nothing left to do. To redo a finished stage,
delete its `results/<stage>/` directory. Downstream stages are not
invalidated automatically, so delete those too.

**Grid chunking.** Grid chunking doesn't affect results: a design point's seed
depends only on (run, donors, mode, depth, iteration). `EPA_GRID_CHUNKS` and
`EPA_AD_GRID_CHUNKS` only trade task length against task count.

**Dates in outputs.** Jobs set `SOURCE_DATE_EPOCH` to the commit time of the
code snapshot, so files that record an export date (the `predict_power` bundle,
its README) are byte-identical across reruns of the same commit.

## Stages and what they produce

| Stage | Depends on | Produces (paper asset) | Rough time |
|---|---|---|---|
| `shard` | - | donor shards, filtered cohort (age >= 40, >= 30 cells) and full cohort, PMI annotated | 4 h |
| `eda` | shard | PCA gallery (Supp. Fig.), sex-label validation (Supp. Table 4), cohort statistics (donor counts; Supp. Data 1) | < 1 h |
| `fgt` | shard | full-cohort sex DE, gene-universe audit (Supp. Data 1), DEG chromosome grid (Supp. Fig.) | 1-2 h |
| `deg_loss` | shard | depth / cell-count equalisation towards microglia, astrocytes and oligodendrocytes (Supp. Table 3) | 4-5 h |
| `gt_manifest` | shard | ground truth for 10 runs, grid job manifests | 4-8 h |
| `grid` (array, 2 x 40) | shard, gt_manifest | per-iteration DE results for the whole design grid | ~6.5 h per task |
| `cartesian_cells`, `cartesian_reads` | grid | aggregated recall/precision, all grid and surface figures | 2-6 h each |
| `crossscore` | grid | cross-method scoring | ~10 h |
| `pcf_aggregate` | shard, grid | power-model training tables; truth-half DESeq2 checkpoints | ~16 h |
| `xgb_search` (array 108) | pcf_aggregate | hyperparameter search of the power model (n_cells >= 100 rows): 4 depths x 3 learning rates x 3 min_child_weights x 3 reweighting exponents, early stopping up to 20,000 rounds | 0.5-6 h per task |
| `xgb_select` | xgb_search | the power model, config ranking with per-fold MAE (Supp. Data 4), CV metrics quoted in text, 5 SHAP refits | ~6 h |
| `xgb_shap` | xgb_select | SHAP figures (Supp. Figs.) | 2-4 h |
| `xgb_post` | xgb_select | Figure 5 (incl. hyperparameter bracketing, per-class precision/recall/F1), exported `predict_power` bundle | ~2 h |
| `gene_class` | shard, xgb_post | Figure 6 | ~2 h |
| `tables` | fgt, pcf_aggregate, cartesian_* | Supp. Table 2, recall grid CSV (numbers in text) | < 1 h |
| `compose` | cartesian_*, crossscore | Figures 2-4, cross-method and per-method Supp. Figs. | < 1 h |
| `supp_data` | several | Supplementary Data files (SD1 and SD4 are cited in the manuscript) | < 1 h |
| `paper_values` | eda, tables, deg_loss, pcf_aggregate, xgb_*, gene_class | `Paper_Source_Data/paper_values.csv`: every data-derived number quoted in the text, with its definition and source file | minutes |
| `knit` | all of the above | `Manuscript/nature_comms_{draft,supplementary}.html` | < 1 h |
| `ad_shard` | shard | AD vs control donor shards: AD as the only diagnosis vs no diagnosis | < 1 h |
| `ad_gt_manifest`, `ad_grid` (array, `EPA_AD_GRID_CHUNKS`, cells mode) | ad_shard | the AD vs control grid, 10 runs x 3 methods, same seeds | ~5 h per task |
| `ad_cartesian` | ad_grid | AD recall/precision grids and surfaces (method-pooled), p-value-quartile and expression panels | 2-4 h |
| `ad_de` | shard, ad_shard | full-cohort DE for both contrasts; AD vs autosomal sex-biased DEG |LFC| and expression distributions | 1-2 h |
| `finalise` | knit, AD stages | `MANIFEST.sha256`, `PROVENANCE.txt` | minutes |

**Critical path.** `shard`, then `gt_manifest`, `grid`, `pcf_aggregate`,
`xgb_search` and `xgb_select`, `xgb_post`, `gene_class` and finally `knit`.
That is about 3-4 days of compute plus queueing. With 20 grid tasks at a time,
the grid itself takes about a day. The AD vs control stages run alongside.

**Space on the cluster filesystem.**
- Shards: 11 GB.
- Grid tars: ~310 GB (sex) plus ~100 GB (AD vs control).
- All other results: about 10 GB.
- Peak node-local scratch is ~330 GB, in the stages that unpack the whole grid
  (`crossscore`, `pcf_aggregate`).

## What differs from the figures in the preprint

These differences are intentional, decided on 2026-10-01 and 2026-10-06, and
listed so that nobody mistakes them for reproduction failures:
- **All figures use all 10 robustness runs.** Before this pipeline, the
  cross-method figures and some cells-mode surfaces had been redrawn from runs
  1-5, while the reads-mode figures used runs 1-10.
- **There is one power model.** The former "baseline" (all grid rows) model is
  gone; the model is trained on the n_cells >= 100 rows and no longer called
  "realistic". Its outputs live in `Models/xgboost/` and `Exports/xgboost/`.
- **The hyperparameter grid is wider** (108 configurations instead of 36): one
  more value beyond the best-performing end of each tree hyperparameter
  (depth 10, learning rate 0.01, min_child_weight 25), so the search brackets
  its optimum. Figure 5 panel B shows that bracketing (best MAE per value, with
  fold-to-fold standard errors) instead of the 36-bar chart.
- **Figure 5 panel C is out-of-fold** and shows precision, recall and F1 of the
  model per class (it previously compared in-sample F1 with the baseline).
- **Figure 5 panels A and D** use the current panel code (10-90% bands, renamed
  features), without the baseline facet.
- **Supplementary Table 3** has oligodendrocyte rows as well as astrocyte rows.
- **The AD vs control analysis is new** (`Datalink/AD_vs_Control/`); where it
  goes in the manuscript is not decided.
- **Supplementary Table 2** pools all 10 truth halves. The earlier version used
  the 6 that existed in April.
- **Figure 6 scripts were reconstructed**, because the originals were lost. The
  annotation and prediction tables match the published ones exactly, except two
  non-coding genes that HGNC has since renamed. The plot matches in every
  number; jittered point positions differ.
- **Scope is PsychAD only.** The other cohorts and the GTEx comparison are no
  longer in the paper and are not run.

## Testing

- `bash pipeline/tests/test_pipeline.sh` tests the orchestration layer in
  seconds, with `qsub`, `qstat` and the container stubbed. It covers dependency
  edges, array sizes, resume and partial-array resubmission, reuse of queued
  jobs, preflight refusals, the stage-in/stage-out helpers, and a local-mode
  job.
- `pipeline/run.sh WORK --dry-run` writes every job script and prints the qsub
  chain without submitting.
- `pipeline/run.sh WORK --local` runs the stages one after another on the
  current host (needs apptainer and a large machine). It is meant for testing
  individual stages with `--only`.
