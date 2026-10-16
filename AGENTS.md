# AGENTS.md

Guide for coding agents (and humans) who need to run, inspect, critique or
change this repository. Read this first, then `pipeline/README.md`.

## What this repository is

The code behind *Single-Cell Study Designs Are Systematically Underpowered for
Small-Effect Genes* (bioRxiv, doi 10.64898/2026.09.17.752063). It measures
detection power for pseudobulk differential expression empirically: the
PsychAD cohort is split into sex-stratified halves ten times, one half gives
the ground truth, the other is downsampled across a grid of donors x cells per
donor (and donors x UMIs per cell), each pseudo-experiment is analysed with
DESeq2, edgeR and limma-voom, and recall/precision are scored against the
truth. An XGBoost model is then fitted to the outcomes.

Scope is PsychAD only (astrocytes, microglia, oligodendrocytes). "Roussos" in
paths and code is the lab's name for PsychAD.

## The one supported way to run it

```bash
pipeline/run.sh WORK            # check inputs, freeze code, submit every stage, return
pipeline/run.sh WORK --status   # progress
pipeline/run.sh WORK --dry-run  # write job scripts, submit nothing
pipeline/run.sh WORK --local --only eda   # run stages on this host (testing)
pipeline/unpack_results.sh WORK/results /path/to/checkout   # afterwards
```

`WORK` needs `raw/` (three PsychAD h5ad files and `media-1.csv`),
`reference/` (from `pipeline/reference/fetch_reference.sh`) and `epa.sif`
(see `docs/container.md`). Settings live in `WORK/epa.env`, copied from
`pipeline/epa.env.defaults`; the defaults reproduce the paper.

Do not run `R/<pipeline>/run_pipeline.R` scripts on a bare host to reproduce
anything: package versions would not be the pinned ones. For development, run
them inside the container (`docs/container.md`, "Running the image by hand").

## Map

| Path | What |
| --- | --- |
| `pipeline/run.sh`, `pipeline/lib/common.sh` | Orchestration: job generation, staging in and out, resume |
| `pipeline/stages.tsv` | The stage graph: name, dependencies, resources, array kind |
| `pipeline/stages/<stage>.sh` | One file per stage: which R/Python scripts it runs, which outputs it requires and tars |
| `pipeline/inputs.sha256`, `pipeline/reference/MANIFEST.sha256` | Checksums of every raw input and reference file |
| `R/<pipeline>/` | Analysis code. `config.R` holds constants; numbered scripts are steps |
| `R/shared/` | Helpers used across pipelines (`contrast.R` switches sex vs AD-vs-control) |
| `Dockerfile`, `docker/`, `renv.lock` | The pinned software environment (`docs/container.md`) |
| `nature_comms_draft.rmd`, `nature_comms_supplementary.rmd` | Manuscript sources, knitted by the `knit` stage |
| `_paper_values.R` | `val("key")`: numbers quoted in the text, read at knit time. Keys not registered there come from `Datalink/Paper_Source_Data/paper_values.csv` (stage `paper_values`, `R/figure_data_generation/09_Paper_Values.R`), which states each number's definition and source |
| `_supp_fig_index.R` | `SF("key")`: supplementary figure numbering for both documents |
| `Datalink/` | Output tree. Git tracks only the files the manuscript renders; a run regenerates all of it |

To find what produces a paper asset: look up the file in `README.md` ("Paper
asset lookup"), then grep `pipeline/stages/*.sh` for its path (`epa_require`
lines name each stage's guaranteed outputs).

## How a stage works

Each stage is one PBS job (or job array). It unpacks the frozen code snapshot
(`git archive` of the commit you launched from) and its upstream result tars
into node-local scratch, runs R/Python inside the container with
`apptainer exec --cleanenv --no-home`, checks its declared outputs exist, and
writes them as checksummed tars to `WORK/results/<stage>/`. Compute nodes have
no internet: anything external comes from the reference bundle.

## Invariants: do not break these

- **Seeds.** Every random step seeds from the run, design point and iteration
  (see `R/cartesian_downsampling/config.R` and `utils.R`). Never seed from
  time, PID, chunk number or array index. Grid chunking must not change results.
- **Truth/test separation.** `GT_SPLIT = 0.5`: no donor appears in both the
  ground-truth half and the downsampled half of a run.
- **Load-bearing thresholds.** `GT_LFC_THRESHOLD = 0.2`, FDR 0.05, the effect
  size bins (tiny 0.2-1, small 1-2, medium 2-3, large > 3), minimum gene counts,
  donor filters (age >= 40, >= 30 cells). Changing any of them changes the
  paper; the text must change with them.
- **Shrinkage.** All three DE methods shrink log fold changes with ashr, so
  effect sizes are comparable across methods.
- **Cell type order** is Astro, Micro, Oligo everywhere (figures depend on it).
- **Software.** Every R package comes from `renv.lock`, every Python package
  from `docker/python-requirements.txt` with hashes, the base image by digest.
  To add or upgrade a package, follow "Changing the software" in
  `docs/container.md`. Never install packages at run time.
- **No hard-coded paths.** Code resolves the project root with `here::here()`
  or `POWER_ANALYSIS_ROOT`. `git grep -nE '/mnt/|/home/|/rds/'` over tracked
  files finds only `/home/rstudio` (a path inside the base image) and this line.
- **Numbers in the text** come from `val()` keys, not literals, wherever a key
  exists. To add one, add an `add()` call to `09_Paper_Values.R` with a full
  definition, rather than typing a number.

## Testing

```bash
bash pipeline/tests/test_pipeline.sh       # orchestration tests, seconds, no container or data
pipeline/run.sh WORK --dry-run             # inspect generated job scripts
./docker/docker-smoke.sh                   # environment check of a built image
```

Six self-contained unit tests need no data and run inside the container from
the repository root:

```bash
for t in R/cartesian_downsampling/tests/test_age_parsing.R \
         R/cartesian_downsampling/tests/test_job_planning.R \
         R/cartesian_downsampling/tests/test_reads_sampling.R \
         R/cartesian_downsampling/tests/test_ad_sex_balance.R \
         R/de_diagnostics/toy/run_toy_test.R \
         R/eda_roussos/tests/test_all_utils.R; do
  apptainer exec --cleanenv --no-home --env HOME=/tmp -B "$PWD":/workspace --pwd /workspace epa.sif Rscript "$t"
done
```

Older module tests that need the full data on disk are not shipped.

## Things that are easy to get wrong

- `.gitignore` lists code that exists in the authors' working tree but is not
  part of the reproduction path (superseded runners, other cohorts, the
  archived logistic-regression model, exploratory analyses). It is not
  shipped; do not recreate it or wire it back in.
- Text inside `<!-- -->` or `eval=FALSE` chunks of the `.rmd` files is not part
  of the paper. Ignore it.
- `Datalink/` files tracked in git are a snapshot of the last run. After a new
  run, `unpack_results.sh` replaces them; differences between runs of the same
  commit and image would be a bug.
- Multithreaded XGBoost (8 threads) is deterministic for a given machine and
  thread count; bitwise identity across different CPUs is not guaranteed.
