# EPA HPC scripts (Imperial CX3)

> **Reference only (superseded 2026-10-01).** These hand-driven PBS scripts were
> used for the earlier HPC runs. To reproduce the paper, use the one-command
> pipeline in [`pipeline/README.md`](../pipeline/README.md), which implements the
> same tar-staging design for every stage.

> **REBUILT 2026-06-29.** The original `hpc/` scripts were lost (deleted on
> the analysis server during a move and were never committed to git, so they
> were unrecoverable). These are a faithful reconstruction of the
> agreed tar-staging design, not byte-recovered originals. Run a `MINI_RUN` smoke
> test before trusting a full submit, and confirm the PBS resource directives
> (`ncpus`/`mem`) against your allocation.

## Golden rule

Tar everything, move big, compute on node-local `$TMPDIR`. RDS only ever sees
big files. Never detar onto RDS, never read loose shards off RDS, never write
microfiles back.

```
RDS  ──(shard tars + 1 control tar)──►  compute node $TMPDIR  ──(compute)──►  1 result tar per (mode,chunk) ──► RDS incoming/
                                                                                                                     │
                                                                          login node / iclneuro  ◄──(1 bulk unpack)──┘
```

## Files

| Script | Role |
|--------|------|
| `cx3-gt-manifest.pbs` | Step 1 (ground truth) + Step 1.5 (per-mode manifests). Stages shards to local scratch, emits **one** tag-scoped `HPC_Control_Tar/epa-control-<tag>.tar`. |
| `cx3-grid.pbs` | Step 2 grid array. Each task stages shard tars + the control tar to local scratch, computes one `(mode, chunk)`, writes **one** result tar to `HPC_Result_Tars/incoming/<tag>/`. |
| `cx3-xscore.pbs` | Cross-method 3×3 truth scoring (limma-voom on DESeq2 truth, edgeR on limma truth, ...). Stages `HPC_XGB_Input/cartesian_runs01-05.tar`, the runs 6-10 chunk tars, and the matching control tar into `$TMPDIR`, runs `R/cartesian_downsampling/10_CrossScore_Methods.R`, writes one small result tar (`Cartesian_Downsampling/Cross_Score/...`) into `HPC_Result_Tars/incoming/xscore-<tag>/`. |
| `harvest-grid-results.sh` | Unpacks result tars into the real Datalink in **one** bulk `rsync`. Optional tag arg harvests just one grid. Run on a login node or iclneuro, never a compute node. |

## Workflow

```bash
cd /rds/general/user/cjc20/home/EPA

# 1. ground truth + manifests (one job). Rebuild whenever EPA_N_CHUNKS changes.
qsub -v EPA_N_CHUNKS=40,EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink \
     hpc/cx3-gt-manifest.pbs

# 2. grid array: split mode (default): reads first (1..40), then cells (41..80).
#    -J upper bound MUST equal 2 * EPA_N_CHUNKS in split mode.
qsub -J 1-80%20 \
     -v EPA_N_CHUNKS=40,EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink \
     hpc/cx3-grid.pbs

# 3. harvest whenever convenient (login node or iclneuro). No tag = all grids
#    (incl. legacy untagged tars); pass a tag to harvest just one grid.
bash hpc/harvest-grid-results.sh $PWD/Datalink
```

## Run range + grid tag (incremental runs & isolated grids)

Two independent levers let you extend the robustness runs incrementally and run
a wholly separate grid (e.g. dreamlet) with **zero clashes**:

- **Run range** (`EPA_RUN_START` / `EPA_RUN_END`, R-side `RUN_START`/`RUN_END`)
  selects which robustness runs this invocation builds. Seeds are deterministic
  (`SEED_BASE + run_id`), so computing runs 6–10 leaves the banked 1–5
  byte-identical. `EPA_N_RUNS` is the project total (R validates `END ≤ total`).
- **Grid tag** (`EPA_GRID_TAG`, default `<methods>_runs<START>-<END>`)
  namespaces the control tar, result tars and done-markers:
  `epa-control-<tag>.tar`, `incoming/<tag>/`, `done/<tag>/`. A different method
  set (e.g. `EPA_DE_METHODS=dreamlet`) gets a different default tag, so two
  grids can even run concurrently.

```bash
# Extend an already-banked 1–5 grid with runs 6–10 (deseq2/edger/limma):
qsub -v EPA_RUN_START=6,EPA_RUN_END=10,EPA_N_RUNS=10,EPA_N_CHUNKS=40,\
EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink hpc/cx3-gt-manifest.pbs
qsub -J 1-80%20 -v EPA_RUN_START=6,EPA_RUN_END=10,EPA_N_RUNS=10,EPA_N_CHUNKS=40,\
EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink hpc/cx3-grid.pbs

# A separate dreamlet-only grid (different tag, no clash with deseq results):
qsub -v EPA_DE_METHODS=dreamlet,EPA_N_CHUNKS=40,\
EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink hpc/cx3-gt-manifest.pbs
qsub -J 1-80%20 -v EPA_DE_METHODS=dreamlet,EPA_N_CHUNKS=40,\
EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink hpc/cx3-grid.pbs

# harvest one grid by tag:
bash hpc/harvest-grid-results.sh $PWD/Datalink deseq2-edger-limma_runs06-10
```

> Always run `cx3-gt-manifest.pbs` with the **same** `EPA_RUN_START/END` +
> `EPA_DE_METHODS` as the grid: both derive the same tag, so the grid finds the
> matching `epa-control-<tag>.tar`. Mismatched tags fall back to reading loose
> control files off RDS (a warning), or fail to find the manifest.

## Cross-method truth scoring (3×3)

`cx3-xscore.pbs` runs `R/cartesian_downsampling/10_CrossScore_Methods.R` over the
completed grid to answer "does limma-voom's subsamples recover DESeq2's truth
genes, and vice-versa?" The outputs are cartesian recall/precision grids:
`Datalink/Cartesian_Downsampling/Cross_Score/cross_score_matrix.{qs,csv}`
(compose the figures locally with `Rscript R/figure_composition/compose_crossscore.R`).

The R script needs the raw `Cartesian_Downsampling/{ct}/Roussos/{mode}/{method}/Run_{run}/`
tree. Runs 1-5 were already harvested and repacked as
`HPC_XGB_Input/cartesian_runs01-05.tar`; runs 6-10 are still staged from the
tagged per-chunk result tars plus the matching control tar. Everything is copied
to `$TMPDIR` first and untarred locally.

```bash
# Default: scores runs 1-10 for deseq2+edger+limma across all 3 cell types x
# 2 modes. Stages runs 1-5 from HPC_XGB_Input/cartesian_runs01-05.tar, plus
# runs 6-10 from HPC_Result_Tars/incoming/deseq2-edger-limma_runs06-10/ and
# HPC_Control_Tar/epa-control-deseq2-edger-limma_runs06-10.tar:
qsub -v EPA_SIF=$PWD/epa-cx3-2026-07-19.sif,EPA_DATALINK=$PWD/Datalink \
     hpc/cx3-xscore.pbs

# Harvest the small result tar (login node or iclneuro):
bash hpc/harvest-grid-results.sh $PWD/Datalink xscore-deseq2-edger-limma_runs01-10
```

Scoping env vars (all optional; defaults reproduce the full grid):

| Var | Default | Effect |
|-----|---------|--------|
| `EPA_XSCORE_METHODS` | `deseq2,edger,limma` | R-side `XSCORE_METHODS` (rows AND cols of the 3×3). |
| `EPA_XSCORE_CELLTYPES` | `Astro,Micro,Oligo` | R-side `XSCORE_CELLTYPES`. |
| `EPA_XSCORE_MODES` | `cells,reads` | R-side `XSCORE_MODES`. |
| `EPA_XSCORE_RUNS` | `1,2,3,4,5,6,7,8,9,10` | R-side `XSCORE_RUNS`. Missing runs are skipped silently. |
| `EPA_XSCORE_MAX_ITERS` | *(unset)* | R-side `XSCORE_MAX_ITERS`: per-grid-point iteration cap for cheap probes. |
| `EPA_XSCORE_RUNS15_TAR` | `Datalink/HPC_XGB_Input/cartesian_runs01-05.tar` | Consolidated runs 1-5 cartesian input tar. |
| `EPA_XSCORE_RUNS610_TAG` | `deseq2-edger-limma_runs06-10` | Tag for the runs 6-10 chunk tars and control tar. |
| `EPA_XSCORE_RUNS610_CONTROL_TAR` | `Datalink/HPC_Control_Tar/epa-control-<runs610_tag>.tar` | Matching control tar for runs 6-10 truth files. |
| `EPA_XSCORE_TAG` | `xscore-<methods>_runs<lo>-<hi>` | Result tag; the output tar lands at `HPC_Result_Tars/incoming/<tag>/xscore.tar`. |

Read-only guarantee: the script never writes into the per-run grid tree under
`Cartesian_Downsampling/{ct}/Roussos/{mode}/{method}/Run_*/`; only
`Cartesian_Downsampling/Cross_Score/` is new.

## Split mode & resume

In `EPA_CARTESIAN_MODE=split` (default) a single 80-task array runs **reads
first, cells second**:

| array index | mode | chunk |
|-------------|------|-------|
| `1 .. N`     | reads | index |
| `N+1 .. 2N`  | cells | index − N  (N = `EPA_N_CHUNKS`) |

Each task runs **one mode** (≈ one wall's worth, not two), so at 40 chunks a
~50 h/mode run fits inside the 72 h wall with ~22 h margin. `chunk_archive_exists`
keys on `(mode, chunk)` via `HPC_Result_Tars/done/<tag>/<mode>/chunk_<n>.done`,
so already-banked work is skipped no matter which array index it maps to.

**Resume after timeouts:** just resubmit the same `-J 1-80%20`. Finished
`(mode, chunk)` pairs skip instantly; only missing ones run.

**Sizing if chunks still time out:** raise `EPA_N_CHUNKS` (and rebuild the
manifest with `cx3-gt-manifest.pbs`, since the chunk count is in the manifest
filename and changing it re-partitions the grid). More chunks = smaller units =
more margin, but does not shorten total wall-clock (that is capped by your CPU
allocation, not the chunk count). Note: changing the chunk count orphans
result tars from the old partition.

## Concurrency / queue

`%20` × `ncpus=24` = 480 CPUs, under the 500-CPU cap. Queue `v1_medium72a`
(72 h wall cap). Do not raise `%20`.
