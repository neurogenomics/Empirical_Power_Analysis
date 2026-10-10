# Architecture Decision Records (ADRs)

Each ADR captures a significant design choice, its context, and the reasoning behind
it, so future contributors understand why the code looks the way it does.

Format: [Michael Nygard's ADR template](https://github.com/joelparkerhenderson/architecture-decision-record).

---

## ADR-001: Use Nextflow for pipeline orchestration

**Status**: SUPERSEDED by ADR-001b (2026-04-28). Original decision retained
below for historical context.

**Context**: The Cartesian downsampling step is the hot loop, up to 120,000 DESeq2
runs per cell type × mode × robustness. This exceeds single-machine capacity when
extended (e.g. to neurons or a second cohort). We need a workflow manager that can
(a) dispatch to SLURM with per-step resource declarations, (b) resume from failure,
(c) track provenance, and (d) integrate with containers.

**Options considered**:
- **GNU Make**: rejected. Cannot submit SLURM jobs or declare per-process resources.
- **Snakemake**: viable. Python-based DAG, SLURM profile support. Rejected mainly
  for lab consistency: the Skene/Neurogenomics lab already maintains multiple
  Nextflow repos (`NextFlowIntro`, `nf-imperial-canary`, `wdl-scflow`,
  `scRNAseq-pipeline`, `LigerFlow`). A second paradigm would fragment institutional
  knowledge.
- **R `targets`**: viable for R-native workflows. Rejected because (a) HPC dispatch
  via `future.batchtools::batchtools_slurm` is finicky and (b) lab is Nextflow-first.
- **Nextflow DSL2**: accepted. Native SLURM, container integration, matches lab
  conventions, reviewer-credible at the journal tier we target.

**Consequences**: Maintaining both `R/<pipeline>/run_pipeline.R` (inner parallelism
via `future.apply`) and `modules/local/<name>/main.nf` (outer orchestration). These
must stay in sync when pipelines are added or renamed. CI stub-run enforces this.

---

## ADR-001b: Replace Nextflow with bash drivers + Docker image

**Status**: Accepted, 2026-04-28

**Context**: ADR-001 chose Nextflow for SLURM dispatch and lab consistency.
In practice we never deployed to SLURM for this paper, all production runs
were on a single multi-core box, with the long-pole stages (cartesian
downsampling and power-curve fitting) consuming whole-host resources for
days at a time. The dual-source-of-truth requirement (R `run_pipeline.R`
files plus `modules/local/<name>/main.nf` Nextflow processes) was a
maintenance burden that produced silent drift between the two, and the
Nextflow staging machinery added friction without providing value when the
DAG is run strictly serially on a single host. Reproducibility for external
users is better served by a Docker image than by asking them to install
Nextflow + Java.

**Options considered**:
- **Stay on Nextflow, fix the drift**: viable but does not address the
  external reproducibility friction.
- **Strict-serial bash driver (`run_all.sh`) + Docker image**: accepted.
  Single source of truth (the R pipelines), one shell script for full DAG,
  one for power-curve-only, and a set of small Docker scripts
  (`docker-build.sh` / `docker-smoke.sh` / `docker-run.sh`, chained by
  `docker-all.sh`) that build + smoke-test + run the pipeline inside a
  container. External reproducers run one command and walk away.

**Consequences**: `main.nf`, `nextflow.config`, `iclneuro.config`, `conf/`
and `modules/local/` are removed. The R pipelines and their `config.R`
honour `POWER_ANALYSIS_ROOT` and `POWER_ANALYSIS_FORCE_RECALC` environment
variables so the bash driver can flip behaviour without touching code. ADR-004
is also superseded, the paper is rendered via the bash driver / Docker
image, not via a `RENDER_PAPER` Nextflow process.

---

## ADR-002: Pseudobulk with three DE methods, all shrunk with ashr

**Status**: Accepted. Originally DESeq2 + apeglm (2024, inherited from Salman
Fawad's thesis); extended to three methods in June 2026; uniform ashr
shrinkage adopted 2026-09-04.

**Context**: The grid needs a DE method that scales to millions of runs, gives
calibrated p-values at small n, and produces effect sizes comparable across
methods. A single-method result invites the objection that the power estimates
are specific to that tool.

**Options considered**:
- **Mixed models on cells (NEBULA, glmmTMB)**: rejected. ~10x slower per run,
  and cell-level models are prone to pseudoreplication.
- **DESeq2 only, with apeglm shrinkage**: the original design. Rejected because
  results would be DESeq2-specific and apeglm is DESeq2-only, so effect sizes
  could not be compared with other methods.
- **DESeq2, edgeR and limma-voom on identical pseudobulk samples, all LFCs
  shrunk with ashr**: accepted. ashr works from any method's estimate and
  standard error, so one shrinkage estimator applies to all three.

**Consequences**: Every grid point is run three times. The paper's primary grids
pool outcomes across methods (`all_methods/`); per-method and cross-method
(truth-method × downsampling-method) grids are in the supplement. `de_method` is
a feature of the power model. The ashr settings live in
`R/cartesian_downsampling/utils.R` and `R/full_ground_truth/utils.R`.

---

## ADR-003: Sex-biased DE as ground truth

**Status**: Accepted, 2024 (Salman's thesis)

**Context**: Need a reproducible ground truth for DE detection benchmarking.

**Options considered**:
- **Simulation**: standard approach (scPower, powsimR, muscat). Rejected because
  real empirical ground truth is the paper's contribution.
- **ERCC spike-ins**: doesn't reflect biological complexity.
- **Sex-biased DE**: accepted. Near-50% case-control split in any cohort; XIST and
  Y-linked genes provide a validated biological reference, and the small-effect
  sex-biased genes are predominantly autosomal.

**Consequences**: The paper's headline generalisation claim (to disease contrasts)
depends on whether the (dispersion, baseMean, |LFC|) joint distribution of
sex-biased DE is representative of disease contrasts. This is the paper's single
biggest reviewer-vulnerability. Addressed via the chromosome-stratified recall
heatmaps and the autosomal vs sex-chromosome truth-gene comparison in the
supplement.

---

## ADR-004: Ship paper in repo as Rmd

**Status**: Accepted. The original "render via Nextflow" implementation was
superseded by ADR-001b and then ADR-007.

**Context**: For the paper to be truly reproducible, the manuscript source must be
in the repo alongside the code and data.

**Options considered**:
- **Paper in separate repo**: rejected. Reviewers should reach the manuscript via
  a single clone.
- **Convert Rmd → Quarto**: viable but non-trivial (bookdown cross-references,
  bibliography tooling). Deferred.
- **Keep Rmd, include in repo**: accepted.

**Consequences**: `nature_comms_draft.rmd` and `nature_comms_supplementary.rmd`
are the sources; paths within are `Datalink/`-relative. Quoted numbers are read
at knit time through `val()` (`_paper_values.R`) and supplementary figure numbers
through `SF()` (`_supp_fig_index.R`), so a rerun updates the text as well as the
figures. The `knit` stage of `pipeline/` renders both documents.

---

## ADR-005: Use `bioconductor/bioconductor_docker` as container base

**Status**: Accepted, 2026-04-21. **Amended 2026-06-05**: base moved to
`RELEASE_3_22` (was `RELEASE_3_21`) to match `renv.lock`, which pins
Bioconductor 3.22 package versions throughout. The 3.21 base installed 3.22
package versions anyway (renv resolves from the lockfile), so this only aligns
the base layer with what was already being installed; verified the effective
package versions are unchanged.

**Context**: Need a reproducible R+Bioconductor environment with pinned package
versions.

**Options considered**:
- **Custom image from scratch**: overkill; re-implements Bioconductor tooling.
- **`rocker/r-ver`**: no Bioconductor packages pre-installed; slow build.
- **`bioconductor/bioconductor_docker`**: accepted. Pre-built with Bioconductor
  (R 4.5); release tag tracks `renv.lock` (currently `RELEASE_3_22`).

**Consequences**: Container size ~3 GB (large). The base is pinned by digest
in `Dockerfile`, not only by the mutable release tag; apt packages added on top
are pinned in `docker/apt-packages.txt`, R packages come from `renv.lock`, and
Python packages are installed with `pip --require-hashes`. The image is built
locally with `docker/docker-build.sh` (there is no CI build) and converted to
Apptainer with `docker/build-sif.sh` for cluster use.

---

## ADR-006: AI-agent-first documentation

**Status**: Accepted, 2026-04-21. **Amended**: the agent guide (`AGENTS.md`) is
kept out of the public repository (it is listed in `.gitignore`).

**Context**: The UKDRI AI Initiative is developing standards for lab software in
the agentic-AI era. This repo is a pilot for what those standards should look like.

**Decision**: Keep machine-readable metadata (`CITATION.cff`, `Dockerfile`
LABELs, `pipeline/stages.tsv`) in the repository so tools can discover how the
project is built and run, and keep README and per-pipeline READMEs precise
enough for both humans and agents.

**Consequences**: Significant architecture changes must be reflected in this
file and in the relevant README.

---

## ADR-007: Staged PBS pipeline for full reproduction

**Status**: Accepted, 2026-10-01. Supersedes the single-host route of ADR-001b
for reproducing the paper; the bash drivers remain for single-pipeline reruns.

**Context**: Three DE methods × 10 robustness runs × two downsampling modes made
the grid far too large for one host, and the Imperial CX3 cluster has no internet
on compute nodes and penalises many small files on its shared filesystem.

**Decision**: `pipeline/run.sh` submits each stage in `pipeline/stages.tsv` as a
dependency-chained PBS job (the grid as an array job). Every stage runs in the
pinned Apptainer image, stages inputs to node-local scratch as a few large tars,
and writes its outputs as tars. Raw inputs and an offline `reference/` bundle
(biomaRt caches, annotation tables, citation style) are checked against
checksums. Unpacking results is a separate script (`pipeline/unpack_results.sh`)
that is not part of the pipeline.

**Consequences**: One command reproduces every figure, table, quoted number and
Supplementary Data file and knits the manuscript. The pipeline is written for
PBS Pro; other schedulers need a port of `pipeline/lib/`. The older hand-driven
scripts in `hpc/` are kept for reference only.
