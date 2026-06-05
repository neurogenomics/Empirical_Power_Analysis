# Architecture Decision Records (ADRs)

Each ADR captures a significant design choice, its context, and the reasoning behind
it, so future contributors understand why the code looks the way it does.

Format: [Michael Nygard's ADR template](https://github.com/joelparkerhenderson/architecture-decision-record).

---

## ADR-001 — Use Nextflow for pipeline orchestration

**Status**: SUPERSEDED by ADR-001b (2026-04-28). Original decision retained
below for historical context.

**Context**: The Cartesian downsampling step is the hot loop — up to 120,000 DESeq2
runs per cell type × mode × robustness. This exceeds single-machine capacity when
extended (e.g. to neurons or a second cohort). We need a workflow manager that can
(a) dispatch to SLURM with per-step resource declarations, (b) resume from failure,
(c) track provenance, and (d) integrate with containers.

**Options considered**:
- **GNU Make** — rejected. Cannot submit SLURM jobs or declare per-process resources.
- **Snakemake** — viable. Python-based DAG, SLURM profile support. Rejected mainly
  for lab consistency: the Skene/Neurogenomics lab already maintains multiple
  Nextflow repos (`NextFlowIntro`, `nf-imperial-canary`, `wdl-scflow`,
  `scRNAseq-pipeline`, `LigerFlow`). A second paradigm would fragment institutional
  knowledge.
- **R `targets`** — viable for R-native workflows. Rejected because (a) HPC dispatch
  via `future.batchtools::batchtools_slurm` is finicky and (b) lab is Nextflow-first.
- **Nextflow DSL2** — accepted. Native SLURM, container integration, matches lab
  conventions, reviewer-credible at the journal tier we target.

**Consequences**: Maintaining both `R/<pipeline>/run_pipeline.R` (inner parallelism
via `future.apply`) and `modules/local/<name>/main.nf` (outer orchestration). These
must stay in sync when pipelines are added or renamed. CI stub-run enforces this.

---

## ADR-001b — Replace Nextflow with bash drivers + Docker image

**Status**: Accepted, 2026-04-28

**Context**: ADR-001 chose Nextflow for SLURM dispatch and lab consistency.
In practice we never deployed to SLURM for this paper — all production runs
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
- **Stay on Nextflow, fix the drift** — viable but does not address the
  external reproducibility friction.
- **Strict-serial bash driver (`run_all.sh`) + Docker image** — accepted.
  Single source of truth (the R pipelines), one shell script for full DAG,
  one for power-curve-only, and a set of small Docker scripts
  (`docker-build.sh` / `docker-smoke.sh` / `docker-run.sh`, chained by
  `docker-all.sh`) that build + smoke-test + run the pipeline inside a
  container. External reproducers run one command and walk away.

**Consequences**: `main.nf`, `nextflow.config`, `iclneuro.config`, `conf/`
and `modules/local/` are removed. The R pipelines and their `config.R`
honour `POWER_ANALYSIS_ROOT` and `POWER_ANALYSIS_FORCE_RECALC` environment
variables so the bash driver can flip behaviour without touching code. ADR-004
is also superseded — the paper is rendered via the bash driver / Docker
image, not via a `RENDER_PAPER` Nextflow process.

---

## ADR-002 — Pseudobulk with DESeq2 + apeglm

**Status**: Accepted, 2024 (inherited from Salman Fawad's thesis), 2026-04-21 reaffirmed

**Context**: Differential expression for 22M outcomes requires a DE method that
scales, gives calibrated p-values at small n, and has shrinkage-corrected effect
sizes.

**Options considered**:
- **Mixed models on cells (NEBULA, glmmTMB)** — rejected. ~10x slower per run,
  pseudoreplication debate (Zimmerman 2021 vs Murphy & Skene 2022).
- **edgeR pseudobulk** — viable. Squair et al. 2021 recommends it; marginally more
  sensitive at large n. Rejected for this paper: DESeq2's conservatism is a feature
  here because we are reporting on what the field *will* detect, not what it *could*.
- **DESeq2 pseudobulk + apeglm shrinkage** — accepted.

**Consequences**: Paper results are DESeq2-specific. eLife/Genome Biology revision
will likely require an edgeR sensitivity analysis. This is documented in
Limitations in the paper.

---

## ADR-003 — Sex-biased DE as ground truth

**Status**: Accepted, 2024 (Salman's thesis)

**Context**: Need a reproducible ground truth for DE detection benchmarking
across independent datasets.

**Options considered**:
- **Simulation** — standard approach (scPower, powsimR, muscat). Rejected because
  real empirical ground truth is the paper's contribution.
- **ERCC spike-ins** — doesn't reflect biological complexity.
- **Sex-biased DE** — accepted. Near-50% case-control split in any cohort; XIST and
  Y-linked genes provide a validated biological reference.

**Consequences**: The paper's headline generalisation claim (to disease contrasts)
depends on whether the (dispersion, baseMean, |LFC|) joint distribution of
sex-biased DE is representative of disease contrasts. This is the paper's single
biggest reviewer-vulnerability. Addressed via the chromosome-stratified recall
surfaces (Fig 2 suppl 3).

---

## ADR-004 — Ship paper in repo as Rmd, render via Nextflow

**Status**: PARTIALLY SUPERSEDED by ADR-001b (2026-04-28). The "ship paper as
Rmd in the repo" decision still stands. The "render via Nextflow"
implementation is replaced by rendering inside the Docker image.

**Context**: For the paper to be truly reproducible, the manuscript source must be
in the repo alongside the code and data.

**Options considered**:
- **Paper in separate repo** — rejected. Reviewers should reach the manuscript via
  a single clone.
- **Convert Rmd → Quarto** — viable but non-trivial (Bioconductor R Markdown
  extensions, bibliography tooling). Deferred to post-bioRxiv.
- **Keep Rmd, include in repo, render via Nextflow** — accepted.

**Consequences**: `paper/paper.Rmd` is the source; paths within must be
`Datalink/`-relative. The `RENDER_PAPER` Nextflow process produces `paper.html`
(and `paper.pdf` if TeX available).

---

## ADR-005 — Use `bioconductor/bioconductor_docker` as container base

**Status**: Accepted, 2026-04-21. **Amended 2026-06-05**: base pinned to
`RELEASE_3_22` (was `RELEASE_3_21`) to match `renv.lock`, which pins
Bioconductor 3.22 package versions throughout. The 3.21 base installed 3.22
package versions anyway (renv resolves from the lockfile), so this only aligns
the base layer with what was already being installed; verified the effective
package versions are unchanged.

**Context**: Need a reproducible R+Bioconductor environment with pinned package
versions.

**Options considered**:
- **Custom image from scratch** — overkill; re-implements Bioconductor tooling.
- **`rocker/r-ver`** — no Bioconductor packages pre-installed; slow build.
- **`bioconductor/bioconductor_docker`** — accepted. Pre-built with Bioconductor
  (R 4.5); release tag tracks `renv.lock` (currently `RELEASE_3_22`).

**Consequences**: Container size ~3 GB (large). CI builds and pushes to GHCR on
every main branch push.

---

## ADR-006 — AI-agent-first documentation

**Status**: Accepted, 2026-04-21

**Context**: The UKDRI AI Initiative is developing standards for lab software in
the agentic-AI era. This repo is a pilot for what those standards should look like.

**Decision**: Maintain an `AGENTS.md` at repo root as the canonical agent-facing
entry point. Include machine-readable metadata (`CITATION.cff`, `Dockerfile`
LABELs) to enable tool-independent discovery.

**Consequences**: Every significant architecture change must be reflected in both
`AGENTS.md` (agent guide) and `docs/architecture.md` (this file). PR template
enforces this.
