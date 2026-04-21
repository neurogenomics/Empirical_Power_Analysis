# Architecture Decision Records (ADRs)

Each ADR captures a significant design choice, its context, and the reasoning behind
it, so future contributors understand why the code looks the way it does.

Format: [Michael Nygard's ADR template](https://github.com/joelparkerhenderson/architecture-decision-record).

---

## ADR-001 — Use Nextflow for pipeline orchestration

**Status**: Accepted, 2026-04-21

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

**Status**: Accepted, 2026-04-21

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

## ADR-005 — Use `bioconductor/bioconductor_docker:RELEASE_3_21` as container base

**Status**: Accepted, 2026-04-21

**Context**: Need a reproducible R+Bioconductor environment with pinned package
versions.

**Options considered**:
- **Custom image from scratch** — overkill; re-implements Bioconductor tooling.
- **`rocker/r-ver`** — no Bioconductor packages pre-installed; slow build.
- **`bioconductor/bioconductor_docker`** — accepted. Pre-built with Bioconductor
  3.21 (R 4.5.1), matches what Chris used for the v8 paper.

**Consequences**: Container size ~3 GB (large). Build time ~2 min once cached. CI
builds and pushes to GHCR on every main branch push.

---

## ADR-006 — AI-agent-first documentation

**Status**: Accepted, 2026-04-21

**Context**: The UKDRI AI Initiative is developing standards for lab software in
the agentic-AI era. This repo is a pilot for what those standards should look like.

**Decision**: Maintain an `AGENTS.md` at repo root as the canonical agent-facing
entry point. Include machine-readable metadata (CITATION.cff, nextflow.config
manifest) to enable tool-independent discovery.

**Consequences**: Every significant architecture change must be reflected in both
`AGENTS.md` (agent guide) and `docs/architecture.md` (this file). PR template
enforces this.
