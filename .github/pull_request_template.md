<!-- Thank you for contributing! Please fill in the sections below. -->

## What does this PR do?

<!-- One or two sentences describing the intent -->

## Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Scientific methodology change (results may differ)
- [ ] Refactor (no behaviour change)
- [ ] Documentation
- [ ] CI / build / infrastructure

## Checklist

- [ ] I've read `AGENTS.md` (if working with an AI agent) / `CONTRIBUTING.md` (human)
- [ ] Unit tests added or updated in `R/<pipeline>/tests/`
- [ ] Toy fixtures updated in `R/<pipeline>/toy/` if behaviour changed
- [ ] `nextflow run main.nf -profile test,docker -stub-run` passes locally
- [ ] R unit tests pass locally
- [ ] `lintr::lint_dir('R')` is clean (or rationalised in the PR)
- [ ] Documentation updated (`README.md`, `docs/`, `AGENTS.md` as appropriate)
- [ ] Scientific code safety rails (AGENTS.md) respected — deterministic seeding,
      split integrity, filter thresholds, cell-type ordering
- [ ] If this changes results materially, I've flagged that in the PR description
      and noted paper-text implications
- [ ] If this changes the fitted power model, `docs/model-card-predict_power.md`
      updated
- [ ] `CHANGELOG.md` updated under `## Unreleased`

## Scientific impact

<!-- If results change: which figure/table/claim is affected? -->

## Attribution

<!-- If AI tools assisted with this PR, note which and for which parts.
     See AGENTS.md § Transparency about AI-assisted development. -->
