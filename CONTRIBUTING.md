# Contributing

Thank you for considering contributing to the Empirical Power Analysis pipeline.

## Development loop

1. **Fork** the repo and create a feature branch from `main`:
   `git checkout -b feat/your-change main`
2. **Build and smoke-test the container** to get the pinned environment:
   `./docker/docker-build.sh --image empirical-power-analysis:<tag>`, then
   `./docker/docker-smoke.sh --image empirical-power-analysis:<tag>` (see
   `docs/container.md`)
3. **Make your change**, respecting the scientific safety rails below
4. **Add tests**: every scientific-logic change needs a matching `tests/` update
5. **Run tests locally** before pushing:
   ```bash
   # R unit tests
   for d in R/*/tests; do [ -d "$d" ] && Rscript -e "testthat::test_dir('$d')"; done

   # Reproduction pipeline tests (stubbed scheduler, no cluster needed)
   bash pipeline/tests/test_pipeline.sh
   ```
6. **Open a PR** against `main`. There is no CI at present, so say in the PR
   which tests you ran.

## Scientific safety rails

- Every randomised step is seeded deterministically from the robustness run,
  grid point and iteration; keep it that way.
- The ground-truth and downsampling halves must never share donors.
- Filter and significance thresholds (`GT_LFC_THRESHOLD`, `FDR_THRESHOLD`,
  `MIN_GENE_COUNT`, LFC bins) are load-bearing: changing one changes the paper.
- Report cell types in the order Astro, Micro, Oligo.
- Do not bump DESeq2, edgeR, limma, ashr or xgboost in `renv.lock` without
  re-running the full grid.

## Pull request expectations

- **Scope**: one logical change per PR. Big refactors get their own RFC issue first.
- **Commit messages**: imperative mood ("Add depth-equalised microglia analysis"),
  referencing issue numbers where relevant.
- **Tests**: every scientific-logic change ships with tests
- **Docs**: every public-facing change updates `README.md` and/or `docs/`
- **Model card**: if the fitted power model changes, update
  `docs/model-card-predict_power.md`
- **Paper impact**: if results change, flag it explicitly; the manuscript needs
  regeneration

## Release process

1. Bump `version` in `CITATION.cff`
2. Update `CHANGELOG.md`
3. Tag with semver: `git tag -a v1.0.0 -m "Release v1.0.0"` and push the tag
4. Create a GitHub Release from the tag
5. Add the release date to `CITATION.cff`. The persistent identifier is the
   paper DOI (currently the bioRxiv preprint); no separate archive DOI is minted,
   because every input is already public and the pipeline regenerates the rest.
   Update `preferred-citation` when the journal version is published.

## Reporting bugs

Open an issue using the bug report template. Include:
- How you ran the pipeline: `pipeline/run.sh` on a PBS cluster or with
  `--local`, and any `--only` stages
- The container image you used (tag and sha256; both are recorded in
  `results/PROVENANCE.txt`) and your Apptainer version (`apptainer --version`)
- Log output: the failing stage's PBS log (or the console output with
  `--local`) and `results/PROVENANCE.txt`
- Expected vs actual behaviour

## Scientific questions

For scientific questions about the results, methods, or interpretation, please read
the paper first (`nature_comms_draft.rmd` and `nature_comms_supplementary.rmd`)
and then open an issue tagged `science`. We respond to methodological critiques in
good faith: disagreement is how this gets better.

## Code of conduct

Be excellent to each other. We follow the
[Contributor Covenant v2.1](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).

## License

By contributing, you agree your contributions will be licensed under the MIT License.
