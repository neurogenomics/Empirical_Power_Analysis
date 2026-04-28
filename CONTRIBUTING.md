# Contributing

Thank you for considering contributing to the Empirical Power Analysis pipeline. This
document covers the workflow for humans. Agents should start with `AGENTS.md`.

## Development loop

1. **Fork** the repo and create a feature branch from `dev`:
   `git checkout -b feat/your-change dev`
2. **Install the dev environment**: `docker compose up -d` or pull the container
   directly: `docker pull ghcr.io/neurogenomics/empirical_power_analysis:dev`
3. **Build + smoke-test the container** before making changes to confirm baseline:
   `./docker-overnight.sh --build-only`
4. **Make your change** — respect the "scientific code safety rails" listed in `AGENTS.md`
5. **Add tests** — every scientific-logic change needs a matching `tests/` update
6. **Run tests locally** before pushing:
   ```bash
   # R unit tests
   for d in R/*/tests; do [ -d "$d" ] && Rscript -e "testthat::test_dir('$d')"; done

   # Wrapper sanity check
   ./run_all.sh --dry-run

   # Lint
   Rscript -e 'lintr::lint_dir("R")'
   ```
7. **Open a PR** against `dev`. Fill in the PR template. CI will run.

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

1. Bump `version` in `CITATION.cff` and the `Dockerfile` LABEL block
2. Update `CHANGELOG.md`
3. Merge `dev` → `main`
4. Tag with semver: `git tag -a v1.0.0 -m "Release v1.0.0"`
5. Push tag: `git push origin v1.0.0`
6. GitHub Actions builds the Docker image → `ghcr.io/neurogenomics/empirical_power_analysis:1.0.0`
7. Create a GitHub Release; Zenodo auto-mints a DOI
8. Update `CITATION.cff` with the new DOI and open a follow-up PR

## Reporting bugs

Open an issue using the bug report template. Include:
- How you ran the pipeline (host-native `./run_all.sh`, container via
  `./docker-overnight.sh`, or per-pipeline `Rscript R/<pipeline>/run_pipeline.R`)
- Docker version (`docker --version`) if containerised
- R version (`R --version`) if host-native
- Log output: attach the contents of
  `Datalink/Run_Output_Logs/run_all/latest/` (host-native) or
  `logs/docker-overnight/latest/` (container)
- Expected vs actual behaviour

## Scientific questions

For scientific questions about the results, methods, or interpretation, please read
the paper first (`paper/paper.Rmd`) and then open an issue tagged `science`. We
respond to methodological critiques in good faith — disagreement is how this gets
better.

## Code of conduct

Be excellent to each other. See `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1).

## License

By contributing, you agree your contributions will be licensed under the MIT License.
