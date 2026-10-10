# The container image

Everything in this repository runs inside one container image. The image holds
software only, never data: the pipeline binds the code snapshot, raw inputs and
reference bundle into it at run time. `pipeline/run.sh` runs the Apptainer
(`.sif`) form of the image; Docker is only needed to build it.

Current version: **`empirical-power-analysis:v1.0.0`**.
Registry: `ghcr.io/neurogenomics/empirical_power_analysis` (not yet published;
until it is, build from source as below).

## What is pinned

| Layer | Pinned by | File |
| --- | --- | --- |
| Base OS, R 4.5.2, Bioconductor 3.22 | image digest | `Dockerfile` (`FROM ...@sha256:...`) |
| System packages | exact apt versions | `docker/apt-packages.txt` |
| renv itself | version and tarball SHA-256 | `Dockerfile` (`RENV_VERSION`, `RENV_TARBALL_SHA256`) |
| R packages | exact versions and sources | `renv.lock`, restored by `docker/install_packages.R` |
| Python packages (anndata, matplotlib, ...) | exact versions and wheel hashes | `docker/python-requirements.txt` (`pip --require-hashes --only-binary=:all:`) |

The build stops if any R package ends up at a version other than its
`renv.lock` pin, or any pinned Python package is missing or at another version.

## Smoke test

`docker/smoke_test.R` checks, with no data:
1. every package in `renv.lock` is installed at exactly the locked version;
2. the packages the pipeline uses (DESeq2, edgeR, limma, ashr, xgboost, ...)
   load;
3. every Python pin in `python-requirements.txt` is installed at exactly that
   version, and reticulate can import anndata;
4. the fonts the figure composers need are present.

It runs automatically at the end of `docker-build.sh` (as the last build
step), inside every `.sif` that `build-sif.sh` produces, and on demand:

```bash
./docker/docker-smoke.sh --image empirical-power-analysis:v1.0.0
apptainer exec --cleanenv --no-home epa.sif Rscript /usr/local/share/epa/smoke_test.R
```

It prints `SMOKE TEST PASSED` or lists every mismatch and exits non-zero.

## Getting the `.sif`

Pick one route. Each ends with a `.sif` and its `.sha256` next to it.

### A. Pull the published image (no Docker needed)

Once the image is published, pull it **by digest**, never by tag:

```bash
./docker/build-sif.sh \
  --from-registry ghcr.io/neurogenomics/empirical_power_analysis@sha256:<digest> \
  --out empirical-power-analysis.sif
```

### B. Build from source (Docker + Apptainer on one machine)

```bash
./docker/docker-build.sh                 # ~1 h on a cold cache; ends with the smoke test
./docker/build-sif.sh --out empirical-power-analysis.sif
```

Both default to the tag `empirical-power-analysis:v1.0.0`. The build needs
internet access (Bioconductor, CRAN, PyPI) but no data. Docker needs a running
daemon and either root or membership of the `docker` group. BuildKit is used
if the `buildx` plugin exists; the legacy builder works too.

### C. Build on a workstation, run on a cluster without Docker

Do route B on any machine with Docker and Apptainer, then copy
`empirical-power-analysis.sif` and its `.sha256` to the cluster and check it:

```bash
sha256sum -c empirical-power-analysis.sif.sha256
```

## Using the `.sif` in the pipeline

Put it at `WORK/epa.sif` (or set `EPA_SIF=` in `WORK/epa.env`). To make a run
refuse any other image, copy the checksum into `WORK/epa.env`:

```bash
EPA_SIF_SHA256=$(cut -d' ' -f1 empirical-power-analysis.sif.sha256)
```

Whether or not it is pinned, every run records the image checksum in
`results/PROVENANCE.txt`.

## Running the image by hand

```bash
docker run --rm empirical-power-analysis:v1.0.0             # smoke test (default)
docker run --rm -it empirical-power-analysis:v1.0.0 shell   # interactive bash
docker run --rm -v "$PWD":/workspace empirical-power-analysis:v1.0.0 Rscript -e 'sessionInfo()'

apptainer exec --cleanenv --no-home -B "$PWD":/workspace --pwd /workspace epa.sif Rscript ...
```

With Apptainer, always pass `--cleanenv --no-home`, as `pipeline/run.sh` does.
Without `--no-home`, your home directory is mounted, and an `R_LIBS_USER` set
in `~/.Renviron` puts your own R library ahead of the pinned one.

## Changing the software

Any change to the image can change results, so it gets a new version tag.

1. Add or update the package in an R 4.5 / Bioconductor 3.22 session, then run
   `renv::snapshot()` and commit the new `renv.lock`. For Bioconductor
   packages, check the record has `Repository`, `git_branch` and
   `git_last_commit`. Python: add the exact version and its hashes to
   `docker/python-requirements.txt`.
2. Never add an unpinned `install.packages()`, `BiocManager::install()` or
   `pip install` to the Dockerfile.
3. Rebuild with a new tag (`--image empirical-power-analysis:v1.1.0`), build a
   new `.sif`, and rerun every stage whose code uses the changed package.

## Publishing a release (maintainers)

```bash
echo "$GHCR_TOKEN" | docker login ghcr.io -u <github-user> --password-stdin   # token needs write:packages
docker tag  empirical-power-analysis:v1.0.0 ghcr.io/neurogenomics/empirical_power_analysis:v1.0.0
docker push ghcr.io/neurogenomics/empirical_power_analysis:v1.0.0
docker image inspect --format '{{index .RepoDigests 0}}' ghcr.io/neurogenomics/empirical_power_analysis:v1.0.0
```

Record the printed `@sha256:` digest in this file and in `README.md`, and make
the package public in the GitHub organisation settings. Never push a different
image under a tag that has already been published.
