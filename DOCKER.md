# Running the pipeline in Docker

The Dockerfile bakes the entire compute environment (R 4.5 + Bioconductor 3.22 +
Python/anndata + every pinned package) into a reproducible image, so you can
build and run the pipeline without installing anything on the host except Docker.

The workflow is four small, single-purpose scripts:

| Script | Does | Typical runtime |
|--------|------|-----------------|
| `docker-build.sh` | Build the image | ~30–45 min cold |
| `docker-smoke.sh` | Verify the image (loads every critical R package + reticulate/anndata) | ~30 s |
| `docker-run.sh` | Run the pipeline inside the image | minutes → days |
| `docker-all.sh` | Convenience wrapper: build → smoke → run, with per-phase skips | — |

> These replace the old single `docker-overnight.sh` orchestrator. Same job,
> split into composable pieces with less hidden machinery.

---

## Prerequisites

- **Docker** with the daemon running and reachable by your user
  (`docker info` must succeed; if you see a permissions error, add yourself to
  the `docker` group: `sudo usermod -aG docker $USER`, then start a new shell).
- **buildx** (BuildKit) is used automatically if present, otherwise the build
  falls back to the legacy builder — both produce the same image. BuildKit is
  the default in Docker 23+; if it's missing, install the `docker-buildx`
  plugin (it's the more reliable path on hosts using the containerd snapshotter).
- **Disk**: ~13 GB for the image; ≥50 GB free wherever `Datalink/` lives if you
  run the full pipeline (the Cartesian sweep writes a lot).

---

## Quick start

```bash
# Fresh: build + smoke-test + run power-curve fitting (the default pipeline)
./docker-all.sh

# Just verify the image builds and the environment is correct (no pipeline)
./docker-all.sh --skip-run          # == build + smoke

# Already built the image? Skip ahead and just run:
./docker-all.sh --skip-build --skip-smoke --pipeline all --max-cores 24

# Long unattended run:
nohup ./docker-all.sh --pipeline all --max-cores 24 > nohup.docker-all.out 2>&1 & disown
```

You can always call the phases individually:

```bash
./docker-build.sh
./docker-smoke.sh
./docker-run.sh --pipeline pcf --max-cores 24
```

---

## The scripts

### `docker-build.sh`
Builds `empirical-power-analysis:overnight` from the repo root. Auto-detects
buildx (BuildKit) vs the legacy builder. R/Bioconductor packages are compiled
from source, so a cold build takes a while; this is a one-time cost.

```
--image NAME:TAG    image tag to build (default: empirical-power-analysis:overnight)
-h | --help
```

### `docker-smoke.sh`
Runs the image with `PIPELINE=smoke` using the **same runtime contract as a real
run** (repo + `Datalink/` bind-mounted, dropped to your host UID:GID). Loads every
critical R package plus reticulate/anndata and prints `SMOKE TEST PASSED`. No raw
data needed — this checks the environment, not the science.

```
--image NAME:TAG    image to test
--datalink PATH     host Datalink dir to mount (default: ./Datalink)
-h | --help
```

### `docker-run.sh`
Runs the pipeline inside the image. Foreground `docker run --rm` (PID-1 tini
forwards signals; the exit code is the pipeline's own). Writes a provenance
`manifest.json`, tees output to a log, and translates the container exit code.

```
--pipeline pcf|all       which driver to run inside the container (default: pcf)
--max-cores N            cap each Rscript step (default: min(nproc, 24))
--no-overwrite | --force resume vs regenerate (default: --no-overwrite)
--mode cells|reads|both  Cartesian sweep mode (default: both)
--image NAME:TAG         image to run
--datalink PATH          host Datalink dir (default: ./Datalink)
--extra-arg "STR"        forward an extra arg to the in-container driver (repeatable)
--dry-run                print the docker command + write a manifest, run nothing
-h | --help
```

Exit codes: `0` success · `1` pipeline ran with partial failures · `2` pre-flight
failed · `6` container crashed (OOM / signal / other).

### `docker-all.sh`
Chains build → smoke → run, stopping at the first failure. Skips are orthogonal:

```
--skip-build    reuse the existing image
--skip-smoke    skip the smoke test
--skip-run      build + smoke only (verify the image)
--dry-run       print each phase's command, run nothing
```
All other flags (`--image`, `--datalink`, `--pipeline`, `--max-cores`,
`--no-overwrite`/`--force`, `--mode`, `--extra-arg`) are routed to the right phase.

---

## Pipelines

- `pcf` (default) — power-curve fitting only: the minimum DAG
  (`run_power_curve_fitting.sh`), Roussos shards only.
- `all` — the full DAG (`run_all.sh`).

---

## Logs & provenance

`docker-run.sh` writes to `logs/docker-run/<UTC-timestamp>/` (git-ignored), with a
`latest` symlink:

- `run.log` — full container stdout/stderr
- `manifest.json` — git SHA, image digest, pipeline, args, timings, exit status

The pipeline's own per-stage logs are written under `Datalink/Run_Output_Logs/`.

---

## How the image works

- **Reproducible library**: packages are restored from `renv.lock` into
  `/opt/renv-library` (outside `/workspace` so the runtime bind-mount can't shadow
  them). They're made visible to every R session via `R_LIBS_SITE`, baked into the
  image's `Renviron.site` (see `docker/install_packages.R`).
- **renv autoloader is disabled at runtime** (`RENV_CONFIG_AUTOLOADER_ENABLED=FALSE`
  in `docker/entrypoint.sh`) so the bind-mounted host `.Rprofile` can't re-activate
  renv and shadow the baked library.
- **Python/anndata** is pre-baked into a venv at `/opt/python-env` and used via
  `reticulate` (`RETICULATE_PYTHON`), so no Miniconda download at runtime.
- **Runtime user**: the container starts as root to fix up `/etc/passwd`, then
  drops to your host UID:GID via `gosu`, so files written into `Datalink/` are
  owned by you, not root.
- **Base image**: `bioconductor/bioconductor_docker:RELEASE_3_22` (R 4.5), matching
  the Bioconductor 3.22 versions pinned in `renv.lock`.

---

## Troubleshooting

- **`permission denied … /var/run/docker.sock`** — your user isn't in the `docker`
  group. `sudo usermod -aG docker $USER` and start a new shell.
- **`BuildKit is enabled but the buildx component is missing`** — install the
  `docker-buildx` plugin, or the build will use the legacy builder automatically.
- **Build is slow** — expected; R/Bioconductor packages compile from source
  (`renv::restore()` does not use prebuilt binaries). It's a one-time cost; the
  image caches.
- **Smoke fails on a package** — check the build log; usually a missing system
  `-dev` dependency (add it to the `apt-get install` block in the `Dockerfile`).
- **Pipeline fails immediately on missing inputs** — make sure `Datalink/` (and its
  `raw_data`) is present and passed via `--datalink` if it isn't `./Datalink`.

---

## See also

- `REPRODUCE.md` — full reproduction walkthrough, including host-native (no Docker).
- `Dockerfile`, `docker/entrypoint.sh`, `docker/install_packages.R` — image internals.
