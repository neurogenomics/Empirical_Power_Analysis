# Running the pipeline in Docker

The Dockerfile bakes the entire compute environment (R 4.5 + Bioconductor 3.22 +
Python/anndata + every pinned package) into a reproducible image, so you can
build and run the pipeline without installing anything on the host except Docker.
The base image is pinned by digest, R packages restore from `renv.lock`, apt
additions install from the exact package closure in `docker/apt-packages.txt`,
and Python wheels install from `docker/python-requirements.txt` with
`pip --require-hashes`.

The workflow is four small, single-purpose scripts, which live in `docker/` and
are run from the repo root:

| Script | Does | Typical runtime |
|--------|------|-----------------|
| `docker/docker-build.sh` | Build the image | ~30–45 min cold |
| `docker/docker-smoke.sh` | Verify the image (loads every critical R package + reticulate/anndata) | ~30 s |
| `docker/docker-run.sh` | Run the pipeline inside the image | minutes → days |
| `docker/docker-all.sh` | Convenience wrapper: build → smoke → run, with per-phase skips | — |

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
./docker/docker-all.sh

# Just verify the image builds and the environment is correct (no pipeline)
./docker/docker-all.sh --skip-run          # == build + smoke

# Already built the image? Skip ahead and just run:
./docker/docker-all.sh --skip-build --skip-smoke --pipeline all --max-cores 24

# Long unattended run:
nohup ./docker/docker-all.sh --pipeline all --max-cores 24 > nohup.docker-all.out 2>&1 & disown
```

You can always call the phases individually:

```bash
./docker/docker-build.sh
./docker/docker-smoke.sh
./docker/docker-run.sh --pipeline pcf --max-cores 24
```

For paper-producing runs, use an immutable image tag rather than the mutable
default `empirical-power-analysis:overnight`:

```bash
IMAGE=empirical-power-analysis:2026-06-25-3de-frozen
./docker/docker-build.sh --image "$IMAGE"
./docker/docker-smoke.sh --image "$IMAGE"
./docker/docker-run.sh --image "$IMAGE" --mini-run
```

`docker-run.sh` records both the image tag and image ID in
`logs/docker-run/<timestamp>/manifest.json`; keep that manifest with any result
set that may feed the paper.

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

The `PIPELINE` env var selects what the container runs:

- `pcf` (default) — power-curve fitting only: the minimum DAG
  (`run_power_curve_fitting.sh`), Roussos shards only.
- `all` — the full DAG (`run_all.sh`).
- `smoke` — environment check only (loads every critical package), no pipeline;
  this is what `docker-smoke.sh` drives.
- `shell` — drop into `bash` inside the fully-configured image. Use this to run
  an individual pipeline step (not the whole DAG) or to debug interactively. See
  [Running specific pipeline steps](#running-specific-pipeline-steps-in-the-container).

---

## Running specific pipeline steps in the container

`docker-run.sh` only knows the whole-DAG drivers (`pcf`, `all`). Sometimes you
want to run **a subset of steps** inside the reproducible image — most commonly
to **aggregate a freshly harvested Cartesian grid and redraw the heatmaps**
without re-running the expensive DE fits. Each `R/<pipeline>/run_pipeline.R`
exposes the stable `Rscript run_pipeline.R [steps] [n_cores]` CLI (see
[`AGENTS.md`](AGENTS.md) § Interface contracts), so you drive it directly via the
`shell` pipeline.

**When to use this:** you pulled HPC result tars back to the analysis server,
unpacked them into `Datalink/Cartesian_Downsampling/` (see
[`hpc/README.md`](hpc/README.md) and `hpc/harvest-grid-results.sh`), and now want
the summary checkpoints + heatmaps the visualisation and paper steps consume.
For the Cartesian pipeline: **step 3 = aggregate**, **step 4 = power heatmaps**,
**step 7 = continuous-grid surfaces**.

**Runtime contract** (identical to `docker-run.sh`, just with a custom command):

- Bind-mount the repo to `/workspace` and `Datalink/` under it.
- `POWER_ANALYSIS_ROOT=/workspace`, `HOME=/tmp` (writable),
  `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` (use the baked `/opt/renv-library`, not
  the host `.Rprofile`).
- Pass `TARGET_UID`/`TARGET_GID` so output PNGs are owned by you, not root.
- **`DOWNSAMPLE_MODE` must be set per mode** (`cells` or `reads`) — it is
  config-only and cannot be passed as a positional arg, so run once per mode.

```bash
# On the analysis server (iclneuro), in a tmux — aggregation reads tens of
# thousands of per-fit .qs and takes 20-40 min per mode.
cd /path/to/Empirical_Power_Analysis
tmux new -s agg
for MODE in cells reads; do
  docker run --rm \
    -e PIPELINE=shell \
    -e POWER_ANALYSIS_ROOT=/workspace \
    -e HOME=/tmp \
    -e RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    -e DOWNSAMPLE_MODE="$MODE" \
    -e TARGET_UID=$(id -u) -e TARGET_GID=$(id -g) \
    -v "$PWD:/workspace" \
    -v "$PWD/Datalink:/workspace/Datalink" \
    empirical-power-analysis:overnight \
    -c 'Rscript R/cartesian_downsampling/run_pipeline.R 3,4 8'
done
```

Step 3 writes the aggregated summary checkpoints; step 4 writes the heatmap PNGs
into each method's `…/Plots/` directory. The same pattern runs any pipeline's
steps — swap the `R/<pipeline>/run_pipeline.R` path and step numbers.

**Notes**

- Run it on the analysis server, **not** an HPC compute node — this is
  minutes-of-work post-processing, not the DE grid. The heavy DE fits already ran
  on HPC; post-processing stays where the harvested `Datalink/` lives.
- Use `tmux`/`nohup` for the aggregate: it is long enough that a dropped SSH
  session (or a synchronous, no-backgrounding gateway) would kill it mid-run, and
  the checkpoint is only written at the end.
- The harvested grid must already be in `Datalink/Cartesian_Downsampling/` before
  you aggregate. If a mode is only partially harvested, its surfaces will be
  statistically thinner but still valid; aggregation is glob-based over whatever
  `Run_*` dirs exist.

---

## Logs & provenance

`docker-run.sh` writes to `logs/docker-run/<UTC-timestamp>/` (git-ignored), with a
`latest` symlink:

- `run.log` — full container stdout/stderr
- `manifest.json` — git SHA, image digest, pipeline, args, timings, exit status

The pipeline's own per-stage logs are written under `Datalink/Run_Output_Logs/`.

---

## How the image works

- **Base image**: `bioconductor/bioconductor_docker:RELEASE_3_22` pinned by
  digest in `Dockerfile`, so the tag cannot silently move under future rebuilds.
- **System packages**: apt additions are installed from `docker/apt-packages.txt`,
  an exact dpkg diff against the pinned base image. Packages already present in
  the pinned base image are not re-requested, to avoid silent upgrades.
- **Reproducible R library**: packages are restored from `renv.lock` into
  `/opt/renv-library` (outside `/workspace` so the runtime bind-mount can't shadow
  them). They're made visible to every R session via `R_LIBS_SITE`, baked into the
  image's `Renviron.site` (see `docker/install_packages.R`).
- **renv autoloader is disabled at runtime** (`RENV_CONFIG_AUTOLOADER_ENABLED=FALSE`
  in `docker/entrypoint.sh`) so the bind-mounted host `.Rprofile` can't re-activate
  renv and shadow the baked library.
- **Python/anndata** is pre-baked into a venv at `/opt/python-env` and used via
  `reticulate` (`RETICULATE_PYTHON`). Wheels are installed with
  `--only-binary=:all:` and `--require-hashes` from
  `docker/python-requirements.txt`, so pip cannot choose different versions or
  different files on rebuild.
- **Runtime user**: the container starts as root to fix up `/etc/passwd`, then
  drops to your host UID:GID via `gosu`, so files written into `Datalink/` are
  owned by you, not root.

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

- `reproducing-the-paper.md` — full reproduction walkthrough, including host-native (no Docker).
- `Dockerfile`, `docker/entrypoint.sh`, `docker/install_packages.R` — image internals.
