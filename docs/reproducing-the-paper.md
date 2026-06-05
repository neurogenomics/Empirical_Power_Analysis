# Reproducing the paper

See [`REPRODUCE.md`](../REPRODUCE.md) for the complete cold-start reproduction
guide (Sections 0-3 setup + Sections A-I per-pipeline walkthrough).

For one-command containerised reproduction, run `./docker-all.sh` from the repo
root (see [`DOCKER.md`](../DOCKER.md) and `REPRODUCE.md` §2.6). For a host-native
bash driver that does not use Docker, see `./run_all.sh` (REPRODUCE.md §2.5).

For the design rationale behind pipeline choices, see
[`architecture.md`](architecture.md).
