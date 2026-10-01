#!/usr/bin/env bash
# ==============================================================================
# Harvest staged grid result tars into the real Datalink in ONE bulk unpack.
# ==============================================================================
# REBUILT 2026-06-29 (faithful reconstruction; originals lost, see cx3-grid.pbs).
# Run on a LOGIN NODE or on iclneuro — NOT on a compute node, and NOT as part of
# a grid job. This is the "unpack elsewhere" step of the tar-staging discipline:
# the array writes one tar per (mode, chunk) into HPC_Result_Tars/incoming/; this
# script merges them into Datalink.
#
# DISK-SAFE MODE: tars are processed ONE AT A TIME (extract -> merge -> verify ->
# delete) rather than all-extracted-then-one-rsync. Each source tar is deleted to
# reclaim space ONLY after every one of its members is confirmed present in
# Datalink. A truncated/corrupt tar, a failed extract, or an incomplete merge
# leaves the tar in place and flags it. Re-running is safe and idempotent: an
# already-harvested tar verifies immediately and is deleted, cleaning up leftovers.
# WARNING: on success the source tars in incoming/ are PERMANENTLY removed (no
# harvested/ archive) — by design, to free space on a full disk.
#
# Usage:
#   bash hpc/harvest-grid-results.sh [DATALINK] [TAG]
#   bash hpc/harvest-grid-results.sh                       # ./Datalink, ALL tags
#   bash hpc/harvest-grid-results.sh Datalink              # ALL tags + legacy
#   bash hpc/harvest-grid-results.sh Datalink deseq2-edger-limma_runs06-10
#                                                          # one tag only
# Grids are namespaced by tag under HPC_Result_Tars/incoming/<tag>/. With no TAG
# this harvests every tag's tars AND any legacy untagged tars sitting directly in
# incoming/ (the original 5-run grid). With a TAG it harvests only that grid.
# Merging is idempotent (rsync --ignore-existing) and the in-tar paths are keyed
# by Run_<id>/<method>, so multiple grids coexist in one Datalink with no clash.
# ==============================================================================
set -euo pipefail

DATALINK="${1:-./Datalink}"
TAG="${2:-}"
[ -d "$DATALINK" ] || { echo "FATAL: Datalink not found: $DATALINK" >&2; exit 1; }

INCOMING="$DATALINK/HPC_Result_Tars/incoming"

# Collect result tars. A TAG restricts to incoming/<tag>/; otherwise recurse so
# both tagged subdirs and legacy top-level tars are caught. `-name '*.tar'`
# excludes in-progress hidden `.epa-control...partial.tar` / `.epa_..._partial`.
SCAN_ROOT="$INCOMING"
[ -n "$TAG" ] && SCAN_ROOT="$INCOMING/$TAG"
[ -d "$SCAN_ROOT" ] || { echo "no incoming dir to harvest: $SCAN_ROOT"; exit 0; }

mapfile -d '' tars < <(find "$SCAN_ROOT" -type f -name '*.tar' ! -name '.*' -print0)
if [ ${#tars[@]} -eq 0 ]; then
  echo "no result tars to harvest in $SCAN_ROOT"
  exit 0
fi

echo "harvesting ${#tars[@]} result tar(s) from $SCAN_ROOT — extract, merge, verify, then delete each on success"

harvested=0
kept=0
for t in "${tars[@]}"; do
  # Integrity + manifest in one step: `tar tf` fails on a truncated/corrupt tar,
  # so a bad tar is KEPT (never deleted) and flagged for attention.
  if ! members="$(tar tf "$t" 2>/dev/null)"; then
    echo "WARN: cannot read $t (truncated/corrupt?) — keeping, skipping" >&2
    kept=$((kept + 1))
    continue
  fi
  # Reject tars with absolute or parent-traversal member paths.
  if printf '%s\n' "$members" | grep -qE '(^/|(^|/)\.\.(/|$))'; then
    echo "WARN: unsafe path in $t — keeping, skipping" >&2
    kept=$((kept + 1))
    continue
  fi

  # Fast path: if every member is ALREADY in Datalink, this tar is fully
  # harvested. Delete it with zero extra space (no extract) — this is how a
  # full disk reclaims room before processing the not-yet-merged tars.
  already=1
  while IFS= read -r m; do
    m="${m#./}"
    [ -z "$m" ] && continue
    if [ ! -e "$DATALINK/$m" ]; then already=0; break; fi
  done <<< "$members"
  if [ "$already" -eq 1 ]; then
    rm -f "$t"
    harvested=$((harvested + 1))
    echo "  ok already-merged + deleted $(basename "$t")"
    continue
  fi

  # Extract THIS tar to its own scratch dir, then merge into Datalink. Scratch
  # holds only one tar's worth at a time, so a full disk never has to host all
  # tars' contents simultaneously (the old "extract all then one rsync" path).
  work="$(mktemp -d)"
  if ! tar xf "$t" -C "$work"; then
    echo "WARN: extract failed for $t — keeping, skipping" >&2
    rm -rf "$work"
    kept=$((kept + 1))
    continue
  fi
  # Idempotent merge: --ignore-existing means already-harvested tars are no-ops
  # here and will simply verify + delete below (auto-cleanup of completed tars).
  rsync -a --ignore-existing "$work"/ "$DATALINK"/

  # Verify every member now exists under Datalink BEFORE freeing the source.
  ok=1
  while IFS= read -r m; do
    m="${m#./}"
    [ -z "$m" ] && continue
    if [ ! -e "$DATALINK/$m" ]; then
      echo "WARN: member missing after merge: $m (from $t)" >&2
      ok=0
      break
    fi
  done <<< "$members"
  rm -rf "$work"

  if [ "$ok" -eq 1 ]; then
    rm -f "$t"
    harvested=$((harvested + 1))
    echo "  ok merged + deleted $(basename "$t")"
  else
    echo "  KEPT $(basename "$t") — merge incomplete, NOT deleting" >&2
    kept=$((kept + 1))
  fi
done

echo "harvest complete: $harvested tar(s) merged and deleted, $kept kept (unverified — investigate before re-running)"
