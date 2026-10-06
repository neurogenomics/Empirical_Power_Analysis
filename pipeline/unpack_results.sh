#!/usr/bin/env bash
# =============================================================================
# pipeline/unpack_results.sh -- turn a finished run's result tars into a
# Datalink/ tree. NOT part of the pipeline: run it wherever you want the files
# (e.g. after copying WORKDIR/results/ off the cluster with Globus).
# -----------------------------------------------------------------------------
# Usage: pipeline/unpack_results.sh RESULTS_DIR DEST [--with-grid] [--no-verify]
#
#   RESULTS_DIR  WORKDIR/results (or a copy of it)
#   DEST         directory to create Datalink/ in (e.g. a clone of the repository)
#   --with-grid  also unpack the raw per-iteration grid results of the sex and
#                AD grids (~400 GB, over a million small files); not needed for
#                any figure or table
# The per-task XGBoost search tars (xgb_search/) are never unpacked: the same
# files ship in xgb_select/search.tar.
#   --no-verify  skip the sha256 check against RESULTS_DIR/MANIFEST.sha256
#
# Tars are unpacked in stage order; no two stages write the same file, so the
# order only matters for readability of the log.
# =============================================================================
set -euo pipefail
RES="${1:?usage: unpack_results.sh RESULTS_DIR DEST [--with-grid] [--no-verify]}"
DEST="${2:?usage: unpack_results.sh RESULTS_DIR DEST [--with-grid] [--no-verify]}"
shift 2
WITH_GRID=false VERIFY=true
for a in "$@"; do
  case "$a" in
    --with-grid) WITH_GRID=true ;;
    --no-verify) VERIFY=false ;;
    *) echo "unknown option $a" >&2; exit 2 ;;
  esac
done
RES="$(cd -- "$RES" && pwd)"
[ -f "$RES/MANIFEST.sha256" ] || { echo "no MANIFEST.sha256 in $RES (did the finalise stage run?)" >&2; exit 1; }
mkdir -p "$DEST/Datalink"

if $VERIFY; then
  echo "verifying checksums ..."
  ( cd "$RES" && if $WITH_GRID; then sha256sum -c --quiet MANIFEST.sha256;
                 else grep -vE ' (grid|ad_grid|xgb_search)/' MANIFEST.sha256 | sha256sum -c --quiet; fi )
fi

order="shard eda fgt deg_loss gt_manifest cartesian_cells cartesian_reads crossscore pcf_aggregate
       xgb_select xgb_shap xgb_post gene_class tables compose supp_data knit
       ad_shard ad_gt_manifest ad_cartesian ad_pcf ad_de ad_validate"
$WITH_GRID && order="$order grid ad_grid"
for s in $order; do
  for t in "$RES/$s"/*.tar; do
    [ -e "$t" ] || continue
    echo "unpack $s/$(basename "$t")"
    tar -xf "$t" -C "$DEST/Datalink"
  done
done
cp "$RES/PROVENANCE.txt" "$RES/MANIFEST.sha256" "$DEST/Datalink/" 2>/dev/null || true
echo "done: $DEST/Datalink"
