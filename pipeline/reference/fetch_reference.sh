#!/usr/bin/env bash
# =============================================================================
# pipeline/reference/fetch_reference.sh -- build the offline reference bundle.
# -----------------------------------------------------------------------------
# Usage: pipeline/reference/fetch_reference.sh WORKDIR/reference
#
# Cluster compute nodes have no internet access, so every external annotation
# the pipeline uses is fetched once, here, on a machine that does (a login
# node or a laptop), and then shipped into the work directory. The bundle is:
#
#   biomart/gene_audit_biomart_cache.csv   Ensembl gene -> chromosome, biotype, band
#   biomart/gene_biotype_cache.csv         Ensembl gene -> biotype, chromosome
#       (frozen biomaRt query results the paper was built with; vendored in
#        this repository under static/ because a live query would return a
#        newer Ensembl release)
#   nature.csl                             citation style (CC BY-SA 3.0, CSL project)
#   hgnc_complete_set_2026-09-04.txt       HGNC complete set, monthly archive
#   PTHR19.0_human                         PANTHER 19.0 sequence classification
#   DatabaseExtract_v_1.01.csv             Lambert et al. 2018 human TFs, v1.01
#   gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz  gnomAD v2.1.1 LoF constraint
#
# Every file is checked against MANIFEST.sha256; the script fails on any
# mismatch, so a changed upstream file is caught here rather than silently
# changing Figure 6. run.sh repeats the check before submitting anything.
# =============================================================================
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:?usage: fetch_reference.sh DEST_DIR}"
mkdir -p "$DEST/biomart"

cp "$HERE/static/biomart/gene_audit_biomart_cache.csv" "$DEST/biomart/"
cp "$HERE/static/biomart/gene_biotype_cache.csv"       "$DEST/biomart/"
cp "$HERE/static/nature.csl"                           "$DEST/"

fetch() {  # url file
  if [ -s "$DEST/$2" ]; then echo "have   $2"; return; fi
  echo "fetch  $2"
  curl -fsSL --retry 3 -o "$DEST/$2.part" "$1"
  mv "$DEST/$2.part" "$DEST/$2"
}
fetch https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-09-04.txt \
      hgnc_complete_set_2026-09-04.txt
fetch http://data.pantherdb.org/ftp/sequence_classifications/19.0/PANTHER_Sequence_Classification_files/PTHR19.0_human \
      PTHR19.0_human
fetch http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv \
      DatabaseExtract_v_1.01.csv
fetch https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

( cd "$DEST" && sha256sum -c "$HERE/MANIFEST.sha256" )
echo "reference bundle ready in $DEST"
