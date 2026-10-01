# shellcheck shell=bash
# Stage finalise: write results/MANIFEST.sha256 (every result tar) and
# results/PROVENANCE.txt (code commit, container image, settings, inputs, and
# when each stage finished). Nothing here is recomputed.
stage_main() {
  local res="$EPA_WORK/results"
  ( cd "$res" && find . -name '*.tar' -print | LC_ALL=C sort | sed 's#^\./##' \
      | while read -r t; do cut -d' ' -f1 < "$t.sha256" | sed "s#\$#  $t#"; done ) > "$res/MANIFEST.sha256.tmp"
  mv -f "$res/MANIFEST.sha256.tmp" "$res/MANIFEST.sha256"
  {
    echo "Empirical Power Analysis: reproduction provenance"
    echo "generated : $(date -u +%FT%TZ) on $(hostname)"
    echo "code      : $(cat "$EPA_WORK/code/COMMIT")"
    echo "image     : $EPA_SIF"
    echo "image sha : $(cat "$EPA_WORK/code/SIF_SHA256")"
    echo
    echo "== settings (epa.env at first run) =="
    cat "$EPA_WORK/code/epa.env.at-snapshot"
    echo
    echo "== raw inputs (verified by stage shard) =="
    tar -xOf "$EPA_WORK/code/code.tar" pipeline/inputs.sha256
    echo
    echo "== stages =="
    for c in "$res"/*/.complete; do
      printf '%-18s %s\n' "$(basename "$(dirname "$c")")" "$(grep '^finished=' "$c" | cut -d= -f2)"
    done
  } > "$res/PROVENANCE.txt"
  epa_log "wrote $res/MANIFEST.sha256 and $res/PROVENANCE.txt"
  mkdir -p "$res/finalise"
  EPA_STAGE=finalise epa_complete
}
