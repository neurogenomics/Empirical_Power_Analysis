# shellcheck shell=bash
# Stage compose: multi-panel figures assembled from the cartesian outputs.
#   compose_figures.py --only figure_cartesian_oligo   main-text recall / precision figures
#   compose_recall_panels.R                            main-text recall panel figure
#   compose_crossscore.R (cells, reads)                Supp. cross-method grids
#   compose_supplementary.py (cells, reads)            Supp. per-method surfaces and panels
stage_main() {
  epa_init_ws
  epa_in cartesian_cells
  epa_in cartesian_reads
  epa_in crossscore
  local m
  epa_step "main-text cartesian figures"
  epa_py R/figure_composition/compose_figures.py --only figure_cartesian_oligo
  epa_R R/figure_composition/compose_recall_panels.R --cell Oligo --mode cells
  for m in cells reads; do
    epa_step "cross-score figures ($m)"
    epa_R R/figure_composition/compose_crossscore.R --mode "$m" --figset mega --style discrete
    epa_step "per-method supplementary figures ($m)"
    epa_py R/figure_composition/compose_supplementary.py --mode "$m"
  done
  local d=Composed_Figures n
  epa_require "$d/figure_cartesian_oligo_recall.png" "$d/figure_cartesian_oligo_precision.png" \
              "$d/figure_recall_panels_oligo_cells.png"
  for m in cells reads; do
    for ct in Astro Micro Oligo; do
      epa_require "$d/Cross_Score/crossscore_recall_mega_${ct}_${m}_discrete_raw.png" \
                  "$d/Cross_Score/crossscore_precision_mega_${ct}_${m}_discrete_raw.png"
    done
  done
  # compose_supplementary.py skips a figure (exit 0) when a source is missing.
  n="$(find "$WS/Datalink/$d/Supplementary" -name 'supp_*_bymethod_*.png' | wc -l)"
  [ "$n" -eq 36 ] || epa_die "expected 36 per-method supplementary figures, found $n"
  epa_out compose "$d"
  epa_complete
}
