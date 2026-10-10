# shellcheck shell=bash
# Stage knit: render the manuscript and the supplement from the reproduced
# Datalink/ (HTML, as authored). The citation style is the shipped
# reference/nature.csl, passed to pandoc so no network access is needed; it
# overrides the URL in the documents' YAML.
stage_main() {
  epa_init_ws
  local s
  for s in eda fgt deg_loss cartesian_cells cartesian_reads compose xgb_shap xgb_post gene_class tables supp_data paper_values; do
    epa_in "$s"
  done
  epa_in xgb_select model
  epa_ref nature.csl ../nature.csl
  local doc
  for doc in nature_comms_draft nature_comms_supplementary; do
    epa_step "render $doc.rmd"
    epa_exec Rscript -e "rmarkdown::render('$doc.rmd', output_format = 'bookdown::html_document2', output_options = list(pandoc_args = c('--csl', 'nature.csl')), quiet = TRUE)"
    [ -f "$WS/$doc.html" ] || epa_die "$doc.html was not produced"
  done
  mkdir -p "$WS/Datalink/Manuscript"
  cp "$WS/nature_comms_draft.html" "$WS/nature_comms_supplementary.html" "$WS/Datalink/Manuscript/"
  epa_out manuscript Manuscript
  epa_complete
}
