# shellcheck shell=bash
# Stage gene_class: Figure 6, predicted detectability of every oligodendrocyte gene.
#   classify_oligo_genes.R          full-cohort DESeq2 fit -> per-gene baseMean and dispersion
#   annotate_oligo_genes.R          HGNC symbols, PANTHER classes, Lambert TF families
#   classify_oligo_multi_lfc.R      power at |LFC| 0.2 / 1 / 2 (50 donors, 1,000 cells) + gnomAD pLI
#   plot_oligo_expr_disp_panels.py  the six-panel figure
stage_main() {
  epa_init_ws
  epa_in shard shards_Roussos
  epa_in xgb_post
  epa_ref_dir
  epa_env EPA_REFERENCE_DIR /workspace/reference
  epa_env CLASSIFY_CORES "$EPA_NCPUS"
  local g=Gene_Classification
  epa_step "per-gene features and classification"
  epa_R R/gene_classification/classify_oligo_genes.R
  epa_step "annotation"
  epa_R R/gene_classification/annotate_oligo_genes.R
  epa_step "multi-LFC power"
  epa_R R/gene_classification/classify_oligo_multi_lfc.R
  epa_step "Figure 6"
  epa_py R/gene_classification/plot_oligo_expr_disp_panels.py
  epa_require "$g/oligo_paper_expr_disp_panels.png" "$g/oligo_gene_power_multiLFC.csv" \
              "$g/oligo_gene_group_power_summary.csv"
  epa_out gene_class "$g/*.csv" "$g/oligo_paper_expr_disp_panels.png" "$g/oligo_paper_expr_disp_panels.pdf"
  epa_complete
}
