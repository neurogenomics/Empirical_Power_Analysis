# Reference bundle

Cluster compute nodes have no internet, so all external annotation is shipped
into the work directory as `WORK/reference/`, built once with

```bash
pipeline/reference/fetch_reference.sh WORK/reference
```

| File | Source | Used by |
|---|---|---|
| `biomart/gene_audit_biomart_cache.csv` | frozen biomaRt query (Ensembl gene -> chromosome, biotype, band), vendored in `static/` | gene-universe audit (Supp. Data 1), Supp. Table 2 |
| `biomart/gene_biotype_cache.csv` | frozen biomaRt query (gene -> biotype, chromosome), vendored in `static/` | DEG chromosome grid (Supp. Fig.) |
| `nature.csl` | Citation Style Language project, `zotero.org/styles/nature` (CC BY-SA 3.0), vendored in `static/` | knitting |
| `hgnc_complete_set_2026-09-04.txt` | HGNC monthly archive | Figure 6 (gene symbols, PANTHER join) |
| `PTHR19.0_human` | PANTHER 19.0 sequence classification | Figure 6 (protein classes) |
| `DatabaseExtract_v_1.01.csv` | Lambert et al. 2018, humantfs.ccbr.utoronto.ca, v1.01 | Figure 6 (TF families) |
| `gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` | gnomAD v2.1.1 constraint | Figure 6 (LoF-intolerant genes) |

`MANIFEST.sha256` pins every file; `run.sh` refuses to run if the bundle does
not match. The two biomaRt caches are the exact annotation the paper was built
with (a live query would return a newer Ensembl release), so the R code is run
with `BIOMART_OFFLINE=true` and never queries Ensembl.
