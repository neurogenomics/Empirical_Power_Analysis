# Oligodendrocyte gene detectability by functional group — explainer

**Question this answers:** at a fixed, realistic single-cell design, which kinds of genes are detectable, and are transcription factors (TFs) really lower-expressed than most other groups?

**Setup (fixed for every gene / group):**
- Cell type: oligodendrocytes (Oligo/Roussos, PsychAD), full transcriptome (34,811 genes).
- Design: **100 donors, 250 cells/donor, 5,742 reads/cell**, effect size **|LFC| = 1.5** for all.
- Model: the realistic XGBoost power model → predicted power → **skewed trinary**: power ≤ 0.4 = **Undetectable**, 0.4–0.8 = **Unsure**, > 0.8 = **Detectable**.
- Per-gene inputs: raw pseudobulk **basemean** + **DESeq2 dispersion** (full-cohort sex-design fit — exactly the model's training features).

---

## Part 1 — the per-gene TF result, and why the headline was an artefact

The first pass compared the curated **Lambert TF list vs everything else** and found TFs looked *more* detectable (≈79% vs ≈54%). That gap was **not real biology** — it came from *what* the "everything else" bucket contains:

| non-TF biotype | n | % detectable | median basemean |
|---|---|---|---|
| protein-coding | 15,370 | 78.9 | 90.0 |
| non-coding RNA | 5,199 | 44.7 | 7.2 |
| unannotated (no HGNC) | 9,168 | 30.3 | 3.4 |
| pseudogene | 3,517 | 16.8 | 1.6 |

The Lambert list is 1,400 clean protein-coding regulators; the background was half-full of pseudogenes / lncRNA / unannotated loci that were never going to be detected. **Restricting to protein-coding, the TF "advantage" vanishes: 79.2% (TF) vs 78.9% (non-TF).** So the raw gene-count comparison was never like-for-like.

---

## Part 2 — the fair test: compare GROUP medians (junk cannot dilute a median)

To bypass junk without inventing an arbitrary expression cutoff, we group the **expressed** genes by function, take the **median basemean + median dispersion** per group, and push that single representative entry through the model. Central tendencies, not gene counts.

### TF vs non-TF, on group medians

| Group | n | median basemean | median dispersion | representative power |
|---|---|---|---|---|
| Protein-coding non-TF | 15,370 | **90.0** | 0.081 | 0.991 (Detectable) |
| Lambert TF | 1,402 | **79.0** | 0.056 | 0.993 (Detectable) |
| Other / non-coding | 18,038 | 3.6 | 0.327 | 0.464 (Unsure) |

**Your prior is correct, but the effect is tiny.** TFs sit marginally *below* protein-coding non-TFs on expression (79 vs 90), but they also have *lower* dispersion (0.056 vs 0.081), which cancels the expression gap out — so representative power is effectively identical (0.993 vs 0.991). At |LFC| = 1.5 both groups are far above the ~15-basemean floor where power collapses, so "a bit lower" doesn't translate into "less detectable".

### The real story: the "TF" label is heterogeneous

Splitting TFs by DNA-binding domain (Lambert DBD family) shows the aggregate hides a huge spread:

| TF DBD family | n | median basemean | representative power | class |
|---|---|---|---|---|
| C2H2 ZF | 674 | 103.6 | 0.997 | Detectable |
| HMG/Sox | 26 | 96.0 | 0.989 | Detectable |
| bZIP | 51 | 62.2 | 0.984 | Detectable |
| Nuclear receptor | 46 | 57.1 | 0.978 | Detectable |
| bHLH | 88 | 38.5 | 0.972 | Detectable |
| Ets | 25 | 31.8 | 0.964 | Detectable |
| **Homeodomain** | 103 | **3.8** | **0.576** | **Unsure** |
| **Forkhead** | 37 | **2.4** | **0.246** | **Undetectable** |
| **T-box** | 15 | **1.9** | **0.258** | **Undetectable** |

The **C2H2 zinc-fingers (674 genes — nearly half the TF set)** are broadly expressed, housekeeping-like (KRAB-ZNFs, transposon silencing) and pull the aggregate TF rate up. The **classic low-abundance master regulators — homeodomain, forkhead, T-box — are exactly the poorly-detectable tail you expected.** Both intuitions coexist: master regulators *are* lowly expressed and hard to detect; they're just numerically outnumbered by ZNFs under the flat "TF" label.

### Where genes actually sit (PANTHER protein class, 118 classes)

Ordered by median basemean, the transcriptome runs from housekeeping machinery at the top to secreted/signalling at the bottom:

**Most detectable (high expression):** mRNA polyadenylation factor (382), damaged DNA-binding (343), **zinc-finger TF (316)**, aminoacyl-tRNA synthetase (285), vesicle coat (280), ATP synthase (267), translation initiation (243), chromatin regulator (234), ubiquitin ligase (223), ribosomal protein (198), RNA splicing factor (173).

**Least detectable (low expression, the true low tail):** intermediate filament (1.8, Undetectable), **cytokine (2.9, Unsure)**, MHC protein (3.6), peptide hormone (3.7), immunoglobulin receptor (3.8), **GPCR (8.0)**, growth factor (8.3), **homeodomain TF (8.4)**, intercellular signal molecule (10.7), transmembrane signal receptor (11.7).

**Takeaway:** the genes that are actually hard to detect are the **secreted and cell-surface signalling classes** — cytokines, growth factors, peptide hormones, GPCRs, MHC, immunoglobulin receptors — plus the developmental TF subtypes (homeodomain/forkhead/T-box). These sit *below* the bulk of TFs. So "TFs are lowly expressed" is true only for the master-regulator subtypes; as a whole, TFs are mid-pack, slightly below general protein-coding genes, and dominated by broadly-expressed zinc-fingers.

---

## Files (`Datalink/Gene_Classification/`)
- `oligo_gene_classification_annotated.csv` — per-gene power + trinary class + PANTHER + Lambert.
- `oligo_group_comparison_panther.csv` — 118 PANTHER classes: median/IQR basemean & dispersion, representative power + class, per-gene % detectable.
- `oligo_group_comparison_dbd.csv` — Lambert TF DBD families (the table above).
- `oligo_group_comparison_tf.csv` — TF vs protein-coding non-TF vs other.

**Caveat:** group medians are marginal (median basemean and median dispersion taken independently), so the representative row is a nominal central gene, not a real one. Good for ranking groups; not for a specific gene (use the per-gene table for that).
