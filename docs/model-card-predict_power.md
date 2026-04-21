# Model Card: `predict_power.R`

This model card documents the fitted logistic regression model shipped as
`Datalink/Power_Curve_Fitting/Exports/predict_power.R`. It follows the
[Model Cards for Model Reporting](https://arxiv.org/abs/1810.03993) framework
adapted for biomedical predictive models.

---

## Model details

- **Model type**: Logistic regression (generalised linear model, binomial family,
  logit link)
- **Formula**: `detection ~ log_donors + log_n_cells + log_reads_per_cell +
  log_abs_gt_lfc + log_gt_basemean + log_gt_dispersion`
- **Coefficients**: six main effects; no interactions (see Limitations)
- **Training code**: `R/power_curve_fitting/02_Fit_Models.R`
- **Fitted coefficients**: `Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv`
- **Pre-computed grid**: `Datalink/Power_Curve_Fitting/Exports/power_lookup_table.csv`
- **Version**: 1.0.0 (bioRxiv release)
- **Citation**: Crampton et al. 2026, bioRxiv (DOI on release)

## Intended use

### Primary use
- **Power planning** for pseudobulk DE analyses in snRNA-seq studies of human
  brain tissue, glial cell types, with sex-biased or similarly-distributed case-
  control contrasts.
- **Benchmark comparator** against parametric power calculators
  (scPower, POWSC, powsimR, hierarchicell).

### Out-of-scope uses
- **Neurons** — no neuronal cell types in the training grid
- **Non-brain tissues** — peripheral immune, liver, gut, tumour not evaluated
- **Case-control contrasts with effect-size distributions very different from sex-biased DE** (see Limitations)
- **Single-cell (non-pseudobulk) DE methods** (e.g. NEBULA, dreamlet, glmmTMB)
- **Effect sizes far outside training range**: model is unreliable for
  |LFC| > 10 or very rare genes (baseMean < 10)

---

## Training data

- **Source**: PsychAD cohort (Fullard et al. 2025, *Scientific Data*)
- **Size**: 1,214–1,303 donors depending on cell type (astrocytes, microglia,
  oligodendrocytes after 30-cell minimum filter)
- **Number of per-gene binary detection outcomes**: ≈22 million
- **Grid dimensions**:
  - Donors: {5, 10, 15, 25, 50, 75, 100, 150, 200, 300, 500, 600}
  - Cells per donor: {5, 10, 20, 50, 75, 100, 250, 500, 750, 1000, 1500, 2000, 3000, Max}
  - Reads per donor: {25k, 50k, 100k, 250k, ..., 15M, Max}
  - Robustness runs: 6 (each with an independent ground-truth/validation split)
  - Subsamples per grid point: 20
- **DE tool**: DESeq2 v1.40 with apeglm shrinkage; BH FDR 0.05
- **Ground truth**: sex-biased DE genes in held-out 50% of donors

## Evaluation data

- **Internal hold-out**: 50/50 donor split per robustness run
- **Cross-dataset validation**: IoU against Allen/SEA-AD, DeJager/ROSMAP,
  Rexach, Johansen
- **External anchor**: GTEx v8 bulk sex-biased DE concordance

## Quantitative performance

- **Pseudo-R² (McFadden)**: 0.655 – 0.733 across per-cell-type model variants
- **Pseudo-R² caveat**: At N = 22M observations, pseudo-R² is inflated because
  test statistics are dominated by sample size rather than effect strength.
  **Not a reliable indicator of calibration.** See `Datalink/Power_Curve_Fitting/Plots/Cross_Validation/`.
- **Cross-validation**: leave-one-cell-type-out and per-iteration CV available
  in pipeline output
- **Known calibration gap**: model slightly underestimates power for tiny-effect
  genes at very high cell counts (n_cells > 2000) due to the log-linear
  specification

## Ethical considerations

- **Sex-biased training data**: The model's training ground truth is dominated
  by sex-chromosome genes (XIST, DDX3Y, KDM6C/UTX). Users applying the model
  to contrasts without this dispersion/expression structure (e.g. autosomal
  disease contrasts) should expect meaningful miscalibration.
- **Demographic scope**: PsychAD is US-based, predominantly European-ancestry;
  generalisation to other populations is untested.
- **Not a substitute for statistical power analysis tailored to a specific
  study design.** Treat as a decision-support tool, not an authoritative
  power calculator.

## Limitations

1. **No interaction terms.** The logistic fit assumes additive log-linear
   effects. Reality is known to have interactions (e.g. more cells compensate
   for fewer donors up to a threshold).
2. **Negative reads-per-cell coefficient in some fits.** This is an
   under-fitting artefact of the additive specification on a correlated grid;
   do not interpret it causally.
3. **Trained on one dataset.** Cross-dataset IoU validates the sex-biased
   DEG set; it does not validate the power surfaces themselves on non-PsychAD
   data.
4. **Depth-equalisation one-directional.** Astrocytes were downsampled to
   microglial depth; the reverse was not tested.
5. **DESeq2-specific.** Results may differ under edgeR, glmGamPoi, or
   mixed-model methods. See `docs/architecture.md` for design rationale.

## Caveats for end users

If you use `predict_power()` for grant planning:

- Report the training context (glial, brain cortex, sex contrast, DESeq2)
  when citing the number in a methods section.
- Prefer the lookup table over the closed-form function near grid boundaries.
- For contrasts with very different effect-size distributions, consider
  re-fitting the pipeline on your own pilot data — the framework is
  designed to support this.

## Updates & provenance

This model card is versioned alongside the pipeline in git.  Material changes
(refits, re-training, specification changes) require a corresponding update
here.  See `CHANGELOG.md` for release history.
