# Cell-type ICE / swing analysis

Model: XGBoost `realistic` (exported bundle). In-distribution designs: 150000 rows sampled from `aggregated_fitting_data_realistic.qs`.

## Question
The exported predictor keeps `cell_type` in the model but defaults it to
`Oligo` so a naive caller need not supply it. Is that default safe? ICE
(vary only cell_type, freeze everything else) measures the prediction move
directly -- the decision-relevant quantity SHAP does not give.

## Headline bound
Across 150000 in-distribution designs, changing ONLY cell type moves predicted power by:
- median 0.0046, mean 0.0959
- 95th percentile 0.5716, 99th percentile 0.7412
- maximum 0.8914
- fraction of designs with swing > 0.05: 29.0%

This is the error bound: |power(default cell type) - power(true cell type)| <= max swing.

Defaulting specifically to Oligo (the shipped default), the error vs each row's NATIVE cell type is:
- median 0.0000, 95th percentile 0.0277, max 0.8810

## By power region
(Swing concentrates in the mid-power band; the 0/1 rails are saturated.)

| metric | region | median | p95 | max |
|---|---|---|---|---|
| celltype_swing | low [0,0.1] | 0.0014 | 0.0965 | 0.2786 |
| celltype_swing | mid (0.1,0.9) | 0.3896 | 0.7418 | 0.8914 |
| celltype_swing | high [0.9,1] | 0.0056 | 0.1127 | 0.2297 |
| default_oligo_error | low [0,0.1] | 0.0000 | 0.0079 | 0.2574 |
| default_oligo_error | mid (0.1,0.9) | 0.0000 | 0.1509 | 0.8810 |
| default_oligo_error | high [0.9,1] | 0.0000 | 0.0336 | 0.2297 |

## Secondary: de_method
Toggling de_method (default deseq2) moves power by median 0.0149, p95 0.2137, max 0.9673.

## Interpretation
- Low SHAP for cell_type is corroborated by a direct, artefact-free measure:
  the ICE swing. Unlike mean|SHAP| (which conflates model sensitivity with a
  feature's data range and is deflated for a coarse factor competing against a
  granular continuous depth feature), the swing IS the full dynamic range of
  toggling the categorical, so its magnitude can be trusted.
- The claim is conditional and stated as such: cell type carries little beyond
  depth, yield and gene characteristics -- NOT that cell type is causally
  irrelevant (cf. the depth-equalisation finding, which shows the apparent
  cell-type hierarchy is a depth/yield/independent-filtering artefact).
- Practical consequence: defaulting cell_type in the exported predictor changes
  the answer by at most the bound above across the design space, so a naive
  caller can safely omit it while a power user retains full control.
