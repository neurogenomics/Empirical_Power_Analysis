# Cell-type swing vs OOF model-error floor

Model: XGBoost `realistic` (shipped bundle) for the swing; winner config 25 OOF predictions for the error floor.
Real grid points: 3433289.

## Question
04b showed the cell-type ICE swing is large in the mid-power band (median ~0.39).
Is that a real cell-type effect, or just the model's own resolution -- in the
steep mid-band every input moves the prediction, so is cell type any more than
the model's generic wobble? Compare the swing against the model's OOF error at
the SAME points, paired and split by power region.

## Paired comparison

| region | n | err median | err p95 | swing median | swing p95 | swing/err median | % swing <= err |
|---|---|---|---|---|---|---|---|
| overall | 3433289 | 0.0023 | 0.1752 | 0.0047 | 0.5701 | 1.021 | 27.5% |
| low [0,0.1] | 2867941 | 0.0013 | 0.0722 | 0.0026 | 0.4983 | 0.945 | 25.9% |
| mid (0.1,0.9) | 343264 | 0.1403 | 0.4923 | 0.3066 | 0.7367 | 2.007 | 25.9% |
| high [0.9,1] | 222084 | 0.0039 | 0.0803 | 0.0047 | 0.1208 | 0.744 | 49.9% |

## Reading
- Mid band (the region you plan in): median swing 0.307 vs median OOF error 0.140, ratio 2.01. 25.9% of mid-band points have swing <= OOF error.
- If swing/err >> 1 and few points have swing <= err, the cell-type effect is
  RESOLVABLE: it moves the prediction by more than the model's typical error,
  so it is not just mid-band steepness. Cell type should be kept, not dropped.
- If swing/err ~ 1 with a large fraction swing <= err, the swing sits inside
  the model's own error floor: committing to a specific cell type argues below
  the model's resolution, and defaulting it (to Oligo, the conservative
  reference) is the honest move.

## Caveat (floor definition)
The OOF error e_i mixes model misspecification with irreducible binomial label
noise (obs is a detection rate over finite n_trials, max SE ~0.11 at 20 trials).
Per Chris's choice the floor here is the FULL OOF error (model error only, no
separate label-noise term), which is slightly generous to the 'within the
floor' side because label noise inflates e_i. A model-independent floor (the
detection-rate binomial SE) would be a stricter test if wanted.
