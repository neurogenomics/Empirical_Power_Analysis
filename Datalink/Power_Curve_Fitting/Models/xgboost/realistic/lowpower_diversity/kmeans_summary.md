# Low-power mass k-means (K=4)

Low-power mass = rows with observed detection rate <= 0.10
N low-power rows: 2910179 (of 3433289 filtered)

## Archetypes (auto-labelled from standardised cluster centres)
- Cluster 1 (1057998 rows, 36.4%): high expression
- Cluster 2 (320984 rows, 11.0%): many donors, few cells
- Cluster 3 (59828 rows, 2.1%): high |LFC|, low expression, high dispersion
- Cluster 4 (1471369 rows, 50.6%): low expression

Interpretation: if archetypes span several failure modes (small effect,
low expression, few donors) the low-power mass is DIVERSE and reweighting
has more room to help. If all archetypes concentrate on one axis it is
HOMOGENEOUS and reweighting is fighting an easy battle.
