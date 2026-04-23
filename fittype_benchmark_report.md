# fitType Equivalence + Speedup Benchmark
Date: 2026-04-23 13:50:31 UTC
R: R version 4.5.1 (2025-06-13) | DESeq2: 1.50.2 | glmGamPoi: 1.22.0 | apeglm: 1.32.0

## Timing
| Config | Wall (s) |
|---|---|
| fitType = parametric  | 3822.69 |
| fitType = glmGamPoi   | 3091.45 |
| **Speedup**           | **1.24x** |

## Equivalence
| Metric | Value | Threshold | Pass |
|---|---|---|---|
| Pearson LFC        | 1.0000 | >= 0.99  | PASS |
| Spearman -log10(padj) | 0.9971 | >= 0.99 | PASS |
| Jaccard sig (padj<0.05) | 0.9393 | >= 0.95 | FAIL |
| Sign agreement on sig | 1.0000 | >= 0.99 | PASS |
| Kendall top-100 | 0.9867 | >= 0.95 | PASS |
| No sig-set direction flips | TRUE | TRUE | PASS |

## Sig gene counts
- Parametric: 2186 significant genes
- glmGamPoi: 2315 significant genes
- Shared significant: 2180

