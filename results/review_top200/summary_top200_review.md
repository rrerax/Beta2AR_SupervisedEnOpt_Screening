# Top-200 chemical review

Ready-for-review summary of the Stage-3 leaderboard's top 200
(molecules that the model was never given as actives unless flagged).

| metric | value |
|---|--:|
| compounds analysed | 200 |
| Butina clusters (ECFP4, Tanimoto 0.5) | 143 |
| unique Bemis-Murcko scaffolds | 166 |
| close analogs of a training active (maxT >= 0.6) | 5 |
| structurally novel vs actives (maxT < 0.4) | 194 |
| PAINS/Brenk flagged (aggregator risk) | 90 |

Largest clusters: {0: 8, 1: 6, 2: 5, 4: 5, 5: 4}

Full per-compound data: `top200_chemical_review.csv`. Property figure:
`fig_top200_properties_vs_actives.png`.

Reading: if top-200 collapses into very few clusters, "200 hits" really means
a few chemotypes -- diversity (cluster count) should be inspected before
picking purchase/assay candidates.
