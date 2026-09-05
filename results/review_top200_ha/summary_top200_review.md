# Top-200 chemical review

Ready-for-review summary of the enopt_ha_rank leaderboard's top 200
(molecules that the model was never given as actives unless flagged).

| metric | value |
|---|--:|
| compounds analysed | 200 |
| Butina clusters (ECFP4, Tanimoto 0.5) | 167 |
| unique Bemis-Murcko scaffolds | 170 |
| close analogs of a training active (maxT >= 0.6) | 6 |
| structurally novel vs actives (maxT < 0.4) | 192 |
| PAINS/Brenk flagged (aggregator risk) | 75 |

Largest clusters: {0: 5, 1: 4, 2: 3, 3: 3, 4: 3}

Full per-compound data: `top200_chemical_review.csv`. Property figure:
`fig_top200_properties_vs_actives.png`.

Reading: if top-200 collapses into very few clusters, "200 hits" really means
a few chemotypes -- diversity (cluster count) should be inspected before
picking purchase/assay candidates.
