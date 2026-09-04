# Stage 4 notes: retrospective validation attempt + per-heavy-atom features

Two follow-up tracks were requested on top of the Stage-3 decoy-validated model
(`scripts/07`).  This note records what was tried, what was learned, and what
remains.

## Track A - retrospective top-200 validation against ChEMBL

Script: `scripts/08_retrospective_top200.py` (+ `--only-full`).
Outputs: `results/retrospective_top200/`.

### Finding 1: the chemreps library is almost free of beta-2 ground truth

- Pulled **all 5,686** ADRB2 (CHEMBL210) potency rows from ChEMBL (types
  IC50/Ki/EC50/KD/Kd), 3,638 distinct molecules; 1,608 with pChEMBL >= 6.
- Intersected with the 29,865-compound docking library: only **50** molecules
  are ChEMBL-documented potent beta-2 ligands, of which **48 are the training
  actives**; the remaining 29,817 "no-label" molecules contain just **2**
  documented strong actives (CHEMBL24 rank 17,364; CHEMBL776 rank 27,659) and
  24 weak-only records.
- Consequence: a top-200 (or even whole-library) ChEMBL hit-rate test has
  essentially no statistical power on this library.  Zero ChEMBL records in
  the top-200 means "these molecules were not previously measured on beta-2
  in ChEMBL" - compatible with genuinely new actives *or* false positives.
  It is **not** evidence for ranking quality, and equally not evidence
  against it.

### Finding 2: a concrete failure case (worth reporting honestly)

- CHEMBL776 (orciprenaline, classic beta-agonist, Kd ~ 500 nM, pChEMBL 6.30)
  is ranked 27,659/29,865, while its near-isostere CHEMBL434 (isoprenaline)
  is ranked 14.  Small flexible phenylethanolamines score with very high
  variance across the five conformations - one bad pose sinks a true active.
  This motivates Track C items (pose/redocking checks, higher exhaustiveness).

### Track A pivot: top-200 chemical review (still zero-cost, useful)

Script: `scripts/09_review_top200.py`.  Outputs: `results/review_top200/`.

- 200 compounds -> **143 Butina clusters** (ECFP4, Tanimoto 0.5) and **166
  unique Bemis-Murcko scaffolds**: the top-200 is chemically diverse, not a
  few chemotypes repeated.
- Novelty vs the 206 training actives: only **5** compounds have
  max-Tanimoto >= 0.6 to an active; **194** have max-Tanimoto < 0.4.  The
  leaderboard is mostly new chemistry, not trivial analogs.
- **90/200 flagged by PAINS/Brenk** - do not auto-discard.  Brenk flags
  phenols/catechols that are common genuine beta-agonist pharmacophores
  (isoprenaline itself is flagged).  Flagged hits need human review, not an
  automatic filter.
- Property figure: `fig_top200_properties_vs_actives.png`.

Practical use: pick purchase/assay candidates from diverse clusters, drop or
re-check PAINS-flagged compounds case by case, keep novel scaffolds.

## Track B - per-heavy-atom features (broke the 0.70 wall)

Script: `scripts/10_feature_heavy_atom_experiment.py`.
Outputs: `results/feature_experiment_heavy_atom/`.

Features: the five raw Vina scores plus the same five divided by the number
of heavy atoms (per-conformation), plus optionally the heavy-atom count.
Same 3-fold OOF protocol and XGBoost settings as `scripts/07`.

### Clean-label hard test (206 actives vs 2,978 DUD-E decoys)

| feature set            | AUROC | 95% CI      | EF 1% | EF 5% |
|------------------------|------:|:------------|------:|------:|
| raw only (5)           | 0.696 | 0.655-0.737 | 3.88  | 2.72  |
| **raw + per-ha (10)**  | **0.751** | **0.712-0.791** | **7.28** | **4.76** |
| raw + per-ha + HAC(11) | 0.749 | 0.709-0.788 | 6.80  | 4.18  |
| per-ha only (5)        | 0.637 | 0.595-0.679 | 2.91  | 1.85  |
| simple mean (raw)      | 0.719 | 0.679-0.760 | 1.94  | 3.30  |
| simple mean (per-ha)   | 0.535 | 0.494-0.577 | 0.97  | 0.58  |

- Adding per-heavy-atom scores as **extra** features raises the hard-test
  AUROC from 0.696 to **0.751** and more than doubles top-1% enrichment
  (EF1 3.9 -> 7.3).  CIs overlap only marginally.
- Using per-ha scores *instead of* raw scores (or a plain per-ha mean) is near
  chance: the model needs both absolute and size-normalized views.
- Feature importance is balanced across raw and per-ha features, led by
  `active_bi167107_4LDE` (raw).

### Merged-mirror accounting and leaderboard

| feature set (merged)    | AUROC | 95% CI      | EF 1% | EF 5% |
|-------------------------|------:|:------------|------:|------:|
| raw only                | 0.718 | 0.678-0.758 | 4.37  | 2.72  |
| **raw + per-ha**        | **0.735** | **0.696-0.774** | **4.85** | **3.79** |
| simple mean (raw)       | 0.671 | 0.631-0.712 | 0.97  | 2.33  |

- The raw+ha model re-ranks the library strongly: top-200 overlap with the
  Stage-3 leaderboard is only **24/200** (`library_ranking_heavy_atom_model.csv`).
- Size-bias check: Stage-3 top-200 had median MW 393 (heavier than the 206
  actives' 362); the raw+ha top-200 median MW is **358**, HAC 26 - inside the
  active window.  Per-heavy-atom features correct the "big molecules rank too
  high" artefact while *improving* decoy discrimination.
- Known-active recovery (sanity only): 0/50, 3/200, 6/500, 19/1000, 31/2000
  (Stage-3 was 2/50, 5/200, 8/500, 16/1000, 23/2000): slightly worse at the
  very top, better by rank 1000-2000.

### Honest reading

- The per-heavy-atom feature set is a clear, reproducible win and fixes a real
  artefact.  Recommended as the next leaderboard for candidate picking.
- Remaining Track-C items still need re-docking and cannot be evaluated from
  the current score tables: redocking/pose stability, higher exhaustiveness,
  flexible side chains, protein-ligand interaction fingerprints, and
  cross-checking the CHEMBL776-type small-flexible-molecule failures.

## Files produced

- `scripts/08_retrospective_top200.py`, `scripts/09_review_top200.py`,
  `scripts/10_feature_heavy_atom_experiment.py`
- `results/retrospective_top200/`, `results/review_top200/`,
  `results/feature_experiment_heavy_atom/`
