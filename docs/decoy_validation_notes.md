# DUD-E Decoy Validation (Stage 3)

Date: 2026-09-04
Script: `scripts/07_train_validated_enopt.py`
Deliverables: `results/supervised_enopt_decoy/`

## Why this stage exists

The Stage-2 model card (`results/supervised_enopt/`) had two weaknesses:

1. **Negatives were assumed, not measured.** The model labelled the 29,817
   library molecules that lack reported β2-AR activity as inactive. Most of
   them probably are, but the label is an assumption, and many are large
   molecules that are easy to separate from the (small) known actives by
   molecular size alone — which flatters the AUROC.
2. **Only 48 positives.** With 48 positives, an AUROC estimate of ~0.6 has a
   95% CI spanning roughly 0.47–0.75, so "above chance" could not be
   established reliably.

## What was done

| Step | Detail |
|:-----|:-------|
| True negatives | 3,000 DUD-E decoys generated from the β2-AR active list (property-matched non-binders; 2,978 unique after tautomer-form collapsing) and docked against the same five receptor conformations (15,000 runs). |
| Positive expansion | 158 additional ChEMBL literature actives (pChEMBL ≥ 6, human, direct assays) were selected with MW restricted to the decoy 5–95% window (250–516 Da) so the decoys remain a fair, hard negative set, then docked against the same five conformations. |
| Training set | 206 actives (48 in-library + 158 new) vs 2,978 decoys. |
| Validation | 3-fold stratified out-of-fold cross-validation (no molecule scored by a model trained on it); Hanley–McNeil 95% CIs for AUROC. |

Score tables: `results/tables/docking_scores_decoy_set.csv`
(15,790 rows = 15,000 decoy + 790 expanded-active docking runs).

## Results

### Clean-label hard test (206 actives vs 2,978 matched decoys)

| method                   | AUROC | 95% CI      | EF 1% | EF 5% |
|:-------------------------|------:|:------------|------:|------:|
| EnOpt (supervised OOF)   | 0.696 | 0.655–0.737 | 3.88  | 2.72  |
| ensemble mean (oriented) | 0.719 | 0.679–0.760 | 1.94  | 3.30  |
| ensemble best (oriented) | 0.707 | 0.666–0.748 | 0.97  | 3.01  |

### Screening-like accounting (library + decoys + new actives)

| method                   | AUROC | 95% CI      | EF 1% | EF 5% |
|:-------------------------|------:|:------------|------:|------:|
| EnOpt (supervised OOF)   | 0.718 | 0.678–0.758 | 4.37  | 2.72  |
| ensemble mean (oriented) | 0.671 | 0.631–0.712 | 0.97  | 2.33  |
| ensemble best (oriented) | 0.636 | 0.595–0.677 | 0.49  | 0.97  |

Baselines are sign-corrected (docking score more negative = better), which is
the comparison the Stage-2 README table did not make.

## Interpretation

- **The docking scores carry real signal.** At 206 positives the AUROC CI is
  ±0.02 and clearly above 0.5 (0.655–0.737). The earlier 0.61 at 48 positives
  was a small-sample artefact.
- **The learned model is statistically tied with a simple average on the hard
  test** (CIs overlap: 0.655–0.737 vs 0.679–0.760). What distinguishes EnOpt:
  - top-1% enrichment is ~2× the average baseline (3.9 vs 1.9) — the part of
    the ranking a screening campaign actually buys;
  - in the screening-like accounting EnOpt clearly beats averaging
    (0.718 vs 0.671).
- **The old 0.6647 overstates the model.** It was computed against easy,
  non-matched library negatives and against an unoriented baseline. Against
  matched decoys the honest number is ~0.70 — real, reproducible, but not a
  dramatic model advantage over a correctly used simple ensemble score.
- **Where the project goes from here:** the binding-level bottleneck is the
  information content of five rigid-receptor Vina scores, not the number of
  positives. Improvements would come from better features/protocol
  (pose-quality filters, protein–ligand interaction fingerprints, higher
  exhaustiveness, more receptor states), not from adding more molecules.

## Reproduce

```bash
python scripts/07_train_validated_enopt.py
```

Writes to `results/supervised_enopt_decoy/`:
`metrics_clean_label.csv`, `metrics_merged_mirror.csv`,
`feature_importance_clean.csv`, `library_ranking_decoy_model.csv`,
`enopt_model_decoy_validated.pkl`, `model_card_decoy_validated.json`,
`fig_decoy_validation.png`.

## Inputs

- `results/tables/docking_scores.csv` — 29,865-compound library docking scores.
- `results/tables/docking_scores_decoy_set.csv` — decoy + expanded-active scores.
- `data/training/beta2ar_actives_in_library.csv` — 48 in-library actives.
- `data/training/beta2ar_actives_expand_docked.csv` — 158 new actives (selected + docked).
- `configs/decoy_screen.yml`, `configs/actives_expand.yml` — docking configs used
  to produce the decoy and expanded-active score sets.
