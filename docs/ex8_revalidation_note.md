# Stage 5 note: exhaustiveness-8 re-docking validation (ex8) and leaderboard robustness

Status: **complete (2026-09-04/05)**. This is the follow-up to the "remaining
Track-C" items listed in `stage4_retrospective_and_heavy_atom_notes.md`:
re-dock the validation set and the top-200 shortlist at higher exhaustiveness
to test whether the ex4 results were limited by docking sampling quality.

## What was run

Two independent ex8 (exhaustiveness = 8, num_modes = 20, seed-fixed) re-dock jobs:

1. **Validation set re-dock (cloud, 24 cores)**
   - 3,184 ligands (206 actives + 2,978 DUD-E decoys) x 5 pockets = 15,920 docks.
   - Engine: AutoDock Vina 1.2.5; Open Babel 3.1.0; rdkit (2026.x); venv Python.
   - Score table: `results/ex8_validation/docking_scores_ex8_validation.csv`.
   - Evaluation report: `results/ex8_validation/eval_ex8_report.txt`.

2. **Shortlist re-dock (local, 12 cores)**
   - Top-200 of the raw+ha leaderboard + CHEMBL776 (known failure case) = 201
     ligands x 5 pockets = 1,005 docks, same engine/version.
   - Score table: `results/ex8_rerank/scores_ex8_top200_plus776.csv`.

Rerun entry points live in `revalidation_ex8/`
(`run_stage5_redock.sh`, `make_validation_inputs.py`, `eval_ex8.py`); the
shortlist ex8 re-dock driver is `scripts/11_ex4_ex8_rerank.py` (analysis only —
the dock itself was run with the repository workflow at ex8).

## Results

### 1) Validation-set discrimination (clean-label 3-fold OOF, XGBoost)

| model | ex4 (reference) | ex8 | interpretation |
|---|---:|---:|---|
| raw (5 scores) | 0.696 (0.655-0.737) | 0.654 (0.612-0.695) | CIs overlap -> no change |
| raw + per-heavy-atom (10) | 0.751 (0.712-0.791) | 0.733 (0.692-0.772) | CIs overlap -> no change |

- EF1% ex8: raw 3.40, raw+ha 6.31 (ex4: 3.88 / 7.28).
- Reading: higher sampling quality neither rescued nor degraded
  discrimination. The ex4 conclusion is not an exhaustiveness artefact, and
  increasing exhaustiveness is not the path to higher AUROC.

### 2) Shortlist ranking robustness (ex4 vs ex8, 201 ligands)

| metric | value |
|---|---:|
| Spearman r, best score | 0.874 |
| Spearman r, best / heavy-atom | 0.891 |
| Spearman r, mean score | 0.905 |
| median delta best score (kcal/mol) | +0.00 |
| ligands moving > 1 kcal/mol | 3 / 201 |
| ex4 top-100 still in ex8 top-100 (best) | 89/100 |
| per-pocket Spearman r (5 pockets) | 0.77-0.86 |

- The leaderboard is stable to sampling quality at the top; the ~11% mid-list
  churn argues for re-scoring finalists at ex8 before purchase/assay.
- CHEMBL776 (orciprenaline) remains rank 201/201 at ex8: its failure is a
  systematic rigid-receptor under-scoring problem for small flexible
  phenylethanolamines, not a low-exhaustiveness artefact. This motivates
  interaction-fingerprint and flexible-sidechain features over more sampling.

## Files

- `results/ex8_validation/` - ex8 validation score table + eval report.
- `results/ex8_rerank/` - per-ligand ex4-vs-ex8 comparison, report, figure.
- `results/review_top200_ha/` - chemical review of the **raw+ha** leaderboard
  top-200 (the recommended candidate list; 167 Butina clusters, 170 scaffolds,
  75 PAINS/Brenk flagged - inspect before purchase).
- `data/training/validation_redock_ligands.csv`, `configs/validation_redock.yml`
  - merged validation inputs / ex8 config.
- `revalidation_ex8/` - rerun scripts for the validation-set ex8 job.
- `scripts/11_ex4_ex8_rerank.py` - reruns the shortlist robustness analysis.

## Caveats

- Validation ex8 numbers are single-seed runs; CIs overlap with ex4, so the
  conclusion is "no significant change", not "equivalence proven".
- The top-200 ex8 re-dock used the repository workflow with num_modes 20; pose
  ensemble / IFP analysis of these poses is not yet in the repo.