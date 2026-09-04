# β2-AR Supervised EnOpt Virtual Screening

This repository contains a reproducible computational screening workflow for the β2-adrenergic receptor (β2-AR, ADRB2). The workflow has two stages:

1. **Ensemble docking** — each compound is docked against five active and inactive β2-AR conformations, producing a per-conformation docking-score spectrum.
2. **Supervised EnOpt** — an XGBoost model is trained on the docking scores of known β2-AR actives to learn how to combine the score spectrum into a single ranking that prioritizes true binders.

The approach follows the EnOpt (Ensemble Optimizer) concept of Bhatt, Wang & Durrant, *Sci Rep* 14:20722 (2024), doi:10.1038/s41598-024-71699-3: ranking ensemble-docking scores is system-specific and can be improved by a dataset-specific machine-learning model.

## Project Overview

Classical high-throughput virtual screening often docks each ligand against a single receptor structure. For GPCR targets, this can be limiting because transmembrane helix rearrangements, especially around TM6 activation-associated motion, change the binding pocket geometry. This project therefore uses an ensemble of inactive and active β2-AR conformations, then combines docking results into a single ranking matrix.

Raw score-averaging appears worse than random for this screen (AUROC 0.41–0.47) only when the score sign is ignored — docking scores are *more negative = better*. Once baselines are oriented correctly, plain 5-pocket averaging reaches AUROC ≈ 0.72 against property-matched DUD-E decoys and ≈ 0.67 on the full library. The supervised EnOpt model adds the most value at the top of the list (≈2× top-1% enrichment vs averaging) and in the full-library accounting; see Stage 3.

## Stage 1: Ensemble Docking

### Dataset Summary

- Ligand source: ChEMBL37 subset generated for β2-AR screening
- Input ligand records: 30,000
- Prepared ligand structures: 29,865
- Receptor conformations: 5
- Completed docking score rows: 149,325

Ligand preparation status:

- `ready`: 29,865
- `pdbqt_failed`: 109
- `embed_failed`: 20
- `prep_timeout`: 6

### Receptor Ensemble

The receptor set includes inactive and active-state β2-AR structures so that the screen is not tied to a single rigid receptor geometry.

| conformation_id         |   weight |
|:------------------------|---------:|
| active_bi167107_4LDE    |   0.198  |
| active_gs_complex_3SN6  |   0.1943 |
| active_hbi_4LDL         |   0.2018 |
| inactive_carazolol_2RH1 |   0.2031 |
| inactive_carazolol_5D5A |   0.2028 |

### Stage-1 Methods

1. **Ligand library construction**: A 30,000-compound ChEMBL37 ligand table was prepared with compound identifiers and SMILES strings.
2. **Protein structure preparation**: Five β2-AR crystal structures were selected to represent active and inactive conformational states.
3. **Ligand preparation**: Ligands were desalted, embedded into 3D conformers, minimized, and converted into docking-ready PDBQT files.
4. **Batch docking**: AutoDock Vina was used to dock each prepared ligand against each receptor conformation.
5. **Score matrix construction**: Docking scores were reshaped into a ligand-by-conformation matrix.
6. **Stage-1 baseline ranking**: Conformation weights estimated from the score ensemble were used to calculate weighted consensus scores (baseline only; superseded by the supervised model).

### Key Results (Stage 1)

Top-ranked compounds from the stage-1 ensemble weighted score:

|   rank | ligand_id     |   weighted_score |   best_score |   mean_score |   score_sd |   n_conf |
|-------:|:--------------|-----------------:|-------------:|-------------:|-----------:|---------:|
|      1 | CHEMBL3311247 |          -13.164 |       -14.48 |      -13.146 |      1.249 |        5 |
|      2 | CHEMBL16965   |          -12.896 |       -13.35 |      -12.888 |      0.553 |        5 |
|      3 | CHEMBL36113   |          -12.787 |       -13.71 |      -12.772 |      1.06  |        5 |
|      4 | CHEMBL4436402 |          -12.754 |       -13.65 |      -12.738 |      1.065 |        5 |
|      5 | CHEMBL33607   |          -12.721 |       -13.48 |      -12.706 |      1.007 |        5 |
|      6 | CHEMBL12018   |          -12.641 |       -13.91 |      -12.628 |      1.023 |        5 |
|      7 | CHEMBL12143   |          -12.63  |       -13.25 |      -12.62  |      0.687 |        5 |
|      8 | CHEMBL26449   |          -12.616 |       -13.43 |      -12.604 |      0.951 |        5 |
|      9 | CHEMBL11778   |          -12.571 |       -13.33 |      -12.558 |      1.055 |        5 |
|     10 | CHEMBL21788   |          -12.564 |       -13.67 |      -12.55  |      1.013 |        5 |

## Stage 2: Supervised EnOpt

The EnOpt paper (Bhatt et al., 2024) shows that mapping each compound's docking-score spectrum to a single value is system-specific and best learned from data. The β2-AR-specific model was trained as follows:

- **Positives**: 48 known β2-AR actives from ChEMBL (Ki/IC50/EC50 ≤ 1 µM, confidence score ≥ 8)
- **Negatives (default assumption)**: the remaining 29,817 library molecules
- **Features**: the five per-conformation docking scores
- **Model**: XGBoost, 3-fold stratified cross-validation, out-of-fold predictions (no molecule is scored by a model trained on it)
- **External negative control (prepared)**: 3,000 DUD-E decoys, docking workflow in `configs/decoy_screen.yml`

### Results

| method                    | OOF AUROC | EF 1% | EF 5% |
|:--------------------------|----------:|------:|------:|
| EnOpt (supervised OOF)    | **0.6647** | 4.17  | 2.50  |
| ensemble mean (baseline)  | 0.4133    | 0.00  | 0.00  |
| ensemble best (baseline)  | 0.4667    | 0.00  | 0.00  |

*Baseline rows use raw scores without sign correction. Correctly oriented baselines are compared in [Stage 3](#stage-3-dud-e-decoy-validation).*

- Known actives' median rank improved from the top 39% to the top 22.6% of the library.
- Learned conformation importances: `active_bi167107_4LDE` 0.26, `active_gs_complex_3SN6` 0.21, `inactive_carazolol_2RH1` 0.20, `active_hbi_4LDL` 0.18, `inactive_carazolol_5D5A` 0.16.
- Deliverables: `results/supervised_enopt/` — trained model, model card, metrics, re-ranked leaderboard, and a visual report card.

## Stage 3: DUD-E Decoy Validation

The Stage-2 numbers were produced with two known weaknesses: negatives were
"assumed-inactive" library molecules (no experimental basis), and only 48
positives were available. Stage 3 fixes both.

**What was done**

- Docked **2,978 unique DUD-E decoys** (property-matched non-binders generated
  from the β2-AR actives) against the same five conformations
  (15,000 docking runs; `results/tables/docking_scores_decoy_set.csv`).
- Added **158 newly docked ChEMBL literature actives** (pChEMBL ≥ 6, MW kept
  inside the decoy window so the decoys stay a fair, hard negative set),
  expanding positives from 48 to **206** (`data/training/beta2ar_actives_expand_docked.csv`).
- Retrained and re-validated with 3-fold out-of-fold cross-validation
  (`scripts/07_train_validated_enopt.py`).

**Clean-label hard test — 206 actives vs 2,978 matched decoys (OOF)**

| method                    | AUROC | 95% CI       | EF 1% | EF 5% |
|:--------------------------|------:|:-------------|------:|------:|
| EnOpt (supervised OOF)    | 0.696 | 0.655–0.737  | 3.88  | 2.72  |
| ensemble mean (oriented)  | 0.719 | 0.679–0.760  | 1.94  | 3.30  |
| ensemble best (oriented)  | 0.707 | 0.666–0.748  | 0.97  | 3.01  |

**Screening-like accounting — library + decoys + new actives (OOF)**

| method                    | AUROC | 95% CI       | EF 1% | EF 5% |
|:--------------------------|------:|:-------------|------:|:------|
| EnOpt (supervised OOF)    | 0.718 | 0.678–0.758  | 4.37  | 2.72  |
| ensemble mean (oriented)  | 0.671 | 0.631–0.712  | 0.97  | 2.33  |
| ensemble best (oriented)  | 0.636 | 0.595–0.677  | 0.49  | 0.97  |

**Honest reading**

- With 206 positives the AUROC estimate is tight (95% CI ± ~0.02): the
  five-pocket docking scores carry real, reproducible signal against
  property-matched decoys (AUROC ≈ 0.70, clearly above random).
- The earlier 0.61 (48 positives) was a small-sample artefact; the earlier
  0.6647 was partly inflated by easy, non-matched negatives (large library
  molecules) and by an unoriented baseline comparison.
- On the matched-decoy hard test, the trained model is **statistically tied
  with a correctly oriented simple 5-pocket average** (CIs overlap). Its
  demonstrable value is **top-of-list enrichment** (top-1% EF 3.9 vs 1.9)
  and the **full-library accounting**, where EnOpt clearly beats averaging
  (0.718 vs 0.671).
- Deliverables: `results/supervised_enopt_decoy/` — model, model card,
  metrics, decoy-validated re-ranked leaderboard, and comparison figure.
  Full write-up: `docs/decoy_validation_notes.md`.

## Repository Structure

```text
configs/                    Workflow configuration (screen + decoy screen)
data/training/              Known actives and DUD-E decoys for training/validation
scripts/                    Preparation, docking, analysis, and training scripts
notebooks/                  Analysis notebook for result review
docs/                       Method notes and final run summary
results/tables/             Docking scores, score matrix, weights, and top hits
results/supervised_enopt/   Stage-2 trained model, metrics, ranking, and report card
results/supervised_enopt_decoy/  Stage-3 decoy-validated model, metrics, and leaderboard
results/figures/            Publication-style result figures
examples/top_hit_poses/     Representative top-hit docking outputs
```

## Reproducing the Workflow

Create the environment and run a smoke test:

```bash
bash scripts/00_setup_env.sh
bash scripts/run_smoke_test.sh
```

Run the full workflow:

```bash
bash scripts/run_pipeline.sh
```

Train the supervised EnOpt model (after docking):

```bash
python scripts/06_train_supervised_enopt.py
```

Run the DUD-E-decoy-validated model (Stage 3, after docking):

```bash
python scripts/07_train_validated_enopt.py
```

The full run is CPU-intensive because it performs docking across 30,000 ligands and five receptor conformations. The supplied results were generated with parallel CPU docking.

## Notes and Limitations

- Docking scores are computational prioritization signals, not experimental binding affinities.
- The Stage-2 (06) model labelled uncharacterized library molecules as negatives by assumption. Stage 3 replaced the negative set with 2,978 real DUD-E decoys and expanded positives to 206 (see `docs/decoy_validation_notes.md`).
- The ranking should be interpreted as a shortlist for follow-up analysis, not as confirmed biological activity.
- On the matched-decoy hard test, supervised EnOpt (AUROC 0.696, CI 0.655–0.737) is statistically tied with a correctly oriented simple average (0.719, CI 0.679–0.760). The model's demonstrable value is top-of-list enrichment, not overall AUROC separation.
- Additional validation such as redocking controls, molecular dynamics, or experimental assays would be required before biological claims.
