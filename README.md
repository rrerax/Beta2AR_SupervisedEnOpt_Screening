# β2-AR Supervised EnOpt Virtual Screening

This repository contains a reproducible computational screening workflow for the β2-adrenergic receptor (β2-AR, ADRB2). The workflow has two stages:

1. **Ensemble docking** — each compound is docked against five active and inactive β2-AR conformations, producing a per-conformation docking-score spectrum.
2. **Supervised EnOpt** — an XGBoost model is trained on the docking scores of known β2-AR actives to learn how to combine the score spectrum into a single ranking that prioritizes true binders.

The approach follows the EnOpt (Ensemble Optimizer) concept of Bhatt, Wang & Durrant, *Sci Rep* 14:20722 (2024), doi:10.1038/s41598-024-71699-3: ranking ensemble-docking scores is system-specific and can be improved by a dataset-specific machine-learning model.

## Project Overview

Classical high-throughput virtual screening often docks each ligand against a single receptor structure. For GPCR targets, this can be limiting because transmembrane helix rearrangements, especially around TM6 activation-associated motion, change the binding pocket geometry. This project therefore uses an ensemble of inactive and active β2-AR conformations, then combines docking results into a single ranking matrix.

Naive ways of combining the scores (simple mean or best score) are dominated by molecular size and perform worse than random for this screen (AUROC 0.41–0.47). The supervised EnOpt model instead learns which conformations matter for β2-AR binding and substantially improves ranking.

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

- Known actives' median rank improved from the top 39% to the top 22.6% of the library.
- Learned conformation importances: `active_bi167107_4LDE` 0.26, `active_gs_complex_3SN6` 0.21, `inactive_carazolol_2RH1` 0.20, `active_hbi_4LDL` 0.18, `inactive_carazolol_5D5A` 0.16.
- Deliverables: `results/supervised_enopt/` — trained model, model card, metrics, re-ranked leaderboard, and a visual report card.

## Repository Structure

```text
configs/                    Workflow configuration (screen + decoy screen)
data/training/              Known actives and DUD-E decoys for training/validation
scripts/                    Preparation, docking, analysis, and training scripts
notebooks/                  Analysis notebook for result review
docs/                       Method notes and final run summary
results/tables/             Docking scores, score matrix, weights, and top hits
results/supervised_enopt/   Trained model, metrics, ranking, and report card
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

The full run is CPU-intensive because it performs docking across 30,000 ligands and five receptor conformations. The supplied results were generated with parallel CPU docking.

## Notes and Limitations

- Docking scores are computational prioritization signals, not experimental binding affinities.
- The default negative set assumes that library molecules without reported β2-AR activity are inactive; a DUD-E decoy screen is prepared as an external negative control.
- The ranking should be interpreted as a shortlist for follow-up analysis, not as confirmed biological activity.
- Additional validation such as decoy screening, redocking controls, molecular dynamics, or experimental assays would be required before biological claims.
