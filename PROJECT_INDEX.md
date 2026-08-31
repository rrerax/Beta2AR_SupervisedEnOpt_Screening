# Project Index

This index lists the files needed to review or rerun the β2-AR supervised EnOpt screening project.

## Start Here

- `README.md` — project overview, methods, results, and reproduction commands.
- `notebooks/beta2ar_ensemble_screening_analysis.ipynb` — analysis notebook for the completed screen.
- `docs/final_run_summary.md` — concise final run summary and result counts.
- `docs/method_notes.md` — method rationale and interpretation notes.

## Main Results

- `results/tables/docking_scores.csv` — all completed docking scores (149,325 rows).
- `results/tables/enopt_style_score_matrix.csv` — ligand-level ensemble matrix (29,865 rows).
- `results/tables/enopt_style_top_hits.csv` — top 200 ranked compounds.
- `results/tables/conformation_weights.csv` — receptor conformation weights used in ranking.
- `results/figures/enopt_weighted_top_hits.png` — top-hit weighted score figure.
- `results/figures/score_distributions.png` — docking score distributions across receptor conformations.
- - `results/supervised_enopt/enopt_model.pkl` — trained supervised EnOpt model (XGBoost).
- `results/supervised_enopt/enopt_model_card.json` — model card: training data, features, validation metrics.
- `results/supervised_enopt/enopt_supervised_metrics.csv` — OOF AUROC and enrichment metrics.
- `results/supervised_enopt/enopt_supervised_ranking.csv` — full 30,000-compound re-ranked leaderboard.
- `results/supervised_enopt/enopt_report_card.html` — visual report card (raw vs supervised EnOpt).

## Input and Metadata

- `data/raw/chembl37_ligand_library_30000.csv` — input ligand library.
- `results/tables/ligand_manifest.csv` — ligand preparation status table.
- `results/tables/receptor_manifest.csv` — receptor preparation metadata.
- `configs/beta2ar_screen.yml` — workflow parameters.

## Scripts

- `scripts/make_ligand_library.py` — build the input ligand library.
- `scripts/01_check_inputs.py` — validate project inputs.
- `scripts/02_fetch_receptors.py` — fetch and prepare receptor structures.
- `scripts/03_prepare_ligands.py` — prepare ligand structures.
- `scripts/04_run_vina.py` — run batch AutoDock Vina docking.
- `scripts/05_analyze_enopt.py` — build the ensemble matrix, weights, ranking table, and figures.
- `scripts/06_train_supervised_enopt.py
- `scripts/run_smoke_test.sh` — quick validation run.
- `scripts/run_pipeline.sh` — full pipeline entry point.
