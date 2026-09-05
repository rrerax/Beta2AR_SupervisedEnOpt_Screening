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
- `results/supervised_enopt/enopt_model.pkl` — trained supervised EnOpt model (XGBoost).
- `results/supervised_enopt/enopt_model_card.json` — model card: training data, features, validation metrics.
- `results/supervised_enopt/enopt_supervised_metrics.csv` — OOF AUROC and enrichment metrics.
- `results/supervised_enopt/enopt_supervised_ranking.csv` — full 30,000-compound re-ranked leaderboard.
- `results/supervised_enopt/enopt_report_card.html` — visual report card (raw vs supervised EnOpt).
- `results/supervised_enopt_decoy/` — Stage-3 DUD-E-decoy-validated model, metrics (with 95% CI), re-ranked leaderboard, and comparison figure.
- `results/tables/docking_scores_decoy_set.csv` — DUD-E decoy + expanded-active docking scores (15,790 rows).
- `results/retrospective_top200/` — Stage-4a ChEMBL retrospective: whole-library documented-active flags, top-200 lookup, reports.
- `results/review_top200/` — Stage-4a top-200 chemical review (clusters, scaffolds, novelty, PAINS flags).
- `results/feature_experiment_heavy_atom/` — Stage-4c per-heavy-atom feature experiment (AUROC 0.696 -> 0.751) and heavy-atom re-ranked leaderboard.
- `results/ex8_validation/` — Stage-5 exhaustiveness-8 re-dock of the DUD-E validation set (15,920 rows) + clean-label OOF eval report.
- `results/ex8_rerank/` — Stage-5 ex4-vs-ex8 ranking-robustness analysis of the raw+ha top-200 shortlist (+ CHEMBL776).
- `results/review_top200_ha/` — chemical review of the raw+ha leaderboard top-200 (clusters, scaffolds, novelty, PAINS flags).


## Input and Metadata

- `data/raw/chembl37_ligand_library_30000.csv` — input ligand library.
- `results/tables/ligand_manifest.csv` — ligand preparation status table.
- `results/tables/receptor_manifest.csv` — receptor preparation metadata.
- `configs/beta2ar_screen.yml` — workflow parameters.
- `configs/decoy_screen.yml` — DUD-E decoy docking workflow parameters.
- `data/training/beta2ar_actives_in_library.csv` — 48 known β2-AR actives used as positive labels.
- `data/training/dude_decoys_for_dock.csv` — 3,000 DUD-E decoys for negative-control screening.
- `data/training/beta2ar_actives_expand_docked.csv` — 158 newly docked ChEMBL literature actives (Stage-3 positive expansion).
- `docs/decoy_validation_notes.md` — Stage-3 decoy validation write-up.
- `docs/stage4_retrospective_and_heavy_atom_notes.md` — Stage-4 write-up: ChEMBL retrospective limits, top-200 review, per-heavy-atom experiment.
- `docs/ex8_revalidation_note.md` — Stage-5 write-up: ex8 validation-set re-dock and shortlist ranking robustness.
- `configs/validation_redock.yml` — ex8 validation docking config (exhaustiveness 8, num_modes 20).
- `data/training/validation_redock_ligands.csv` — merged ex8 validation ligand set (206 actives + 2,978 decoys).
- `revalidation_ex8/` — rerun scripts for the ex8 validation job (run_stage5_redock.sh, make_validation_inputs.py, eval_ex8.py).

## Scripts

- `scripts/make_ligand_library.py` — build the input ligand library.
- `scripts/01_check_inputs.py` — validate project inputs.
- `scripts/02_fetch_receptors.py` — fetch and prepare receptor structures.
- `scripts/03_prepare_ligands.py` — prepare ligand structures.
- `scripts/04_run_vina.py` — run batch AutoDock Vina docking.
- `scripts/05_analyze_enopt.py` — build the ensemble matrix, weights, ranking table, and figures.
- `scripts/06_train_supervised_enopt.py` — train the supervised EnOpt model and produce the re-ranked leaderboard (Stage 2).
- `scripts/07_train_validated_enopt.py` — train the DUD-E-decoy-validated EnOpt model (Stage 3).
- `scripts/08_retrospective_top200.py` — Stage-4a ChEMBL retrospective (top-200; `--only-full` whole-library rank test).
- `scripts/09_review_top200.py` — Stage-4a top-200 chemical review (clusters/scaffolds/PAINS).
  (accepts `--ranking` / `--rank-col` so any leaderboard can be reviewed, e.g. the raw+ha one).
- `scripts/10_feature_heavy_atom_experiment.py` — Stage-4c per-heavy-atom feature experiment.
- `scripts/11_ex4_ex8_rerank.py` — Stage-5 ex4-vs-ex8 ranking-robustness analysis (raw+ha top-200 shortlist).
- `scripts/run_smoke_test.sh` — quick validation run.
- `scripts/run_pipeline.sh` — full pipeline entry point.
