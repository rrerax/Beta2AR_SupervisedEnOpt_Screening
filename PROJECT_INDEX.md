# Project Index

This index lists the files needed to review or rerun the β2-AR supervised EnOpt screening project.

Conclusion: after logP is controlled, the five-score model does not stably beat molecular weight, logP, and heavy-atom count. If a docking score must be given, use the oriented five-conformation mean. That mean follows logP and misses the classic phenylethanolamines.

The five-score ranker is not trained again. Equilibration and free energy are not started from the current Ser203–Ser207 pose or from `popc_3sn6_deposit_min.pdb`. The heavy-atom top 200 listed below is an archived table, not a purchase list. Numbers such as AUROC 0.696 and 0.751 are the record from before the logP window. They are not this conclusion.

The matched exam was run on 22 October 2026 and did not pass. Same-scaffold library negatives cover 40 actives and 138 negatives, and the property-matched DUD-E pool covers 141 actives and 561 decoys. Both top 1% counts are below the pre-registered minimum of 10, so representation docking was not started. The record is `followup_2026-10-16/exam_2026-10-22`. Without a property, PAINS, and purchasability filter before docking, plus a small prospective assay of rule-passing molecules together with property-matched negatives, this repository supports a method boundary, not an activity list.

On 26 October 2026 the same Vina protocol was repeated with basic aliphatic amines protonated at pH 7.4. CHEMBL776 ranks 14,360 and CHEMBL434 ranks 14,362 of 14,586. The oriented five-conformation mean still follows Crippen logP (Spearman 0.445). The neutral-amine input does not account for that back rank. No model was trained. The public score table is `followup_2026-10-16/gpr53_2026-10-26/amine_five_scores.csv`, and the read is `followup_2026-10-16/gpr53_2026-10-26/amine_screen_read.json`. Supplier structures stay outside this repository. The five receptor files are `data/raw/receptors` and `data/processed/receptors`.

Allowing Ser203 and Ser207 to move during Vina search on 3SN6 does not give CHEMBL434 or CHEMBL776 all three hydrogen bonds in one pose. Twelve property-matched library molecules were the negatives. The calipers were not widened. The record is `followup_2026-10-16/flex_ser203_207_2026-10-26`.

Scoring the manual three-contact orciprenaline pose, without a new search, gives −4.376 kcal/mol. The best-energy flex modes of the same two drugs score −6.702 and −7.109 on the same call. Their phenolic oxygens reach the serines, and the hydrogens do not. Deposited 2Y03 isoprenaline has all three hydrogen bonds and misses the 2.6 Å line by 0.032 Å at the benzylic oxygen. The record is `followup_2026-10-16/verify_2026-10-26`.

## Start Here

- `README.md` — project overview, methods, results, and reproduction commands.
- `notebooks/beta2ar_ensemble_screening_analysis.ipynb` — analysis notebook for the completed screen.
- `docs/final_run_summary.md` — concise final run summary and result counts.
- `docs/method_notes.md` — method rationale and interpretation notes.
- `docs/preprint_outline.md` — preprint outline bridging the β2-AR mainline with the FtsZ case study.

## Main Results

- `results/tables/docking_scores.csv` — all completed docking scores (149,325 rows).
- `results/tables/enopt_style_score_matrix.csv` — ligand-level ensemble matrix (29,865 rows).
- `results/tables/enopt_style_top_hits.csv` — archived top 200 table, not a purchase list.
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
- `results/feature_experiment_heavy_atom/` — Stage-4c per-heavy-atom archive. The 0.696 and 0.751 figures predate the logP window and are not the conclusion.
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
- `followup_2026-10-16/exam_2026-10-22/exam_matched.py` — 22 October matched, scaffold-held-out exam. The decision is `exam_decision.json` in that folder. The gate did not pass, and the script does not save a model.
