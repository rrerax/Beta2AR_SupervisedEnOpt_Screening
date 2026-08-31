# Method Notes

## Rationale

β2-adrenergic receptor is a GPCR target with well-characterized inactive and active-state structures. A single rigid receptor conformation can miss ligands that fit better into an alternative binding-pocket geometry. The project therefore uses an ensemble docking strategy to reduce dependence on one static receptor model.

## Receptor Selection

Five β2-AR structures were used to represent inactive and active conformational states:

| conformation_id         | pdb_id   | state_label   |
|:------------------------|:---------|:--------------|
| inactive_carazolol_2RH1 | 2RH1     | inactive      |
| inactive_carazolol_5D5A | 5D5A     | inactive      |
| active_gs_complex_3SN6  | 3SN6     | active        |
| active_bi167107_4LDE    | 4LDE     | active        |
| active_hbi_4LDL         | 4LDL     | active        |

## Ligand Preparation

The screen starts from a 30,000-compound ChEMBL37 ligand table. Structures that failed desalting, 3D embedding, or PDBQT conversion were kept in the ligand manifest instead of silently dropping them. This makes the screen auditable and keeps the prepared-compound count separate from the initial library size.

## Docking and Ranking

Each prepared ligand was docked against each receptor conformation. Docking scores were merged into a single ligand-by-conformation matrix. A stage-1 EnOpt-style weighted consensus score provides a baseline ranking. The final ranking comes from a supervised EnOpt model (XGBoost) trained on the receptor's own experimental data: 48 known β2-AR actives from ChEMBL as positives and the remaining library molecules as assumed negatives. The model learns how to combine the per-conformation scores for β2-AR specifically; out-of-fold predictions reach an AUROC of 0.66, and known actives' median rank improves from the top 39% to the top 22.6% of the library.

The final ranking table reports:

- `enopt_style_weighted_score_kcal_mol`: ensemble weighted docking score used for ranking.
- `ensemble_best_score_kcal_mol`: best single-conformation docking score.
- `ensemble_mean_score_kcal_mol`: average docking score across conformations.
- `ensemble_score_std_kcal_mol`: score variability across conformations.

## Interpretation

Both the stage-1 and supervised results should be interpreted as a computational shortlist. The supervised EnOpt model improves ranking over the stage-1 baselines (OOF AUROC 0.66; actives median rank 39% → 22.6%), and a 3,000-compound DUD-E decoy screen (`configs/decoy_screen.yml`) is prepared as an external negative control. Strong docking scores can suggest promising candidates, but they do not prove binding or biological activity. Follow-up work could include binding-site visual inspection, decoy benchmarking, molecular dynamics, MM/GBSA rescoring, or experimental validation.
