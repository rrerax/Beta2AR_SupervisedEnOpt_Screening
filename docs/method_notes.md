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

Each prepared ligand was docked against each receptor conformation. Docking scores were merged into a single ligand-by-conformation matrix.

Conclusion: after logP is controlled, the five-score model does not stably beat molecular weight, logP, and heavy-atom count. If a docking score must be given, use the oriented five-conformation mean. That mean follows logP and misses the classic phenylethanolamines.

Three things stay stopped. The five-score ranker is not trained again. Equilibration and free energy are not started from the current Ser203–Ser207 pose or from `popc_3sn6_deposit_min.pdb`. The heavy-atom model's top 200 is not a purchase list.

The ranking table still reports:

- `enopt_style_weighted_score_kcal_mol`: ensemble weighted docking score used for ranking.
- `ensemble_best_score_kcal_mol`: best single-conformation docking score.
- `ensemble_mean_score_kcal_mol`: average docking score across conformations.
- `ensemble_score_std_kcal_mol`: score variability across conformations.

## Matched exam, 22 October 2026

The exam in `followup_2026-10-16/exam_2026-10-22` did not pass. Same-scaffold library negatives inside the frozen property calipers cover 40 actives and 138 negatives, so the top 1% is one molecule. Property-matched DUD-E decoys cover 141 actives and 561 decoys, so the top 1% is seven molecules. The pre-registered minimum is ten. The five-score ranker stays unused, and representation docking was not started.

The stored AUROC values are not that gate. On both pools the five-score model is behind molecular weight, Crippen logP, and heavy-atom count, and behind 4LDE. The write-up is in `followup_2026-10-16/README.md`.

## Protonated amine screen, 26 October 2026

The earlier ChEMBL ranks used neutral amines. GPR-53 keeps the five receptors, the original boxes, and exhaustiveness 4, and docks 14,586 basic aliphatic amines after the pH 7.4 rule. CHEMBL776 is 14,360 and CHEMBL434 is 14,362 of 14,586 on the oriented five-conformation mean. Both means are −7.28 kcal/mol. Spearman with Crippen logP is 0.445. Protonation does not move the two phenylethanolamines forward, and the score still follows logP. No model was trained.
