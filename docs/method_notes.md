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

## If calculation continues

Change the exam first, then the representation. This exam has not been run. If it does not pass, stop.

Rebuild the negatives from the same chemotypes as the 206 actives. Match logP, molecular weight, hydrogen-bond donor and acceptor counts, and charge. Validate by Bemis–Murcko scaffold groups. Report both baselines in the same exam: molecular weight, logP, and heavy-atom count alone, and 4LDE alone.

Prepare each ligand as the protonated amine. Before docking, require the Asp113 salt bridge and the Ser203/Ser207 hydrogen bonds. Relax only a serine side chain that blocks a contact. The success line is fixed: CHEMBL776 and CHEMBL434 each keep all three hydrogen bonds at once, every other heavy-atom contact is at least 2.6 Å, and the score beats the logP baseline on that matched, scaffold-grouped exam. Raising Vina exhaustiveness does not change the present conclusion.

Filter the screening library for properties, PAINS, and purchasability before docking. Otherwise the ranking still rewards high logP. Without that filter and a small prospective assay of the molecules that pass the rule, tested together with property-matched negatives, this repository supports a method boundary, not an activity list.
