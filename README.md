# β2-AR ensemble docking screen

Reproducible AutoDock Vina screen of the β2-adrenergic receptor across five crystal conformations, with a pre-specified test of whether a learned combination of those five scores beats a simple average once molecular properties are controlled.

**Current conclusion.** After Crippen logP is controlled, the five-score XGBoost ranker does not stably beat a descriptor model that sees only molecular weight, logP, and heavy-atom count. Reported docking scores use the oriented five-conformation mean (more negative Vina score ranks first). That mean is a reporting choice: ranking each conformation alone shows that one structure (4LDE) already matches or beats the mean on the unmatched labeled set. This repository is a computational priority ranking and a method boundary. It is not an activity list and not a purchase list.

Full dated checks live under [`followup_2026-10-16`](followup_2026-10-16/README.md). File map: [`PROJECT_INDEX.md`](PROJECT_INDEX.md).

## What was docked

| Item | Count |
|---|---:|
| Input ChEMBL ligands | 30,000 |
| Prepared ligands | 29,865 |
| Receptor conformations | 5 |
| Completed docking rows | 149,325 |

Receptors (inactive and active):

| conformation_id | PDB |
|---|---|
| `inactive_carazolol_2RH1` | 2RH1 |
| `inactive_carazolol_5D5A` | 5D5A |
| `active_gs_complex_3SN6` | 3SN6 |
| `active_bi167107_4LDE` | 4LDE |
| `active_hbi_4LDL` | 4LDL |

Engine: AutoDock Vina 1.2, exhaustiveness 4 for the main library. Labels for the hard test: 206 literature actives and 2,978 DUD-E decoys. Score matrix: [`results/tables/docking_scores.csv`](results/tables/docking_scores.csv).

The workflow follows the EnOpt idea of Bhatt, Wang & Durrant (*Sci Rep* 14:20722, 2024): how to combine a ligand’s per-conformation scores is system-specific. This repository tests that claim on β2-AR and records where it fails the pre-specified gate.

## Decision record

Hyperparameters and the pass rule were fixed before the logP-controlled fit: XGBoost with 15 trees, learning rate 0.3, max depth 6; ten seeds (42–51); three-fold stratified out-of-fold prediction. The descriptor control used only molecular weight, Crippen logP, and heavy-atom count.

1. **Unmatched set (206 actives / 2,978 decoys).** Score-only AUROC 0.676–0.715 and top-1% enrichment 2.4–5.3× across seeds. Oriented five-conformation mean: AUROC 0.719, top-1% enrichment 1.9×. Descriptor-only model: AUROC 0.773–0.819, top-1% enrichment 5.3–9.7×. Active and decoy median logP are 3.00 and 1.91; the oriented mean correlates with logP (Spearman 0.37).

2. **Frozen logP window.** Intersection of the two Crippen logP IQRs (2.40–2.97) leaves 48 actives and 379 decoys. Inside that window the score-only top-1% enrichment beats the descriptor model on 2 of 10 seeds. The top 1% has four molecules. The gate does not pass. The learned ranker stays unused.

3. **Matched exam (22 October 2026).** Same-scaffold library negatives: 40 actives / 138 negatives (top 1% = 1 molecule). Property-matched DUD-E rematch: 141 / 561 (top 1% = 7). Pre-registered minimum is 10 molecules in the top 1%. Representation docking was not started. Record: [`exam_2026-10-22`](followup_2026-10-16/exam_2026-10-22/exam_decision.json).

4. **One conformation at a time (29 October 2026).** No model was fit. On 206 / 2,978, 4LDE alone reaches AUROC 0.725 and top-1% enrichment 2.43×, against 0.719 and 1.94× for the oriented mean. Inside the logP window, 3SN6 and 4LDL each meet or beat the mean. Table: [`conformation_2026-10-29`](followup_2026-10-16/conformation_2026-10-29/conformation_contribution.csv).

5. **Protonated amines (pH 7.4).** Same five receptors and boxes; 14,586 basic aliphatic amines. CHEMBL776 ranks 14,360 and CHEMBL434 ranks 14,362. Unrounded means −7.2834 and −7.2802 kcal/mol (both −7.28 to two decimals). Spearman with logP is 0.445. Protonation does not move these two phenylethanolamines forward. Table: [`gpr53_2026-10-26/amine_five_scores.csv`](followup_2026-10-16/gpr53_2026-10-26/amine_five_scores.csv). On the earlier neutral-amine ChEMBL library the oriented-mean ranks were 25,527 and 26,697 of 29,865.

6. **Flexible Ser203 / Ser207 on 3SN6.** Nine Vina modes; energy is not the pass rule. Neither CHEMBL434 nor CHEMBL776 places both phenolic hydrogens on the two serines in one pose. Record: [`flex_ser203_207_2026-10-26`](followup_2026-10-16/flex_ser203_207_2026-10-26/flex_decision.json).

7. **Manual pose and crystal geometry.** The accepted three-contact orciprenaline pose scores −4.376 kcal/mol by Vina `score_only` on the screen receptor. Best-energy flex modes of the same two drugs score −6.702 and −7.109. Deposited 2Y03 isoprenaline has the three hydrogen bonds and misses the frozen 2.6 Å heavy-contact line by 0.032 Å at the benzylic oxygen. Record: [`verify_2026-10-26`](followup_2026-10-16/verify_2026-10-26/).

8. **2Y03 membrane MD.** Three 50 ns replicates from deposited isoprenaline (turkey β1-AR). A rule written before production called a contact persistent at occupancy ≥ 0.20 in at least two replicates. The amine–aspartate hydrogen bond and both phenolic hydrogen bonds meet that rule (one phenolic contact is essentially absent in the third replicate). Bridging waters do not. This is not an affinity and does not compare CHEMBL776 with CHEMBL434. Trajectories stay outside git. Record: [`y03_md_2026-10-26`](followup_2026-10-16/y03_md_2026-10-26/).

## Earlier unmatched record (Stage 3–5)

These numbers are the archive from before the logP window. They are not the decision above.

| Comparison (206 actives / 2,978 decoys, OOF) | AUROC | Top-1% EF |
|---|---:|---:|
| Five-score model (raw) | 0.696 | 3.88× |
| Oriented five-conformation mean | 0.719 | 1.94× |
| Raw + per-heavy-atom features | 0.751 | 7.28× |

On that unmatched set the five-score model and the oriented mean have overlapping AUROC confidence intervals. Top-1% enrichment looked higher for the model before property control. Appending per-heavy-atom scores moves the unmatched AUROC to 0.751; that move is not treated as a docking-feature gain over molecular weight, logP, and heavy-atom count. Re-docking the validation set at exhaustiveness 8 leaves those unmatched intervals overlapping. Details: `docs/decoy_validation_notes.md`, `docs/stage4_retrospective_and_heavy_atom_notes.md`, `docs/ex8_revalidation_note.md`.

## Repository layout

```text
configs/                 Screen and decoy docking parameters
data/training/           Actives, expanded actives, DUD-E decoys
scripts/                 Preparation, docking, training, Stage 4–5 analyses
results/tables/          Docking scores and ensemble matrices
results/supervised_enopt*/  Archived Stage-2/3 models and leaderboards
results/feature_experiment_heavy_atom/  Archived Stage-4 raw+ha table
followup_2026-10-16/     Pre-specified checks after Stage 3–4 (decision record)
docs/                    Method notes and Stage write-ups
```

## Reproduce

```bash
bash scripts/00_setup_env.sh
bash scripts/run_smoke_test.sh
bash scripts/run_pipeline.sh
```

Archived Stage-3/4/5 scripts (they do not change the October decision):

```bash
python scripts/07_train_validated_enopt.py
python scripts/10_feature_heavy_atom_experiment.py
python scripts/11_ex4_ex8_rerank.py
```

Per-conformation ranking on the stored matrix (no training):

```bash
python followup_2026-10-16/conformation_2026-10-29/conformation_contribution.py
```

## Limits

- Docking scores are computational priorities, not measured affinities.
- The learned ranker stays unused. Where a docking score is reported, use the oriented five-conformation mean.
- The five-conformation mean is not evidence that every conformation is required.
- Archived top-200 and raw+ha tables are not purchase lists.
- Supplier catalog structures and MD trajectories are outside this repository.
- Without a property / PAINS / purchasability filter before docking, and a small prospective assay of rule-passing molecules with property-matched negatives, this repository supports a method boundary, not an activity claim.
