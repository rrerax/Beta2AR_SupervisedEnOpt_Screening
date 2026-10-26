# Follow-up, 16 October 2026

These folders record checks run after the Stage 3 and Stage 4 numbers in the
main README. No new docking of the 30,000-compound library was done. The
learned ranker stays unused. Docking stops at the oriented five-conformation
ranking below. No further calculation starts from the Ser203–Ser207 pair or
from `popc_3sn6_deposit_min.pdb`. The 22 October matched exam did not pass.
Its record is at the end of this file.

The week-folder scripts keep the absolute paths of the machine on which they
were run. `exam_2026-10-22/exam_matched.py` uses paths relative to this
repository. The CSV and JSON files in each folder are the archived outputs.

## Decision

On the unmatched labeled set (206 actives, 2,978 DUD-E decoys), a model that
uses only the five Vina scores reaches AUROC 0.676–0.715 and top-1% enrichment
2.4–5.3× across seeds 42–51. That enrichment is above the oriented
five-conformation mean (top-1% enrichment 1.9×). A model that uses only
molecular weight, Crippen logP, and heavy-atom count is higher (AUROC
0.773–0.819; top-1% enrichment 5.3–9.7×). Active and decoy median Crippen logP
values are 3.00 and 1.91. The oriented mean correlates with logP (Spearman
0.37 on all 3,184 molecules, 0.54 within the actives, 0.33 within the decoys).

The logP window was frozen before refitting as the overlap of the two
interquartile ranges: Crippen logP 2.3971–2.9703. It leaves 48 actives and
379 decoys, with medians 2.75 and 2.71. Inside that window the score-only
top-1% enrichment exceeds the descriptor model on 2 of 10 seeds (47 and 51).
The top 1% contains 4 molecules, so the enrichment is unstable. `pass_gate`
is false. That failed exam is the decision. The unmatched AUROC range
0.676–0.715, and the Stage-4 figures 0.696 and 0.751, are the record from
before this window. The unmatched set already had 206 actives, so the failed
gate is not a claim that a larger training set would pass.

Where a docking score is required, use the oriented five-conformation mean.
On the unmatched set the strongest single structure is 4LDE (AUROC 0.725).
The five-structure mean is not a claim that every conformation is required,
and it does not beat that single structure on global AUROC.

Stage 4's per-heavy-atom AUROC 0.751 and top-1% enrichment 7.3× describe the
unmatched comparison. They are not a docking-feature gain over molecular
properties. The heavy-atom top 200 (median logP 3.81, 3 known actives) and
the learned full-library top 200 (median logP 3.97, 5 known actives) are
shortlist chemistry. Decoy median logP is 1.91. Neither list is a purchase
list.

## CHEMBL434 and orciprenaline

CHEMBL434 (3,4-catechol) and CHEMBL776 (3,5-resorcinol, orciprenaline) share
molecular weight 211.26, Crippen logP 1.13, four hydrogen-bond donors, four
acceptors, and 15 heavy atoms. The largest absolute Vina gap between them is
0.546 kcal/mol. Oriented-mean ranks on the 29,865-ligand library are 26,697
and 25,527. The learned rank of 14 for CHEMBL434 is training-set membership.
Held out, that rank is 19,910–29,587, and it leads CHEMBL776 on 2 of 10 seeds.
CHEMBL13 repeats the pattern: in-library rank 23, held-out ranks
27,045–29,695.

That ranking is the docking result. CHEMBL776 is 25,527 and CHEMBL434 is
26,697 in the 29,865-ligand library, both toward the back. The oriented mean
follows Crippen logP (Spearman 0.37 on the 3,184 labeled molecules). Both
ligands have Crippen logP 1.13.

The geometry and membrane numbers are appended after this rank. They do not
move it. On the fixed Ser203–Ser207 side chains, the hydrogens written for
protonated orciprenaline are not all inside the rule at once. Before
minimization, H15 is 4.161 Å from Ser203 OG at 17.2°. After minimization
only the Ser207 hydrogen remains inside the rule. The other five
regioisomers do not pass the same rule. No phospholipid atom lies within
4 Å of the ligand. This low-ranked ligand is not sent to equilibration or
to a relative free energy.

The conclusion that stays: the ranking follows logP, and these two known
ligands rank toward the back. The contacts placed afterward are not the
hydrogens left in the coordinate file.

`week3_2026-10-16/week3_summary.json` fields `training_catechols` (7) and
`training_resorcinols` (2) come from an ether-matching SMARTS and are not
used. The rank fields in that file are the ones cited above. Phenolic
contacts for the pair do not differ, so side chains stay fixed.

## Regioisomers

`generation_2026-10-16` enumerates two phenolic hydroxy groups on
`CC(C)NCC(O)c1ccccc1`. The property window was written before filtering.
Ten substitution patterns collapse to six constitutions. All six stay inside
the window. Two are CHEMBL434 and CHEMBL776. Four are new. They are stored
as SMILES and properties. Synthesizability was not computed, and these poses
were not docked.

## 4LDE pose gate

`pose_gate_2026-10-17` applies one rule, written before the angles were
read, to every saved Vina mode of CHEMBL434 and CHEMBL776 on 4LDE. A contact
passes when the hydrogen–acceptor distance is at most 2.5 Å and the
donor–hydrogen–acceptor angle is at least 120°. The amine acceptor is ASP
A1113. The phenol acceptors are SER A1203, SER A1204, and SER A1207. Heavy
atoms stay at the docked coordinates. The docked amine is neutral; the
second ammonium hydrogen is placed along the lone pair. Each phenolic
hydrogen is also tested after a 180° flip about the C–O bond.

Model 1 passes the amine contact for both ligands and fails both phenolic
contacts. Across 20 modes per ligand, none pass. No free-energy sampling
was started from these poses.

## Orthosteric serine check, 19 October 2026

`orthosteric_2026-10-19` asks whether protonated orciprenaline can meet the
three orthosteric contacts on the 4LDE docking receptor, and whether serine
side-chain rotamers can bring one phenolic oxygen pair to 4.7–4.9 Å.
File numbers are the crystal numbers plus 1000: ASP A1113 is Asp113, and
SER A1203, A1204, and A1207 are Ser203, Ser204, and Ser207.

With the crystal side chains fixed, conformer phenolic O–O distances are
4.72–4.87 Å. Crystal OG–OG distances are 4.121 Å (1203–1204), 5.513 Å
(1203–1207), and 6.937 Å (1204–1207). Rigid placement onto the amine site
and a serine pair either misses the three-point target (RMSD 1.89 Å) or
penetrates TYR A1316 (closest heavy-atom contact 1.0 Å). No pose was written.

Cβ–Cβ distances are 5.415 Å, 5.471 Å, and 6.063 Å. The Ser204–Ser207 pair
was kept because 6.063 Å is not far above 6 Å. Only χ1 of these three
serines was changed, to 63°, 178°, and −65°. Backbone atoms and Cβ were not
moved. Crystal χ1 values are −66.6°, 56.0°, and 71.9°.

One rotamer pair falls in 4.7–4.9 Å: Ser203 gauche− with Ser207 trans,
OG–OG 4.735 Å. That Ser207 oxygen is 2.08 Å from the Ser203 backbone
carbonyl, so the pair is not clear of the pocket. No other rotamer pair
falls in the window. A free χ1 scan puts the Ser204–Ser207 minimum at
4.69 Å. Orciprenaline was not placed on a rotated receptor. The other four
crystals were not measured, because one rotamer distance does fall inside
4.7–4.9 Å.

The hydrogen-bond rule is unchanged: hydrogen–acceptor distance at most
2.5 Å and donor–hydrogen–acceptor angle at least 120°. The 20 October
section continues this scan.

## 3SN6 serine pair, 20 October 2026

`orthosteric_2026-10-20` keeps the backbone and Cβ fixed. A serine pair is
accepted when the hydroxyl-oxygen distance is 4.7–4.9 Å and each oxygen is
at least 2.6 Å from backbone heavy atoms of every other residue. A ligand
pose passes only when the amine and both phenolic hydrogens meet the
2.5 Å / 120° rule and every other ligand–receptor heavy-atom contact is at
least 2.6 Å. Serines are taken from the sequence AIASSIVSF, and Asp113 from
EFWTSIDVLCVTA.

On 4LDE the continuous χ1 scan does enter 4.7–4.9 Å. The 4.69 Å value in
the 19 October note is the Ser204–Ser207 minimum, and the distance also
runs above it. The χ1 7°–16° arrangement meets the three hydrogen bonds,
but another heavy atom approaches to 2.24 Å, so that pose fails the 2.6 Å
rule.

On 3SN6 the representative angle was fixed before the ligand was moved:
the clear angle closest to 4.80 Å. Chain R uses residue numbers 113, 203,
204, and 207. Ser203–Ser204 at χ1 93° and 2° fails because Phe290 CE2 is
2.273 Å away. Ser203–Ser207 at χ1 −167° and 49° passes. The amine hydrogen
is 2.298 Å from Asp113 at 156.6°. The phenolic hydrogens are 2.437 Å at
131.6° and 2.482 Å at 127.8°. The nearest other heavy atom is Val114 CG2
at 2.734 Å. This side chain was built by the rule. It is not a crystal
structure of orciprenaline. 2RH1, 4LDL, and 5D5A were not used after this
pass.

The other five regioisomers were measured as protonated conformers.
Phenolic O–O distances are 2.632–2.719 Å (2,3), 4.730–4.800 Å (2,4),
5.485–5.505 Å (2,5), 4.764–4.793 Å (2,6), and 2.707–2.792 Å (3,4,
CHEMBL434). Only 2,4 and 2,6 fall in 4.7–4.9 Å, so only those two were
placed on the same 3SN6 χ1 pair. The 2,4 isomer meets the hydrogen bonds
and comes within 1.639 Å of Asn312 ND2. The 2,6 isomer does not meet all
three hydrogen bonds at once: with both phenols passing, the amine
hydrogen is 3.099 Å from Asp113 at 76.8°. No second isomer passes, and no
relative free-energy calculation was started.

## Written hydrogens, 22 October 2026

These measurements are the numbers appended after the ranking. The 20 October
paragraph records the search. The saved file
`3sn6_orciprenaline.pdb` contains one hydrogen on each oxygen. The Ser203
figure in that paragraph, and the 2.441 Å / 131.9° line in
`3sn6_pose_export.txt`, count a 180° phenolic flip that was not written
into the file. The hydrogen that is in the file, H15 on the oxygen aimed
at Ser203, is 4.161 Å from Ser203 OG at 17.2°. The amine hydrogen aimed at
the nearer Asp113 oxygen is 2.432 Å at 149.4°. The other phenolic hydrogen
is 2.479 Å from Ser207 OG at 128.0°. The three contacts are not present
together. No side chain was moved after this measurement, and no failed
isomer was given another serine pair.

The same hydrogens were copied into a POPC membrane. Minimization used
only heavy atoms present in deposited 3SN6. Phe193 stops at CB. The two
chain breaks stay open, with C–N distances 10.776 Å and 8.564 Å. After
minimization the three hydrogen distances and angles are 2.838 Å at 76.0°,
5.036 Å at 18.4°, and 1.861 Å at 152.2°. Only the Ser207 hydrogen remains
inside the rule. No other ligand–receptor heavy contact is shorter than
2.6 Å. No phospholipid atom lies within 4 Å of the ligand. The nearest
phospholipid heavy atom is 10.613 Å from a ligand heavy atom.

Equilibration was not started: the three hydrogens are not all inside the
rule after minimization. Relative free energy was not started: one ligand
was placed, and the three contacts do not remain together.

## What this follow-up does not contain

No free-energy number or activity assay was produced here. The membrane
file is a minimized geometry check. Feature importances are not a result.
The six regioisomers are not a purchase list.

## Matched exam, 22 October 2026

`exam_2026-10-22/exam_matched.py` reads this repository through relative
paths. The calipers and the pass rule were fixed before the metrics. No
model file was saved, and no docking was started.

The calipers are absolute: Crippen logP 0.5, molecular weight 50 Da,
hydrogen-bond donors 1, hydrogen-bond acceptors 1, and equal formal charge
in the input SMILES. Each active could take at most five negatives. The
gate pool is docked library molecules that share a Bemis–Murcko scaffold
with one of the 206 actives, lie inside those calipers, and are not
themselves in the 206. The split is five-fold scaffold holdout, seed 42.
The model is the earlier XGBoost setting: 15 trees, learning rate 0.3,
depth 6. A pool passes only when it has at least 30 covered actives, 100
negatives, 10 active scaffolds, and 10 molecules in the top 1%, the
five-score out-of-fold AUROC exceeds each baseline by more than its own
AUROC standard error, and its top-1% enrichment is higher. The baselines
are molecular weight, Crippen logP, and heavy-atom count, and the oriented
4LDE score. The DUD-E rematch is archived beside the gate and cannot pass
it.

The same-chemotype pool has 40 actives, 138 negatives, 22 scaffolds, and
178 molecules. Its top 1% is one molecule. The DUD-E rematch has 141
actives, 561 decoys, 106 scaffolds, and 702 molecules. Its top 1% is seven
molecules. `exam_pass_gate` is false. The reason is
`top1_percent_below_10_molecules`. The learned ranker stays unused.
Representation docking was not started.

AUROC was stored after that size failure. It is not the gate. On the
same-chemotype pool the five-score AUROC is 0.488 (0.386–0.589), the
descriptor model is 0.500 (0.398–0.602), 4LDE is 0.471 (0.370–0.571), and
the oriented five-conformation mean is 0.504 (0.402–0.606). On the DUD-E
rematch the five-score AUROC is 0.549 (0.495–0.603), the descriptor model
is 0.568 (0.514–0.623), 4LDE is 0.628 (0.574–0.682), and the five-conformation
mean is 0.615 (0.561–0.669). `numerical_comparison_passes` is false in both
pools. The oriented mean still follows Crippen logP: Spearman 0.388 on the
same-chemotype pool and 0.432 on the DUD-E rematch.

The unmatched decoy table remains shifted. The 206 actives have median
Crippen logP 3.00; the 2,978 decoys have median 1.91. Formal charge in the
stored SMILES is 0 for 204 of the 206 actives. Among the decoys, 1,848 are
written as +1 and 1,048 as 0. The charge caliper kept the 561 neutral
matched decoys. Their medians are Crippen logP 3.16 and molecular weight
446, against 3.07 and 427 for the 141 covered actives. The same-chemotype
medians are 2.50 against 2.46, and 318 against 307.

CHEMBL434 is inside the 206, with Crippen logP 1.13 and scaffold
`c1ccccc1`. Both matches cover it. This recount places its oriented-mean
rank at 26,698 of 29,865. CHEMBL776 is not inside the 206. This recount
places it at 25,529. Both stay toward the back of the library.

## GPR-53 preparation, 22 October 2026

The plated Enamine GPCR library is outside this repository, at
`D:\repositories\Enamine_GPCR_Library_plated`. `prepare_gpr53.py` read the
4 October 2026 CSV of 53,440 catalog IDs. The pH 7.4 rule was fixed before
that file was parsed: aliphatic amines, amidines, and guanidines are
protonated; carboxylic and sulfonic acids are deprotonated; anilines,
amides, sulfonamides, pyridines, imidazoles, and phenols stay neutral.
PAINS and Brenk are columns. They do not delete rows.

The prepared table has 53,442 rows. CHEMBL434 and CHEMBL776 were not in
the catalog and were appended. CHEMBL434 is `CC(C)[NH2+]CC(O)c1ccc(O)c(O)c1`
and keeps the catechol flags. CHEMBL776 is `CC(C)[NH2+]CC(O)c1cc(O)cc(O)c1`.
14,586 rows have at least one basic aliphatic amine. 541 rows carry a PAINS
flag and 2,847 carry a Brenk flag. No model was trained. The docking list
is those 14,586 rows, with CHEMBL434 and CHEMBL776 first, in
`prepared\gpr53_amine_ligands.csv`. The learned ranker stays unused.

The five receptor PDBQT files and AutoDock Vina 1.2.5 were added on the
same day. The original box centers and sizes were unchanged. A two-ligand
smoke test wrote ten ready scores in
`Enamine_GPCR_Library_plated\docking\smoke_scores.csv`. The full amine
screen writes its manifest, ligand files, and scores under that docking
folder, not into the original ChEMBL tables. Vendor files stay out of this
Git repository.

## Protonated amine screen, 26 October 2026

The ChEMBL-library ranks CHEMBL776 25,527 and CHEMBL434 26,697 of 29,865
were computed with neutral amines. The β2 orthosteric site recognizes a
protonated amine, so those ranks did not separate logP collinearity from
that input. This screen holds the receptor set, the original box centers
and sizes, and exhaustiveness 4 fixed. The pH 7.4 protonation rule was
fixed before the Enamine file was read. 14,586 molecules with a basic
aliphatic amine were docked, with protonated CHEMBL434 and CHEMBL776
appended. No model was trained. All 72,930 scores are ready. The score
table stays in the local docking folder. The read is
`gpr53_2026-10-26/amine_screen_read.json`.

On the oriented five-conformation mean, CHEMBL776 is 14,360 and CHEMBL434
is 14,362 of 14,586. Both means are −7.28 kcal/mol, and the absolute gap
is 0.003 kcal/mol. Both Crippen logP values, taken on the neutral parent,
are 1.129. The oriented mean still follows logP: Spearman 0.445 on these
14,586 molecules. The library median logP is 2.664. The top 1% is 146
molecules and has median logP 3.812. Oriented 4LDE versus logP is Spearman
0.412. Protonating the amine does not move the two phenylethanolamines
forward, and the score still follows logP. The neutral-amine input does
not account for the back rank. The learned ranker stays unused. This table
is not an activity list.
