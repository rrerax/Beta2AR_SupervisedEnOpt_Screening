# Follow-up, 16 October 2026

These folders record checks run after the Stage 3 and Stage 4 numbers in the
main README. No new docking of the 30,000-compound library was done. The
learned ranker stays unused. Docking stops at the oriented five-conformation
ranking below. No further calculation starts from the Ser203–Ser207 pair or
from `popc_3sn6_deposit_min.pdb`.

The scripts keep the absolute paths of the machine on which they were run.
The CSV and JSON files in each folder are the archived outputs.

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
