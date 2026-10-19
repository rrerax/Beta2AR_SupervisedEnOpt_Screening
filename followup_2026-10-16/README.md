# Follow-up, 16 October 2026

These folders record checks run after the Stage 3 and Stage 4 numbers in the
main README. No new docking of the 30,000-compound library was done. The
learned ranker stays unused.

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
2.5 Å and donor–hydrogen–acceptor angle at least 120°.

## What this follow-up does not contain

No free-energy number, membrane system, protonation coordinate file, or
activity assay was produced here. Feature importances are not a result.
The six regioisomers were not docked and are not a purchase list.
