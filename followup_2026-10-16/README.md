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
top-1% enrichment exceeds the descriptor model on 2 of 10 seeds. The top 1%
contains 4 molecules, so the enrichment is unstable. `pass_gate` is false.

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

## What this follow-up does not contain

No free-energy result, membrane system, protonation file, or activity assay
was produced here. Feature importances are not a result.
