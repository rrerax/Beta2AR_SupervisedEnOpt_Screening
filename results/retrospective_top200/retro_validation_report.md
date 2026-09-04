# Retrospective top-200 validation (ChEMBL lookup)

Lookup target: **CHEMBL210** (human beta-2 adrenergic receptor), potency types
IC50,Ki,EC50,KD,Kd, "documented strong" = at least one record with pChEMBL >= 6
(<= 1 uM). Library molecules never received a literature label at training
time, so these top-ranked molecules were **unknown to the model**.

| set | n | documented strong | hit rate |
|---|--:|--:|--:|
| no-label top | 195 | 0 | 0.0% |
| no-label control (random) | 500 | 0 | 0.0% |

- Enrichment ratio (top vs control): **nanx**, Fisher exact one-sided
  p = **1**.
- Sanity check only (these were training positives): known in-library actives
  in top 200 = 5/48 vs chance 0.32
  (15.6x).

Files: `top200_no_label_with_chembl_hits.csv` (ranked list with best pChEMBL),
`control_sample_with_chembl_hits.csv`, `chembl_lookup_raw_rows.csv`,
`retro_summary.csv`, `fig_top200_vs_base_rate.png`.

Caveats: retrospective literature recovery is a lower bound (not everything
measured is deposited in ChEMBL; BindingDB / patents not yet merged). It is
still a direct check of ranking power, at zero experimental cost.
