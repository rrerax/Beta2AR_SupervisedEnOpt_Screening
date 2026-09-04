# Whole-library retrospective validation (ChEMBL ADRB2)

Pulled all 5686 ADRB2 (CHEMBL210) potency rows (types IC50,Ki,EC50,KD,Kd), flagged library molecules with pChEMBL >= 6 (<= 1 uM), excluded the training actives, and tested whether the remaining documented actives rank high in the Stage-3 leaderboard.

- Documented strong actives among library unknowns: **2** / 29817
- Rank-AUC (documented vs rest): **0.2466** (0.5 = random)
- Mann-Whitney one-sided p: **0.878**
- Median rank of documented actives: **22512** / 29817 (**75.5%**; random = 50%)

Top-N enrichment:

| top_n | documented_in_top | chance_expect | enrichment |
|---|---|---|---|
| 50.0 | 0.0 | 0.0 | 0.0 |
| 200.0 | 0.0 | 0.01 | 0.0 |
| 500.0 | 0.0 | 0.03 | 0.0 |
| 1000.0 | 0.0 | 0.07 | 0.0 |
| 2000.0 | 0.0 | 0.13 | 0.0 |
| 5000.0 | 0.0 | 0.34 | 0.0 |

Files: `full_library_documented_flags.csv`, `chembl_adrb2_all_rows.csv`, `full_library_retro_summary.csv`.
