# ex4 vs ex8 re-dock robustness (heavy-atom top-200 shortlist + CHEMBL776)

Re-docked 201 ligands (top-200 of the raw+ha leaderboard + CHEMBL776, the known failure case) at exhaustiveness 8, num_modes 20, seed-fixed; compared against their ex4 (exhaustiveness 4) scores from the library run.

| metric | value |
|---|---:|
| ligands scored in both runs | 201 |
| Spearman r, best score ex4 vs ex8 | 0.874 |
| Pearson r, best score ex4 vs ex8 | 0.874 |
| Spearman r, best score / heavy-atom | 0.891 |
| Spearman r, mean score ex4 vs ex8 | 0.905 |
| median delta best score (ex8-ex4, kcal/mol) | +0.00 |
| ligands moving > 1 kcal/mol (abs) | 3 / 201 |

## Top-N stability (subset ranks among these 201 ligands)

### best score
| ex4 top-N | in ex8 top-50 | in ex8 top-100 | in ex8 top-150 | in ex8 top-200 |
|---|--:|--:|--:|--:|
| **50** | 41 | 48 | 50 | 50 |
| **100** | 50 | 89 | 99 | 100 |
| **150** | 51 | 96 | 139 | 150 |
| **200** | 51 | 100 | 150 | 200 |

### best/HAC
| ex4 top-N | in ex8 top-50 | in ex8 top-100 | in ex8 top-150 | in ex8 top-200 |
|---|--:|--:|--:|--:|
| **50** | 41 | 50 | 50 | 50 |
| **100** | 46 | 85 | 99 | 100 |
| **150** | 49 | 99 | 141 | 150 |
| **200** | 50 | 100 | 150 | 199 |

## Per-pocket score stability (Spearman, ex4 vs ex8)
| conformation | Spearman r | n |
|---|---:|---:|
| active_bi167107_4LDE | 0.837 | 201 |
| active_gs_complex_3SN6 | 0.856 | 201 |
| active_hbi_4LDL | 0.772 | 201 |
| inactive_carazolol_2RH1 | 0.846 | 201 |
| inactive_carazolol_5D5A | 0.834 | 201 |

## Known failure case CHEMBL776 (orciprenaline)

- ex4 best -8.37 (rank 201.0/201) -> ex8 best -8.07 (rank 201.0/201); best/HAC ex4 -0.558 -> ex8 -0.538
- ex4 leaderboard absolute rank (raw+ha model): 19209 (old decoy-model rank 27659)

## Biggest rank movers (best score)

| ligand_id | ex4 best | ex8 best | rank ex4 | rank ex8 | d_rank | enopt_ha_rank | known |
|---|---:|---:|---:|---:|---:|---:|---:|
| CHEMBL25334 | -10.49 | -11.95 | 129 | 7 | -122 | 141 | 0 |
| CHEMBL20721 | -10.20 | -11.33 | 171 | 55 | -116 | 158 | 0 |
| CHEMBL36360 | -11.70 | -10.74 | 17 | 105 | +88 | 155 | 0 |
| CHEMBL3309727 | -10.18 | -10.93 | 178 | 93 | -85 | 161 | 0 |
| CHEMBL31134 | -10.59 | -10.02 | 113 | 194 | +81 | 143 | 0 |
| CHEMBL25623 | -11.50 | -10.64 | 44 | 116 | +72 | 46 | 0 |
| CHEMBL19028 | -10.30 | -11.08 | 155 | 84 | -71 | 131 | 0 |
| CHEMBL4436519 | -10.21 | -10.83 | 168 | 98 | -70 | 49 | 0 |
| CHEMBL34287 | -11.31 | -10.51 | 62 | 131 | +69 | 93 | 0 |
| CHEMBL3310238 | -10.14 | -10.62 | 186 | 118 | -68 | 153 | 0 |
| CHEMBL17703 | -11.06 | -10.33 | 84 | 151 | +67 | 16 | 0 |
| CHEMBL10815 | -10.19 | -10.72 | 174 | 107 | -67 | 198 | 0 |

