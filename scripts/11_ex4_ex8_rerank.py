#!/usr/bin/env python3
"""ex4 vs ex8 leaderboard robustness on the re-docked shortlist (201 ligands)."""
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

import argparse
REPO = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(description="ex4 vs ex8 leaderboard robustness "
                                         "(heavy-atom top-200 shortlist + CHEMBL776)")
ap.add_argument("--ex8", default=str(REPO / "results/ex8_rerank/scores_ex8_top200_plus776.csv"),
                help="ex8 re-dock score file (long format)")
ap.add_argument("--ex4", default=str(REPO / "results/tables/docking_scores.csv"))
ap.add_argument("--ranking", default=str(REPO / "results/feature_experiment_heavy_atom/library_ranking_heavy_atom_model.csv"))
ap.add_argument("--library", default=str(REPO / "data/raw/chembl37_ligand_library_30000.csv"))
ap.add_argument("--out", default=str(REPO / "results/ex8_rerank"))
args = ap.parse_args()
OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)

# --- load ex4 library scores (long) ---
ex4 = pd.read_csv(args.ex4)
ex4 = ex4[ex4["status"] == "ready"].copy()
ex4["best_score_kcal_mol"] = pd.to_numeric(ex4["best_score_kcal_mol"], errors="coerce")
ex4 = ex4.dropna(subset=["best_score_kcal_mol"])

# --- load ex8 re-dock scores (long) ---
ex8 = pd.read_csv(args.ex8)
ex8["best_score_kcal_mol"] = pd.to_numeric(ex8["best_score_kcal_mol"], errors="coerce")
ex8 = ex8.dropna(subset=["best_score_kcal_mol"])

# --- leaderboard (heavy-atom model, ex4-derived) ---
rank = pd.read_csv(args.ranking)

# --- properties (heavy atoms) ---
lib = pd.read_csv(args.library)

def best_wide(df):
    g = df.groupby(["ligand_id", "conformation_id"])["best_score_kcal_mol"].min()
    w = g.unstack("conformation_id")
    w["best"] = w.min(axis=1)
    w["mean"] = w.mean(axis=1)
    return w

w4 = best_wide(ex4)
w8 = best_wide(ex8)
common = sorted(set(w4.index) & set(w8.index))
print("ex4 ligands:", len(w4), "| ex8 ligands:", len(w8), "| common:", len(common))

d = pd.DataFrame(index=common)
d["ex4_best"] = w4.loc[common, "best"]
d["ex8_best"] = w8.loc[common, "best"]
d["ex4_mean"] = w4.loc[common, "mean"]
d["ex8_mean"] = w8.loc[common, "mean"]
d["n_conf_ex8"] = w8.loc[common].notna().sum(axis=1)

hac = lib.set_index("ligand_id")["heavy_atoms"].reindex(common).astype(float)
d["heavy_atoms"] = hac
d["ex4_ha"] = d["ex4_best"] / hac
d["ex8_ha"] = d["ex8_best"] / hac

# subset ranks (1 = best/most negative)
d["rank_ex4_best"] = d["ex4_best"].rank(ascending=True).astype(int)
d["rank_ex8_best"] = d["ex8_best"].rank(ascending=True).astype(int)
d["rank_ex4_ha"] = d["ex4_ha"].rank(ascending=True).astype(int)
d["rank_ex8_ha"] = d["ex8_ha"].rank(ascending=True).astype(int)

# leaderboard absolute rank from heavy-atom model (ex4 features)
lb = rank.set_index("ligand_id")[["enopt_ha_rank", "old_decoy_rank", "is_known_active"]]
d = d.join(lb)
d["delta_best"] = d["ex8_best"] - d["ex4_best"]  # + means worse at ex8
d["rank_delta_best"] = d["rank_ex8_best"] - d["rank_ex4_best"]
d = d.sort_values("rank_ex4_best")

n = len(d)
def sp(a, b):
    v = d[a].corr(d[b], method="spearman")
    return v
def pear(a, b):
    return d[a].corr(d[b], method="pearson")

report = []
report.append("# ex4 vs ex8 re-dock robustness (heavy-atom top-200 shortlist + CHEMBL776)")
report.append("")
report.append(f"Re-docked {n} ligands (top-200 of the raw+ha leaderboard + CHEMBL776, "
              "the known failure case) at exhaustiveness 8, num_modes 20, seed-fixed; "
              "compared against their ex4 (exhaustiveness 4) scores from the library run.")
report.append("")
report.append("| metric | value |")
report.append("|---|---:|")
report.append(f"| ligands scored in both runs | {n} |")
report.append(f"| Spearman r, best score ex4 vs ex8 | {sp('ex4_best','ex8_best'):.3f} |")
report.append(f"| Pearson r, best score ex4 vs ex8 | {pear('ex4_best','ex8_best'):.3f} |")
report.append(f"| Spearman r, best score / heavy-atom | {sp('ex4_ha','ex8_ha'):.3f} |")
report.append(f"| Spearman r, mean score ex4 vs ex8 | {sp('ex4_mean','ex8_mean'):.3f} |")
md = d["delta_best"].median()
report.append(f"| median delta best score (ex8-ex4, kcal/mol) | {md:+.2f} |")
big = int((d["delta_best"].abs() > 1.0).sum())
report.append(f"| ligands moving > 1 kcal/mol (abs) | {big} / {n} |")
report.append("")
report.append("## Top-N stability (subset ranks among these %d ligands)" % n)
for name, r4, r8 in [("best score", "rank_ex4_best", "rank_ex8_best"),
                     ("best/HAC", "rank_ex4_ha", "rank_ex8_ha")]:
    report.append(f"\n### {name}")
    report.append("| ex4 top-N | in ex8 top-50 | in ex8 top-100 | in ex8 top-150 | in ex8 top-200 |")
    report.append("|---|--:|--:|--:|--:|")
    for N in (50, 100, 150, 200):
        top4 = set(d[d[r4] <= N].index)
        row = [f"**{N}**"]
        for M in (50, 100, 150, 200):
            top8 = set(d[d[r8] <= M].index)
            row.append(f"{len(top4 & top8)}")
        report.append("| " + " | ".join(row) + " |")
report.append("")
report.append("## Per-pocket score stability (Spearman, ex4 vs ex8)")
conf_rows = []
for c in w4.columns.drop(["best", "mean"]):
    if c not in w8.columns:
        continue
    s = pd.concat([w4[c], w8[c]], axis=1, keys=["a", "b"]).dropna()
    r = s["a"].corr(s["b"], method="spearman")
    conf_rows.append((c, r, len(s)))
report.append("| conformation | Spearman r | n |")
report.append("|---|---:|---:|")
for c, r, nn in conf_rows:
    report.append(f"| {c} | {r:.3f} | {nn} |")
report.append("")

# known failure case
if "CHEMBL776" in d.index:
    r = d.loc["CHEMBL776"]
    report.append("## Known failure case CHEMBL776 (orciprenaline)")
    report.append("")
    report.append(f"- ex4 best {r['ex4_best']:.2f} (rank {r['rank_ex4_best']}/{n}) -> "
                  f"ex8 best {r['ex8_best']:.2f} (rank {r['rank_ex8_best']}/{n}); "
                  f"best/HAC ex4 {r['ex4_ha']:.3f} -> ex8 {r['ex8_ha']:.3f}")
    report.append(f"- ex4 leaderboard absolute rank (raw+ha model): {r['enopt_ha_rank']:.0f} "
                  f"(old decoy-model rank {r['old_decoy_rank']:.0f})")
    report.append("")

# biggest movers
mv = d.reindex(d["rank_delta_best"].abs().sort_values(ascending=False).index).head(12)
report.append("## Biggest rank movers (best score)")
report.append("")
report.append("| ligand_id | ex4 best | ex8 best | rank ex4 | rank ex8 | d_rank | enopt_ha_rank | known |")
report.append("|---|---:|---:|---:|---:|---:|---:|---:|")
for lid, r in mv.iterrows():
    ha_rank = f"{r['enopt_ha_rank']:.0f}" if pd.notna(r['enopt_ha_rank']) else "-"
    report.append(f"| {lid} | {r['ex4_best']:.2f} | {r['ex8_best']:.2f} | {int(r['rank_ex4_best'])} | "
                  f"{int(r['rank_ex8_best'])} | {int(r['rank_delta_best']):+d} | {ha_rank} | {int(r['is_known_active'])} |")
report.append("")

txt = "\n".join(report) + "\n"
(OUT / "ex4_ex8_report.md").write_text(txt, encoding="utf-8")

# per-ligand csv
d.to_csv(OUT / "ex4_ex8_per_ligand.csv", encoding="utf-8")

# figure
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax, (a, b, t) in zip(axes,
                         [("ex4_best", "ex8_best", "best score (kcal/mol)"),
                          ("ex4_ha", "ex8_ha", "best score / heavy atom")]):
    ax.scatter(d[a], d[b], s=14, alpha=0.6)
    lo = min(d[a].min(), d[b].min()) - 0.5
    hi = max(d[a].max(), d[b].max()) + 0.5
    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_xlabel(f"ex4 {t}"); ax.set_ylabel(f"ex8 {t}")
    ax.set_title(f"{t}\nSpearman r = {d[a].corr(d[b], method='spearman'):.3f}")
if "CHEMBL776" in d.index:
    for ax, a, b in zip(axes, ["ex4_best", "ex4_ha"], ["ex8_best", "ex8_ha"]):
        ax.scatter([d.loc["CHEMBL776", a]], [d.loc["CHEMBL776", b]], marker="*", s=120, c="red")
fig.tight_layout()
fig.savefig(OUT / "fig_ex4_vs_ex8.png", dpi=150)
print("\n".join(report[:40]))
print("saved ->", OUT)