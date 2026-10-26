"""Read the finished amine screen. No model is trained.

The public table keeps catalog ids, Vina scores, and descriptors. Supplier
SMILES and docked coordinates stay in the local library folder.
"""
import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
DOCK = Path(r"D:\repositories\Enamine_GPCR_Library_plated\docking")
PROPS = Path(r"D:\repositories\Enamine_GPCR_Library_plated\prepared\gpr53_pH74_table.csv")
PUBLIC = HERE / "amine_five_scores.csv"
COLS = [
    "inactive_carazolol_2RH1",
    "inactive_carazolol_5D5A",
    "active_gs_complex_3SN6",
    "active_bi167107_4LDE",
    "active_hbi_4LDL",
]
CONTROLS = ["CHEMBL434", "CHEMBL776"]

scores = pd.read_csv(DOCK / "docking_scores.csv", usecols=["ligand_id", "conformation_id", "best_score_kcal_mol", "status"])
if not (scores["status"] == "ready").all():
    raise SystemExit("score table has rows that are not ready")
wide = scores.pivot(index="ligand_id", columns="conformation_id", values="best_score_kcal_mol")
missing = [column for column in COLS if column not in wide.columns]
if missing:
    raise SystemExit(f"missing conformations: {missing}")
wide = wide[COLS].apply(pd.to_numeric, errors="coerce").dropna()
props = pd.read_csv(
    PROPS,
    usecols=["ligand_id", "crippen_logp", "mw", "hbd", "hba", "hac", "formal_charge_pH74", "pains", "brenk"],
)
table = wide.join(props.set_index("ligand_id"), how="left")
if table["crippen_logp"].isna().any():
    raise SystemExit(f"ligands without properties: {int(table['crippen_logp'].isna().sum())}")

table["five_mean"] = table[COLS].mean(axis=1)
table["oriented"] = -table["five_mean"]
table["rank"] = table["five_mean"].rank(method="min", ascending=True).astype(int)
table["inactive_mean"] = table[COLS[:2]].mean(axis=1)
table["active_mean"] = table[COLS[2:]].mean(axis=1)
n = len(table)
rho, _ = spearmanr(table["oriented"], table["crippen_logp"])
rho_4lde, _ = spearmanr(-table["active_bi167107_4LDE"], table["crippen_logp"])
top = table.nsmallest(max(int(round(n * 0.01)), 1), "five_mean")

def ligand_record(ligand_id: str) -> dict:
    row = table.loc[ligand_id]
    return {
        "ligand_id": ligand_id,
        "rank": int(row["rank"]),
        "library_size": n,
        "five_mean": round(float(row["five_mean"]), 3),
        "crippen_logp": round(float(row["crippen_logp"]), 3),
        "mw": round(float(row["mw"]), 3),
        "hbd": int(row["hbd"]),
        "hba": int(row["hba"]),
        "hac": int(row["hac"]),
        "formal_charge_pH74": int(row["formal_charge_pH74"]),
        "pains": "" if pd.isna(row["pains"]) else str(row["pains"]),
        "brenk": "" if pd.isna(row["brenk"]) else str(row["brenk"]),
        "inactive_mean": round(float(row["inactive_mean"]), 3),
        "active_mean": round(float(row["active_mean"]), 3),
        "scores": {column: round(float(row[column]), 3) for column in COLS},
    }

summary = {
    "ligands": n,
    "scores": int(len(scores)),
    "model_trained": False,
    "score": "oriented five-conformation mean; more negative raw mean ranks first",
    "spearman_oriented_mean_vs_crippen_logp": round(float(rho), 3),
    "spearman_oriented_4lde_vs_crippen_logp": round(float(rho_4lde), 3),
    "library_logp_median": round(float(table["crippen_logp"].median()), 3),
    "top_1pct_n": int(len(top)),
    "top_1pct_logp_median": round(float(top["crippen_logp"].median()), 3),
    "controls": [ligand_record(ligand_id) for ligand_id in CONTROLS],
}
gap = abs(table.loc["CHEMBL434", "five_mean"] - table.loc["CHEMBL776", "five_mean"])
summary["control_five_mean_abs_gap"] = round(float(gap), 3)
(DOCK / "amine_screen_read.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
(HERE / "amine_screen_read.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
public = table.reset_index()[
    ["ligand_id", *COLS, "five_mean", "rank", "crippen_logp", "mw", "hbd", "hba", "hac", "formal_charge_pH74", "pains", "brenk"]
].sort_values(["rank", "ligand_id"])
public.to_csv(PUBLIC, index=False)
print(json.dumps(summary, indent=2))
