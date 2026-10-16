#!/usr/bin/env python3
"""Week-1 controls on the existing beta2 score matrix.

Single conformations, oriented consensus, learned ranker, descriptor-only
model, and ten repeated 3-fold splits. Hyperparameters stay at the EnOpt
defaults. No new docking.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier

WS = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs")
REPO = Path(
    r"D:\repositories\Codex\2026-08-04\rrerax-beta2ar-enopt-gpcr-screening-https"
    r"\work\Beta2AR_EnOpt_GPCR_Screening-main"
)
OUT = WS / "outputs" / "week1_2026-10-16"
OUT.mkdir(parents=True, exist_ok=True)

N_EST, LR, DEPTH = 15, 0.3, 6
SEEDS = list(range(42, 52))


def load_wide(path):
    df = pd.read_csv(path)
    return df.pivot_table(
        index="ligand_id",
        columns="conformation_id",
        values="best_score_kcal_mol",
        aggfunc="first",
    ).dropna()


def auc_ci(auc, n_pos, n_neg):
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc ** 2 / (1.0 + auc)
    se = math.sqrt(
        (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2) + (n_neg - 1) * (q2 - auc ** 2))
        / (n_pos * n_neg)
    )
    return se, auc - 1.96 * se, auc + 1.96 * se


def ef(scores, y, percent):
    scores = np.asarray(scores, dtype=float)
    y = np.asarray(y)
    k = max(int(len(scores) * percent), 1)
    order = np.argsort(-scores)
    return float(y[order[:k]].sum() / max(int(y.sum()), 1)) / percent


def metrics(name, scores, y):
    y = np.asarray(y)
    auc = float(roc_auc_score(y, scores))
    se, lo, hi = auc_ci(auc, int(y.sum()), int(len(y) - y.sum()))
    return {
        "method": name,
        "auroc": round(auc, 4),
        "auroc_se": round(se, 4),
        "auroc_ci95_lo": round(lo, 3),
        "auroc_ci95_hi": round(hi, 3),
        "ef_1pct": round(ef(scores, y, 0.01), 3),
        "ef_5pct": round(ef(scores, y, 0.05), 3),
    }


def oof_predict(X, y, seed):
    y = pd.Series(np.asarray(y), index=X.index)
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    oof = np.full(len(X), np.nan)
    pos = int(y.sum())
    neg = int(len(y) - pos)
    for tr, te in skf.split(X, y):
        model = XGBClassifier(
            n_estimators=N_EST,
            learning_rate=LR,
            max_depth=DEPTH,
            scale_pos_weight=neg / pos,
            eval_metric="auc",
            random_state=seed,
            n_jobs=1,
        )
        model.fit(X.iloc[tr], y.iloc[tr])
        oof[te] = model.predict_proba(X.iloc[te])[:, 1]
    return oof


def mol_from_smiles(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles)


lib = load_wide(REPO / "results/tables/docking_scores.csv")
exp = load_wide(WS / "outputs/data_prep/docking_scores_decoy_plus_expand.csv")
cols = list(lib.columns)
exp = exp[cols]
lib = lib[cols]
Xdec = exp[exp.index.str.startswith("decoy_")]
Xnew = exp[~exp.index.str.startswith("decoy_")]
lib48 = set(pd.read_csv(WS / "outputs/beta2ar_actives_in_library.csv")["molecule_chembl_id"])
act48 = set(lib.index) & lib48
act_all = act48 | set(Xnew.index)
Xc = pd.concat([lib[lib.index.isin(act_all)], Xnew, Xdec])
Xc = Xc[~Xc.index.duplicated(keep="first")]
y = pd.Series(Xc.index.isin(act_all).astype(int), index=Xc.index)

smiles = {}
for _, row in pd.read_csv(WS / "outputs/beta2ar_actives_in_library.csv").iterrows():
    smiles[row["molecule_chembl_id"]] = row["smiles"]
for _, row in pd.read_csv(WS / "outputs/data_prep/beta2ar_actives_expand_selected.csv").iterrows():
    smiles[row["ligand_id"]] = row["smiles"]
manifest = pd.read_csv(REPO / "results/tables/ligand_manifest.csv", usecols=["ligand_id", "smiles"])
for _, row in manifest.iterrows():
    smiles.setdefault(row["ligand_id"], row["smiles"])
for _, row in pd.read_csv(WS / "outputs/data_prep/dude_decoys_for_dock.csv").iterrows():
    smiles[row["ligand_id"]] = row["smiles"]

records = []
desc_rows = []
missing = []
for lig in Xc.index:
    mol = mol_from_smiles(smiles.get(lig))
    if mol is None:
        missing.append(lig)
        desc_rows.append((lig, np.nan, np.nan, np.nan, None))
        continue
    desc_rows.append(
        (
            lig,
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            float(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
            MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False),
        )
    )
desc = pd.DataFrame(desc_rows, columns=["ligand_id", "mw", "logp", "hac", "scaffold"]).set_index("ligand_id")
desc = desc.loc[Xc.index]

# Consensus and single conformations. More negative Vina score is better.
records.append(metrics("oriented_mean", -Xc.mean(axis=1), y))
records.append(metrics("oriented_best", -Xc.min(axis=1), y))
for col in cols:
    records.append(metrics(f"single_{col}", -Xc[col], y))

raw_oof = {}
desc_oof = {}
Xraw = Xc[cols]
Xdesc = desc[["mw", "logp", "hac"]]
usable = Xdesc.dropna().index
print(f"molecules {len(Xc)} actives {int(y.sum())} decoys {int((1 - y).sum())} missing_smiles {len(missing)}")
print(f"descriptor rows {len(usable)}")

for seed in SEEDS:
    raw_oof[seed] = oof_predict(Xraw, y, seed)
    desc_oof[seed] = oof_predict(Xdesc.loc[usable], y.loc[usable], seed)
    records.append(metrics(f"learned_raw_seed{seed}", raw_oof[seed], y))
    records.append(metrics(f"descriptor_only_seed{seed}", desc_oof[seed], y.loc[usable]))

# Scaffolds among actives only.
active_scaffolds = desc.loc[y[y == 1].index, "scaffold"].dropna()
n_scaffolds = int(active_scaffolds.nunique())
sizes = active_scaffolds.value_counts()
can_fill = n_scaffolds >= 3 and int((sizes > 0).sum()) >= 3

out = pd.DataFrame(records)
out.to_csv(OUT / "week1_metrics.csv", index=False)

learned = out[out.method.str.startswith("learned_raw")]
descriptor = out[out.method.str.startswith("descriptor_only")]
mean_row = out[out.method == "oriented_mean"].iloc[0]
summary = {
    "n_molecules": int(len(Xc)),
    "n_actives": int(y.sum()),
    "n_decoys": int((1 - y).sum()),
    "missing_smiles": missing,
    "hyperparameters": {"n_estimators": N_EST, "learning_rate": LR, "max_depth": DEPTH},
    "seeds": SEEDS,
    "oriented_mean_ef1": float(mean_row.ef_1pct),
    "learned_auroc_min": float(learned.auroc.min()),
    "learned_auroc_max": float(learned.auroc.max()),
    "learned_ef1_min": float(learned.ef_1pct.min()),
    "learned_ef1_max": float(learned.ef_1pct.max()),
    "learned_ef1_above_mean_every_seed": bool((learned.ef_1pct > mean_row.ef_1pct).all()),
    "descriptor_auroc_min": float(descriptor.auroc.min()),
    "descriptor_auroc_max": float(descriptor.auroc.max()),
    "descriptor_ef1_min": float(descriptor.ef_1pct.min()),
    "descriptor_ef1_max": float(descriptor.ef_1pct.max()),
    "descriptor_ef1_below_learned_every_seed": bool(
        descriptor.ef_1pct.max() < learned.ef_1pct.min()
    ),
    "active_scaffold_count": n_scaffolds,
    "active_scaffold_singleton_count": int((sizes == 1).sum()),
    "largest_active_scaffold": int(sizes.max()) if len(sizes) else 0,
    "scaffold_split_decision": "可以做" if can_fill else "不做",
}
(OUT / "week1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(out.to_string(index=False))
print(json.dumps(summary, indent=2))
