#!/usr/bin/env python3
"""logP-window comparison on the existing beta2 score matrix.

The window is the intersection of the active and decoy Crippen-logP
interquartile ranges. It is written before any model is fit. Hyperparameters
stay at the week-1 defaults. No new docking.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier

WS = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs")
REPO = Path(
    r"D:\repositories\Codex\2026-08-04\rrerax-beta2ar-enopt-gpcr-screening-https"
    r"\work\Beta2AR_EnOpt_GPCR_Screening-main"
)
OUT = WS / "outputs" / "logp_match_2026-10-16"
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
        "top_1pct_k": int(max(int(len(scores) * 0.01), 1)),
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

desc_rows = []
for lig in Xc.index:
    mol = Chem.MolFromSmiles(smiles.get(lig)) if isinstance(smiles.get(lig), str) else None
    if mol is None:
        desc_rows.append((lig, np.nan, np.nan, np.nan))
        continue
    desc_rows.append(
        (
            lig,
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            float(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
        )
    )
desc = pd.DataFrame(desc_rows, columns=["ligand_id", "mw", "logp", "hac"]).set_index("ligand_id")
desc = desc.loc[Xc.index]

active_logp = desc.loc[y == 1, "logp"].dropna()
decoy_logp = desc.loc[y == 0, "logp"].dropna()
active_q25 = float(np.quantile(active_logp, 0.25))
active_q75 = float(np.quantile(active_logp, 0.75))
decoy_q25 = float(np.quantile(decoy_logp, 0.25))
decoy_q75 = float(np.quantile(decoy_logp, 0.75))
lo = max(active_q25, decoy_q25)
hi = min(active_q75, decoy_q75)
if lo > hi:
    raise SystemExit("IQR intersection is empty; stop before fitting.")

keep = desc["logp"].between(lo, hi, inclusive="both")
Xk = Xc.loc[keep]
yk = y.loc[keep]
dk = desc.loc[keep]
window = {
    "rule": "intersection of active and decoy Crippen logP interquartile ranges",
    "quantile": "numpy.quantile linear",
    "active_q25": round(active_q25, 4),
    "active_q75": round(active_q75, 4),
    "decoy_q25": round(decoy_q25, 4),
    "decoy_q75": round(decoy_q75, 4),
    "logp_low_inclusive": round(lo, 4),
    "logp_high_inclusive": round(hi, 4),
    "n_before": int(len(Xc)),
    "n_actives_before": int(y.sum()),
    "n_decoys_before": int((1 - y).sum()),
    "active_logp_median_before": round(float(active_logp.median()), 4),
    "decoy_logp_median_before": round(float(decoy_logp.median()), 4),
    "n_after": int(len(Xk)),
    "n_actives_after": int(yk.sum()),
    "n_decoys_after": int((1 - yk).sum()),
    "active_logp_median_after": round(float(dk.loc[yk == 1, "logp"].median()), 4),
    "decoy_logp_median_after": round(float(dk.loc[yk == 0, "logp"].median()), 4),
    "top_1pct_k": int(max(int(len(Xk) * 0.01), 1)),
}
(OUT / "logp_match_window.json").write_text(json.dumps(window, indent=2), encoding="utf-8")
print(json.dumps(window, indent=2))

kept = dk.copy()
kept["y"] = yk
kept.reset_index().to_csv(OUT / "logp_match_kept.csv", index=False)

records = []
records.append(metrics("oriented_mean", -Xk.mean(axis=1), yk))
Xraw = Xk[cols]
Xdesc = dk[["mw", "logp", "hac"]]
for seed in SEEDS:
    raw_oof = oof_predict(Xraw, yk, seed)
    desc_oof = oof_predict(Xdesc, yk, seed)
    learned = metrics(f"learned_raw_seed{seed}", raw_oof, yk)
    descriptor = metrics(f"descriptor_only_seed{seed}", desc_oof, yk)
    learned["learned_ef1_above_descriptor"] = bool(learned["ef_1pct"] > descriptor["ef_1pct"])
    descriptor["learned_ef1_above_descriptor"] = learned["learned_ef1_above_descriptor"]
    records.append(learned)
    records.append(descriptor)

out = pd.DataFrame(records)
out.to_csv(OUT / "logp_match_metrics.csv", index=False)
learned = out[out.method.str.startswith("learned_raw")].reset_index(drop=True)
descriptor = out[out.method.str.startswith("descriptor_only")].reset_index(drop=True)
summary = {
    "window": window,
    "hyperparameters": {"n_estimators": N_EST, "learning_rate": LR, "max_depth": DEPTH},
    "seeds": SEEDS,
    "oriented_mean_auroc": float(out.loc[out.method == "oriented_mean", "auroc"].iloc[0]),
    "oriented_mean_ef1": float(out.loc[out.method == "oriented_mean", "ef_1pct"].iloc[0]),
    "learned_auroc_min": float(learned.auroc.min()),
    "learned_auroc_max": float(learned.auroc.max()),
    "learned_ef1_min": float(learned.ef_1pct.min()),
    "learned_ef1_max": float(learned.ef_1pct.max()),
    "descriptor_auroc_min": float(descriptor.auroc.min()),
    "descriptor_auroc_max": float(descriptor.auroc.max()),
    "descriptor_ef1_min": float(descriptor.ef_1pct.min()),
    "descriptor_ef1_max": float(descriptor.ef_1pct.max()),
    "learned_ef1_above_descriptor_every_seed": bool(learned.learned_ef1_above_descriptor.all()),
    "seeds_learned_ef1_above_descriptor": int(learned.learned_ef1_above_descriptor.sum()),
    "pass_gate": bool(learned.learned_ef1_above_descriptor.all()),
    "ef1_note": "top 1 percent contains only a handful of molecules, so this gate is unstable",
}
(OUT / "logp_match_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(out.to_string(index=False))
print(json.dumps(summary, indent=2))
