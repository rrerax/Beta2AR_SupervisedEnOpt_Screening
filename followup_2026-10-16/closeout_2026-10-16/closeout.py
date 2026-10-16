#!/usr/bin/env python3
"""Close-out after the logP window failed its gate.

Four calculations, no new docking:
1. Spearman of the five Vina scores with Crippen logP.
2. Descriptor-only models on the 48 library actives and the 158 added actives.
3. logP of the existing top-200 lists, plus a descriptor top-200 inside the labeled set.
4. Hold CHEMBL13 out of training and rank the full library.
"""
import json
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
OUT = WS / "outputs" / "closeout_2026-10-16"
OUT.mkdir(parents=True, exist_ok=True)
N_EST, LR, DEPTH = 15, 0.3, 6
SEEDS = list(range(42, 52))
HA_TOP = WS / "work" / "review_ha" / "out" / "top200_chemical_review.csv"
DECOY_RANK = WS / "outputs" / "decoy_validation" / "library_ranking_with_decoy_model_v2.csv"


def load_wide(path):
    df = pd.read_csv(path)
    return df.pivot_table(
        index="ligand_id",
        columns="conformation_id",
        values="best_score_kcal_mol",
        aggfunc="first",
    ).dropna()


def spearman(a, b):
    pair = pd.DataFrame({"a": np.asarray(a, dtype=float), "b": np.asarray(b, dtype=float)}).dropna()
    if len(pair) < 3:
        return np.nan, 0
    rho = float(pair.a.rank().corr(pair.b.rank()))
    return rho, int(len(pair))


def ef(scores, y, percent):
    scores = np.asarray(scores, dtype=float)
    y = np.asarray(y)
    k = max(int(len(scores) * percent), 1)
    order = np.argsort(-scores)
    return float(y[order[:k]].sum() / max(int(y.sum()), 1)) / percent


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


def props_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return np.nan, np.nan, np.nan
    return (
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        float(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
    )


lib = load_wide(REPO / "results/tables/docking_scores.csv")
exp = load_wide(WS / "outputs/data_prep/docking_scores_decoy_plus_expand.csv")
cols = list(lib.columns)
exp = exp[cols]
lib = lib[cols]
Xdec = exp[exp.index.str.startswith("decoy_")]
Xnew = exp[~exp.index.str.startswith("decoy_")]
lib48_ids = set(pd.read_csv(WS / "outputs/beta2ar_actives_in_library.csv")["molecule_chembl_id"])
act48 = set(lib.index) & lib48_ids
act158 = set(Xnew.index) - act48
act_all = act48 | set(Xnew.index)
Xc = pd.concat([lib[lib.index.isin(act_all)], Xnew, Xdec])
Xc = Xc[~Xc.index.duplicated(keep="first")][cols]
y = pd.Series(Xc.index.isin(act_all).astype(int), index=Xc.index)
assert len(act48) == 48
assert len(act158) == 158
assert int(y.sum()) == 206
assert "CHEMBL13" in act48 and y.loc["CHEMBL13"] == 1

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

prop_rows = []
for lig in Xc.index:
    mw, logp, hac = props_from_smiles(smiles.get(lig))
    prop_rows.append((lig, mw, logp, hac))
desc = pd.DataFrame(prop_rows, columns=["ligand_id", "mw", "logp", "hac"]).set_index("ligand_id").loc[Xc.index]
assert int(desc.logp.isna().sum()) == 0

# 1. Spearman
rho_rows = []
groups = {
    "all_3184": Xc.index,
    "actives_206": y[y == 1].index,
    "decoys_2978": y[y == 0].index,
}
for group_name, idx in groups.items():
    logp = desc.loc[idx, "logp"]
    for col in cols:
        rho, n = spearman(Xc.loc[idx, col], logp)
        rho_rows.append({"group": group_name, "feature": f"vina_{col}", "spearman_vs_logp": round(rho, 4), "n": n})
        rho, n = spearman(-Xc.loc[idx, col], logp)
        rho_rows.append({"group": group_name, "feature": f"oriented_{col}", "spearman_vs_logp": round(rho, 4), "n": n})
    rho, n = spearman(Xc.loc[idx, cols].mean(axis=1), logp)
    rho_rows.append({"group": group_name, "feature": "vina_mean", "spearman_vs_logp": round(rho, 4), "n": n})
    rho, n = spearman(-Xc.loc[idx, cols].mean(axis=1), logp)
    rho_rows.append({"group": group_name, "feature": "oriented_mean", "spearman_vs_logp": round(rho, 4), "n": n})
rho_df = pd.DataFrame(rho_rows)
rho_df.to_csv(OUT / "spearman_vina_logp.csv", index=False)

# 2. 48 vs 158, descriptor only
split_rows = []
Xdesc = desc[["mw", "logp", "hac"]]
for split_name, active_ids in (("library_48", act48), ("expanded_158", act158)):
    idx = Xc.index.isin(active_ids) | Xc.index.isin(Xdec.index)
    Xs = Xdesc.loc[idx]
    ys = pd.Series(Xc.index.isin(active_ids).astype(int), index=Xc.index).loc[idx]
    for seed in SEEDS:
        scores = oof_predict(Xs, ys, seed)
        split_rows.append(
            {
                "split": split_name,
                "seed": seed,
                "n": int(len(ys)),
                "n_actives": int(ys.sum()),
                "n_decoys": int((1 - ys).sum()),
                "auroc": round(float(roc_auc_score(ys, scores)), 4),
                "ef_1pct": round(ef(scores, ys, 0.01), 3),
                "ef_5pct": round(ef(scores, ys, 0.05), 3),
                "active_logp_median": round(float(desc.loc[list(active_ids), "logp"].median()), 4),
                "decoy_logp_median": round(float(desc.loc[Xdec.index, "logp"].median()), 4),
            }
        )
split_df = pd.DataFrame(split_rows)
split_df.to_csv(OUT / "split_48_158_metrics.csv", index=False)

# 3. Top-200 logP
rank = pd.read_csv(DECOY_RANK).set_index("ligand_id")
list_rows = []


def summarize_list(name, ids, score_note):
    ids = list(ids)
    missing = [i for i in ids if i not in smiles]
    vals = [props_from_smiles(smiles.get(i)) for i in ids]
    frame = pd.DataFrame(vals, index=ids, columns=["mw", "logp", "hac"])
    known = int(sum(i in act_all for i in ids))
    return {
        "list": name,
        "n": int(len(ids)),
        "missing_smiles": missing,
        "n_known_actives": known,
        "logp_median": round(float(frame.logp.median()), 4),
        "logp_q25": round(float(frame.logp.quantile(0.25)), 4),
        "logp_q75": round(float(frame.logp.quantile(0.75)), 4),
        "mw_median": round(float(frame.mw.median()), 2),
        "note": score_note,
    }, frame


active_ref = {
    "list": "labeled_actives_206",
    "n": 206,
    "missing_smiles": [],
    "n_known_actives": 206,
    "logp_median": round(float(desc.loc[y == 1, "logp"].median()), 4),
    "logp_q25": round(float(desc.loc[y == 1, "logp"].quantile(0.25)), 4),
    "logp_q75": round(float(desc.loc[y == 1, "logp"].quantile(0.75)), 4),
    "mw_median": round(float(desc.loc[y == 1, "mw"].median()), 2),
    "note": "Crippen logP on the 206 docked actives",
}
decoy_ref = {
    "list": "labeled_decoys_2978",
    "n": 2978,
    "missing_smiles": [],
    "n_known_actives": 0,
    "logp_median": round(float(desc.loc[y == 0, "logp"].median()), 4),
    "logp_q25": round(float(desc.loc[y == 0, "logp"].quantile(0.25)), 4),
    "logp_q75": round(float(desc.loc[y == 0, "logp"].quantile(0.75)), 4),
    "mw_median": round(float(desc.loc[y == 0, "mw"].median()), 2),
    "note": "Crippen logP on the 2978 docked decoys",
}
list_rows.extend([active_ref, decoy_ref])

ha = pd.read_csv(HA_TOP)
ha_summary, _ = summarize_list(
    "heavy_atom_model_top200",
    ha.ligand_id,
    "existing per-heavy-atom shortlist in work/review_ha/out/top200_chemical_review.csv",
)
list_rows.append(ha_summary)
decoy_top = rank.sort_values("enopt_decoy_rank").head(200).index
decoy_summary, _ = summarize_list(
    "learned_full_library_top200",
    decoy_top,
    "library_ranking_with_decoy_model_v2.csv ranks 1-200",
)
list_rows.append(decoy_summary)

seed42_scores = None
desc_medians = []
for seed in SEEDS:
    scores = pd.Series(oof_predict(Xdesc, y, seed), index=Xc.index)
    top = scores.sort_values(ascending=False).head(200)
    desc_medians.append(
        {
            "seed": seed,
            "logp_median": round(float(desc.loc[top.index, "logp"].median()), 4),
            "mw_median": round(float(desc.loc[top.index, "mw"].median()), 2),
            "n_actives": int(y.loc[top.index].sum()),
        }
    )
    if seed == 42:
        seed42_scores = top
desc_top_summary, _ = summarize_list(
    "descriptor_oof_top200_seed42_within_3184",
    seed42_scores.index,
    "OOF probabilities from MW, logP, and heavy-atom count; labeled 3184 only",
)
list_rows.append(desc_top_summary)
pd.DataFrame(desc_medians).to_csv(OUT / "descriptor_top200_by_seed.csv", index=False)
pd.DataFrame(list_rows).to_csv(OUT / "top200_logp.csv", index=False)
seed42_scores.rename("oof_score").reset_index().to_csv(OUT / "descriptor_top200_seed42.csv", index=False)

# 4. CHEMBL13 holdout
held_rows = []
for seed in SEEDS:
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    for tr, te in skf.split(Xc, y):
        if "CHEMBL13" not in set(Xc.index[te]):
            continue
        pos = int(y.iloc[tr].sum())
        neg = int(len(tr) - pos)
        model = XGBClassifier(
            n_estimators=N_EST,
            learning_rate=LR,
            max_depth=DEPTH,
            scale_pos_weight=neg / pos,
            eval_metric="auc",
            random_state=seed,
            n_jobs=1,
        )
        model.fit(Xc.iloc[tr], y.iloc[tr])
        p_lib = pd.Series(model.predict_proba(lib[cols])[:, 1], index=lib.index)
        order = p_lib.rank(ascending=False, method="min")
        held_rows.append(
            {
                "seed": seed,
                "p_13_held_out": round(float(p_lib.loc["CHEMBL13"]), 4),
                "rank_13_held_out": int(order.loc["CHEMBL13"]),
                "n_train_actives": pos,
            }
        )
held = pd.DataFrame(held_rows)
held.to_csv(OUT / "chembl13_holdout_ranks.csv", index=False)
oriented_order = (-lib[cols].mean(axis=1)).rank(ascending=False, method="min")

summary = {
    "spearman_oriented_mean_vs_logp": rho_df[
        (rho_df.group == "all_3184") & (rho_df.feature == "oriented_mean")
    ].iloc[0].to_dict(),
    "spearman_oriented_mean_within_actives": rho_df[
        (rho_df.group == "actives_206") & (rho_df.feature == "oriented_mean")
    ].iloc[0].to_dict(),
    "spearman_oriented_mean_within_decoys": rho_df[
        (rho_df.group == "decoys_2978") & (rho_df.feature == "oriented_mean")
    ].iloc[0].to_dict(),
    "split_48": {
        "auroc_min": float(split_df.loc[split_df.split == "library_48", "auroc"].min()),
        "auroc_max": float(split_df.loc[split_df.split == "library_48", "auroc"].max()),
        "ef1_min": float(split_df.loc[split_df.split == "library_48", "ef_1pct"].min()),
        "ef1_max": float(split_df.loc[split_df.split == "library_48", "ef_1pct"].max()),
    },
    "split_158": {
        "auroc_min": float(split_df.loc[split_df.split == "expanded_158", "auroc"].min()),
        "auroc_max": float(split_df.loc[split_df.split == "expanded_158", "auroc"].max()),
        "ef1_min": float(split_df.loc[split_df.split == "expanded_158", "ef_1pct"].min()),
        "ef1_max": float(split_df.loc[split_df.split == "expanded_158", "ef_1pct"].max()),
    },
    "top200": list_rows,
    "descriptor_top200_logp_median_across_seeds": {
        "min": min(r["logp_median"] for r in desc_medians),
        "max": max(r["logp_median"] for r in desc_medians),
        "n_actives_min": min(r["n_actives"] for r in desc_medians),
        "n_actives_max": max(r["n_actives"] for r in desc_medians),
    },
    "chembl13": {
        "in_sample_learned_rank": int(rank.loc["CHEMBL13", "enopt_decoy_rank"]),
        "oriented_mean_rank": int(oriented_order.loc["CHEMBL13"]),
        "heldout_rank_min": int(held.rank_13_held_out.min()),
        "heldout_rank_max": int(held.rank_13_held_out.max()),
        "heldout_probability_min": float(held.p_13_held_out.min()),
        "heldout_probability_max": float(held.p_13_held_out.max()),
    },
}
(OUT / "closeout_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(rho_df.to_string(index=False))
print(split_df.groupby("split")[["auroc", "ef_1pct"]].agg(["min", "max"]).to_string())
print(pd.DataFrame(list_rows).to_string(index=False))
print(held.to_string(index=False))
print(json.dumps(summary, indent=2))
