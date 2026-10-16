#!/usr/bin/env python3
"""Week 3: is CHEMBL434's rank 14 an in-sample score?

Side-chain relaxation is not run. The same XGBoost settings are kept.
For each seed, the one cross-validation fold that holds CHEMBL434 out
scores both CHEMBL434 and CHEMBL776, then ranks the full library.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier

WS = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs")
REPO = Path(
    r"D:\repositories\Codex\2026-08-04\rrerax-beta2ar-enopt-gpcr-screening-https"
    r"\work\Beta2AR_EnOpt_GPCR_Screening-main"
)
OUT = WS / "outputs" / "week3_2026-10-16"
OUT.mkdir(parents=True, exist_ok=True)
N_EST, LR, DEPTH = 15, 0.3, 6
SEEDS = list(range(42, 52))
PAIR = ("CHEMBL434", "CHEMBL776")
CATECHOL = Chem.MolFromSmarts("c1ccc(O)c(O)c1")
RESORCINOL = Chem.MolFromSmarts("c1cc(O)cc(O)c1")


def load_wide(path):
    df = pd.read_csv(path)
    return df.pivot_table(
        index="ligand_id",
        columns="conformation_id",
        values="best_score_kcal_mol",
        aggfunc="first",
    ).dropna()


def pattern(smiles, query):
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return bool(mol and mol.HasSubstructMatch(query))


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
Xc = Xc[~Xc.index.duplicated(keep="first")][cols]
y = pd.Series(Xc.index.isin(act_all).astype(int), index=Xc.index)
assert "CHEMBL434" in Xc.index and y.loc["CHEMBL434"] == 1
assert "CHEMBL776" in lib.index and "CHEMBL776" not in act_all

manifest = pd.read_csv(REPO / "results/tables/ligand_manifest.csv").set_index("ligand_id")
rank = pd.read_csv(
    WS / "outputs/decoy_validation/library_ranking_with_decoy_model_v2.csv"
).set_index("ligand_id")

score_rows = []
for seed in SEEDS:
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    for tr, te in skf.split(Xc, y):
        test_ids = set(Xc.index[te])
        if "CHEMBL434" not in test_ids:
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
        p_lib = model.predict_proba(lib[cols])[:, 1]
        order = pd.Series(p_lib, index=lib.index).rank(ascending=False, method="min")
        p_434 = float(model.predict_proba(Xc.loc[["CHEMBL434"], cols])[:, 1][0])
        p_776 = float(model.predict_proba(lib.loc[["CHEMBL776"], cols])[:, 1][0])
        score_rows.append(
            {
                "seed": seed,
                "p_434_held_out": round(p_434, 4),
                "p_776": round(p_776, 4),
                "rank_434_held_out": int(order.loc["CHEMBL434"]),
                "rank_776": int(order.loc["CHEMBL776"]),
                "n_train_actives": pos,
            }
        )

held = pd.DataFrame(score_rows)
held.to_csv(OUT / "week3_holdout_ranks.csv", index=False)

actives = manifest.loc[manifest.index.isin(act_all), ["smiles"]].copy()
actives["catechol"] = actives.smiles.map(lambda s: pattern(s, CATECHOL))
actives["resorcinol"] = actives.smiles.map(lambda s: pattern(s, RESORCINOL))
actives["learned_rank"] = rank["enopt_decoy_rank"]
actives["in_library"] = actives.index.isin(lib.index)

pair = lib.loc[list(PAIR), cols].copy()
pair.insert(0, "in_training_206", [True, False])
pair.insert(1, "full_model_rank", [int(rank.loc[i, "enopt_decoy_rank"]) for i in PAIR])
pair.to_csv(OUT / "week3_pair_scores.csv")

summary = {
    "CHEMBL434_in_206": True,
    "CHEMBL776_in_206": False,
    "full_model_rank_434": int(rank.loc["CHEMBL434", "enopt_decoy_rank"]),
    "full_model_rank_776": int(rank.loc["CHEMBL776", "enopt_decoy_rank"]),
    "heldout_rank_434_min": int(held.rank_434_held_out.min()),
    "heldout_rank_434_max": int(held.rank_434_held_out.max()),
    "heldout_rank_776_min": int(held.rank_776.min()),
    "heldout_rank_776_max": int(held.rank_776.max()),
    "heldout_434_ahead_of_776_every_seed": bool(
        (held.rank_434_held_out < held.rank_776).all()
    ),
    "score_vector_max_abs_diff": round(
        float((lib.loc["CHEMBL434", cols] - lib.loc["CHEMBL776", cols]).abs().max()), 3
    ),
    "training_catechols": int(actives.catechol.sum()),
    "training_resorcinols": int(actives.resorcinol.sum()),
    "library_catechol_actives_best_learned_rank": (
        None
        if actives.loc[actives.catechol & actives.in_library, "learned_rank"].dropna().empty
        else int(actives.loc[actives.catechol & actives.in_library, "learned_rank"].min())
    ),
}
(OUT / "week3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(pair.round(3).to_string())
print(held.to_string(index=False))
print(json.dumps(summary, indent=2))
print("catechol actives in library:")
print(
    actives.loc[actives.catechol & actives.in_library, ["learned_rank", "smiles"]]
    .sort_values("learned_rank")
    .head(15)
    .to_string()
)
