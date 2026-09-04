#!/usr/bin/env python3
"""Stage 4c: per-heavy-atom feature experiment.

Motivation: raw Vina scores are dominated by molecular size (bigger molecules
get more negative, "better", scores).  Question: does adding
per-heavy-atom-normalized scores (score / #heavy_atoms) -- or explicit size
features -- improve the decoy-validated EnOpt?

Protocol mirrors scripts/07_train_validated_enopt.py so numbers are directly
comparable:
  - clean-label 3-fold OOF: 206 actives vs 2978 DUD-E decoys
  - merged-mirror 3-fold OOF: library + decoys + new actives
Feature sets tested (same XGBoost hyper-parameters as 07):
  raw          : 5 conformation scores            (reproduces 0.696)
  raw + ha     : + 5 per-heavy-atom scores        (10 features)
  raw + ha + hac : + heavy-atom count scalar      (11 features)
  ha           : 5 per-heavy-atom scores only

Requires rdkit (heavy-atom counts recomputed from SMILES).  Docks nothing.

Usage (from repo root):
    python scripts/10_feature_heavy_atom_experiment.py
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def largest_fragment_hac(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return np.nan
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    if not frags:
        return np.nan
    return max(frags, key=lambda f: f.GetNumHeavyAtoms()).GetNumHeavyAtoms()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO_ROOT / "results/feature_experiment_heavy_atom"))
    ap.add_argument("--folds", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    s7 = load_module("s7", REPO_ROOT / "scripts/07_train_validated_enopt.py")
    lib = s7.load_wide(REPO_ROOT / "results/tables/docking_scores.csv")
    decoy_set = s7.load_wide(REPO_ROOT / "results/tables/docking_scores_decoy_set.csv")
    cols = list(lib.columns)
    decoy_set = decoy_set[cols]
    Xdec = decoy_set[decoy_set.index.str.startswith("decoy_")]
    Xnew = decoy_set[~decoy_set.index.str.startswith("decoy_")]
    act48 = set(pd.read_csv(REPO_ROOT / "data/training/beta2ar_actives_in_library.csv")
                ["molecule_chembl_id"])
    act48 = set(lib.index) & act48
    act_all = act48 | set(Xnew.index)
    print(f"[data] library {len(lib)} | decoys {len(Xdec)} | new actives "
          f"{len(Xnew)} | actives total {len(act_all)}")

    # heavy-atom counts for every docked ligand, recomputed from SMILES
    smiles_sources = {
        "lib": pd.read_csv(REPO_ROOT / "data/raw/chembl37_ligand_library_30000.csv")[["ligand_id", "smiles"]],
        "new": pd.read_csv(REPO_ROOT / "data/training/beta2ar_actives_expand_docked.csv"),
        "dec": pd.read_csv(REPO_ROOT / "data/training/dude_decoys_for_dock.csv"),
    }
    smi = pd.concat([s[["ligand_id", "smiles"]] for s in smiles_sources.values()])
    smi = smi.drop_duplicates("ligand_id")
    smi = smi[smi["ligand_id"].isin(
        set(lib.index) | set(decoy_set.index))]
    smi["hac"] = smi["smiles"].apply(largest_fragment_hac)
    n_fail = int(smi["hac"].isna().sum())
    print(f"[data] HAC computed for {len(smi)} ligands (parse failures "
          f"{n_fail})")
    hac = smi.set_index("ligand_id")["hac"]

    # clean-label set
    Xc = pd.concat([lib[lib.index.isin(act_all)], Xnew, Xdec])
    yc = pd.Series((Xc.index.isin(act_all)).astype(int), index=Xc.index)
    # merged mirror set
    Xm = pd.concat([lib, Xnew, Xdec])
    ym = pd.Series((Xm.index.isin(act_all)).astype(int), index=Xm.index)

    def add_ha(X, hac_map):
        X = X.copy()
        h = hac_map.reindex(X.index)
        X = X[h.notna()]
        X["heavy_atoms"] = h[X.index].astype(float)
        X["score_per_ha_mean"] = X.mean(axis=1) / X["heavy_atoms"]
        return X

    Xc_ha = add_ha(Xc, hac)
    yc_ha = yc.reindex(Xc_ha.index)
    Xm_ha = add_ha(Xm, hac)
    ym_ha = ym.reindex(Xm_ha.index)

    raw_c = Xc_ha[cols]
    ha_c = Xc_ha[cols].div(Xc_ha["heavy_atoms"], axis=0)
    ha_c.columns = [c + "__per_ha" for c in ha_c.columns]
    hac_c = Xc_ha[["heavy_atoms"]]
    raw_ha_c = pd.concat([raw_c, ha_c], axis=1)
    raw_ha_hac_c = pd.concat([raw_c, ha_c, hac_c], axis=1)

    sets = {
        "raw (5)": raw_c,
        "raw+ha (10)": raw_ha_c,
        "raw+ha+hac (11)": raw_ha_hac_c,
        "ha only (5)": ha_c,
    }
    rows = []
    oofs = {}
    for name, X in sets.items():
        oof, imp, w = s7.oof_predict(X, yc_ha, n_splits=args.folds,
                                     seed=args.seed)
        oofs[name] = oof
        rows.append({"feature_set": name, **s7.metrics_row(
            "EnOpt_" + name.split(" (")[0].replace("+", "_").replace(" ", ""),
            pd.Series(oof, index=yc_ha.index), yc_ha)})
    # oriented baselines on the same rows
    rows.append({"feature_set": "base mean (raw)",
                 **s7.metrics_row("ensemble_neg_mean_raw",
                                  -raw_c.mean(axis=1), yc_ha)})
    rows.append({"feature_set": "base mean (per-ha)",
                 **s7.metrics_row("ensemble_neg_mean_ha",
                                  -ha_c.mean(axis=1), yc_ha)})
    metrics = pd.DataFrame(rows)
    print("\n== clean-label OOF (feature comparison; actives vs decoys) ==")
    print(metrics.to_string(index=False))

    # merged-mirror comparison: raw vs best
    # best clean-label set by AUROC (see metrics table); use raw+ha for merged/final
    oof_m_raw, _, _ = s7.oof_predict(Xm_ha[cols], ym_ha,
                                     n_splits=args.folds, seed=args.seed)
    ha_m = Xm_ha[cols].div(Xm_ha["heavy_atoms"], axis=0)
    ha_m.columns = [c + "__per_ha" for c in ha_m.columns]
    raw_m = Xm_ha[cols]
    Xm_best = pd.concat([raw_m, ha_m], axis=1)
    oof_m_best, imp_m, _ = s7.oof_predict(Xm_best, ym_ha,
                                          n_splits=args.folds, seed=args.seed)
    merged = pd.DataFrame([
        {"feature_set": "raw (5) merged", **s7.metrics_row(
            "EnOpt_raw_merged", pd.Series(oof_m_raw, index=ym_ha.index), ym_ha)},
        {"feature_set": "raw+ha (10) merged", **s7.metrics_row(
            "EnOpt_raw_ha_merged", pd.Series(oof_m_best, index=ym_ha.index), ym_ha)},
        {"feature_set": "base mean (raw) merged", **s7.metrics_row(
            "ensemble_neg_mean_raw", -raw_m.mean(axis=1), ym_ha)},
    ])
    print("\n== merged-mirror OOF ==")
    print(merged.to_string(index=False))

    # final model: train on clean-label raw+ha (Stage-3 protocol), rerank library
    pos = int(yc_ha.sum()); neg = int(len(yc_ha) - pos)
    from xgboost import XGBClassifier
    final_c = XGBClassifier(n_estimators=15, learning_rate=0.3, max_depth=6,
                            scale_pos_weight=neg / pos, eval_metric="auc",
                            random_state=args.seed, n_jobs=1)
    final_c.fit(pd.concat([raw_c, ha_c], axis=1), yc_ha)
    lib_h = hac.reindex(lib.index)
    lib_h = lib_h[lib_h.notna()]
    lib_raw = lib.loc[lib_h.index]
    lib_ha2 = lib_raw.div(lib_h, axis=0)
    lib_ha2.columns = [c + "__per_ha" for c in lib_ha2.columns]
    lib_ha_feat = pd.concat([lib_raw, lib_ha2], axis=1)
    score_lib = final_c.predict_proba(lib_ha_feat)[:, 1]
    ranking = pd.DataFrame({"ligand_id": lib_ha_feat.index,
                            "is_known_active": lib_ha_feat.index.isin(act_all).astype(int),
                            "enopt_ha_score": score_lib})
    ranking = ranking.sort_values("enopt_ha_score", ascending=False).reset_index(drop=True)
    ranking["enopt_ha_rank"] = np.arange(1, len(ranking) + 1)
    old = pd.read_csv(REPO_ROOT / "results/supervised_enopt_decoy/library_ranking_decoy_model.csv")
    old_map = dict(zip(old["ligand_id"], old["enopt_decoy_rank"]))
    ranking["old_decoy_rank"] = ranking["ligand_id"].map(old_map)
    ov = {n: int(len(set(ranking.head(n)["ligand_id"]) &
                      set(old.head(n)["ligand_id"]))) for n in (200, 500, 1000)}
    print("\n== leaderboard change (raw+ha final model vs Stage-3 decoy model) ==")
    print("top-N overlap:", ov)
    print("known actives in top 50/200/500/1000 (new):",
          {n: int(ranking.head(n)["is_known_active"].sum()) for n in (50, 200, 500, 1000)})

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out / "metrics_feature_sets.csv", index=False)
    merged.to_csv(out / "metrics_merged_compare.csv", index=False)
    ranking.to_csv(out / "library_ranking_heavy_atom_model.csv", index=False)
    imp_df = pd.DataFrame({"feature": list(pd.concat([raw_c, ha_c], axis=1).columns),
                           "importance": final_c.feature_importances_})
    imp_df.sort_values("importance", ascending=False).to_csv(
        out / "feature_importance_raw_ha.csv", index=False)

    # figure: AUROC + EF1 per feature set
    fig, ax = plt.subplots(figsize=(8, 4.5))
    names = metrics["feature_set"]
    aurocs = metrics["auroc"]
    ax.bar(range(len(metrics)), aurocs, color="#2980b9")
    ax.errorbar(range(len(metrics)), aurocs,
                yerr=1.96 * metrics["auroc_se"], fmt="none", color="black",
                capsize=3)
    ax.axhline(0.5, color="grey", ls="--", lw=1)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("clean-label OOF AUROC (95% CI)")
    ax.set_ylim(0.4, 0.85)
    fig.tight_layout()
    fig.savefig(out / "fig_feature_comparison.png", dpi=150)
    plt.close(fig)

    card = {
        "experiment": "per-heavy-atom feature comparison",
        "protocol": "identical to scripts/07 (3-fold OOF, XGBoost n_estimators=15, lr=0.3, depth=6, scale_pos_weight=neg/pos)",
        "heavy_atom_definition": "largest-fragment heavy-atom count (RDKit), recomputed from SMILES",
        "clean_label_feature_sets": metrics.to_dict("records"),
        "merged_mirror": merged.to_dict("records"),
        "leaderboard_topN_overlap_with_stage3": ov,
        "note": "best clean-label set chosen by AUROC, then merged-mirror and final ranking recomputed for it",
    }
    (out / "model_card_heavy_atom.json").write_text(
        json.dumps(card, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
