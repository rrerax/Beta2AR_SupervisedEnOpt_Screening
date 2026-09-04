#!/usr/bin/env python3
"""Stage 3: DUD-E decoy-validated supervised EnOpt.

Fixes the two known weaknesses of the Stage-2 model card (06):
  1) negatives were "assumed inactive" library molecules (no experimental basis)
  2) only 48 positives -> too few for a reliable AUROC estimate

This script docks nothing. It uses already-computed score tables:
  - results/tables/docking_scores.csv            : 29,865-compound library
  - results/tables/docking_scores_decoy_set.csv  : 2,978 DUD-E decoys + 158 newly
                                                   docked ChEMBL actives (5 pockets each)
  - data/training/beta2ar_actives_in_library.csv : 48 known actives in the library

Reported numbers:
  - clean-label OOF: 206 actives vs 2,978 decoys, 3-fold, out-of-fold AUROC with
    Hanley-McNeil 95% CI, EF1%/EF5%, and baselines (simple mean/best, correctly
    oriented: more-negative = better).
  - merged mirror OOF: old 06 accounting on the library + decoys + new actives.
  - final model (fit on actives + decoys) re-ranks the full 29,865 library.

Usage (from repo root):
    python scripts/07_train_validated_enopt.py
"""
import argparse
import json
import math
import numpy as np
import pandas as pd
from pathlib import Path
import joblib
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_wide(csv_path):
    dock = pd.read_csv(csv_path)
    return dock.pivot_table(index="ligand_id", columns="conformation_id",
                            values="best_score_kcal_mol", aggfunc="first").dropna()


def auc_ci(auc, n_pos, n_neg):
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc ** 2 / (1.0 + auc)
    se = math.sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
                    + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg))
    return se, auc - 1.96 * se, auc + 1.96 * se


def oof_predict(X, y, n_splits=3, seed=42):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.full(len(X), np.nan)
    importances = []
    pos = int(y.sum()); neg = int(len(y) - pos)
    for tr, te in skf.split(X, y):
        model = XGBClassifier(n_estimators=15, learning_rate=0.3, max_depth=6,
                              scale_pos_weight=neg / pos, eval_metric="auc",
                              random_state=seed, n_jobs=1)
        model.fit(X.iloc[tr], y.iloc[tr])
        oof[te] = model.predict_proba(X.iloc[te])[:, 1]
        importances.append(model.feature_importances_)
    return oof, np.mean(importances, axis=0), neg / pos


def enrichment_factor(scores, y, percent):
    k = max(int(len(scores) * percent), 1)
    order = np.argsort(-scores)
    top = y.iloc[order[:k]]
    n_act = max(int(y.sum()), 1)
    return float(top.sum() / n_act) / percent


def metrics_row(name, scores, y):
    auc = roc_auc_score(y, scores)
    se, lo, hi = auc_ci(auc, int(y.sum()), int(len(y) - y.sum()))
    return {"method": name, "auroc": round(auc, 4), "auroc_se": round(se, 4),
            "auroc_ci95": f"{lo:.3f}-{hi:.3f}",
            "ef_1pct": round(enrichment_factor(scores, y, 0.01), 3),
            "ef_5pct": round(enrichment_factor(scores, y, 0.05), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib-scores", default=str(REPO_ROOT / "results/tables/docking_scores.csv"))
    ap.add_argument("--decoy-scores", default=str(REPO_ROOT / "results/tables/docking_scores_decoy_set.csv"))
    ap.add_argument("--in-library-actives", default=str(REPO_ROOT / "data/training/beta2ar_actives_in_library.csv"))
    ap.add_argument("--out", default=str(REPO_ROOT / "results/supervised_enopt_decoy"))
    ap.add_argument("--folds", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    lib = load_wide(args.lib_scores)
    decoy_set = load_wide(args.decoy_scores)
    cols = list(lib.columns)
    if list(decoy_set.columns) != cols:
        raise SystemExit(f"conformation columns differ: {list(decoy_set.columns)} vs {cols}")
    lib, decoy_set = lib[cols], decoy_set[cols]

    Xdec = decoy_set[decoy_set.index.str.startswith("decoy_")]
    Xnew = decoy_set[~decoy_set.index.str.startswith("decoy_")]
    lib_48 = set(pd.read_csv(args.in_library_actives)["molecule_chembl_id"])
    act48 = set(lib.index) & lib_48
    act_all = act48 | set(Xnew.index)
    print(f"[data] library {len(lib)} | decoys {len(Xdec)} | new actives {len(Xnew)} "
          f"(lib overlap {len(set(Xnew.index) & set(lib.index))}) | actives total {len(act_all)}")

    # 1) clean-label OOF: all actives vs DUD-E decoys
    Xc = pd.concat([lib[lib.index.isin(act_all)], Xnew, Xdec])
    yc = pd.Series((Xc.index.isin(act_all)).astype(int), index=Xc.index)
    oof_c, imp_c, w_c = oof_predict(Xc, yc, n_splits=args.folds, seed=args.seed)
    metrics_c = pd.DataFrame([
        metrics_row("EnOpt_clean_OOF", pd.Series(oof_c, index=Xc.index), yc),
        metrics_row("ensemble_neg_mean", -Xc.mean(axis=1), yc),
        metrics_row("ensemble_neg_best", -Xc.min(axis=1), yc)])
    n_c = len(Xc); order_c = np.argsort(-oof_c)
    top_c = {"top_1pct_actives": int(yc.iloc[order_c[:int(n_c * 0.01)]].sum()),
             "top_5pct_actives": int(yc.iloc[order_c[:int(n_c * 0.05)]].sum())}

    # 2) merged mirror OOF (old accounting, library + decoys + new actives)
    Xm = pd.concat([lib, Xnew, Xdec])
    ym = pd.Series((Xm.index.isin(act_all)).astype(int), index=Xm.index)
    oof_m, _, w_m = oof_predict(Xm, ym, n_splits=args.folds, seed=args.seed)
    metrics_m = pd.DataFrame([
        metrics_row("EnOpt_merged_OOF", pd.Series(oof_m, index=Xm.index), ym),
        metrics_row("ensemble_neg_mean", -Xm.mean(axis=1), ym),
        metrics_row("ensemble_neg_best", -Xm.min(axis=1), ym)])

    # 3) final model on clean labels -> re-rank library
    pos = int(yc.sum()); neg = int(len(yc) - pos)
    final = XGBClassifier(n_estimators=15, learning_rate=0.3, max_depth=6,
                          scale_pos_weight=neg / pos, eval_metric="auc",
                          random_state=args.seed, n_jobs=1)
    final.fit(Xc, yc)
    score_lib = final.predict_proba(lib)[:, 1]
    ranking = pd.DataFrame({"ligand_id": lib.index,
                            "is_known_active": lib.index.isin(act_all).astype(int),
                            "enopt_decoy_score": score_lib})
    ranking = ranking.sort_values("enopt_decoy_score", ascending=False).reset_index(drop=True)
    ranking["enopt_decoy_rank"] = np.arange(1, len(ranking) + 1)

    imp_df = pd.DataFrame({"conformation_id": cols,
                           "feature_importance": np.round(imp_c, 4)})
    imp_df = imp_df.sort_values("feature_importance", ascending=False)

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    metrics_c.to_csv(out / "metrics_clean_label.csv", index=False)
    metrics_m.to_csv(out / "metrics_merged_mirror.csv", index=False)
    imp_df.to_csv(out / "feature_importance_clean.csv", index=False)
    ranking.to_csv(out / "library_ranking_decoy_model.csv", index=False)
    joblib.dump(final, out / "enopt_model_decoy_validated.pkl")

    card = {
        "模型名称": "Beta2AR 监督式 EnOpt（DUD-E decoy 验证版, XGBoost）",
        "训练时间": pd.Timestamp.now().isoformat(timespec="minutes"),
        "相比 06 的改动": ["负样本从'假设无效的库分子'换成 2978 个真实 DUD-E decoy",
                          "正例从 48 扩到 206（新增 158 个文献活性分子并按 decoy 分子量窗口筛选后对接）"],
        "训练数据": {"活性(正例)": int(pos), "DUD-E decoy(负例)": int(neg),
                    "正负平衡权重 scale_pos_weight": round(w_c, 3)},
        "输入特征(顺序)": list(cols),
        "超参数": {"n_estimators": 15, "learning_rate": 0.3, "max_depth": 6,
                   "eval_metric": "auc", "random_state": args.seed},
        "clean_label_OOF_metrics": metrics_c.to_dict("records"),
        "clean_top_enrichment": top_c,
        "merged_mirror_OOF_metrics": metrics_m.to_dict("records"),
        "旧 06 成绩单": {"auroc": 0.6647, "ef_1pct": 4.167, "ef_5pct": 2.5},
        "注意": "clean-label 硬测试下 EnOpt 与'简单平均分(正确方向)'的 AUROC 区间重叠;"
                "EnOpt 的优势在 top1% 富集与全库口径下相对平均分的领先,不是绝对的 AUROC 胜利",
    }
    (out / "model_card_decoy_validated.json").write_text(json.dumps(card, ensure_ascii=False, indent=2))

    print("\n== clean-label OOF ({} actives vs {} decoys) ==".format(pos, neg))
    print(metrics_c.to_string(index=False))
    print("actives in top1%/top5%:", top_c)
    print("\n== merged mirror OOF ==")
    print(metrics_m.to_string(index=False))
    print("\n== library re-rank ==")
    print("known actives in top 50/200/500/1000:", {n: int(ranking.head(n)["is_known_active"].sum()) for n in (50, 200, 500, 1000)})
    print("feature importance:\n", imp_df.to_string(index=False))
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
