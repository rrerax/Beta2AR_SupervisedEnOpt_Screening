#!/usr/bin/env python3
"""Supervised EnOpt training for the Beta2AR ensemble screen.

Trains a gradient-boosted-tree model (XGBoost) on the 5-conformation docking
scores of the 30,000-compound library, using known beta2-AR actives as
positive labels (and all other library molecules as negatives, per the
EnOpt paper default assumption). Predictions are made out-of-fold so no
molecule is ever scored by a model that saw it during training.

Usage (from the repo root):
    python scripts/06_train_supervised_enopt.py

Outputs (written under results/tables/):
    enopt_supervised_metrics.csv      AUROC + enrichment factors vs baselines
    enopt_supervised_ranking.csv      per-molecule EnOpt score + new ranking
    enopt_supervised_feature_importance.csv  per-conformation importance
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]  # repo root when script sits in scripts/


def load_feature_matrix(docking_scores_csv):
    """One row per ligand, one column per receptor conformation."""
    dock = pd.read_csv(docking_scores_csv)
    pivot = dock.pivot_table(
        index="ligand_id", columns="conformation_id",
        values="best_score_kcal_mol", aggfunc="first")
    pivot = pivot.dropna()
    return pivot


def load_labels(active_csv):
    actives = pd.read_csv(active_csv)
    return set(actives["molecule_chembl_id"])


def out_of_fold_predict(X, y, n_splits=3, seed=42):
    """3-fold CV; every molecule is scored only by a model trained without it."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.full(len(X), np.nan)
    importances = []
    pos = int(y.sum()); neg = int(len(y) - y.sum())
    for train_idx, test_idx in skf.split(X, y):
        model = XGBClassifier(
            n_estimators=15,          # low tree count: avoid overfitting with few positives
            learning_rate=0.3,
            max_depth=6,
            scale_pos_weight=neg / pos,  # correct extreme class imbalance (1 : ~624)
            eval_metric="auc",
            random_state=seed,
        )
        model.fit(X.iloc[train_idx], y.iloc[train_idx])
        oof[test_idx] = model.predict_proba(X.iloc[test_idx])[:, 1]
        importances.append(model.feature_importances_)
    return oof, np.mean(importances, axis=0)


def enrichment_factor(scores, y, percent):
    """EF(x%) = (fraction of actives in top x%) / x%."""
    k = max(int(len(scores) * percent), 1)
    order = np.argsort(-scores)
    top = y.iloc[order[:k]]
    n_act = max(int(y.sum()), 1)
    return float(top.sum() / n_act) / percent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=str(REPO_ROOT / "results/tables/docking_scores.csv"))
    ap.add_argument("--actives", default=str(REPO_ROOT / "data/training/beta2ar_actives_in_library.csv"))
    ap.add_argument("--out", default=str(REPO_ROOT / "results/tables"))
    ap.add_argument("--folds", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    X = load_feature_matrix(args.scores)
    active_ids = load_labels(args.actives)
    y = pd.Series(X.index.isin(active_ids).astype(int), index=X.index)

    n_pos = int(y.sum())
    print(f"[data] {len(X):,} ligands x {X.shape[1]} conformations; {n_pos} known actives; {int(len(y)-n_pos):,} assumed negatives")

    oof, imp = out_of_fold_predict(X, y, n_splits=args.folds, seed=args.seed)
    result = pd.DataFrame({
        "ligand_id": X.index,
        "is_known_active": y.values,
        "enopt_score": oof,
    })
    result = result.sort_values("enopt_score", ascending=False).reset_index(drop=True)
    result["enopt_rank"] = np.arange(1, len(result) + 1)

    # Baselines: simple ensemble average and ensemble best
    mean_scores = X.mean(axis=1)
    best_scores = X.min(axis=1)
    rows = [];
    for name, s in [("EnOpt (supervised OOF)", pd.Series(oof, index=X.index)),
                    ("ensemble_mean", mean_scores),
                    ("ensemble_best", best_scores)]:
        auroc = roc_auc_score(y, s)
        ef1 = enrichment_factor(s, y, 0.01)
        ef5 = enrichment_factor(s, y, 0.05)
        rows.append({"method": name, "auroc": round(auroc, 4), "ef_1pct": round(ef1, 3), "ef_5pct": round(ef5, 3)})
    metrics = pd.DataFrame(rows)

    imp_df = pd.DataFrame({"conformation_id": X.columns, "feature_importance": np.round(imp, 4)})
    imp_df = imp_df.sort_values("feature_importance", ascending=False)

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out / "enopt_supervised_metrics.csv", index=False)
    result.to_csv(out / "enopt_supervised_ranking.csv", index=False)
    imp_df.to_csv(out / "enopt_supervised_feature_importance.csv", index=False)

    print(); print(metrics.to_string(index=False)); print();
    print("Top 10 by EnOpt score:"); print(result.head(10).to_string(index=False));
    print(); print("Conformation importance:"); print(imp_df.to_string(index=False));
    print(f"\nSaved metrics/ranking/importance to {out}")


if __name__ == "__main__":
    from pathlib import Path
    main()