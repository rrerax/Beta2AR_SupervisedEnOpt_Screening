#!/usr/bin/env python3
"""Does each of the five conformations contribute, or does one match the mean?

Uses the stored 29,865-ligand score matrix and the stored decoy-set scores.
No model is fit. The learned-ranker rows are the 16 October 2026 tables.
The logP window membership is the frozen file logp_match_kept.csv.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def roc_auc_score(labels: np.ndarray, scores: np.ndarray) -> float:
    """Mann-Whitney AUROC with average ranks for ties. No scikit-learn."""
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    order = np.argsort(scores, kind="mergesort")
    ordered = scores[order]
    ranks = np.empty(len(scores), dtype=float)
    start = 0
    while start < len(scores):
        stop = start
        while stop + 1 < len(scores) and ordered[stop + 1] == ordered[start]:
            stop += 1
        ranks[order[start : stop + 1]] = (start + 1 + stop + 1) / 2.0
        start = stop + 1
    n_pos = int(labels.sum())
    n_neg = int(len(labels) - n_pos)
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
FOLLOW = REPO / "followup_2026-10-16"
SCORE_COLS = [
    "inactive_carazolol_2RH1",
    "inactive_carazolol_5D5A",
    "active_gs_complex_3SN6",
    "active_bi167107_4LDE",
    "active_hbi_4LDL",
]
SHORT = {
    "inactive_carazolol_2RH1": "2RH1",
    "inactive_carazolol_5D5A": "5D5A",
    "active_gs_complex_3SN6": "3SN6",
    "active_bi167107_4LDE": "4LDE",
    "active_hbi_4LDL": "4LDL",
}


def load_wide(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, usecols=["ligand_id", "conformation_id", "best_score_kcal_mol"])
    frame["ligand_id"] = frame["ligand_id"].astype(str)
    wide = frame.pivot_table(
        index="ligand_id",
        columns="conformation_id",
        values="best_score_kcal_mol",
        aggfunc="first",
    )
    missing = [column for column in SCORE_COLS if column not in wide.columns]
    if missing:
        raise SystemExit(f"{path} is missing {missing}")
    return wide[SCORE_COLS]


def enrichment(scores: np.ndarray, labels: np.ndarray) -> tuple[float, int]:
    count = max(int(len(scores) * 0.01), 1)
    order = np.argsort(-np.asarray(scores, dtype=float))
    positives = max(int(np.asarray(labels).sum()), 1)
    value = float(np.asarray(labels)[order[:count]].sum() / positives) / 0.01
    return value, count


def measure(pool: str, name: str, scores: np.ndarray, labels: np.ndarray, source: str) -> dict:
    labels = np.asarray(labels).astype(int)
    auc = float(roc_auc_score(labels, scores))
    n_pos = int(labels.sum())
    n_neg = int(len(labels) - n_pos)
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc ** 2 / (1.0 + auc)
    se = math.sqrt(
        (auc * (1.0 - auc) + (n_pos - 1) * (q1 - auc ** 2) + (n_neg - 1) * (q2 - auc ** 2))
        / (n_pos * n_neg)
    )
    ef1, top1 = enrichment(scores, labels)
    return {
        "pool": pool,
        "method": name,
        "source": source,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "top1_count": top1,
        "auroc": round(auc, 4),
        "auroc_se": round(se, 4),
        "ef_1pct": round(ef1, 3),
    }


def ranking_rows(pool: str, scores: pd.DataFrame, labels: pd.Series) -> list[dict]:
    labels = labels.loc[scores.index]
    rows = []
    for column in SCORE_COLS:
        rows.append(
            measure(
                pool,
                f"single_{SHORT[column]}",
                -scores[column].to_numpy(dtype=float),
                labels.to_numpy(),
                "recomputed from stored Vina scores",
            )
        )
    rows.append(
        measure(
            pool,
            "oriented_five_mean",
            -scores[SCORE_COLS].mean(axis=1).to_numpy(dtype=float),
            labels.to_numpy(),
            "recomputed from stored Vina scores",
        )
    )
    return rows


def stored_ranker(pool: str, path: Path) -> dict:
    frame = pd.read_csv(path)
    learned = frame[frame["method"].str.startswith("learned_raw_seed")]
    return {
        "pool": pool,
        "method": "learned_ranker_ten_seeds",
        "source": str(path.relative_to(REPO)).replace("\\", "/"),
        "n_pos": None,
        "n_neg": None,
        "top1_count": int(frame["top_1pct_k"].iloc[0]) if "top_1pct_k" in frame.columns else None,
        "auroc": None,
        "auroc_min": round(float(learned["auroc"].min()), 4),
        "auroc_max": round(float(learned["auroc"].max()), 4),
        "ef_1pct": None,
        "ef_1pct_min": round(float(learned["ef_1pct"].min()), 3),
        "ef_1pct_max": round(float(learned["ef_1pct"].max()), 3),
        "seeds": 10,
        "refit": False,
    }


def main() -> None:
    library = load_wide(REPO / "results/tables/docking_scores.csv").dropna()
    decoy_scores = load_wide(REPO / "results/tables/docking_scores_decoy_set.csv").dropna()
    in_library = set(
        pd.read_csv(REPO / "data/training/beta2ar_actives_in_library.csv")["molecule_chembl_id"].astype(str)
    )
    expanded = set(
        pd.read_csv(REPO / "data/training/beta2ar_actives_expand_docked.csv")["ligand_id"].astype(str)
    )
    scores = pd.concat([decoy_scores, library])
    scores = scores[~scores.index.duplicated(keep="first")]
    actives = (in_library & set(library.index)) | (expanded & set(scores.index))
    decoys = {ligand for ligand in decoy_scores.index if ligand.startswith("decoy_")}
    labeled = scores.loc[sorted(actives | decoys)]
    labels = pd.Series(labeled.index.isin(actives).astype(int), index=labeled.index)
    if int(labels.sum()) != 206 or int((1 - labels).sum()) != 2978:
        raise SystemExit(f"expected 206/2978, got {int(labels.sum())}/{(1 - labels).sum()}")

    kept = pd.read_csv(FOLLOW / "logp_match_2026-10-16/logp_match_kept.csv")
    kept["ligand_id"] = kept["ligand_id"].astype(str)
    missing = sorted(set(kept["ligand_id"]) - set(labeled.index))
    if missing:
        raise SystemExit(f"frozen logP window has {len(missing)} ids without scores")
    window_scores = labeled.loc[kept["ligand_id"]]
    window_labels = kept.set_index("ligand_id").loc[window_scores.index, "y"].astype(int)
    if int(window_labels.sum()) != 48 or int((1 - window_labels).sum()) != 379:
        raise SystemExit(
            f"expected 48/379, got {int(window_labels.sum())}/{int((1 - window_labels).sum())}"
        )

    rows = ranking_rows("unmatched_206_2978", labeled, labels)
    rows.extend(ranking_rows("logp_window_48_379", window_scores, window_labels))
    table = pd.DataFrame(rows)
    rankers = [
        stored_ranker("unmatched_206_2978", FOLLOW / "week1_2026-10-16/week1_metrics.csv"),
        stored_ranker("logp_window_48_379", FOLLOW / "logp_match_2026-10-16/logp_match_metrics.csv"),
    ]
    table.to_csv(OUT / "conformation_contribution.csv", index=False)
    pd.DataFrame(rankers).to_csv(OUT / "learned_ranker_stored.csv", index=False)

    def lead(pool: str) -> dict:
        part = table[table["pool"] == pool]
        mean = part[part["method"] == "oriented_five_mean"].iloc[0]
        singles = part[part["method"].str.startswith("single_")]
        both = singles[
            (singles["auroc"] >= mean["auroc"]) & (singles["ef_1pct"] >= mean["ef_1pct"])
        ]
        return {
            "pool": pool,
            "mean_auroc": float(mean["auroc"]),
            "mean_ef_1pct": float(mean["ef_1pct"]),
            "best_single_auroc": singles.loc[singles["auroc"].idxmax(), "method"],
            "best_single_auroc_value": float(singles["auroc"].max()),
            "best_single_ef": singles.loc[singles["ef_1pct"].idxmax(), "method"],
            "best_single_ef_value": float(singles["ef_1pct"].max()),
            "singles_at_least_the_mean_on_both": both["method"].tolist(),
        }

    decision = {
        "question": "Do all five conformations contribute, or does one alone match the mean?",
        "rule": (
            "A single conformation matches the mean when its AUROC and its top-1% "
            "enrichment are both at least the oriented five-conformation mean. "
            "More negative Vina score ranks first. No model was fit."
        ),
        "pools": [lead("unmatched_206_2978"), lead("logp_window_48_379")],
        "learned_ranker_not_refit": rankers,
    }
    (OUT / "conformation_contribution.json").write_text(
        json.dumps(decision, indent=2), encoding="utf-8"
    )
    print(table.to_string(index=False))
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
