#!/usr/bin/env python3
"""Property-matched, scaffold-held-out exam of the existing five Vina scores.

The constants below were fixed before any metric in this folder was read.
No model is written. No library is re-ranked. No docking is started.

Question for the gate. Among docked library molecules that share a
Bemis-Murcko scaffold with one of the 206 actives, and that sit inside the
property calipers of those actives, does an out-of-fold five-score model beat
both baselines under scaffold holdout? The baselines are molecular weight,
Crippen logP, and heavy-atom count, and the oriented 4LDE score alone.

The same calipers are also applied to the docked DUD-E decoys. That pool is a
different-chemotype control. It is archived and cannot flip the gate.

A pool can pass only when it has at least 30 covered actives, 100 negatives,
10 active scaffolds, and at least 10 molecules in the top 1%. AUROC may be
stored when the first three counts are met and the top 1% is smaller than 10.
Those stored numbers are not the gate.

Passing a baseline means the five-score AUROC exceeds that baseline by more
than the five-score AUROC standard error, and the five-score top-1%
enrichment is higher. The gate passes only when the same-chemotype pool is
evaluable and passes both baselines. Ties do not pass.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBClassifier

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

CALIPER = {"logp": 0.5, "mw": 50.0, "hbd": 1.0, "hba": 1.0, "charge": 0.0}
MATCH_FEATURES = ["logp", "mw", "hbd", "hba", "charge"]
DESCRIPTOR_FEATURES = ["mw", "logp", "hac"]
SCORE_COLS = [
    "inactive_carazolol_2RH1",
    "inactive_carazolol_5D5A",
    "active_gs_complex_3SN6",
    "active_bi167107_4LDE",
    "active_hbi_4LDL",
]
CONF_4LDE = "active_bi167107_4LDE"
K_NEGATIVES = 5
MIN_ACTIVES = 30
MIN_NEGATIVES = 100
MIN_ACTIVE_SCAFFOLDS = 10
MIN_TOP1 = 10
N_SPLITS = 5
SEED = 42
N_ESTIMATORS = 15
LEARNING_RATE = 0.3
MAX_DEPTH = 6


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
        raise SystemExit(f"{path} is missing conformations: {missing}")
    return wide[SCORE_COLS]


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())


def describe(smiles: str) -> dict | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = largest_fragment(mol)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if not scaffold:
        scaffold = "ACYCLIC:" + Chem.MolToSmiles(mol, canonical=True)
    return {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "hac": float(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
        "hbd": float(Lipinski.NumHDonors(mol)),
        "hba": float(Lipinski.NumHAcceptors(mol)),
        "charge": float(Chem.GetFormalCharge(mol)),
        "scaffold": scaffold,
    }


def auc_se(auc: float, n_pos: int, n_neg: int) -> float:
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc ** 2 / (1.0 + auc)
    return math.sqrt(
        (auc * (1.0 - auc) + (n_pos - 1) * (q1 - auc ** 2) + (n_neg - 1) * (q2 - auc ** 2))
        / (n_pos * n_neg)
    )


def enrichment(scores: np.ndarray, labels: np.ndarray, percent: float) -> tuple[float, int]:
    count = max(int(len(scores) * percent), 1)
    order = np.argsort(-np.asarray(scores, dtype=float))
    positives = max(int(np.asarray(labels).sum()), 1)
    value = float(np.asarray(labels)[order[:count]].sum() / positives) / percent
    return value, int(len(scores) * percent)


def measure(name: str, scores: np.ndarray, labels: np.ndarray) -> dict:
    labels = np.asarray(labels).astype(int)
    auc = float(roc_auc_score(labels, scores))
    n_pos = int(labels.sum())
    n_neg = int(len(labels) - n_pos)
    se = auc_se(auc, n_pos, n_neg)
    ef1, top1 = enrichment(scores, labels, 0.01)
    ef5, _ = enrichment(scores, labels, 0.05)
    return {
        "method": name,
        "n_molecules": int(len(labels)),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "top1_count": top1,
        "auroc": round(auc, 6),
        "auroc_se": round(se, 6),
        "auroc_ci95_lo": round(auc - 1.96 * se, 6),
        "auroc_ci95_hi": round(auc + 1.96 * se, 6),
        "ef_1pct": round(ef1, 6),
        "ef_5pct": round(ef5, 6),
    }


def beats(challenger: dict, baseline: dict) -> bool:
    return (
        challenger["auroc"] > baseline["auroc"] + challenger["auroc_se"]
        and challenger["ef_1pct"] > baseline["ef_1pct"]
    )


def oof_predict(features: pd.DataFrame, labels: np.ndarray, groups: np.ndarray) -> np.ndarray:
    splitter = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    scores = np.full(len(features), np.nan)
    frame = features.reset_index(drop=True)
    for train_idx, test_idx in splitter.split(frame, labels, groups):
        train_labels = labels[train_idx]
        n_pos = int(train_labels.sum())
        n_neg = int(len(train_labels) - n_pos)
        if n_pos == 0 or n_neg == 0:
            raise RuntimeError("a training fold has only one class")
        model = XGBClassifier(
            n_estimators=N_ESTIMATORS,
            learning_rate=LEARNING_RATE,
            max_depth=MAX_DEPTH,
            scale_pos_weight=n_neg / n_pos,
            random_state=SEED,
            n_jobs=1,
            verbosity=0,
        )
        model.fit(frame.iloc[train_idx], train_labels)
        scores[test_idx] = model.predict_proba(frame.iloc[test_idx])[:, 1]
    if np.isnan(scores).any():
        raise RuntimeError("out-of-fold scores are incomplete")
    return scores


def match_pairs(actives: pd.DataFrame, negatives: pd.DataFrame) -> pd.DataFrame:
    """Greedy nearest match inside the calipers. Uses properties only."""
    if actives.empty or negatives.empty:
        return pd.DataFrame(columns=["active_id", "negative_id", "distance", "active_scaffold"])
    means = actives[MATCH_FEATURES].mean()
    stds = actives[MATCH_FEATURES].std(ddof=0).replace(0.0, np.nan)
    active_z = ((actives[MATCH_FEATURES] - means) / stds).fillna(0.0).to_numpy(dtype=float)
    negative_z = ((negatives[MATCH_FEATURES] - means) / stds).fillna(0.0).to_numpy(dtype=float)
    active_raw = actives[MATCH_FEATURES].to_numpy(dtype=float)
    negative_raw = negatives[MATCH_FEATURES].to_numpy(dtype=float)
    delta = np.abs(active_raw[:, None, :] - negative_raw[None, :, :])
    allowed = (
        (delta[:, :, 0] <= CALIPER["logp"])
        & (delta[:, :, 1] <= CALIPER["mw"])
        & (delta[:, :, 2] <= CALIPER["hbd"])
        & (delta[:, :, 3] <= CALIPER["hba"])
        & (delta[:, :, 4] <= CALIPER["charge"])
    )
    distance = np.sqrt(np.sum((active_z[:, None, :] - negative_z[None, :, :]) ** 2, axis=2))
    active_idx, negative_idx = np.where(allowed)
    if len(active_idx) == 0:
        return pd.DataFrame(columns=["active_id", "negative_id", "distance", "active_scaffold"])
    pair_distance = distance[active_idx, negative_idx]
    order = np.lexsort((negative_idx, active_idx, pair_distance))
    active_ids = actives.index.to_numpy()
    negative_ids = negatives.index.to_numpy()
    used: set[str] = set()
    taken: dict[str, int] = {}
    rows = []
    for position in order:
        aid = str(active_ids[active_idx[position]])
        nid = str(negative_ids[negative_idx[position]])
        if taken.get(aid, 0) >= K_NEGATIVES or nid in used:
            continue
        used.add(nid)
        taken[aid] = taken.get(aid, 0) + 1
        rows.append(
            {
                "active_id": aid,
                "negative_id": nid,
                "distance": round(float(pair_distance[position]), 6),
                "active_scaffold": actives.at[aid, "scaffold"],
            }
        )
    return pd.DataFrame(rows)


def match_same_scaffold(actives: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for scaffold, active_group in actives.groupby("scaffold", sort=True):
        negative_group = library[library["scaffold"] == scaffold]
        if negative_group.empty:
            continue
        frames.append(match_pairs(active_group, negative_group))
    if not frames:
        return pd.DataFrame(columns=["active_id", "negative_id", "distance", "active_scaffold"])
    return pd.concat(frames, ignore_index=True)


def median_block(frame: pd.DataFrame, label: str) -> dict:
    columns = ["logp", "mw", "hac", "hbd", "hba", "charge"]
    out = {"set": label, "n": int(len(frame))}
    for column in columns:
        out[f"median_{column}"] = None if frame.empty else round(float(frame[column].median()), 4)
    return out


def evaluate_pool(name: str, pairs: pd.DataFrame, properties: pd.DataFrame, scores: pd.DataFrame) -> dict:
    covered = sorted(set(pairs["active_id"])) if not pairs.empty else []
    negatives = sorted(set(pairs["negative_id"])) if not pairs.empty else []
    active_scaffolds = properties.loc[covered, "scaffold"].nunique() if covered else 0
    molecules = covered + negatives
    top1_possible = int(len(molecules) * 0.01)
    counts_ok = (
        len(covered) >= MIN_ACTIVES
        and len(negatives) >= MIN_NEGATIVES
        and int(active_scaffolds) >= MIN_ACTIVE_SCAFFOLDS
    )
    enrichment_estimable = top1_possible >= MIN_TOP1
    evaluable = counts_ok and enrichment_estimable
    if not counts_ok:
        reason = "matched_set_below_minimum"
    elif not enrichment_estimable:
        reason = "top1_percent_below_10_molecules"
    else:
        reason = "metrics_not_computed"
    result = {
        "pool": name,
        "n_pairs": int(len(pairs)),
        "n_covered_actives": int(len(covered)),
        "n_negatives": int(len(negatives)),
        "n_molecules": int(len(molecules)),
        "n_active_scaffolds": int(active_scaffolds),
        "top1_count": int(top1_possible),
        "evaluable": bool(evaluable),
        "pass_pool": False,
        "reason": reason,
        "metrics_used_for_gate": bool(evaluable),
        "metrics": [],
    }
    if not counts_ok:
        return result
    chosen = properties.loc[molecules]
    labels = np.array([1 if ligand in set(covered) else 0 for ligand in molecules], dtype=int)
    groups = chosen["scaffold"].to_numpy()
    try:
        five_scores = oof_predict(scores.loc[molecules, SCORE_COLS], labels, groups)
        descriptor_scores = oof_predict(chosen[DESCRIPTOR_FEATURES], labels, groups)
    except (RuntimeError, ValueError) as exc:
        result["reason"] = f"scaffold_split_unavailable: {exc}"
        return result
    oriented_4lde = -scores.loc[molecules, CONF_4LDE].to_numpy(dtype=float)
    oriented_mean = -scores.loc[molecules, SCORE_COLS].mean(axis=1).to_numpy(dtype=float)
    metrics = [
        measure("five_score_oof", five_scores, labels),
        measure("descriptor_oof", descriptor_scores, labels),
        measure("oriented_4lde", oriented_4lde, labels),
        measure("oriented_five_mean", oriented_mean, labels),
    ]
    by_name = {row["method"]: row for row in metrics}
    numerical_pass = beats(by_name["five_score_oof"], by_name["descriptor_oof"]) and beats(
        by_name["five_score_oof"], by_name["oriented_4lde"]
    )
    result["metrics"] = metrics
    result["numerical_comparison_passes"] = bool(numerical_pass)
    result["spearman_oriented_mean_logp"] = round(float(spearman(
        oriented_mean, chosen["logp"].to_numpy(dtype=float)
    )), 6)
    if not evaluable:
        return result
    result["pass_pool"] = bool(numerical_pass)
    result["reason"] = (
        "five_score_beats_both_baselines" if numerical_pass else "five_score_does_not_beat_both_baselines"
    )
    return result


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    ranked_left = pd.Series(left).rank().to_numpy(dtype=float)
    ranked_right = pd.Series(right).rank().to_numpy(dtype=float)
    if np.std(ranked_left) == 0.0 or np.std(ranked_right) == 0.0:
        return float("nan")
    return float(np.corrcoef(ranked_left, ranked_right)[0, 1])


def charge_counts(frame: pd.DataFrame) -> dict:
    counts = frame["charge"].value_counts().sort_index()
    return {
        str(int(charge)) if float(charge).is_integer() else str(charge): int(count)
        for charge, count in counts.items()
    }


def lookup_known(
    ligand_id: str,
    actives: pd.DataFrame,
    library_scores: pd.DataFrame,
    same_pairs: pd.DataFrame,
    dude_pairs: pd.DataFrame,
) -> dict:
    record = {"ligand_id": ligand_id, "in_206_actives": ligand_id in actives.index}
    if ligand_id in actives.index:
        row = actives.loc[ligand_id]
        record.update({
            "logp": round(float(row.logp), 4),
            "mw": round(float(row.mw), 4),
            "formal_charge": int(row.charge),
            "scaffold": row.scaffold,
            "covered_by_same_chemotype_match": bool(
                not same_pairs.empty and (same_pairs["active_id"] == ligand_id).any()
            ),
            "covered_by_dude_match": bool(
                not dude_pairs.empty and (dude_pairs["active_id"] == ligand_id).any()
            ),
        })
    if ligand_id in library_scores.index:
        oriented = -library_scores[SCORE_COLS].mean(axis=1)
        order = oriented.sort_values(ascending=False)
        record["library_oriented_mean_rank"] = int(order.index.get_loc(ligand_id)) + 1
        record["library_size"] = int(len(order))
    return record


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    library_scores = load_wide(REPO / "results/tables/docking_scores.csv").dropna()
    decoy_scores = load_wide(REPO / "results/tables/docking_scores_decoy_set.csv").dropna()
    scores = pd.concat([decoy_scores, library_scores])
    scores = scores[~scores.index.duplicated(keep="first")]

    actives_in_library = pd.read_csv(
        REPO / "data/training/beta2ar_actives_in_library.csv",
        usecols=["molecule_chembl_id", "smiles"],
    ).rename(columns={"molecule_chembl_id": "ligand_id"})
    actives_expanded = pd.read_csv(
        REPO / "data/training/beta2ar_actives_expand_docked.csv",
        usecols=["ligand_id", "smiles"],
    )
    decoy_table = pd.read_csv(
        REPO / "data/training/dude_decoys_for_dock.csv",
        usecols=["ligand_id", "smiles"],
    )
    manifest = pd.read_csv(
        REPO / "results/tables/ligand_manifest.csv",
        usecols=["ligand_id", "smiles"],
    )
    active_table = pd.concat([actives_in_library, actives_expanded], ignore_index=True)
    active_table = active_table.drop_duplicates("ligand_id", keep="first")

    smiles = {}
    for frame in (active_table, decoy_table, manifest):
        for ligand_id, value in frame[["ligand_id", "smiles"]].itertuples(index=False):
            smiles.setdefault(str(ligand_id), value)

    described = {}
    missing = []
    needed = set(active_table["ligand_id"].astype(str))
    needed.update(str(ligand) for ligand in decoy_scores.index if str(ligand).startswith("decoy_"))
    needed.update(str(ligand) for ligand in library_scores.index)
    for ligand_id in sorted(needed):
        record = describe(smiles.get(ligand_id, ""))
        if record is None:
            missing.append(ligand_id)
            continue
        described[ligand_id] = record
    properties = pd.DataFrame.from_dict(described, orient="index")
    properties.index.name = "ligand_id"

    active_ids = [ligand for ligand in active_table["ligand_id"].astype(str) if ligand in properties.index]
    actives = properties.loc[active_ids]
    actives = actives[actives.index.isin(scores.index)]
    active_canonical = set(actives["canonical_smiles"])

    decoy_ids = [
        ligand for ligand in properties.index
        if ligand.startswith("decoy_")
        and ligand in decoy_scores.index
        and properties.at[ligand, "canonical_smiles"] not in active_canonical
    ]
    decoys = properties.loc[decoy_ids]

    library_ids = [
        ligand for ligand in library_scores.index.astype(str)
        if ligand in properties.index
        and ligand not in actives.index
        and not ligand.startswith("decoy_")
        and properties.at[ligand, "canonical_smiles"] not in active_canonical
    ]
    library = properties.loc[library_ids]

    same_pairs = match_same_scaffold(actives, library)
    dude_pairs = match_pairs(actives, decoys)
    same_result = evaluate_pool("same_chemotype_library", same_pairs, properties, scores)
    dude_result = evaluate_pool("dude_property_matched", dude_pairs, properties, scores)

    gate = bool(same_result["pass_pool"])
    decision = {
        "gate_constants_fixed_before_metrics": True,
        "auroc_below_top1_minimum_is_not_the_gate": True,
        "question": (
            "On property-matched, same-Bemis-Murcko-scaffold library molecules, "
            "does the five-score model beat molecular weight, Crippen logP, and "
            "heavy-atom count, and beat 4LDE alone, under scaffold holdout?"
        ),
        "calipers": CALIPER,
        "negatives_per_active_cap": K_NEGATIVES,
        "minimums": {
            "covered_actives": MIN_ACTIVES,
            "negatives": MIN_NEGATIVES,
            "active_scaffolds": MIN_ACTIVE_SCAFFOLDS,
            "top1_count": MIN_TOP1,
        },
        "model": {
            "n_estimators": N_ESTIMATORS,
            "learning_rate": LEARNING_RATE,
            "max_depth": MAX_DEPTH,
            "n_splits": N_SPLITS,
            "seed": SEED,
            "saved": False,
        },
        "n_active_ids_in_tables": int(len(active_table)),
        "n_actives_with_descriptors_and_scores": int(len(actives)),
        "n_dude_decoys_complete": int(len(decoys)),
        "n_library_non_actives_complete": int(len(library)),
        "n_missing_descriptors": int(len(missing)),
        "same_chemotype": same_result,
        "dude_property_matched": dude_result,
        "exam_pass_gate": gate,
        "use_learned_ranker": False,
        "start_representation_docking": False,
        "reason": same_result["reason"] if not gate else "same_chemotype_pool_passed",
        "formal_charge_counts": {
            "actives": charge_counts(actives),
            "dude_decoys": charge_counts(decoys),
            "same_chemotype_negatives": charge_counts(properties.loc[same_pairs["negative_id"].unique()])
            if not same_pairs.empty else {},
            "dude_matched_decoys": charge_counts(properties.loc[dude_pairs["negative_id"].unique()])
            if not dude_pairs.empty else {},
        },
        "chembl434": lookup_known("CHEMBL434", actives, library_scores, same_pairs, dude_pairs),
        "chembl776": lookup_known("CHEMBL776", actives, library_scores, same_pairs, dude_pairs),
        "property_medians": [
            median_block(actives, "actives_206_complete"),
            median_block(decoys, "dude_decoys_complete"),
            median_block(properties.loc[same_pairs["active_id"].unique()], "same_chemotype_covered_actives")
            if not same_pairs.empty else median_block(actives.iloc[0:0], "same_chemotype_covered_actives"),
            median_block(properties.loc[same_pairs["negative_id"].unique()], "same_chemotype_negatives")
            if not same_pairs.empty else median_block(library.iloc[0:0], "same_chemotype_negatives"),
            median_block(properties.loc[dude_pairs["active_id"].unique()], "dude_matched_actives")
            if not dude_pairs.empty else median_block(actives.iloc[0:0], "dude_matched_actives"),
            median_block(properties.loc[dude_pairs["negative_id"].unique()], "dude_matched_decoys")
            if not dude_pairs.empty else median_block(decoys.iloc[0:0], "dude_matched_decoys"),
        ],
    }
    if gate:
        decision["reason"] = (
            "same_chemotype_pool_passed; representation docking was not started by this script"
        )
        decision["start_representation_docking"] = False

    same_pairs.to_csv(OUT / "same_chemotype_pairs.csv", index=False)
    dude_pairs.to_csv(OUT / "dude_matched_pairs.csv", index=False)
    pd.DataFrame(same_result["metrics"]).to_csv(OUT / "same_chemotype_metrics.csv", index=False)
    pd.DataFrame(dude_result["metrics"]).to_csv(OUT / "dude_matched_metrics.csv", index=False)
    pd.DataFrame(decision["property_medians"]).to_csv(OUT / "property_medians.csv", index=False)
    (OUT / "exam_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    print(json.dumps({
        "exam_pass_gate": decision["exam_pass_gate"],
        "reason": decision["reason"],
        "same_chemotype": {
            key: same_result[key]
            for key in (
                "n_covered_actives", "n_negatives", "n_active_scaffolds",
                "top1_count", "evaluable", "pass_pool", "reason",
            )
        },
        "dude_property_matched": {
            key: dude_result[key]
            for key in (
                "n_covered_actives", "n_negatives", "n_active_scaffolds",
                "top1_count", "evaluable", "pass_pool", "reason",
            )
        },
    }, indent=2))


if __name__ == "__main__":
    main()
