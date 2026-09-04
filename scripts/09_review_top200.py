#!/usr/bin/env python3
"""Stage 4b: chemical review of the top-200 leaderboard.

ChEMBL retrospective (08) shows the 30k chemreps library contains almost no
beta-2-documented actives outside the training labels (2 strong / 29817), so
literature lookup cannot score the top-200 statistically.  This script
therefore reviews the top-200 as a *candidate shortlist*:
  - structural diversity (Butina clusters + Bemis-Murcko scaffolds)
  - novelty vs the 206 training actives (max Tanimoto to an active)
  - PAINS / Brenk pan-assay-interference flags (aggregator risk)
  - physicochemical window vs the training actives

Requires rdkit.  Docks nothing.

Usage (from repo root):
    python scripts/09_review_top200.py
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

REPO_ROOT = Path(__file__).resolve().parents[1]


def ecfp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def max_tanimoto_to_actives(fp, active_fps):
    best = 0.0
    for afp in active_fps:
        t = DataStructs.TanimotoSimilarity(fp, afp)
        if t > best:
            best = t
    return round(float(best), 3)


def pains_flags(mol, catalogs):
    out = []
    for name, cat in catalogs.items():
        if cat.HasMatch(mol):
            out.append(name)
    return ";".join(out)


def butina_clusters(fps, cutoff=0.5):
    from rdkit.ML.Cluster import Butina
    n = len(fps)
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranking", default=str(REPO_ROOT / "results/supervised_enopt_decoy/library_ranking_decoy_model.csv"))
    ap.add_argument("--library", default=str(REPO_ROOT / "data/raw/chembl37_ligand_library_30000.csv"))
    ap.add_argument("--actives-48", default=str(REPO_ROOT / "data/training/beta2ar_actives_in_library.csv"))
    ap.add_argument("--actives-158", default=str(REPO_ROOT / "data/training/beta2ar_actives_expand_docked.csv"))
    ap.add_argument("--out", default=str(REPO_ROOT / "results/review_top200"))
    ap.add_argument("--top-n", type=int, default=200)
    args = ap.parse_args()

    rank = pd.read_csv(args.ranking).sort_values("enopt_decoy_rank").head(args.top_n)
    lib = pd.read_csv(args.library)[["ligand_id", "smiles", "mol_weight",
                                     "heavy_atoms", "logp", "tpsa"]]
    df = rank.merge(lib, on="ligand_id", how="left")

    act = pd.concat([
        pd.read_csv(args.actives_48)[["molecule_chembl_id", "smiles"]]
        .rename(columns={"molecule_chembl_id": "ligand_id"}),
        pd.read_csv(args.actives_158),
    ])

    params = FilterCatalogParams()
    for c in (FilterCatalogParams.FilterCatalogs.PAINS,
              FilterCatalogParams.FilterCatalogs.BRENK):
        params.AddCatalog(c)
    catalogs = {"PAINS": FilterCatalog.FilterCatalog(params)}
    params2 = FilterCatalogParams()
    params2.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalogs["Brenk"] = FilterCatalog.FilterCatalog(params2)

    active_fps = []
    for smi in act["smiles"]:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            active_fps.append(ecfp(m))
    print(f"[data] top {len(df)} | actives {len(active_fps)}")

    fps, clusters = [], None
    records = []
    for _, row in df.iterrows():
        m = Chem.MolFromSmiles(row["smiles"])
        if m is None:
            records.append({"ligand_id": row["ligand_id"],
                            "parse_fail": True})
            continue
        fp = ecfp(m)
        fps.append(fp)
        rec = {
            "ligand_id": row["ligand_id"],
            "rank": int(row["enopt_decoy_rank"]),
            "is_known_active": int(row["is_known_active"]),
            "smiles": row["smiles"],
            "mol_weight": row["mol_weight"],
            "heavy_atoms": row["heavy_atoms"],
            "logp": row["logp"],
            "tpsa": row["tpsa"],
            "max_tanimoto_active": max_tanimoto_to_actives(fp, active_fps),
            "flags": pains_flags(m, catalogs),
            "parse_fail": False,
        }
        records.append(rec)
    out_df = pd.DataFrame(records)
    ok = out_df[~out_df["parse_fail"]].copy()

    clusters = butina_clusters(fps, cutoff=0.5)
    cid = np.zeros(len(ok), dtype=int)
    for i, cl in enumerate(clusters):
        for idx in cl:
            cid[idx] = i
    ok["cluster_id"] = cid
    scaffolds = []
    for smi in ok["smiles"]:
        m = Chem.MolFromSmiles(smi)
        sc = MurckoScaffold.MakeScaffoldGeneric(
            MurckoScaffold.GetScaffoldForMol(m))
        scaffolds.append(Chem.MolToSmiles(sc))
    ok["scaffold"] = scaffolds

    n_clusters = ok["cluster_id"].nunique()
    top_clusters = (ok.groupby("cluster_id").size()
                    .sort_values(ascending=False).head(5))
    n_scaffolds = ok["scaffold"].nunique()
    flagged = ok[ok["flags"] != ""]
    n_close = int((ok["max_tanimoto_active"] >= 0.6).sum())
    n_novel = int((ok["max_tanimoto_active"] < 0.4).sum())

    print("\n== top-200 chemical review ==")
    print(f"Butina clusters (Tanimoto 0.5): {n_clusters} for {len(ok)} "
          f"compounds | unique scaffolds: {n_scaffolds}")
    print("largest clusters:", top_clusters.to_dict())
    print(f"close analogs of training actives (maxT >= 0.6): {n_close}")
    print(f"novel scaffolds/signals (maxT < 0.4): {n_novel}")
    print(f"PAINS/Brenk flagged: {len(flagged)}")
    if len(flagged):
        print(flagged[["rank", "ligand_id", "flags"]].head(20).to_string(index=False))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ok.sort_values("rank").to_csv(out / "top200_chemical_review.csv", index=False)

    # property comparison vs actives
    act_mols = []
    for smi in act["smiles"]:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            act_mols.append(m)
    props = {"mol_weight": lambda m: Descriptors.MolWt(m),
             "logp": lambda m: Crippen.MolLogP(m),
             "tpsa": lambda m: rdMolDescriptors.CalcTPSA(m)}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, (name, fn) in zip(axes, props.items()):
        top_vals = [fn(Chem.MolFromSmiles(s)) for s in ok["smiles"]]
        act_vals = [fn(m) for m in act_mols]
        ax.boxplot([top_vals, act_vals])
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["top-200", "206 actives"])
        ax.set_title(name)
        ax.tick_params(axis="x", labelrotation=15, labelsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_top200_properties_vs_actives.png", dpi=150)
    plt.close(fig)

    md = f"""# Top-{args.top_n} chemical review

Ready-for-review summary of the Stage-3 leaderboard's top {args.top_n}
(molecules that the model was never given as actives unless flagged).

| metric | value |
|---|--:|
| compounds analysed | {len(ok)} |
| Butina clusters (ECFP4, Tanimoto 0.5) | {n_clusters} |
| unique Bemis-Murcko scaffolds | {n_scaffolds} |
| close analogs of a training active (maxT >= 0.6) | {n_close} |
| structurally novel vs actives (maxT < 0.4) | {n_novel} |
| PAINS/Brenk flagged (aggregator risk) | {len(flagged)} |

Largest clusters: {top_clusters.to_dict()}

Full per-compound data: `top200_chemical_review.csv`. Property figure:
`fig_top200_properties_vs_actives.png`.

Reading: if top-200 collapses into very few clusters, "200 hits" really means
a few chemotypes -- diversity (cluster count) should be inspected before
picking purchase/assay candidates.
"""
    (out / "summary_top200_review.md").write_text(md, encoding="utf-8")
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
