#!/usr/bin/env python3
"""Stage 4a: retrospective validation of the top-200 leaderboard against ChEMBL.

Question: the Stage-3 decoy-validated model re-ranks the 29,865-compound
library, but the library molecules were never labelled by literature activity
(only the 48 curated in-library actives carry a label).  Of the top 200 ranked
compounds that the model was never told were active, how many are *already
documented* as human beta-2 adrenergic receptor ligands in ChEMBL -- and is
that proportion higher than the base rate in the rest of the library?

This is a free, retrospective proxy for experimental validation:
  1) top-200 "no-label" compounds  -> query ChEMBL activities (CHEMBL210)
  2) random 500-compound control of the no-label remainder -> base rate
  3) enrichment ratio + Fisher exact test

Requires network access to the public ChEMBL REST API.  Docks nothing.

Usage (from repo root):
    python scripts/08_retrospective_top200.py
"""
import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu

REPO_ROOT = Path(__file__).resolve().parents[1]
API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
STANDARD_TYPES = "IC50,Ki,EC50,KD,Kd"
STRONG_PCHEMBL = 6.0


def query_chunk(ids, target="CHEMBL210", retries=4):
    params = {
        "target_chembl_id": target,
        "molecule_chembl_id__in": ",".join(ids),
        "standard_type__in": STANDARD_TYPES,
        "limit": "1000",
    }
    url = API + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "beta2ar-retro-validation/0.1"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = json.load(resp)
            return data.get("activities", [])
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                print(f"  [warn] chunk failed after retries ({exc}); returning empty")
                return []
            time.sleep(2.0 * (attempt + 1))


def fetch_activities(ligand_ids, chunk=50, pause=0.4):
    rows = []
    for i in range(0, len(ligand_ids), chunk):
        acts = query_chunk(ligand_ids[i:i + chunk])
        for a in acts:
            p = a.get("pchembl_value")
            if p is None:
                continue
            rows.append({
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "activity_id": a.get("activity_id"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "pchembl_value": float(p),
                "assay_chembl_id": a.get("assay_chembl_id"),
                "document_chembl_id": a.get("document_chembl_id"),
            })
        time.sleep(pause)
        print(f"  ... chunk {i // chunk + 1} done ({len(rows)} rows so far)", flush=True)
    return rows


def summarise(rows):
    if rows is None or len(rows) == 0:
        return {"best_pchembl": np.nan, "n_measurements": 0, "n_strong": 0,
                "types": ""}
    df = pd.DataFrame(rows)
    strong = df[df["pchembl_value"] >= STRONG_PCHEMBL]
    best = strong["pchembl_value"].max() if len(strong) else \
        df["pchembl_value"].max()
    return {"best_pchembl": round(float(best), 2),
            "n_measurements": int(len(df)),
            "n_strong": int(len(strong)),
            "types": ",".join(sorted(df["standard_type"].unique()))[:80]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranking", default=str(REPO_ROOT / "results/supervised_enopt_decoy/library_ranking_decoy_model.csv"))
    ap.add_argument("--library", default=str(REPO_ROOT / "data/raw/chembl37_ligand_library_30000.csv"))
    ap.add_argument("--out", default=str(REPO_ROOT / "results/retrospective_top200"))
    ap.add_argument("--top-n", type=int, default=200)
    ap.add_argument("--control-n", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--only-full", action="store_true",
                    help="only run the whole-library rank enrichment analysis")
    args = ap.parse_args()
    if args.only_full:
        full_library_retro(args)
        return

    rank = pd.read_csv(args.ranking).sort_values("enopt_decoy_rank")
    lib = pd.read_csv(args.library)
    known = set(rank.loc[rank["is_known_active"] == 1, "ligand_id"])
    top = rank.head(args.top_n)
    unknowns = rank.loc[~rank["ligand_id"].isin(known)]
    top_unknowns = top.loc[~top["ligand_id"].isin(known)]
    rng = np.random.default_rng(args.seed)
    pool = unknowns.loc[~unknowns["ligand_id"].isin(top_unknowns["ligand_id"])]
    control = pool.sample(n=min(args.control_n, len(pool)),
                          random_state=args.seed)
    print(f"[data] library {len(rank)} | known actives {len(known)} | "
          f"top {len(top)} (known-in-top {int(top['is_known_active'].sum())}) | "
          f"top unknowns to look up {len(top_unknowns)} | "
          f"control sample {len(control)}")

    lookup = pd.concat([top_unknowns, control])["ligand_id"].unique().tolist()
    print(f"[chembl] querying {len(lookup)} molecules in chunks ...")
    raw_rows = fetch_activities(lookup)
    raw = pd.DataFrame(raw_rows)
    if raw.empty:
        raw["molecule_chembl_id"] = []
    print(f"[chembl] returned {len(raw)} potency rows for "
          f"{raw['molecule_chembl_id'].nunique() if len(raw) else 0} molecules")

    per_mol = (raw.groupby("molecule_chembl_id")
               .apply(lambda g: pd.Series(summarise(g)), include_groups=False)
               .reset_index())

    def attach(frame):
        frame = frame.merge(per_mol, left_on="ligand_id",
                            right_on="molecule_chembl_id", how="left")
        frame["n_measurements"] = frame["n_measurements"].fillna(0).astype(int)
        frame["n_strong"] = frame["n_strong"].fillna(0).astype(int)
        frame["best_pchembl"] = frame["best_pchembl"].fillna(np.nan)
        frame["documented_strong"] = frame["n_strong"] >= 1
        return frame

    top_unknowns = attach(top_unknowns.copy())
    control = attach(control.copy())
    top_all = attach(top.copy())

    n_top_hit = int(top_unknowns["documented_strong"].sum())
    n_ctrl_hit = int(control["documented_strong"].sum())
    n_top = len(top_unknowns)
    n_ctrl = len(control)
    top_rate = n_top_hit / n_top
    ctrl_rate = n_ctrl_hit / n_ctrl
    ratio = top_rate / ctrl_rate if ctrl_rate > 0 else np.nan
    _, p = fisher_exact([[n_top_hit, n_top - n_top_hit],
                         [n_ctrl_hit, n_ctrl - n_ctrl_hit]],
                        alternative="greater")

    n_known_top = int(top_all["is_known_active"].sum())
    expect_known = args.top_n * len(known) / len(rank)
    print("\n== retrospective result (ChEMBL CHEMBL210, pChEMBL >= 6) ==")
    print(f"known in-library actives recovered in top {args.top_n}: "
          f"{n_known_top}/{len(known)} "
          f"(chance expectation {expect_known:.2f}, "
          f"enrichment {n_known_top / expect_known:.1f}x)  [sanity check]")
    print(f"no-label top-{args.top_n - n_known_top} documented strong: "
          f"{n_top_hit}/{n_top} ({top_rate * 100:.1f}%)")
    print(f"control (no-label remainder, n={n_ctrl}) documented strong: "
          f"{n_ctrl_hit}/{n_ctrl} ({ctrl_rate * 100:.1f}%)")
    print(f"enrichment ratio (top vs control): {ratio:.2f}x  "
          f"(Fisher one-sided p = {p:.4g})")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    desc = lib[["ligand_id", "smiles", "mol_weight", "heavy_atoms", "logp"]]
    top_unknowns.drop(columns=["enopt_decoy_rank"], errors="ignore") \
        .merge(desc, on="ligand_id", how="left").to_csv(
            out / "top200_no_label_with_chembl_hits.csv", index=False)
    control.to_csv(out / "control_sample_with_chembl_hits.csv", index=False)
    raw.to_csv(out / "chembl_lookup_raw_rows.csv", index=False)

    summary = pd.DataFrame([{
        "set": "no-label top", "n": n_top, "documented_strong": n_top_hit,
        "rate": round(top_rate, 4),
    }, {
        "set": "no-label control", "n": n_ctrl,
        "documented_strong": n_ctrl_hit, "rate": round(ctrl_rate, 4),
    }])
    summary.to_csv(out / "retro_summary.csv", index=False)

    # figure
    fig, ax = plt.subplots(figsize=(5, 4))
    labels = ["top-200\n(no-label)", "control\n(no-label, n=%d)" % n_ctrl]
    rates = [top_rate * 100, ctrl_rate * 100]
    bars = ax.bar(labels, rates, color=["#c0392b", "#7f8c8d"])
    ax.bar_label(bars, fmt="%.1f%%")
    ax.set_ylabel("fraction with documented\nhuman beta-2 activity (pChEMBL>=6)")
    ax.set_title(f"Retrospective hit rate: {ratio:.2f}x enrichment "
                 f"(Fisher p={p:.3g})")
    ax.set_ylim(0, max(max(rates) * 1.35, 5))
    fig.tight_layout()
    fig.savefig(out / "fig_top200_vs_base_rate.png", dpi=160)
    plt.close(fig)

    md = f"""# Retrospective top-{args.top_n} validation (ChEMBL lookup)

Lookup target: **CHEMBL210** (human beta-2 adrenergic receptor), potency types
{STANDARD_TYPES}, "documented strong" = at least one record with pChEMBL >= 6
(<= 1 uM). Library molecules never received a literature label at training
time, so these top-ranked molecules were **unknown to the model**.

| set | n | documented strong | hit rate |
|---|--:|--:|--:|
| no-label top | {n_top} | {n_top_hit} | {top_rate * 100:.1f}% |
| no-label control (random) | {n_ctrl} | {n_ctrl_hit} | {ctrl_rate * 100:.1f}% |

- Enrichment ratio (top vs control): **{ratio:.2f}x**, Fisher exact one-sided
  p = **{p:.4g}**.
- Sanity check only (these were training positives): known in-library actives
  in top {args.top_n} = {n_known_top}/{len(known)} vs chance {expect_known:.2f}
  ({n_known_top / expect_known:.1f}x).

Files: `top200_no_label_with_chembl_hits.csv` (ranked list with best pChEMBL),
`control_sample_with_chembl_hits.csv`, `chembl_lookup_raw_rows.csv`,
`retro_summary.csv`, `fig_top200_vs_base_rate.png`.

Caveats: retrospective literature recovery is a lower bound (not everything
measured is deposited in ChEMBL; BindingDB / patents not yet merged). It is
still a direct check of ranking power, at zero experimental cost.
"""
    (out / "retro_validation_report.md").write_text(md, encoding="utf-8")
    print(f"\nSaved to {out}")


def _markdown_table(df):
    lines = ["| " + " | ".join(map(str, df.columns)) + " |",
             "|" + "---|" * len(df.columns)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines)


def full_library_retro(args):
    """Whole-library retrospective.

    Pull every ADRB2 (CHEMBL210) potency row from ChEMBL, flag library
    molecules with pChEMBL >= 6, then ask: among library molecules that were
    NEVER given an activity label at training time, do the ChEMBL-documented
    actives rank higher than random in the Stage-3 leaderboard?
    """
    print("[full-lib] pulling all ADRB2 potency rows from ChEMBL ...",
          flush=True)
    rows, offset, limit = [], 0, 1000
    while True:
        params = {"target_chembl_id": "CHEMBL210",
                  "standard_type__in": STANDARD_TYPES,
                  "limit": limit, "offset": offset}
        url = API + "?" + urllib.parse.urlencode(params)
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "beta2ar-retro-validation/0.1"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = json.load(resp)
        except Exception as exc:  # noqa: BLE001
            print("  [warn] offset %d failed: %s" % (offset, exc))
            break
        acts = data.get("activities", [])
        rows.extend(acts)
        total = int(data["page_meta"]["total_count"])
        offset += limit
        print("  ... offset %d / %d" % (offset, total), flush=True)
        if offset >= total or not acts:
            break
        time.sleep(0.4)
    df = pd.DataFrame(rows)
    print("[full-lib] fetched %d rows for %d molecules"
          % (len(df), df["molecule_chembl_id"].nunique()))
    keep = df[df["pchembl_value"].notna()].copy()
    keep["pchembl_value"] = keep["pchembl_value"].astype(float)
    strong = keep[keep["pchembl_value"] >= STRONG_PCHEMBL]
    weak = keep[keep["pchembl_value"] < STRONG_PCHEMBL]
    doc_strong = set(strong["molecule_chembl_id"].unique())
    doc_weak = set(weak["molecule_chembl_id"].unique()) - doc_strong
    print("[full-lib] molecules pChEMBL>=6: %d | weak-only: %d"
          % (len(doc_strong), len(doc_weak)))

    rank = pd.read_csv(args.ranking).sort_values("enopt_decoy_rank")
    act48 = set(pd.read_csv(
        REPO_ROOT / "data/training/beta2ar_actives_in_library.csv")
        ["molecule_chembl_id"]) & set(rank["ligand_id"])
    unk = rank.loc[~rank["ligand_id"].isin(act48)].copy()
    unk["doc_strong"] = unk["ligand_id"].isin(doc_strong).astype(int)
    unk["doc_weak"] = unk["ligand_id"].isin(doc_weak).astype(int)
    D = int(unk["doc_strong"].sum())
    n = len(unk)
    auc = roc_auc_score(unk["doc_strong"], -unk["enopt_decoy_rank"])
    mann = mannwhitneyu(
        unk.loc[unk["doc_strong"] == 1, "enopt_decoy_rank"],
        unk.loc[unk["doc_strong"] == 0, "enopt_decoy_rank"],
        alternative="less")
    print("\n== whole-library retrospective (documented actives excluded "
          "from training) ==")
    print("documented strong actives among library unknowns: %d / %d" % (D, n))
    print("AUC (documented vs rest): %.4f | Mann-Whitney p = %.3g"
          % (auc, mann.pvalue))
    med = unk.loc[unk["doc_strong"] == 1, "enopt_decoy_rank"].median()
    print("median rank of documented actives: %.0f / %d (%.1f%%; random=50%%)"
          % (med, n, med / n * 100))
    topn_rows = []
    for nn in (50, 200, 500, 1000, 2000, 5000):
        exp = D * nn / n
        obs = int(unk.head(nn)["doc_strong"].sum())
        topn_rows.append({"top_n": nn, "documented_in_top": obs,
                          "chance_expect": round(exp, 2),
                          "enrichment": round(obs / exp, 2)
                          if exp else float("nan")})
    topn = pd.DataFrame(topn_rows)
    print(topn.to_string(index=False))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    unk[["ligand_id", "enopt_decoy_score", "enopt_decoy_rank",
         "doc_strong", "doc_weak"]].to_csv(
        out / "full_library_documented_flags.csv", index=False)
    keep.to_csv(out / "chembl_adrb2_all_rows.csv", index=False)
    pd.DataFrame([
        {"metric": "documented_strong_unknowns", "value": D},
        {"metric": "library_unknowns", "value": n},
        {"metric": "auc_doc_vs_rest", "value": round(auc, 4)},
        {"metric": "mannwhitney_p", "value": mann.pvalue},
        {"metric": "median_rank_pct", "value": round(med / n, 4)},
    ]).to_csv(out / "full_library_retro_summary.csv", index=False)

    md = (
        "# Whole-library retrospective validation (ChEMBL ADRB2)\n\n"
        "Pulled all %d ADRB2 (CHEMBL210) potency rows (types %s), flagged "
        "library molecules with pChEMBL >= 6 (<= 1 uM), excluded the training "
        "actives, and tested whether the remaining documented actives rank "
        "high in the Stage-3 leaderboard.\n\n"
        "- Documented strong actives among library unknowns: **%d** / %d\n"
        "- Rank-AUC (documented vs rest): **%.4f** (0.5 = random)\n"
        "- Mann-Whitney one-sided p: **%.3g**\n"
        "- Median rank of documented actives: **%.0f** / %d (**%.1f%%**; "
        "random = 50%%)\n\n"
        "Top-N enrichment:\n\n"
        "%s\n\n"
        "Files: `full_library_documented_flags.csv`, "
        "`chembl_adrb2_all_rows.csv`, `full_library_retro_summary.csv`.\n"
        % (len(df), STANDARD_TYPES, D, n, auc, mann.pvalue,
           med, n, med / n * 100, _markdown_table(topn))
    )
    (out / "full_library_retro_report.md").write_text(md, encoding="utf-8")
    print("\nSaved full-library analysis to %s" % out)


if __name__ == "__main__":
    main()
