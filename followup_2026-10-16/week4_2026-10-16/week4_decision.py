#!/usr/bin/env python3
"""Week 4: decision page from the week 1-3 tables. No new docking."""
import json
from pathlib import Path

import pandas as pd

WS = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs\outputs")
OUT = WS / "week4_2026-10-16"
OUT.mkdir(parents=True, exist_ok=True)

week1 = json.loads((WS / "week1_2026-10-16" / "week1_summary.json").read_text(encoding="utf-8"))
metrics = pd.read_csv(WS / "week1_2026-10-16" / "week1_metrics.csv")
props = pd.read_csv(WS / "week2_2026-10-16" / "week2_properties.csv").set_index("ligand_id")
held = pd.read_csv(WS / "week3_2026-10-16" / "week3_holdout_ranks.csv")
week3 = json.loads((WS / "week3_2026-10-16" / "week3_summary.json").read_text(encoding="utf-8"))

singles = metrics[metrics.method.str.startswith("single_")].sort_values("auroc", ascending=False)
best_single = singles.iloc[0]
mean_ef = float(week1["oriented_mean_ef1"])
use_learned = bool(
    week1["learned_ef1_above_mean_every_seed"]
    and week1["descriptor_ef1_below_learned_every_seed"]
)
pair = props.loc[["CHEMBL434", "CHEMBL776"]]
same_properties = bool(
    (pair.loc["CHEMBL434", ["mw", "logp", "hbd", "hba", "hac"]].to_numpy()
     == pair.loc["CHEMBL776", ["mw", "logp", "hbd", "hba", "hac"]].to_numpy()).all()
)
held_ahead = int((held.rank_434_held_out < held.rank_776).sum())

decision = {
    "use_learned_ranker": use_learned,
    "global_score": "molecular weight, logP, and heavy-atom count",
    "docking_score_if_required": "oriented five-conformation average",
    "learned_ef1_range": [week1["learned_ef1_min"], week1["learned_ef1_max"]],
    "oriented_mean_ef1": mean_ef,
    "descriptor_ef1_range": [week1["descriptor_ef1_min"], week1["descriptor_ef1_max"]],
    "learned_auroc_range": [week1["learned_auroc_min"], week1["learned_auroc_max"]],
    "descriptor_auroc_range": [week1["descriptor_auroc_min"], week1["descriptor_auroc_max"]],
    "best_single_conformation": best_single.method,
    "best_single_auroc": float(best_single.auroc),
    "best_single_ef1": float(best_single.ef_1pct),
    "keep_4_hydroxy_over_5_hydroxy": False,
    "properties_identical": same_properties,
    "full_model_ranks": {
        "CHEMBL434": int(props.loc["CHEMBL434", "rank_learned"]),
        "CHEMBL776": int(props.loc["CHEMBL776", "rank_learned"]),
    },
    "oriented_mean_ranks": {
        "CHEMBL434": int(props.loc["CHEMBL434", "rank_oriented_mean"]),
        "CHEMBL776": int(props.loc["CHEMBL776", "rank_oriented_mean"]),
    },
    "heldout_434_ahead_count": held_ahead,
    "heldout_seeds": int(len(held)),
    "heldout_rank_434": [week3["heldout_rank_434_min"], week3["heldout_rank_434_max"]],
    "heldout_rank_776": [week3["heldout_rank_776_min"], week3["heldout_rank_776_max"]],
    "side_chain_relaxation": "not indicated",
}
(OUT / "week4_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
print(json.dumps(decision, indent=2))
