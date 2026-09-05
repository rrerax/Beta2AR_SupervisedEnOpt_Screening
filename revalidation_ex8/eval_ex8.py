#!/usr/bin/env python3
"""Evaluate an ex8 re-docked validation score file exactly like the Stage-3/4
protocols, so numbers are directly comparable:

  clean-label 3-fold OOF, 206 actives vs 2,978 DUD-E decoys, XGBoost settings
  identical to scripts/07.  Reports BOTH the 5 raw-score model (ex4 -> 0.696)
  and the raw + per-heavy-atom model (ex4 -> 0.751).

Usage (cloud machine, repo root):
    python eval_ex8.py --repo . --score results/tables/docking_scores_ex8_validation.csv
"""
import argparse
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ap = argparse.ArgumentParser()
ap.add_argument("--repo", default=".")
ap.add_argument("--score", required=True)
args = ap.parse_args()
REPO = Path(args.repo)

spec = importlib.util.spec_from_file_location("s7", REPO / "scripts/07_train_validated_enopt.py")
s7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s7)

wide = s7.load_wide(args.score)
print(f"[data] ex8 score rows -> wide matrix {wide.shape}")
dec = wide[wide.index.str.startswith("decoy_")]
act = wide[~wide.index.str.startswith("decoy_")]
Xc = wide.copy()
yc = pd.Series((~Xc.index.str.startswith("decoy_")).astype(int), index=Xc.index)
print(f"[data] actives {len(act)} | decoys {len(dec)}")

# heavy-atom counts from the merged smiles file
smi = pd.read_csv(REPO / "data/training/validation_redock_ligands.csv")
def hac(s):
    m = Chem.MolFromSmiles(s)
    if m is None:
        return np.nan
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda f: f.GetNumHeavyAtoms()).GetNumHeavyAtoms() if frags else np.nan
smi["hac"] = smi["smiles"].apply(hac)
hmap = smi.set_index("ligand_id")["hac"]

cols = list(wide.columns)
raw = Xc[cols]
X = Xc[hmap.reindex(Xc.index).notna()]
y = yc.reindex(X.index)
h = hmap.reindex(X.index)
ha = X[cols].div(h, axis=0)
ha.columns = [c + "__per_ha" for c in ha.columns]

def report(tag, Xf, y):
    oof, imp, w = s7.oof_predict(Xf, y, n_splits=3, seed=42)
    r = s7.metrics_row(tag, pd.Series(oof, index=y.index), y)
    print(r)
    return r

print("\n== ex8 re-docked validation (clean-label OOF) ==")
r_raw = report("EnOpt_raw_ex8", raw.reindex(X.index), y)
r_ha = report("EnOpt_raw_ha_ex8", pd.concat([raw.reindex(X.index), ha], axis=1), y)
print("\nreference (ex4): raw 0.696 (0.655-0.737), raw+ha 0.751 (0.712-0.791)")
print("\nNOTE: if ex8 CI overlaps ex4 CI, higher exhaustiveness did not change "
      "discrimination; if it is clearly above, the ex4 model was limited by "
      "docking quality. Also check raw+ha vs raw: if the per-ha gap narrows, "
      "part of the size-bias correction was a docking-sampling artefact.")
