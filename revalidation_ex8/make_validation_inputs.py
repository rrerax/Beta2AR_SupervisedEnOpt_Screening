#!/usr/bin/env python3
"""Build the merged validation ligand set for ex8 re-docking.

Merges 48 in-library actives + 158 expanded ChEMBL actives + 3,000 DUD-E
decoys into one CSV and writes a config for scripts/03 + 04 that points at it
(num_modes 20 so pose ensembles can be analysed afterwards).

Usage (on the cloud machine, from repo root):
    python make_validation_inputs.py --repo .
"""
import argparse
from pathlib import Path

import pandas as pd
import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--repo", default=".")
args = ap.parse_args()
REPO = Path(args.repo)

a48 = pd.read_csv(REPO / "data/training/beta2ar_actives_in_library.csv")
a158 = pd.read_csv(REPO / "data/training/beta2ar_actives_expand_docked.csv")
dec = pd.read_csv(REPO / "data/training/dude_decoys_for_dock.csv")

cols = {"ligand_id", "smiles"}
a48 = a48.rename(columns={"molecule_chembl_id": "ligand_id"})[["ligand_id", "smiles"]]
merged = pd.concat([a48, a158[["ligand_id", "smiles"]], dec[["ligand_id", "smiles"]]])
merged = merged.drop_duplicates("ligand_id").reset_index(drop=True)
out_csv = REPO / "data/training/validation_redock_ligands.csv"
merged.to_csv(out_csv, index=False)
n_act = int(merged["ligand_id"].str.startswith("CHEMBL").sum())
n_dec = len(merged) - n_act
print(f"validation ligands written: {len(merged)} total "
      f"({n_act} actives, {n_dec} decoys) -> {out_csv}")

cfg = yaml.safe_load((REPO / "configs/decoy_screen.yml").read_text())
cfg["ligands"]["input_csv"] = "data/training/validation_redock_ligands.csv"
cfg["ligands"]["full_limit"] = 40000
cfg["docking"]["num_modes"] = 20
cfg["docking"]["exhaustiveness"] = 8
out_yml = REPO / "configs/validation_redock.yml"
out_yml.write_text(yaml.safe_dump(cfg, sort_keys=False))
print(f"config written: {out_yml}")
