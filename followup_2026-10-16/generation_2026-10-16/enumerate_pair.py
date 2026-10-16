#!/usr/bin/env python3
"""Enumerate phenolic regioisomers of the CHEMBL434 / orciprenaline pair.

The window is fixed before any structure is kept. The ethanolamine side
chain stays in place. Only the two phenolic hydroxy positions change.
No stereoisomers are added. No docking score is used for selection.
"""
import json
from itertools import combinations
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

OUT = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs\outputs\generation_2026-10-16")
OUT.mkdir(parents=True, exist_ok=True)

WINDOW = {
    "parent_smiles": "CC(C)NCC(O)c1ccccc1",
    "phenolic_oh_count": 2,
    "mw_center": 211.26,
    "mw_tolerance": 0.5,
    "logp_center": 1.13,
    "logp_tolerance": 0.5,
    "hbd": 4,
    "hba": 4,
    "hac": 15,
    "max_molecules": 30,
    "stereoisomers": False,
}
(OUT / "generation_window.json").write_text(json.dumps(WINDOW, indent=2), encoding="utf-8")

KNOWN = {
    "CHEMBL434": "CC(C)NCC(O)c1ccc(O)c(O)c1",
    "CHEMBL776": "CC(C)NCC(O)c1cc(O)cc(O)c1",
}


def props(mol):
    return {
        "mw": round(float(Descriptors.MolWt(mol)), 2),
        "logp": round(float(Crippen.MolLogP(mol)), 2),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hac": int(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
    }


def in_window(row):
    return (
        abs(row["mw"] - WINDOW["mw_center"]) <= WINDOW["mw_tolerance"]
        and abs(row["logp"] - WINDOW["logp_center"]) <= WINDOW["logp_tolerance"]
        and row["hbd"] == WINDOW["hbd"]
        and row["hba"] == WINDOW["hba"]
        and row["hac"] == WINDOW["hac"]
    )


base = Chem.MolFromSmiles(WINDOW["parent_smiles"])
ring = [atom.GetIdx() for atom in base.GetAromaticAtoms()]
attachment = None
for atom in base.GetAtoms():
    if not atom.GetIsAromatic():
        continue
    for neighbor in atom.GetNeighbors():
        if not neighbor.GetIsAromatic():
            attachment = atom.GetIdx()
            break
order = []
seen = {attachment}
current = attachment
while len(order) < len(ring):
    order.append(current)
    nxt = None
    for neighbor in base.GetAtomWithIdx(current).GetNeighbors():
        if neighbor.GetIsAromatic() and neighbor.GetIdx() not in seen:
            nxt = neighbor.GetIdx()
            break
    if nxt is None:
        break
    seen.add(nxt)
    current = nxt
# order[0] is the attachment carbon. Ring positions 2..6 are order[1:].
position_to_idx = {pos: order[pos - 1] for pos in range(1, 7)}

def phenolic_positions(mol):
    attachment = None
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and any(not n.GetIsAromatic() for n in atom.GetNeighbors()):
            attachment = atom.GetIdx()
            break
    order = [attachment]
    seen = {attachment}
    current = attachment
    while len(order) < 6:
        nxt = None
        for neighbor in mol.GetAtomWithIdx(current).GetNeighbors():
            if neighbor.GetIsAromatic() and neighbor.GetIdx() not in seen:
                nxt = neighbor.GetIdx()
                break
        seen.add(nxt)
        order.append(nxt)
        current = nxt
    found = []
    for pos, idx in enumerate(order, start=1):
        for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
            if neighbor.GetSymbol() == "O" and neighbor.GetTotalNumHs() >= 1:
                found.append(pos)
    mirrored = tuple(sorted(8 - pos for pos in found))
    return ",".join(str(pos) for pos in min(tuple(found), mirrored))


rows = []
for combo in combinations(range(2, 7), 2):
    mol = Chem.RWMol(base)
    for pos in combo:
        oxygen = mol.AddAtom(Chem.Atom(8))
        mol.AddBond(position_to_idx[pos], oxygen, Chem.BondType.SINGLE)
    parsed = Chem.MolFromSmiles(Chem.MolToSmiles(mol, canonical=True))
    row = props(parsed)
    row.update(
        {
            "positions": phenolic_positions(parsed),
            "smiles": Chem.MolToSmiles(parsed, canonical=True),
        }
    )
    rows.append(row)

# Collapse ring-flip duplicates.
unique = {}
for row in rows:
    unique.setdefault(row["smiles"], row)
kept = []
for row in unique.values():
    row["in_window"] = bool(in_window(row))
    row["known_id"] = ""
    for chembl_id, raw in KNOWN.items():
        if Chem.MolToSmiles(Chem.MolFromSmiles(raw), canonical=True) == row["smiles"]:
            row["known_id"] = chembl_id
    if row["in_window"]:
        kept.append(row)

kept = sorted(kept, key=lambda row: row["positions"])[: WINDOW["max_molecules"]]
for i, row in enumerate(kept, start=1):
    row["isomer_id"] = f"ISO{i:02d}"

import pandas as pd

out = pd.DataFrame(kept)
out.to_csv(OUT / "regioisomers.csv", index=False)
summary = {
    "window": WINDOW,
    "raw_combinations": len(rows),
    "unique_constitutions": len(unique),
    "kept": len(kept),
    "known_found": sorted(row["known_id"] for row in kept if row["known_id"]),
}
(OUT / "generation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(out.to_string(index=False))
print(json.dumps(summary, indent=2))
