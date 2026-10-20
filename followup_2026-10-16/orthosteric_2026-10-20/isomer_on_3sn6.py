#!/usr/bin/env python3
"""Phenolic O-O distances of the other five regioisomers.

Only an isomer whose phenolic oxygens sit in 4.7-4.9 A is placed on the
already chosen 3SN6 side chains, chi1 -167 and 49. The same hydrogen-bond
and 2.6 A heavy-atom rules are used. No free-energy calculation.
"""

import csv
import importlib.util
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("pr", HERE / "place_representative.py")
pr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pr)
OUT = HERE / "isomer_on_3sn6.json"
CSV = Path(
    r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs"
    r"\outputs\generation_2026-10-16\regioisomers.csv"
)
CHI = {"203": -167, "207": 49}


def protonated(smiles):
    return smiles.replace("CC(C)NCC", "CC(C)[NH2+]CC", 1)


def phenolic_distance(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(protonated(smiles)))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261020
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=30, params=params))
    oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and any(nbr.GetIsAromatic() for nbr in atom.GetNeighbors()):
            oxygens.append(atom.GetIdx())
    distances = []
    for cid in ids:
        AllChem.MMFFOptimizeMolecule(mol, confId=int(cid))
        conf = mol.GetConformer(int(cid))

        def xyz(idx, conf=conf):
            pos = conf.GetAtomPosition(int(idx))
            return (pos.x, pos.y, pos.z)

        distances.append(pr.sp.sc._dist(xyz(oxygens[0]), xyz(oxygens[1])))
    return [round(min(distances), 3), round(max(distances), 3)]


def confs_for(smiles):
    saved = pr.sp.SMILES
    pr.sp.SMILES = protonated(smiles)
    try:
        return pr.ligand_confs()
    finally:
        pr.sp.SMILES = saved


def main():
    rows = list(csv.DictReader(CSV.open(encoding="utf-8")))
    measured = []
    for row in rows:
        if row["isomer_id"] == "ISO06":
            continue
        span = phenolic_distance(row["smiles"])
        measured.append(
            {
                "isomer_id": row["isomer_id"],
                "positions": row["positions"],
                "known_id": row["known_id"],
                "phenolic_O_O": span,
                "inside_4.7_4.9": span[0] >= 4.7 and span[1] <= 4.9,
            }
        )
        print(row["isomer_id"], row["positions"], span, flush=True)

    path = pr.sp.sc.ROOT / pr.FILES["3SN6"]
    serines = pr.sp.motif_serines(path)
    asp = pr.find_asp(path)
    _, protein, _, side, asp_atoms = pr.load_side(path, serines, asp)
    ogs = {crystal: pr.og_at(side, crystal, CHI[crystal]) for crystal in CHI}
    og_distance = round(pr.sp.sc._dist(ogs["203"], ogs["207"]), 3)
    posed, cells = pr.receptor_with_ogs(protein, side, ("203", "207"), ogs)
    placed = []
    for row, source in zip(measured, [item for item in rows if item["isomer_id"] != "ISO06"]):
        if not row["inside_4.7_4.9"]:
            row["placement"] = None
            continue
        found = pr.search(
            confs_for(source["smiles"]),
            posed,
            cells,
            [asp_atoms["OD1"], asp_atoms["OD2"]],
            asp_atoms["CG"],
            ogs["203"],
            ogs["207"],
            side["203"]["file_resi"],
            side["207"]["file_resi"],
            asp["file_resi"],
        )
        row["placement"] = found
        placed.append(row["isomer_id"])
        print(row["isomer_id"], "placement", found, flush=True)
    payload = {
        "fixed_chi1": CHI,
        "serine_OG_OG": og_distance,
        "isomers": measured,
        "placed": placed,
        "passed_same_standard": [row["isomer_id"] for row in measured if row.get("placement") and row["placement"].get("passes")],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("passed", payload["passed_same_standard"], flush=True)


if __name__ == "__main__":
    main()
