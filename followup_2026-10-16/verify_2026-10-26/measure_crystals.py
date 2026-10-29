#!/usr/bin/env python3
"""Measure deposited catecholamine contacts. No docking and no atom moves.

Heavy atoms stay at deposited coordinates. The crystal files do not contain
the ligand hydrogens, so phenolic hydrogens are built in the ring plane at
0.97 A with a 109.5 degree C-O-H angle, both flips, and a secondary amine
receives two tetrahedral hydrogens at 1.01 A. A contact can use either phenol
flip. The pass rule is unchanged: amine hydrogen to Asp OD1/OD2, one phenolic
hydrogen each to the motif Ser203 and Ser207, every other ligand-protein
heavy contact at least 2.6 A. Ser204 and bridging waters are recorded and
cannot make a failing pose pass.
"""

from __future__ import annotations

import json
import math
import urllib.request
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CRYSTAL_DIR = HERE / "crystals"
HBOND_MAX_A = 2.5
HBOND_MIN_DEG = 120.0
CLEARANCE_MIN_A = 2.6
WATER_BRIDGE_A = 3.2
OH_LENGTH = 0.97
NH_LENGTH = 1.01
COH_DEG = 109.5
SER_MOTIFS = ("AIASSIVSF", "AIASSIISF")
ASP_MOTIFS = (("EFWTSIDVLCVTA", 6), ("WTSLDVLCVTAS", 4), ("SLDVLCVTA", 2), ("SIDVLCVTA", 2))
AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
SKIP = {"HOH", "WAT", "NA", "CL", "ZN", "MG", "CA", "SO4", "GOL", "EDO", "PEG", "ACT", "OLC", "OLB", "CLR", "CDL"}
SOURCES = {
    "2Y03": "https://files.rcsb.org/download/2Y03.pdb",
    "4LDO": "https://files.rcsb.org/download/4LDO.pdb",
    "4LDL": None,
}


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, value):
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    return scale(a, 1.0 / length)


def dist(a, b):
    return norm(sub(a, b))


def angle(donor, hydrogen, acceptor):
    ba, bc = sub(donor, hydrogen), sub(acceptor, hydrogen)
    denom = norm(ba) * norm(bc)
    cosv = max(-1.0, min(1.0, dot(ba, bc) / denom))
    return math.degrees(math.acos(cosv))


def contact(donor, hydrogen, acceptor):
    distance = dist(hydrogen, acceptor)
    value = angle(donor, hydrogen, acceptor)
    return {
        "passes": distance <= HBOND_MAX_A and value >= HBOND_MIN_DEG,
        "distance": round(distance, 3),
        "angle": round(value, 1),
    }


def parse_pdb(path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        alt = line[16]
        if alt not in " A":
            continue
        element = line[76:78].strip() if len(line) >= 78 else ""
        if not element:
            element = "".join(ch for ch in line[12:16] if ch.isalpha())[:1]
        atoms.append(
            {
                "record": line[:6].strip(),
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "element": element.upper(),
            }
        )
    return atoms


def download(pdb_id):
    CRYSTAL_DIR.mkdir(parents=True, exist_ok=True)
    destination = CRYSTAL_DIR / f"{pdb_id}.pdb"
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    url = SOURCES[pdb_id]
    request = urllib.request.Request(url, headers={"User-Agent": "beta2ar-verify"})
    with urllib.request.urlopen(request, timeout=60) as response:
        destination.write_bytes(response.read())
    return destination


def motifs(atoms):
    rows = []
    seen = set()
    for atom in atoms:
        if atom["record"] != "ATOM" or atom["name"] != "CA":
            continue
        key = (atom["chain"], atom["resi"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(atom)
    by_chain = defaultdict(list)
    for atom in rows:
        by_chain[atom["chain"]].append(atom)
    found = []
    for chain, items in by_chain.items():
        seq = "".join(AA.get(atom["resn"], "X") for atom in items)
        ser_at = next((seq.find(motif) for motif in SER_MOTIFS if seq.find(motif) >= 0), -1)
        asp_hit = next(((seq.find(motif), offset) for motif, offset in ASP_MOTIFS if seq.find(motif) >= 0), None)
        if ser_at < 0 or asp_hit is None:
            continue
        serine = {
            "203": items[ser_at + 3],
            "204": items[ser_at + 4],
            "207": items[ser_at + 7],
        }
        found.append({"chain": chain, "asp": items[asp_hit[0] + asp_hit[1]], "serine": serine})
    return found


def bonded(atoms, cutoff=1.7):
    heavy = [atom for atom in atoms if atom["element"] != "H"]
    adj = [[] for _ in heavy]
    for i, left in enumerate(heavy):
        for j in range(i + 1, len(heavy)):
            right = heavy[j]
            limit = 1.9 if "S" in {left["element"], right["element"]} else cutoff
            if dist(left["xyz"], right["xyz"]) <= limit:
                adj[i].append(j)
                adj[j].append(i)
    return heavy, adj


def rings6(adj):
    found = []
    n = len(adj)

    def rec(start, current, visited):
        for nxt in adj[current]:
            if nxt == start and len(visited) == 6:
                found.append(list(visited))
            elif nxt > start and nxt not in visited and len(visited) < 6:
                visited.append(nxt)
                rec(start, nxt, visited)
                visited.pop()

    for start in range(n):
        rec(start, start, [start])
    unique = []
    seen = set()
    for ring in found:
        key = tuple(sorted(ring))
        if key not in seen and len(key) == 6:
            seen.add(key)
            unique.append(ring)
    return unique


def build_hydrogens(heavy, adj):
    carbon_rings = set()
    for ring in rings6(adj):
        if all(heavy[index]["element"] == "C" for index in ring):
            carbon_rings.update(ring)
    phenols = []
    amines = []
    for index, atom in enumerate(heavy):
        if atom["element"] == "O" and len(adj[index]) == 1:
            carbon = adj[index][0]
            if carbon not in carbon_rings:
                continue
            others = [nbr for nbr in adj[carbon] if nbr != index and nbr in carbon_rings]
            if not others:
                continue
            plane = unit(cross(sub(atom["xyz"], heavy[carbon]["xyz"]), sub(heavy[others[0]]["xyz"], heavy[carbon]["xyz"])))
            oc = unit(sub(heavy[carbon]["xyz"], atom["xyz"]))
            inplane = unit(cross(plane, oc))
            theta = math.radians(COH_DEG)
            hydrogens = []
            for sign in (1.0, -1.0):
                direction = unit(add(scale(oc, math.cos(theta)), scale(inplane, math.sin(theta) * sign)))
                hydrogens.append(add(atom["xyz"], scale(direction, OH_LENGTH)))
            phenols.append({"oxygen": atom, "hydrogens": hydrogens})
        if atom["element"] == "N":
            carbons = [heavy[nbr]["xyz"] for nbr in adj[index] if heavy[nbr]["element"] == "C"]
            if len(carbons) != 2:
                amines.append({"nitrogen": atom, "heavy_neighbors": len(adj[index]), "hydrogens": []})
                continue
            u = unit(sub(carbons[0], atom["xyz"]))
            v = unit(sub(carbons[1], atom["xyz"]))
            half = scale(sub((0.0, 0.0, 0.0), add(u, v)), 0.5)
            side = unit(cross(u, v))
            length2 = 1.0 - dot(half, half)
            span = math.sqrt(max(0.0, length2))
            hydrogens = []
            for sign in (1.0, -1.0):
                direction = unit(add(half, scale(side, span * sign)))
                hydrogens.append(add(atom["xyz"], scale(direction, NH_LENGTH)))
            amines.append({"nitrogen": atom, "heavy_neighbors": 2, "hydrogens": hydrogens})
    return phenols, amines


def nearest_ligand(atoms, asp_xyz):
    groups = defaultdict(list)
    for atom in atoms:
        if atom["record"] != "HETATM" or atom["resn"] in SKIP or atom["element"] == "H":
            continue
        groups[(atom["chain"], atom["resi"], atom["resn"])].append(atom)
    ranked = []
    for key, group in groups.items():
        carbons = sum(atom["element"] == "C" for atom in group)
        if carbons < 6:
            continue
        closest = min(dist(atom["xyz"], asp_xyz) for atom in group)
        ranked.append((closest, key, group))
    ranked.sort(key=lambda item: item[0])
    return ranked[:1]


def measure_one(path, pdb_id):
    atoms = parse_pdb(path)
    reports = []
    found = motifs(atoms)
    if not found:
        return [{"pdb": pdb_id, "status": "motif_not_found"}]
    for motif in found:
        chain_atoms = [atom for atom in atoms if atom["record"] == "ATOM" and atom["chain"] == motif["chain"] and atom["element"] != "H"]
        asp_atoms = [
            atom for atom in chain_atoms if atom["resi"] == motif["asp"]["resi"] and atom["name"] in {"OD1", "OD2", "CG"}
        ]
        asp_od = [atom["xyz"] for atom in asp_atoms if atom["name"] in {"OD1", "OD2"}]
        cg = next(atom["xyz"] for atom in asp_atoms if atom["name"] == "CG")
        ser = {}
        for label, ca in motif["serine"].items():
            og = next(
                atom["xyz"]
                for atom in chain_atoms
                if atom["resi"] == ca["resi"] and atom["resn"] == "SER" and atom["name"] == "OG"
            )
            ser[label] = {"file_resi": ca["resi"], "og": og}
        chosen = nearest_ligand(atoms, cg)
        if not chosen:
            reports.append({"pdb": pdb_id, "chain": motif["chain"], "status": "no_ligand"})
            continue
        _distance, key, ligand_atoms = chosen[0]
        heavy, adj = bonded(ligand_atoms)
        phenols, amines = build_hydrogens(heavy, adj)
        amine_hits = []
        for amine in amines:
            for hydrogen in amine["hydrogens"]:
                for od in asp_od:
                    amine_hits.append(contact(amine["nitrogen"]["xyz"], hydrogen, od))
        amine_best = max(amine_hits, key=lambda item: (item["passes"], -item["distance"])) if amine_hits else None
        phenol_table = []
        for phenol_index, phenol in enumerate(phenols):
            for label, target in ser.items():
                flips = [contact(phenol["oxygen"]["xyz"], hydrogen, target["og"]) for hydrogen in phenol["hydrogens"]]
                best = max(flips, key=lambda item: (item["passes"], -item["distance"]))
                phenol_table.append(
                    {
                        "phenol_index": phenol_index,
                        "oxygen_name": phenol["oxygen"]["name"],
                        "serine": label,
                        "file_resi": target["file_resi"],
                        "oxygen_og": round(dist(phenol["oxygen"]["xyz"], target["og"]), 3),
                        "best_flip": best,
                    }
                )
        pair_pass = False
        pair_used = None
        pair_index = None
        if len(phenols) >= 2:
            for i in range(len(phenols)):
                for j in range(len(phenols)):
                    if i == j:
                        continue
                    left = next(item for item in phenol_table if item["phenol_index"] == i and item["serine"] == "203")
                    right = next(item for item in phenol_table if item["phenol_index"] == j and item["serine"] == "207")
                    if left["best_flip"]["passes"] and right["best_flip"]["passes"]:
                        pair_pass = True
                        pair_used = [phenols[i]["oxygen"]["name"], phenols[j]["oxygen"]["name"]]
                        pair_index = (i, j)
                        break
                if pair_pass:
                    break
        allowed = []
        if amine_best and amine_best["passes"]:
            allowed.extend((amine["nitrogen"]["xyz"], od) for amine in amines for od in asp_od)
        if pair_index is not None:
            allowed.append((phenols[pair_index[0]]["oxygen"]["xyz"], ser["203"]["og"]))
            allowed.append((phenols[pair_index[1]]["oxygen"]["xyz"], ser["207"]["og"]))

        def excluded(left, right):
            return any(dist(left, donor) < 0.05 and dist(right, acceptor) < 0.05 for donor, acceptor in allowed)

        ligand_heavies = [atom["xyz"] for atom in heavy]
        protein_heavies = [(atom, atom["xyz"]) for atom in chain_atoms]
        closest = min(
            (
                (dist(left, right), heavy[index]["name"], atom["resn"], atom["resi"], atom["name"])
                for index, left in enumerate(ligand_heavies)
                for atom, right in protein_heavies
                if not excluded(left, right)
            ),
            key=lambda item: item[0],
        )
        clearance = closest[0]
        waters = []
        nearest_water = None
        for atom in atoms:
            if atom["resn"] not in {"HOH", "WAT"} or atom["element"] != "O":
                continue
            for phenol in phenols:
                for label, target in ser.items():
                    d_lig = dist(atom["xyz"], phenol["oxygen"]["xyz"])
                    d_ser = dist(atom["xyz"], target["og"])
                    bridge = max(d_lig, d_ser)
                    if nearest_water is None or bridge < nearest_water["bridge"]:
                        nearest_water = {
                            "resi": atom["resi"],
                            "phenol": phenol["oxygen"]["name"],
                            "serine": label,
                            "to_phenol_o": round(d_lig, 3),
                            "to_og": round(d_ser, 3),
                            "bridge": round(bridge, 3),
                        }
                    if d_lig <= WATER_BRIDGE_A and d_ser <= WATER_BRIDGE_A:
                        waters.append(
                            {
                                "resi": atom["resi"],
                                "phenol": phenol["oxygen"]["name"],
                                "serine": label,
                                "to_phenol_o": round(d_lig, 3),
                                "to_og": round(d_ser, 3),
                            }
                        )
        reports.append(
            {
                "pdb": pdb_id,
                "chain": motif["chain"],
                "ligand": {"resn": key[2], "chain": key[0], "resi": key[1], "heavy_atoms": len(heavy)},
                "asp_file_resi": motif["asp"]["resi"],
                "serine_file_resi": {label: item["file_resi"] for label, item in ser.items()},
                "amine_neighbors": [amine["heavy_neighbors"] for amine in amines],
                "phenol_count": len(phenols),
                "amine_best": amine_best,
                "phenol_contacts": phenol_table,
                "pair_203_207": pair_pass,
                "pair_oxygens": pair_used,
                "min_heavy_contact": round(clearance, 3),
                "closest_other": {
                    "distance": round(closest[0], 3),
                    "ligand_atom": closest[1],
                    "residue": closest[2],
                    "file_resi": closest[3],
                    "atom": closest[4],
                },
                "frozen_rule_pass": bool(
                    pair_pass and amine_best and amine_best["passes"] and clearance >= CLEARANCE_MIN_A
                ),
                "bridging_waters": waters,
                "nearest_water": nearest_water,
            }
        )
    return reports


def main():
    paths = {
        "2Y03": download("2Y03"),
        "4LDO": download("4LDO"),
        "4LDL": REPO / "data" / "raw" / "receptors" / "4LDL.pdb",
    }
    reports = []
    for pdb_id, path in paths.items():
        reports.extend(measure_one(path, pdb_id))
    payload = {
        "rules": {
            "heavy_atoms": "deposited coordinates, altLoc blank or A, first model",
            "phenol_h": "both in-plane flips, O-H 0.97 A, C-O-H 109.5 deg",
            "amine_h": "two tetrahedral H at 1.01 A when nitrogen has two carbon neighbors",
            "pass": "amine to Asp OD and both Ser203 and Ser207, clearance >= 2.6 A",
            "ser204_and_water": "recorded, not allowed to turn a failure into a pass",
            "docking": "none",
        },
        "structures": reports,
    }
    (HERE / "crystal_contacts.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    brief = [
        {
            "pdb": item["pdb"],
            "ligand": item.get("ligand", {}).get("resn"),
            "phenols": item.get("phenol_count"),
            "amine": None if item.get("amine_best") is None else item["amine_best"]["passes"],
            "pair_203_207": item.get("pair_203_207"),
            "clearance": item.get("min_heavy_contact"),
            "pass": item.get("frozen_rule_pass"),
            "waters": len(item.get("bridging_waters", [])),
        }
        for item in reports
    ]
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
