#!/usr/bin/env python3
"""H-bond angle gate for the 4LDE Vina model-1 poses.

Written before the angles were inspected. Heavy atoms stay where Vina put
them. The docked amine is neutral (one hydrogen). The second ammonium
hydrogen is placed along the lone pair. Each phenolic hydrogen is also
tested after a 180 degree flip about the C-O bond, because that hydrogen
is a free rotor and is not part of the heavy-atom pose.

Both ligands must pass, or sampling does not start.
A contact passes when H...acceptor <= 2.5 A and donor-H-acceptor >= 120 deg.
Amine acceptor: ASP A1113 OD1/OD2.
Phenol acceptor: SER A1203, SER A1204, SER A1207 OG.
"""

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REC = ROOT / "active_bi167107_4LDE.pdb"
POSES = {
    "CHEMBL434": ROOT / "CHEMBL434_4LDE_model1.pdbqt",
    "CHEMBL776": ROOT / "CHEMBL776_4LDE_model1.pdbqt",
}
ASP = ("ASP", "1113")
SERS = {("SER", "1203"), ("SER", "1204"), ("SER", "1207")}
H_MAX = 2.5
ANGLE_MIN = 120.0


def parse_atoms(path):
    atoms = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "resi": line[22:26].strip(),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "atype": line[77:79].strip() if len(line) >= 79 else "",
                "element": "".join(ch for ch in line[12:16] if ch.isalpha())[:1],
            }
        )
    return atoms


def sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def add(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def scale(a, s):
    return tuple(a[i] * s for i in range(3))


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n < 1e-8:
        raise ValueError("zero vector")
    return scale(a, 1.0 / n)


def dist(a, b):
    return norm(sub(a, b))


def angle_at(donor, hydrogen, acceptor):
    u = sub(donor, hydrogen)
    v = sub(acceptor, hydrogen)
    cosang = max(-1.0, min(1.0, dot(u, v) / (norm(u) * norm(v))))
    return math.degrees(math.acos(cosang))


def nearest_heavy(atom, atoms, cutoff):
    hits = []
    for other in atoms:
        if other is atom or other["element"] == "H":
            continue
        d = dist(atom["xyz"], other["xyz"])
        if d <= cutoff:
            hits.append(other)
    return hits


def flip_oh(carbon, oxygen, hydrogen):
    axis = unit(sub(oxygen, carbon))
    rel = sub(hydrogen, oxygen)
    parallel = scale(axis, dot(rel, axis))
    perp = sub(rel, parallel)
    return add(oxygen, sub(parallel, perp))


def lone_pair_hydrogen(nitrogen, neighbors):
    direction = (0.0, 0.0, 0.0)
    for other in neighbors:
        direction = add(direction, unit(sub(other["xyz"], nitrogen)))
    return add(nitrogen, scale(unit(scale(direction, -1.0)), 1.01))


def best_contact(donor, hydrogens, acceptors):
    best = None
    for hydrogen in hydrogens:
        for acceptor in acceptors:
            row = {
                "h_acceptor": round(dist(hydrogen, acceptor["xyz"]), 3),
                "angle": round(angle_at(donor, hydrogen, acceptor["xyz"]), 2),
                "acceptor": f"{acceptor['resn']}{acceptor['resi']} {acceptor['name']}",
            }
            row["pass"] = row["h_acceptor"] <= H_MAX and row["angle"] >= ANGLE_MIN
            if best is None or (row["pass"], -row["h_acceptor"]) > (best["pass"], -best["h_acceptor"]):
                best = row
    return best


def ligand_report(name, atoms, receptor):
    heavy = [a for a in atoms if a["element"] != "H"]
    aromatic = [a for a in heavy if a["atype"] == "A"]
    phenols = []
    for oxygen in [a for a in heavy if a["element"] == "O"]:
        carbons = [
            a
            for a in nearest_heavy(oxygen, aromatic, 1.7)
            if a["element"] == "C"
        ]
        if not carbons:
            continue
        hydrogens = [a for a in atoms if a["element"] == "H" and dist(a["xyz"], oxygen["xyz"]) < 1.25]
        if len(hydrogens) != 1 or len(carbons) != 1:
            phenols.append({"error": f"expected one H and one C, got {len(hydrogens)} H {len(carbons)} C"})
            continue
        flipped = flip_oh(carbons[0]["xyz"], oxygen["xyz"], hydrogens[0]["xyz"])
        phenols.append({"oxygen": oxygen, "hydrogens": [hydrogens[0]["xyz"], flipped]})
    nitrogens = [a for a in heavy if a["element"] == "N"]
    if len(nitrogens) != 1:
        raise SystemExit(f"{name}: expected one nitrogen, found {len(nitrogens)}")
    nitrogen = nitrogens[0]
    neighbors = nearest_heavy(nitrogen, atoms, 1.7)
    existing_h = [a for a in atoms if a["element"] == "H" and dist(a["xyz"], nitrogen["xyz"]) < 1.25]
    if len(existing_h) != 1:
        raise SystemExit(f"{name}: expected one docked amine hydrogen, found {len(existing_h)}")
    second = lone_pair_hydrogen(nitrogen["xyz"], neighbors + existing_h)
    asp = [a for a in receptor if (a["resn"], a["resi"]) == ASP and a["name"] in {"OD1", "OD2"}]
    ser = [a for a in receptor if (a["resn"], a["resi"]) in SERS and a["name"] == "OG"]
    if len(asp) != 2 or len(ser) != 3:
        raise SystemExit(f"receptor atoms missing: asp {len(asp)} ser {len(ser)}")
    amine = best_contact(nitrogen["xyz"], [existing_h[0]["xyz"], second], asp)
    phenol_hits = []
    for phenol in phenols:
        if "error" in phenol:
            phenol_hits.append(phenol)
            continue
        phenol_hits.append(best_contact(phenol["oxygen"]["xyz"], phenol["hydrogens"], ser))
    amine_pass = bool(amine and amine["pass"])
    phenol_pass = any(hit.get("pass") for hit in phenol_hits)
    return {
        "ligand": name,
        "docked_amine": "neutral secondary amine, one hydrogen in the Vina file",
        "angle_amine": "second hydrogen placed on the lone pair; heavy atoms fixed",
        "n_phenolic_oh": len(phenols),
        "amine_contact": amine,
        "phenol_contacts": phenol_hits,
        "pass": amine_pass and phenol_pass and len(phenols) == 2,
    }


def main():
    receptor = parse_atoms(REC)
    reports = [ligand_report(name, parse_atoms(path), receptor) for name, path in POSES.items()]
    payload = {
        "receptor": "4LDE docking receptor",
        "model": 1,
        "h_acceptor_max_angstrom": H_MAX,
        "donor_h_acceptor_min_degree": ANGLE_MIN,
        "protonation_for_sampling": "both ligands protonated secondary amine, neutral phenols, net charge +1",
        "pka_computed": False,
        "ligands": reports,
        "pass_gate": all(row["pass"] for row in reports),
    }
    out = ROOT / "geometry_gate.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"pass_gate": payload["pass_gate"], "ligands": reports}, indent=2))


if __name__ == "__main__":
    main()
