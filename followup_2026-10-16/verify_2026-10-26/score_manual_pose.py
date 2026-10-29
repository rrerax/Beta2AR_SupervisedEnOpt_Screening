#!/usr/bin/env python3
"""Score the accepted 3SN6 orciprenaline pose. No new search.

The October 20 search accepted a phenolic-hydrogen flip that was not written
into 3sn6_orciprenaline.pdb. This script puts that flip back, checks the
frozen 2.5 A / 120 degree rule, and scores the pose with the same Vina
function as the screen. Ser203 and Ser207 are rotated to the accepted chi1
inside the screen receptor. Vina does not move any atom.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

from rdkit import Chem

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ORTHO = REPO / "followup_2026-10-16" / "orthosteric_2026-10-20"
VINA = REPO / "tools" / "vina.exe"
RECEPTOR = REPO / "data" / "processed" / "receptors" / "active_gs_complex_3SN6.pdbqt"
OBABEL = Path(r"C:\Users\victor_rax\miniconda3\envs\beta2ar-enopt-screen\Library\bin\obabel.exe")
POSE_ROOT = Path(r"D:\repositories\Enamine_GPCR_Library_plated\docking\flex_3sn6_ser203_207")
SMILES = "CC(C)[NH2+]CC(O)c1cc(O)cc(O)c1"
HBOND_MAX_A = 2.5
HBOND_MIN_DEG = 120.0
CLEARANCE_MIN_A = 2.6
BOX = (19.724, 6.018, -11.183, 24.501, 25.011, 25.452)


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, value):
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    if length == 0:
        raise SystemExit("zero vector")
    return scale(a, 1.0 / length)


def dist(a, b):
    return norm(sub(a, b))


def angle(donor, hydrogen, acceptor):
    ba, bc = sub(donor, hydrogen), sub(acceptor, hydrogen)
    denom = norm(ba) * norm(bc)
    if denom == 0:
        return 0.0
    cosv = max(-1.0, min(1.0, dot(ba, bc) / denom))
    return math.degrees(math.acos(cosv))


def hbond(donor, hydrogen, acceptor):
    distance = dist(hydrogen, acceptor)
    value = angle(donor, hydrogen, acceptor)
    return distance <= HBOND_MAX_A and value >= HBOND_MIN_DEG, round(distance, 3), round(value, 1)


def flip_h(carbon, oxygen, hydrogen):
    axis = unit(sub(oxygen, carbon))
    rel = sub(hydrogen, oxygen)
    proj = scale(axis, dot(rel, axis))
    return add(oxygen, sub(scale(proj, 2), rel))


def dihedral(a, b, c, d):
    b0 = unit(sub(a, b))
    b1 = unit(sub(c, b))
    b2 = unit(sub(d, c))
    n1 = cross(b0, b1)
    n2 = cross(b1, b2)
    y = dot(cross(n1, unit(b1)), n2)
    x = dot(n1, n2)
    return math.degrees(math.atan2(y, x))


def rodrigues(point, origin, axis, degrees_angle):
    vector = sub(point, origin)
    k = unit(axis)
    theta = math.radians(degrees_angle)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    rotated = add(
        add(scale(vector, cos_t), scale(cross(k, vector), sin_t)),
        scale(k, dot(k, vector) * (1.0 - cos_t)),
    )
    return add(origin, rotated)


def xyz_of(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def with_xyz(line, xyz):
    return f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def parse_pdb_atoms(path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "xyz": xyz_of(line),
            }
        )
    return atoms


def load_manual_coordinates():
    lines = (ORTHO / "3sn6_orciprenaline.pdb").read_text(encoding="utf-8").splitlines()
    coords = [xyz_of(line) for line in lines if line.startswith("HETATM")]
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    if len(coords) != mol.GetNumAtoms():
        raise SystemExit(f"ligand atom count {len(coords)} != {mol.GetNumAtoms()}")
    phenolic = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetIsAromatic()]
        if not carbons:
            continue
        hydrogen = next(nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
        phenolic.append((atom.GetIdx(), carbons[0].GetIdx(), hydrogen))
    if len(phenolic) != 2:
        raise SystemExit(f"phenolic oxygen count {len(phenolic)}")
    nitrogen = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
    amine_h = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(nitrogen).GetNeighbors() if nbr.GetSymbol() == "H"]
    return mol, coords, phenolic, nitrogen, amine_h


def receptor_acceptors():
    wanted = {203: {}, 207: {}}
    asp = {}
    for atom in parse_pdb_atoms(ORTHO / "3sn6_ser203_207_receptor.pdb"):
        if atom["chain"] != "R":
            continue
        if atom["resn"] == "SER" and atom["resi"] in wanted and atom["name"] == "OG":
            wanted[atom["resi"]] = atom["xyz"]
        if atom["resn"] == "ASP" and atom["resi"] == 113 and atom["name"] in {"OD1", "OD2"}:
            asp[atom["name"]] = atom["xyz"]
    if set(wanted) != {203, 207} or set(asp) != {"OD1", "OD2"}:
        raise SystemExit("manual receptor is missing OG or OD")
    return wanted, [asp["OD1"], asp["OD2"]]


def choose_flips(coords, phenolic, ogs):
    """Use the flip the October 20 search counted, assigned by nearest OG."""
    chosen = {}
    records = []
    unused = set(ogs)
    for oxygen, carbon, hydrogen in phenolic:
        target = min(unused, key=lambda resi: dist(coords[oxygen], ogs[resi]))
        unused.remove(target)
        original = coords[hydrogen]
        flipped = flip_h(coords[carbon], coords[oxygen], original)
        hits = {
            "written": hbond(coords[oxygen], original, ogs[target]),
            "flipped": hbond(coords[oxygen], flipped, ogs[target]),
        }
        use_flip = hits["flipped"][0] and (hits["flipped"][0], -hits["flipped"][1]) >= (
            hits["written"][0],
            -hits["written"][1],
        )
        chosen[hydrogen] = flipped if use_flip else original
        records.append(
            {
                "oxygen_index": oxygen,
                "serine": target,
                "written": {"passes": hits["written"][0], "distance": hits["written"][1], "angle": hits["written"][2]},
                "flipped": {"passes": hits["flipped"][0], "distance": hits["flipped"][1], "angle": hits["flipped"][2]},
                "used": "flipped" if use_flip else "written",
            }
        )
    return chosen, records


def measure(coords, phenolic, nitrogen, amine_h, ogs, ods, receptor_heavies):
    amine_hits = [hbond(coords[nitrogen], coords[hydrogen], od) for hydrogen in amine_h for od in ods]
    amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1]))
    phenols = []
    allowed = []
    for oxygen, _carbon, hydrogen in phenolic:
        target = min(ogs, key=lambda resi: dist(coords[oxygen], ogs[resi]))
        passed, distance, value = hbond(coords[oxygen], coords[hydrogen], ogs[target])
        phenols.append({"serine": target, "passes": passed, "distance": distance, "angle": value})
        allowed.append((coords[oxygen], ogs[target]))
    for od in ods:
        allowed.append((coords[nitrogen], od))
    heavies = []
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            heavies.append(coords[atom.GetIdx()])

    def excluded(left, right):
        return any(dist(left, donor) < 0.05 and dist(right, acceptor) < 0.05 for donor, acceptor in allowed)

    clearance = min(dist(left, right) for left in heavies for right in receptor_heavies if not excluded(left, right))
    return {
        "amine": {"passes": amine[0], "distance": amine[1], "angle": amine[2]},
        "phenols": phenols,
        "hydrogen_bonds": int(amine[0]) + sum(int(item["passes"]) for item in phenols),
        "min_heavy_contact": round(clearance, 3),
    }


def write_ligand(coords):
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    conf = Chem.Conformer(mol.GetNumAtoms())
    for index, xyz in enumerate(coords):
        conf.SetAtomPosition(index, xyz)
    mol.AddConformer(conf, assignId=True)
    sdf = HERE / "manual_orciprenaline.sdf"
    Chem.MolToMolFile(mol, str(sdf))
    pdbqt = HERE / "manual_orciprenaline.pdbqt"
    command = [str(OBABEL), str(sdf), "-O", str(pdbqt), "--partialcharge", "gasteiger"]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0 or not pdbqt.exists():
        raise SystemExit(completed.stderr)
    return pdbqt


def rotate_receptor(manual_ogs):
    lines = RECEPTOR.read_text(encoding="utf-8", errors="ignore").splitlines()
    parsed = []
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parsed.append(
            {
                "index": index,
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "xyz": xyz_of(line),
                "element": line.split()[-1],
            }
        )
    heavies = []
    replacements = {}
    for resi, target in manual_ogs.items():
        group = [atom for atom in parsed if atom["resn"] == "SER" and atom["chain"] == "R" and atom["resi"] == resi]
        by_name = {atom["name"]: atom for atom in group if atom["name"] in {"N", "CA", "CB", "OG"}}
        hydrogens = [atom for atom in group if atom["element"] == "HD"]
        hydroxyl = min(hydrogens, key=lambda atom: dist(atom["xyz"], by_name["OG"]["xyz"]))
        current = dihedral(by_name["N"]["xyz"], by_name["CA"]["xyz"], by_name["CB"]["xyz"], by_name["OG"]["xyz"])
        wanted = dihedral(by_name["N"]["xyz"], by_name["CA"]["xyz"], by_name["CB"]["xyz"], target)
        delta = wanted - current
        origin = by_name["CA"]["xyz"]
        axis = sub(by_name["CB"]["xyz"], origin)
        rotated_og = rodrigues(by_name["OG"]["xyz"], origin, axis, delta)
        if dist(rotated_og, target) > 0.02:
            rotated_og = rodrigues(by_name["OG"]["xyz"], origin, axis, -delta)
            hydroxyl_xyz = rodrigues(hydroxyl["xyz"], origin, axis, -delta)
            if dist(rotated_og, target) > 0.02:
                raise SystemExit(f"Ser{resi} rotation missed the manual OG by {dist(rotated_og, target):.3f} A")
        else:
            hydroxyl_xyz = rodrigues(hydroxyl["xyz"], origin, axis, delta)
        replacements[by_name["OG"]["index"]] = rotated_og
        replacements[hydroxyl["index"]] = hydroxyl_xyz
    out = []
    for index, line in enumerate(lines):
        out.append(with_xyz(line, replacements[index]) if index in replacements else line)
    path = HERE / "3sn6_manual_chi1.pdbqt"
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    for atom in parsed:
        if atom["element"] == "HD":
            continue
        xyz = replacements.get(atom["index"], atom["xyz"])
        heavies.append(xyz)
    return path, heavies


def score_only(receptor, ligand):
    command = [
        str(VINA),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--score_only",
        "--cpu",
        "1",
        "--center_x",
        str(BOX[0]),
        "--center_y",
        str(BOX[1]),
        "--center_z",
        str(BOX[2]),
        "--size_x",
        str(BOX[3]),
        "--size_y",
        str(BOX[4]),
        "--size_z",
        str(BOX[5]),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    log = HERE / f"{ligand.stem}_{receptor.stem}_score.log"
    log.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    affinity = None
    for line in completed.stdout.splitlines():
        if "Estimated Free Energy of Binding" in line:
            affinity = float(line.split(":")[1].split()[0])
    return {"returncode": completed.returncode, "affinity": affinity, "log": str(log)}


def pdbqt_atoms(path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL") and atoms:
            break
        if line.startswith("BEGIN_RES"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        kind = line.split()[-1]
        element = {"HD": "H", "OA": "O", "NA": "N", "A": "C", "C": "C", "N": "N"}.get(kind, kind[:1])
        atoms.append({"element": element, "xyz": xyz_of(line), "kind": kind})
    return atoms


def confirm_pdbqt(path, source_heavies, ogs, ods, receptor_heavies):
    atoms = pdbqt_atoms(path)
    heavies = [atom for atom in atoms if atom["element"] != "H"]
    if len(heavies) != len(source_heavies):
        raise SystemExit(f"pdbqt heavy atoms {len(heavies)} != {len(source_heavies)}")
    used = set()
    for atom in heavies:
        nearest = min(
            (dist(atom["xyz"], xyz), index)
            for index, xyz in enumerate(source_heavies)
            if index not in used
        )
        if nearest[0] > 0.05:
            raise SystemExit(f"obabel moved a heavy atom by {nearest[0]:.3f} A")
        used.add(nearest[1])
    hydrogens = [atom for atom in atoms if atom["element"] == "H"]
    donors = []
    for hydrogen in hydrogens:
        parent = min(heavies, key=lambda atom: dist(atom["xyz"], hydrogen["xyz"]))
        if dist(parent["xyz"], hydrogen["xyz"]) > 1.25 or parent["element"] not in {"N", "O"}:
            continue
        donors.append((parent, hydrogen))
    amine_hits = [
        hbond(parent["xyz"], hydrogen["xyz"], od)
        for parent, hydrogen in donors
        if parent["element"] == "N"
        for od in ods
    ]
    amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1]))
    oxygens = [atom for atom in heavies if atom["element"] == "O"]
    phenols = []
    for resi, og in ogs.items():
        best = None
        for parent, hydrogen in donors:
            if parent["element"] != "O":
                continue
            hit = hbond(parent["xyz"], hydrogen["xyz"], og)
            if best is None or (hit[0], -hit[1]) > (best[0], -best[1]):
                best = hit
        phenols.append({"serine": resi, "passes": best[0], "distance": best[1], "angle": best[2]})
    allowed = [(oxygen["xyz"], og) for oxygen in oxygens for og in ogs.values()]
    nitrogen = next(atom["xyz"] for atom in heavies if atom["element"] == "N")
    allowed.extend((nitrogen, od) for od in ods)

    def excluded(left, right):
        return any(dist(left, donor) < 0.05 and dist(right, acceptor) < 0.05 for donor, acceptor in allowed)

    clearance = min(
        dist(left["xyz"], right)
        for left in heavies
        for right in receptor_heavies
        if not excluded(left["xyz"], right)
    )
    return {
        "amine": {"passes": amine[0], "distance": amine[1], "angle": amine[2]},
        "phenols": phenols,
        "hydrogen_bonds": int(amine[0]) + sum(int(item["passes"]) for item in phenols),
        "min_heavy_contact": round(clearance, 3),
        "polar_hydrogens": len(donors),
    }


def main():
    _mol, coords, phenolic, nitrogen, amine_h = load_manual_coordinates()
    ogs, ods = receptor_acceptors()
    flips, flip_records = choose_flips(coords, phenolic, ogs)
    scored_coords = list(coords)
    for hydrogen, xyz in flips.items():
        scored_coords[hydrogen] = xyz
    receptor, heavies = rotate_receptor(ogs)
    before = measure(coords, phenolic, nitrogen, amine_h, ogs, ods, heavies)
    after = measure(scored_coords, phenolic, nitrogen, amine_h, ogs, ods, heavies)
    if after["hydrogen_bonds"] != 3 or after["min_heavy_contact"] < CLEARANCE_MIN_A:
        raise SystemExit(f"manual pose no longer passes: {after}")
    ligand = write_ligand(scored_coords)
    source_heavies = [scored_coords[atom.GetIdx()] for atom in Chem.AddHs(Chem.MolFromSmiles(SMILES)).GetAtoms() if atom.GetSymbol() != "H"]
    confirmed = confirm_pdbqt(ligand, source_heavies, ogs, ods, heavies)
    if confirmed["hydrogen_bonds"] != 3 or confirmed["min_heavy_contact"] < CLEARANCE_MIN_A:
        raise SystemExit(f"pdbqt lost the three contacts: {confirmed}")
    result = score_only(receptor, ligand)
    payload = {
        "rules": {
            "search": "none",
            "score": "vina score_only, screen box, screen receptor with Ser203/Ser207 chi1 set to the accepted rotamer",
            "hydrogens": "phenolic flip counted by the October 20 search and omitted from the written PDB",
            "hbond": "H...acceptor <= 2.5 A and donor-H-acceptor >= 120 deg",
            "clearance": "ligand-receptor heavy atoms >= 2.6 A",
        },
        "flip": flip_records,
        "written_hydrogens": before,
        "accepted_hydrogens": after,
        "pdbqt_hydrogens": confirmed,
        "vina": result,
    }
    (HERE / "manual_score.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": before, "accepted": after, "pdbqt": confirmed, "vina": result}, indent=2))


if __name__ == "__main__":
    main()
