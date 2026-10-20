#!/usr/bin/env python3
"""One representative serine pair, then protonated orciprenaline.

The representative is fixed before the ligand is moved: among angles whose
OG-OG distance is 4.7-4.9 A and whose hydroxyls stay at least 2.6 A from
other residues' backbone, take the angle closest to 4.80 A. Backbone stays
fixed. A pose passes only when all three hydrogen bonds pass and every other
ligand-receptor heavy-atom contact is at least 2.6 A.
"""

import importlib.util
import json
import math
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("sp", HERE / "ser203_204_place.py")
sp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp)
OUT = HERE / "place_representative.json"
ASP_MOTIF = "EFWTSIDVLCVTA"
FILES = {
    "3SN6": "active_gs_complex_3SN6.pdb",
    "2RH1": "inactive_carazolol_2RH1.pdb",
    "4LDL": "active_hbi_4LDL.pdb",
    "5D5A": "inactive_carazolol_5D5A.pdb",
}
ORDER = ("3SN6", "2RH1", "4LDL", "5D5A")
PAIRS = (("203", "204"), ("203", "207"))


def find_asp(path):
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in " A":
            resn = line[17:20].strip()
            rows.append((line[21].strip(), line[22:26].strip(), sp.AA.get(resn, "X"), resn))
    by_chain = defaultdict(list)
    for chain, resi, letter, resn in rows:
        by_chain[chain].append((resi, letter, resn))
    hits = []
    for chain, items in by_chain.items():
        seq = "".join(item[1] for item in items)
        start = seq.find(ASP_MOTIF)
        if start < 0:
            continue
        resi, letter, resn = items[start + 6]
        if letter != "D":
            raise SystemExit(f"{path.name} motif aspartate is {resn}")
        hits.append({"chain": chain, "file_resi": resi, "resn": resn})
    if len(hits) != 1:
        raise SystemExit(f"{path.name} Asp motif count {len(hits)}")
    return hits[0]


def og_at(side, crystal, chi):
    row = side[crystal]
    return sp.sc._rotate(row["OG"], row["CA"], sp.sc._sub(row["CB"], row["CA"]), chi - row["chi1"])


def load_side(path, serines, asp):
    atoms = sp.sc.load(path)
    chain = serines["203"]["chain"]
    if asp["chain"] != chain:
        raise SystemExit("serine and Asp chains differ")
    wanted = {row["file_resi"]: crystal for crystal, row in serines.items()}
    side = {crystal: {} for crystal in serines}
    backbone = []
    asp_atoms = {}
    protein = []
    for resn, atom_chain, resi, name, xyz in atoms:
        if atom_chain != chain or name[:1] == "H":
            continue
        protein.append((resn, resi, name, xyz))
        if name in sp.sc.BACKBONE_NAMES:
            backbone.append((resn, resi, name, xyz))
        crystal = wanted.get(resi)
        if crystal and name in {"N", "CA", "CB", "OG"}:
            side[crystal][name] = xyz
        if resn == "ASP" and resi == asp["file_resi"] and name in {"OD1", "OD2", "CG"}:
            asp_atoms[name] = xyz
    for crystal, row in side.items():
        row["chi1"] = sp.sc._dihedral(row["N"], row["CA"], row["CB"], row["OG"])
        row["file_resi"] = serines[crystal]["file_resi"]
    return atoms, protein, backbone, side, asp_atoms


def passing_angles(backbone, side, pair):
    a, b = pair
    nearby = {}
    for crystal in pair:
        own = side[crystal]["file_resi"]
        nearby[crystal] = [
            xyz for resn, resi, name, xyz in backbone if resi != own and sp.sc._dist(xyz, side[crystal]["CB"]) <= 4.5
        ]
    clear = []
    for i in range(-180, 180):
        og_a = og_at(side, a, i)
        for j in range(-180, 180):
            og_b = og_at(side, b, j)
            distance = sp.sc._dist(og_a, og_b)
            if not 4.7 <= distance <= 4.9:
                continue
            da = min(sp.sc._dist(og_a, xyz) for xyz in nearby[a])
            db = min(sp.sc._dist(og_b, xyz) for xyz in nearby[b])
            clearance = min(da, db)
            if clearance >= 2.6:
                clear.append((i, j, distance, clearance))
    return clear


def representative(clear):
    return min(clear, key=lambda row: (abs(row[2] - 4.8), -row[3]))


def perp(normal):
    seed = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    x = sp.sc._unit(sp.sc._cross(normal, seed))
    return x, sp.sc._cross(normal, x)


def nitrogen_circle(o1, r1, o2, r2):
    diff = sp.sc._sub(o2, o1)
    dist = sp.sc._norm(diff)
    if dist < 1e-4 or dist > r1 + r2 or dist < abs(r1 - r2):
        return []
    along = (r1 * r1 - r2 * r2 + dist * dist) / (2 * dist)
    height2 = r1 * r1 - along * along
    if height2 < -1e-6:
        return []
    base = sp.sc._add(o1, sp.sc._scale(sp.sc._unit(diff), along))
    x, y = perp(sp.sc._unit(diff))
    height = math.sqrt(max(0.0, height2))
    return [
        sp.sc._add(
            base,
            sp.sc._add(sp.sc._scale(x, height * math.cos(theta)), sp.sc._scale(y, height * math.sin(theta))),
        )
        for theta in (2 * math.pi * k / 12 for k in range(12))
    ]


def sphere(center, radius):
    points = []
    for i in range(8):
        pol = math.pi * (i + 0.5) / 8
        for j in range(16):
            az = 2 * math.pi * j / 16
            points.append(
                sp.sc._add(
                    center,
                    (
                        radius * math.sin(pol) * math.cos(az),
                        radius * math.sin(pol) * math.sin(az),
                        radius * math.cos(pol),
                    ),
                )
            )
    return points


def ligand_confs():
    mol = Chem.AddHs(Chem.MolFromSmiles(sp.SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261020
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=12, params=params))
    n_idx = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
    phenolic, phenol_c, phenol_h = [], {}, {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetIsAromatic()]
        if carbons:
            phenolic.append(atom.GetIdx())
            phenol_c[atom.GetIdx()] = carbons[0].GetIdx()
            phenol_h[atom.GetIdx()] = next(nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
    amine_h = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(n_idx).GetNeighbors() if nbr.GetSymbol() == "H"]
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"]
    confs = []
    for cid in ids:
        conf = mol.GetConformer(int(cid))

        def xyz(idx, conf=conf):
            pos = conf.GetAtomPosition(int(idx))
            return (pos.x, pos.y, pos.z)

        needed = set(heavy + amine_h + list(phenol_h.values()) + list(phenol_c.values()))
        confs.append(
            {
                "n": n_idx,
                "phenolic": phenolic,
                "phenol_c": phenol_c,
                "phenol_h": phenol_h,
                "amine_h": amine_h,
                "heavy": heavy,
                "pts": {idx: xyz(idx) for idx in needed},
            }
        )
    return confs


def receptor_with_ogs(protein, side, pair, ogs):
    replaced = []
    for resn, resi, name, xyz in protein:
        for crystal in pair:
            if resn == "SER" and name == "OG" and resi == side[crystal]["file_resi"]:
                xyz = ogs[crystal]
        replaced.append((resn, resi, name, xyz))
    protein = replaced
    cells = defaultdict(list)
    for row in replaced:
        xyz = row[3]
        cells[(int(math.floor(xyz[0] / 4)), int(math.floor(xyz[1] / 4)), int(math.floor(xyz[2] / 4)))].append(row)
    return protein, cells


def brute_clearance(protein, placed, heavy, allowed):
    worst = 99.0
    worst_atom = None
    for idx in heavy:
        for resn, resi, name, xyz in protein:
            if (idx, resi, name) in allowed:
                continue
            dist = sp.sc._dist(placed[idx], xyz)
            if dist < worst:
                worst = dist
                worst_atom = [resn, resi, name, round(dist, 3)]
    return worst, worst_atom


def other_clearance(cells, placed, heavy, allowed):
    worst = 99.0
    worst_atom = None
    for idx in heavy:
        ix, iy, iz = (int(math.floor(placed[idx][k] / 4)) for k in range(3))
        for a in (ix - 1, ix, ix + 1):
            for b in (iy - 1, iy, iy + 1):
                for c in (iz - 1, iz, iz + 1):
                    for resn, resi, name, xyz in cells.get((a, b, c), ()):
                        if (idx, resi, name) in allowed:
                            continue
                        dist = sp.sc._dist(placed[idx], xyz)
                        if dist < worst:
                            worst = dist
                            worst_atom = [resn, resi, name, round(dist, 3)]
    return worst, worst_atom


def search(confs, protein, cells, ods, cg, og_a, og_b, resi_a, resi_b, asp_resi):
    clouds = {radius: (sphere(og_a, radius), sphere(og_b, radius)) for radius in (2.8, 3.0, 3.2)}
    best = None
    partial = None
    tested = 0
    for radius, (cloud_a, cloud_b) in clouds.items():
        pairs = [(sp.sc._dist(p1, p2), p1, p2) for p1 in cloud_a for p2 in cloud_b]
        for conf in confs:
            pts = conf["pts"]
            oo = sp.sc._dist(pts[conf["phenolic"][0]], pts[conf["phenolic"][1]])
            radii = {
                conf["phenolic"][0]: sp.sc._dist(pts[conf["n"]], pts[conf["phenolic"][0]]),
                conf["phenolic"][1]: sp.sc._dist(pts[conf["n"]], pts[conf["phenolic"][1]]),
            }
            near = [row for row in pairs if abs(row[0] - oo) <= 0.3]
            for _, p1, p2 in near:
                for oa, ob, t1, t2 in (
                    (conf["phenolic"][0], conf["phenolic"][1], p1, p2),
                    (conf["phenolic"][1], conf["phenolic"][0], p1, p2),
                ):
                    for target_n in nitrogen_circle(t1, radii[oa], t2, radii[ob]):
                        if sp.sc._dist(target_n, cg) < 2.4:
                            continue
                        if not 2.6 <= min(sp.sc._dist(target_n, od) for od in ods) <= 3.5:
                            continue
                        tested += 1
                        rot, trans = sp.kabsch((pts[conf["n"]], pts[oa], pts[ob]), (target_n, t1, t2))
                        placed = {idx: sp.apply_rt(rot, trans, point) for idx, point in pts.items()}
                        amine_hits = [sp.hbond(placed[conf["n"]], placed[h], od) for h in conf["amine_h"] for od in ods]
                        amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1], hit[2]))
                        phenols = []
                        for oxy, target in ((oa, og_a), (ob, og_b)):
                            oxygen = placed[oxy]
                            hydrogen = placed[conf["phenol_h"][oxy]]
                            flipped = sp.flip_h(placed[conf["phenol_c"][oxy]], oxygen, hydrogen)
                            phenols.append(
                                max(
                                    (sp.hbond(oxygen, hydrogen, target), sp.hbond(oxygen, flipped, target)),
                                    key=lambda hit: (hit[0], -hit[1], hit[2]),
                                )
                            )
                        n_pass = int(amine[0]) + sum(int(hit[0]) for hit in phenols)
                        phenols_ok = sum(int(hit[0]) for hit in phenols)
                        partial_row = {
                            "passes": False,
                            "passed_contacts": n_pass,
                            "phenols_passed": phenols_ok,
                            "amine_H_to_OD": amine[1],
                            "amine_angle": amine[2],
                            "phenol_H_to_OG": [hit[1] for hit in phenols],
                            "phenol_angles": [hit[2] for hit in phenols],
                            "target_radius": radius,
                        }
                        new_key = (phenols_ok, n_pass, -partial_row["amine_H_to_OD"])
                        old_key = None if partial is None else (partial["phenols_passed"], partial["passed_contacts"], -partial["amine_H_to_OD"])
                        if old_key is None or new_key > old_key:
                            partial = partial_row
                        if n_pass < 3:
                            continue
                        allowed = {
                            (conf["n"], asp_resi, "OD1"),
                            (conf["n"], asp_resi, "OD2"),
                            (oa, resi_a, "OG"),
                            (ob, resi_b, "OG"),
                        }
                        clearance, atom = other_clearance(cells, placed, conf["heavy"], allowed)
                        row = {
                            "passes": clearance >= 2.6,
                            "other_heavy_min": round(clearance, 3),
                            "closest_other": atom,
                            "amine_H_to_OD": amine[1],
                            "amine_angle": amine[2],
                            "phenol_H_to_OG": [hit[1] for hit in phenols],
                            "phenol_angles": [hit[2] for hit in phenols],
                            "target_radius": radius,
                        }
                        if best is None or clearance > best["other_heavy_min"]:
                            best = row
                        if row["passes"]:
                            brute, brute_atom = brute_clearance(protein, placed, conf["heavy"], allowed)
                            row["other_heavy_min"] = round(brute, 3)
                            row["closest_other"] = brute_atom
                            row["passes"] = brute >= 2.6
                            if row["passes"]:
                                best = row
                                best["arrangements_tested"] = tested
                                return best
    if best:
        best["arrangements_tested"] = tested
        return best
    if partial:
        partial["arrangements_tested"] = tested
    return partial


def main():
    confs = ligand_confs()
    records = []
    winner = None
    for pdb_id in ORDER:
        path = sp.sc.ROOT / FILES[pdb_id]
        serines = sp.motif_serines(path)
        asp = find_asp(path)
        _, protein, backbone, side, asp_atoms = load_side(path, serines, asp)
        ods = [asp_atoms["OD1"], asp_atoms["OD2"]]
        for pair in PAIRS:
            clear = passing_angles(backbone, side, pair)
            chosen = representative(clear)
            chi = {pair[0]: chosen[0], pair[1]: chosen[1]}
            ogs = {pair[0]: og_at(side, pair[0], chosen[0]), pair[1]: og_at(side, pair[1], chosen[1])}
            print(pdb_id, pair, "representative", chi, "OG", round(chosen[2], 3), "backbone", round(chosen[3], 3), "n", len(clear), flush=True)
            posed_protein, cells = receptor_with_ogs(protein, side, pair, ogs)
            found = search(
                confs,
                posed_protein,
                cells,
                ods,
                asp_atoms["CG"],
                ogs[pair[0]],
                ogs[pair[1]],
                side[pair[0]]["file_resi"],
                side[pair[1]]["file_resi"],
                asp["file_resi"],
            )
            record = {
                "structure": pdb_id,
                "pair": list(pair),
                "sequence": {"serines": serines, "asp": asp},
                "clear_angles": len(clear),
                "representative_chi1": chi,
                "OG_OG": round(chosen[2], 3),
                "backbone_clearance": round(chosen[3], 3),
                "ligand": found,
            }
            records.append(record)
            print(" ", None if found is None else found, flush=True)
            OUT.write_text(json.dumps({"attempts": records}, indent=2), encoding="utf-8")
            if found and found["passes"]:
                winner = record
                break
        if winner:
            break
    payload = {"winner": winner, "attempts": records}
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("done", None if winner is None else winner["structure"], flush=True)


if __name__ == "__main__":
    main()
