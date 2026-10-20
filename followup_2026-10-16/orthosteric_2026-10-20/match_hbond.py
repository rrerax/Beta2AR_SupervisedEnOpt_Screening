#!/usr/bin/env python3
"""Place protonated orciprenaline on a backbone-clear Ser203/Ser204 pair.

The amine nitrogen is not fixed on the outward carboxylate axis. It may sit
anywhere that can still satisfy the original hydrogen-bond rule.
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
OUT = HERE / "orciprenaline_hbond_match.json"


def perp(normal):
    seed = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    x = sp.sc._unit(sp.sc._cross(normal, seed))
    y = sp.sc._cross(normal, x)
    return x, y


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
    normal = sp.sc._unit(diff)
    x, y = perp(normal)
    height = math.sqrt(max(0.0, height2))
    points = []
    for k in range(16):
        theta = 2 * math.pi * k / 16
        points.append(
            sp.sc._add(base, sp.sc._add(sp.sc._scale(x, height * math.cos(theta)), sp.sc._scale(y, height * math.sin(theta))))
        )
    return points


def prepare_receptor(scanned, og203, og204):
    protein = []
    asp = {}
    for resn, chain, resi, name, xyz in scanned["atoms"]:
        if chain != scanned["chain"] or name[:1] == "H":
            continue
        if resn == "SER" and name == "OG" and resi == "1203":
            xyz = og203
        elif resn == "SER" and name == "OG" and resi == "1204":
            xyz = og204
        protein.append((resn, resi, name, xyz))
        if resn == "ASP" and resi == "1113" and name in {"OD1", "OD2", "CG"}:
            asp[name] = xyz
    cells = defaultdict(list)
    for row in protein:
        xyz = row[3]
        cells[(int(math.floor(xyz[0] / 4)), int(math.floor(xyz[1] / 4)), int(math.floor(xyz[2] / 4)))].append(row)
    return protein, cells, [asp["OD1"], asp["OD2"]], asp["CG"]


def ligand_confs():
    mol = Chem.AddHs(Chem.MolFromSmiles(sp.SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261020
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=16, params=params))
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


def score_pose(cells, ods, item, target_n, t1, t2, og203, og204):
    pts = item["pts"]
    oa, ob = item["assign"]
    rot, trans = sp.kabsch((pts[item["n"]], pts[oa], pts[ob]), (target_n, t1, t2))
    placed = {idx: sp.apply_rt(rot, trans, pts[idx]) for idx in pts}
    n_clash = 0
    closest = 99.0
    examples = []
    for idx in item["heavy"]:
        ix, iy, iz = (int(math.floor(placed[idx][k] / 4)) for k in range(3))
        for a in (ix - 1, ix, ix + 1):
            for b in (iy - 1, iy, iy + 1):
                for c in (iz - 1, iz, iz + 1):
                    for resn, resi, name, rxyz in cells.get((a, b, c), ()):
                        dist = sp.sc._dist(placed[idx], rxyz)
                        if dist < closest:
                            closest = dist
                        if dist < 2.2:
                            n_clash += 1
                            if len(examples) < 5:
                                examples.append([resn, resi, name, round(dist, 2)])
    amine_hits = [sp.hbond(placed[item["n"]], placed[h], od) for h in item["amine_h"] for od in ods]
    amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1], hit[2]))
    phenols = []
    for oxy, target in ((oa, og203), (ob, og204)):
        oxygen, hydrogen = placed[oxy], placed[item["phenol_h"][oxy]]
        flipped = sp.flip_h(placed[item["phenol_c"][oxy]], oxygen, hydrogen)
        phenols.append(max((sp.hbond(oxygen, hydrogen, target), sp.hbond(oxygen, flipped, target)), key=lambda hit: (hit[0], -hit[1], hit[2])))
    passed = int(amine[0]) + sum(int(hit[0]) for hit in phenols)
    return {
        "all_three_hbonds": bool(amine[0] and all(hit[0] for hit in phenols)),
        "passed_contacts": passed,
        "clashes_under_2.2A": n_clash,
        "closest_heavy_contact": round(closest, 3),
        "amine_H_to_OD": amine[1],
        "amine_angle": amine[2],
        "phenol_H_to_OG": [hit[1] for hit in phenols],
        "phenol_angles": [hit[2] for hit in phenols],
        "clash_examples": examples,
        "N_to_OD": round(min(sp.sc._dist(target_n, od) for od in ods), 3),
        "O_to_OG": [round(sp.sc._dist(t1, og203), 3), round(sp.sc._dist(t2, og204), 3)],
    }


def match_one(scanned, confs, chi, og203, og204):
    _, cells, ods, cg = prepare_receptor(scanned, og203, og204)
    clouds = [sp.sphere(og203, radius) + [] for radius in (2.8,)]
    cloud_a = sp.sphere(og203, 2.85)
    cloud_b = sp.sphere(og204, 2.85)
    pairs = [(sp.sc._dist(a, b), a, b) for a in cloud_a for b in cloud_b]
    best = None
    heavy_hits = 0
    for conf in confs:
        oo = sp.sc._dist(conf["pts"][conf["phenolic"][0]], conf["pts"][conf["phenolic"][1]])
        r1 = sp.sc._dist(conf["pts"][conf["n"]], conf["pts"][conf["phenolic"][0]])
        r2 = sp.sc._dist(conf["pts"][conf["n"]], conf["pts"][conf["phenolic"][1]])
        near = [row for row in pairs if abs(row[0] - oo) <= 0.3]
        for _, p1, p2 in near:
            for oa, ob, t1, t2, na, nb in (
                (conf["phenolic"][0], conf["phenolic"][1], p1, p2, r1, r2),
                (conf["phenolic"][1], conf["phenolic"][0], p1, p2, r2, r1),
            ):
                for target_n in nitrogen_circle(t1, na, t2, nb):
                    nod = min(sp.sc._dist(target_n, od) for od in ods)
                    if not 2.6 <= nod <= 3.5 or sp.sc._dist(target_n, cg) < 2.2:
                        continue
                    heavy_hits += 1
                    item = dict(conf)
                    item["assign"] = (oa, ob)
                    row = score_pose(cells, ods, item, target_n, t1, t2, og203, og204)
                    row["chi1"] = chi
                    key = (not row["all_three_hbonds"], -row["passed_contacts"], row["clashes_under_2.2A"], row["amine_H_to_OD"])
                    if best is None or key < best[0]:
                        best = (key, row)
    return heavy_hits, None if best is None else best[1]


def main():
    rec = sp.sc.ROOT / "active_bi167107_4LDE.pdb"
    scanned = sp.scan_pair(rec, sp.motif_serines(rec))
    confs = ligand_confs()
    chosen = sp.choose_pose(scanned["components"][0])
    samples = [chosen]
    samples.extend(row for row in scanned["slice_ser203_chi_m65"] if row["clearance"] >= 2.6)
    ordered = sorted(scanned["components"][0], key=lambda row: (row["chi1"]["203"], row["chi1"]["204"]))
    step = max(1, len(ordered) // 8)
    samples.extend(ordered[::step])
    seen = set()
    unique = []
    for row in samples:
        key = (row["chi1"]["203"], row["chi1"]["204"])
        if key not in seen:
            seen.add(key)
            unique.append(row)
    results = []
    for row in unique:
        og203 = scanned["circles"]["203"][row["chi1"]["203"] + 180]
        og204 = scanned["circles"]["204"][row["chi1"]["204"] + 180]
        hits, best = match_one(scanned, confs, row["chi1"], og203, og204)
        results.append({"chi1": row["chi1"], "OG_OG": row["OG_OG"], "heavy_atom_arrangements": hits, "best": best})
        print(row["chi1"], "hits", hits, "pass", None if best is None else best["passed_contacts"], "clash", None if best is None else best["clashes_under_2.2A"], flush=True)
    passed = [row for row in results if row["best"] and row["best"]["all_three_hbonds"]]
    payload = {
        "serine_angles_tested": len(unique),
        "poses_passing_all_three_hbonds": len(passed),
        "passing": passed[:3],
        "best_attempt": min(results, key=lambda row: (not (row["best"] or {}).get("all_three_hbonds", False), -((row["best"] or {}).get("passed_contacts", -1)), (row["best"] or {}).get("clashes_under_2.2A", 999)))["best"],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("passing", len(passed))


if __name__ == "__main__":
    main()
