#!/usr/bin/env python3
"""4LDE Ser203-Ser204 continuous chi1, then the other four crystals.

Window 4.7-4.9 A. An angle is rejected when either OG is closer than 2.6 A
to a backbone heavy atom of another residue. Backbone and CB stay fixed.
A clear stretch is required before protonated orciprenaline is placed.
Hydrogen bonds: H...acceptor <= 2.5 A and donor-H-acceptor >= 120 degrees.
Serines on every crystal are taken from the TM5 sequence AIASSIVSF.
"""

import importlib.util
import json
import math
from collections import defaultdict, deque
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("serine_continuous", HERE / "serine_continuous.py")
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)

OUT = HERE / "ser203_204_place.json"
SMILES = "CC(C)[NH2+]CC(O)c1cc(O)cc(O)c1"
MOTIF = "AIASSIVSF"
AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def motif_serines(path):
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
            continue
        resn = line[17:20].strip()
        rows.append((line[21].strip(), line[22:26].strip(), AA.get(resn, "X"), resn))
    by_chain = defaultdict(list)
    for chain, resi, letter, resn in rows:
        by_chain[chain].append((chain, resi, letter, resn))
    hits = []
    for chain, items in by_chain.items():
        seq = "".join(item[2] for item in items)
        start = seq.find(MOTIF)
        if start < 0:
            continue
        ser_index = {"203": start + 3, "204": start + 4, "207": start + 7}
        found = {}
        for crystal, index in ser_index.items():
            chain_i, resi, letter, resn = items[index]
            if letter != "S":
                raise SystemExit(f"{path.name} motif serine {crystal} is {resn}")
            found[crystal] = {"chain": chain_i, "file_resi": resi, "resn": resn}
        hits.append(found)
    if len(hits) != 1:
        raise SystemExit(f"{path.name} motif count {len(hits)}")
    return hits[0]


def scan_pair(path, side_ids):
    atoms = sc.load(path)
    chain = side_ids["203"]["chain"]
    wanted = {row["file_resi"]: crystal for crystal, row in side_ids.items()}
    side = {crystal: {} for crystal in side_ids}
    backbone = []
    for resn, atom_chain, resi, name, xyz in atoms:
        if atom_chain != chain:
            continue
        if name in sc.BACKBONE_NAMES:
            backbone.append((resn, resi, name, xyz))
        crystal = wanted.get(resi)
        if crystal and name in {"N", "CA", "CB", "OG"}:
            side[crystal][name] = xyz
    for crystal, row in side.items():
        row["chi1"] = sc._dihedral(row["N"], row["CA"], row["CB"], row["OG"])
        row["file_resi"] = side_ids[crystal]["file_resi"]

    def circle(crystal):
        origin = side[crystal]["CA"]
        axis = sc._sub(side[crystal]["CB"], origin)
        base = side[crystal]["chi1"]
        return [sc._rotate(side[crystal]["OG"], origin, axis, chi - base) for chi in range(-180, 180)]

    circles = {crystal: circle(crystal) for crystal in ("203", "204")}
    nearby = {}
    for crystal, points in circles.items():
        own = side[crystal]["file_resi"]
        nearby[crystal] = [
            (resn, resi, name, xyz)
            for resn, resi, name, xyz in backbone
            if resi != own and sc._dist(xyz, side[crystal]["CB"]) <= 4.5
        ]

    def closest(point, crystal):
        best_d, best_atom = 99.0, None
        for resn, resi, name, xyz in nearby[crystal]:
            d = sc._dist(point, xyz)
            if d < best_d:
                best_d, best_atom = d, [resn, resi, name, round(d, 3)]
        return best_d, best_atom

    clear = {}
    window_n = 0
    min_d, max_d = 99.0, 0.0
    slice_rows = []
    for i, og_a in enumerate(circles["203"]):
        chi_a = i - 180
        for j, og_b in enumerate(circles["204"]):
            distance = sc._dist(og_a, og_b)
            min_d = min(min_d, distance)
            max_d = max(max_d, distance)
            if not (4.7 <= distance <= 4.9):
                continue
            window_n += 1
            da, atom_a = closest(og_a, "203")
            db, atom_b = closest(og_b, "204")
            clearance = min(da, db)
            rec = {
                "chi1": {"203": chi_a, "204": j - 180},
                "OG_OG": round(distance, 3),
                "clearance": round(clearance, 3),
                "closest": {"203": atom_a, "204": atom_b},
            }
            if chi_a == -65:
                slice_rows.append(rec)
            if clearance >= 2.6:
                clear[(i, j)] = rec
    components = []
    seen = set()
    for start in clear:
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        comp = []
        while queue:
            i, j = queue.popleft()
            comp.append(clear[(i, j)])
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    nb = ((i + di) % 360, (j + dj) % 360)
                    if nb in clear and nb not in seen:
                        seen.add(nb)
                        queue.append(nb)
        components.append(comp)
    components.sort(key=len, reverse=True)
    return {
        "side": side,
        "circles": circles,
        "OG_OG_min": round(min_d, 3),
        "OG_OG_max": round(max_d, 3),
        "window_angles": window_n,
        "clear_angles": len(clear),
        "components": components,
        "slice_ser203_chi_m65": slice_rows,
        "atoms": atoms,
        "chain": chain,
    }


def choose_pose(component):
    interior = [row for row in component if 4.75 <= row["OG_OG"] <= 4.85 and row["clearance"] >= 2.8]
    pool = interior or [row for row in component if row["clearance"] >= 2.6]
    pool.sort(key=lambda row: (abs(row["OG_OG"] - 4.8), -row["clearance"]))
    return pool[0]


def grid_of(atoms):
    cells = defaultdict(list)
    for row in atoms:
        xyz = row[4]
        cells[(int(math.floor(xyz[0] / 4)), int(math.floor(xyz[1] / 4)), int(math.floor(xyz[2] / 4)))].append(row)
    return cells


def neighbors(cells, point):
    ix, iy, iz = int(math.floor(point[0] / 4)), int(math.floor(point[1] / 4)), int(math.floor(point[2] / 4))
    for a in (ix - 1, ix, ix + 1):
        for b in (iy - 1, iy, iy + 1):
            for c in (iz - 1, iz, iz + 1):
                yield from cells.get((a, b, c), ())


def sphere(center, radius):
    points = []
    for i in range(10):
        pol = math.pi * (i + 0.5) / 10
        for j in range(20):
            az = 2 * math.pi * j / 20
            points.append(
                sc._add(
                    center,
                    (
                        radius * math.sin(pol) * math.cos(az),
                        radius * math.sin(pol) * math.sin(az),
                        radius * math.cos(pol),
                    ),
                )
            )
    return points


def flip_h(carbon, oxygen, hydrogen):
    axis = sc._unit(sc._sub(oxygen, carbon))
    rel = sc._sub(hydrogen, oxygen)
    proj = sc._scale(axis, sc._dot(rel, axis))
    return sc._add(oxygen, sc._sub(sc._scale(proj, 2), rel))


def hbond(donor, hydrogen, acceptor):
    dist = sc._dist(hydrogen, acceptor)
    ba, bc = sc._sub(donor, hydrogen), sc._sub(acceptor, hydrogen)
    denom = math.sqrt(sc._dot(ba, ba) * sc._dot(bc, bc))
    cosv = max(-1.0, min(1.0, sc._dot(ba, bc) / denom))
    ang = math.degrees(math.acos(cosv))
    return dist <= 2.5 and ang >= 120.0, round(dist, 3), round(ang, 1)


def _mean(pts):
    n = len(pts)
    return (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n, sum(p[2] for p in pts) / n)


def kabsch(src, dst):
    sc_m, dc = _mean(src), _mean(dst)
    s = [sc._sub(p, sc_m) for p in src]
    d = [sc._sub(p, dc) for p in dst]

    def basis(pts):
        x = sc._unit(sc._sub(pts[1], pts[0]))
        y = sc._sub(pts[2], pts[0])
        y = sc._sub(y, sc._scale(x, sc._dot(y, x)))
        y = sc._unit(y)
        z = sc._cross(x, y)
        return (x, y, z)

    sb, db = basis(s), basis(d)
    rot = tuple(tuple(sum(db[k][i] * sb[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    rotated = tuple(sum(rot[i][k] * sc_m[k] for k in range(3)) for i in range(3))
    return rot, sc._sub(dc, rotated)


def apply_rt(rot, trans, p):
    return sc._add(
        (
            rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2],
            rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2],
            rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2],
        ),
        trans,
    )


def place_ligand(atoms, chain, og203, og204):
    protein = []
    asp = {}
    for resn, atom_chain, resi, name, xyz in atoms:
        if atom_chain != chain:
            continue
        element = "".join(ch for ch in name if ch.isalpha())[:1]
        if name in {"OG"} and resi in {None}:
            continue
        protein.append((resn, resi, name, xyz))
        if resn == "ASP" and resi == "1113" and name in {"OD1", "OD2", "CG"}:
            asp[name] = xyz
    # Replace the two serine OG coordinates with the rotated positions.
    replaced = []
    for resn, resi, name, xyz in protein:
        if resn == "SER" and name == "OG" and resi == "1203":
            xyz = og203
        elif resn == "SER" and name == "OG" and resi == "1204":
            xyz = og204
        replaced.append((resn, resi, name, xyz))
    protein = [row for row in replaced if row[2][0] != "H"]
    ods = [asp["OD1"], asp["OD2"]]
    mid = sc._scale(sc._add(ods[0], ods[1]), 0.5)
    n_tgt = sc._add(mid, sc._scale(sc._unit(sc._sub(mid, asp["CG"])), 2.9))
    cells = defaultdict(list)
    for row in protein:
        xyz = row[3]
        cells[(int(math.floor(xyz[0] / 4)), int(math.floor(xyz[1] / 4)), int(math.floor(xyz[2] / 4)))].append(row)

    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261020
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
    n_idx = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
    phenolic = []
    phenol_c, phenol_h = {}, {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetIsAromatic()]
        if not carbons:
            continue
        phenolic.append(atom.GetIdx())
        phenol_c[atom.GetIdx()] = carbons[0].GetIdx()
        phenol_h[atom.GetIdx()] = next(nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
    amine_h = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(n_idx).GetNeighbors() if nbr.GetSymbol() == "H"]
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"]
    clouds = []
    for radius in (2.7, 2.9, 3.1):
        c203, c204 = sphere(og203, radius), sphere(og204, radius)
        pairs = []
        for p1 in c203:
            for p2 in c204:
                pairs.append((sc._dist(p1, p2), p1, p2))
        clouds.append(pairs)

    best = None
    tested = 0
    for cid in ids:
        conf = mol.GetConformer(int(cid))

        def xyz(idx, conf=conf):
            pos = conf.GetAtomPosition(int(idx))
            return (pos.x, pos.y, pos.z)

        pts = {idx: xyz(idx) for idx in list(heavy) + amine_h + list(phenol_h.values())}
        oo = sc._dist(pts[phenolic[0]], pts[phenolic[1]])
        for pairs in clouds:
            near = [row for row in pairs if abs(row[0] - oo) <= 0.25]
            near.sort(key=lambda row: abs(row[0] - oo))
            for _, p1, p2 in near[:25]:
                for oa, ob, t1, t2 in (
                    (phenolic[0], phenolic[1], p1, p2),
                    (phenolic[1], phenolic[0], p1, p2),
                ):
                    rot, trans = kabsch((pts[n_idx], pts[oa], pts[ob]), (n_tgt, t1, t2))
                    mapped = tuple(apply_rt(rot, trans, p) for p in (pts[n_idx], pts[oa], pts[ob]))
                    rmsd = math.sqrt(sum(sc._dist(a, b) ** 2 for a, b in zip(mapped, (n_tgt, t1, t2))) / 3)
                    if rmsd > 0.5:
                        continue
                    tested += 1
                    placed = {idx: apply_rt(rot, trans, pts[idx]) for idx in pts}
                    clashes = []
                    closest_d = 99.0
                    for idx in heavy:
                        for resn, resi, name, rxyz in neighbors(cells, placed[idx]):
                            d = sc._dist(placed[idx], rxyz)
                            if d < closest_d:
                                closest_d = d
                            if d < 2.2 and len(clashes) < 6:
                                clashes.append([resn, resi, name, round(d, 2)])
                            elif d < 2.2:
                                clashes.append(None)
                    n_clash = len(clashes)
                    amine_hits = []
                    for h in amine_h:
                        for od in ods:
                            amine_hits.append(hbond(placed[n_idx], placed[h], od))
                    amine = max(amine_hits, key=lambda item: (item[0], -item[1], item[2]))
                    phenols = []
                    for oxy, target in ((oa, og203), (ob, og204)):
                        oxygen, hydrogen = placed[oxy], placed[phenol_h[oxy]]
                        flipped = flip_h(placed[phenol_c[oxy]], oxygen, hydrogen)
                        choices = (
                            hbond(oxygen, hydrogen, target),
                            hbond(oxygen, flipped, target),
                        )
                        phenols.append(max(choices, key=lambda item: (item[0], -item[1], item[2])))
                    n_pass = int(amine[0]) + sum(int(item[0]) for item in phenols)
                    row = {
                        "all_three_hbonds": bool(amine[0] and all(item[0] for item in phenols)),
                        "passed_contacts": n_pass,
                        "three_point_rmsd": round(rmsd, 3),
                        "clashes_under_2.2A": n_clash,
                        "closest_heavy_contact": round(closest_d, 3),
                        "amine_H_to_OD": amine[1],
                        "amine_angle": amine[2],
                        "phenol_H_to_OG": [item[1] for item in phenols],
                        "phenol_angles": [item[2] for item in phenols],
                        "clash_examples": [item for item in clashes if item is not None][:6],
                    }
                    key = (not row["all_three_hbonds"], -n_pass, n_clash, rmsd)
                    if best is None or key < best[0]:
                        best = (key, row)
    return {"trials_under_0.5A_rmsd": tested, "best": None if best is None else best[1]}


def summarize_component(comp):
    return {
        "angles": len(comp),
        "chi1_203": [min(row["chi1"]["203"] for row in comp), max(row["chi1"]["203"] for row in comp)],
        "chi1_204": [min(row["chi1"]["204"] for row in comp), max(row["chi1"]["204"] for row in comp)],
        "OG_OG": [min(row["OG_OG"] for row in comp), max(row["OG_OG"] for row in comp)],
        "clearance": [round(min(row["clearance"] for row in comp), 3), round(max(row["clearance"] for row in comp), 3)],
    }


def main():
    four = []
    files = {
        "4LDE": "active_bi167107_4LDE.pdb",
        "4LDL": "active_hbi_4LDL.pdb",
        "3SN6": "active_gs_complex_3SN6.pdb",
        "5D5A": "inactive_carazolol_5D5A.pdb",
        "2RH1": "inactive_carazolol_2RH1.pdb",
    }
    path = sc.ROOT / files["4LDE"]
    ids = motif_serines(path)
    scanned = scan_pair(path, ids)
    largest = scanned["components"][0] if scanned["components"] else []
    chosen = choose_pose(largest) if largest else None
    ligand = None
    if chosen:
        i = chosen["chi1"]["203"] + 180
        j = chosen["chi1"]["204"] + 180
        ligand = place_ligand(
            scanned["atoms"],
            scanned["chain"],
            scanned["circles"]["203"][i],
            scanned["circles"]["204"][j],
        )
    slice_clear = [row for row in scanned["slice_ser203_chi_m65"] if row["clearance"] >= 2.6]
    payload = {
        "criteria": {
            "window": [4.7, 4.9],
            "backbone_cutoff": 2.6,
            "hbond": "H...acceptor <= 2.5 A and angle >= 120 deg",
            "motif": MOTIF,
        },
        "4LDE": {
            "sequence_match": ids,
            "chi1_crystal": {k: round(scanned["side"][k]["chi1"], 1) for k in ("203", "204")},
            "OG_OG_range": [scanned["OG_OG_min"], scanned["OG_OG_max"]],
            "window_angles": scanned["window_angles"],
            "clear_angles": scanned["clear_angles"],
            "n_stretches": len(scanned["components"]),
            "largest_stretch": None if not largest else summarize_component(largest),
            "ser203_fixed_at_m65": {
                "window_angles": len(scanned["slice_ser203_chi_m65"]),
                "clear_angles": len(slice_clear),
            },
            "chosen_for_ligand": chosen,
            "ligand": ligand,
        },
    }
    for pdb_id in ("4LDL", "3SN6", "5D5A", "2RH1"):
        rec_path = sc.ROOT / files[pdb_id]
        found = motif_serines(rec_path)
        result = sc.scan_receptor(pdb_id, files[pdb_id])
        mapped = result["file_residue_numbers"]
        if {k: found[k]["file_resi"] for k in found} != mapped:
            raise SystemExit(f"{pdb_id} motif residues {found} do not match {mapped}")
        payload[pdb_id] = {
            "sequence_match": found,
            "chain": result["chain"],
            "pairs": {
                pair: {
                    "OG_OG_min": row["OG_OG_min"],
                    "window_angles": row["window_angles"],
                    "clear_angles": row["window_angles_clear_of_backbone"],
                    "best_clearance": None if row["best_in_window"] is None else row["best_in_window"]["backbone_clearance"],
                    "usable": row["usable"],
                }
                for pair, row in result["pairs"].items()
            },
            "any_usable_pair": result["any_usable_pair"],
        }
        print(pdb_id, "usable" if result["any_usable_pair"] else "none", flush=True)
    usable = ["4LDE"] if scanned["clear_angles"] else []
    usable += [pdb_id for pdb_id in ("4LDL", "3SN6", "5D5A", "2RH1") if payload[pdb_id]["any_usable_pair"]]
    payload["structures_with_a_clear_pair"] = usable
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lig = payload["4LDE"]["ligand"]
    print("4LDE window", scanned["window_angles"], "clear", scanned["clear_angles"], "stretches", len(scanned["components"]))
    print("ligand", None if lig is None else json.dumps(lig["best"]))


if __name__ == "__main__":
    main()
