#!/usr/bin/env python3
"""Place protonated orciprenaline onto the 4LDE orthosteric contacts.

Targets, fixed side chains: amine N toward ASP A1113, the two phenolic
oxygens toward a pair of SER A1203 / A1204 / A1207. No side-chain motion.
"""

import json
import math
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

REC = Path(
    r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs"
    r"\work\local_redock\receptors\active_bi167107_4LDE.pdb"
)
OUT = Path(__file__).resolve().parent / "orciprenaline_pocket_placement.json"
SMILES = "CC(C)[NH2+]CC(O)c1cc(O)cc(O)c1"


def load_receptor():
    atoms = []
    for line in REC.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20].strip()
        resi = line[22:26].strip()
        name = line[12:16].strip()
        element = (line[76:78].strip() or name[:1]).upper()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        atoms.append((resn, resi, name, element, xyz))
    return atoms


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _unit(a):
    n = math.sqrt(_dot(a, a))
    return _scale(a, 1.0 / n)


def _mean(pts):
    n = len(pts)
    return (
        sum(p[0] for p in pts) / n,
        sum(p[1] for p in pts) / n,
        sum(p[2] for p in pts) / n,
    )


def kabsch(src, dst):
    """Rigid map of three points, using only plain arithmetic."""
    sc, dc = _mean(src), _mean(dst)
    s = [_sub(p, sc) for p in src]
    d = [_sub(p, dc) for p in dst]

    def basis(pts):
        x = _unit(_sub(pts[1], pts[0]))
        y = _sub(pts[2], pts[0])
        y = _sub(y, _scale(x, _dot(y, x)))
        y = _unit(y)
        z = _cross(x, y)
        return (x, y, z)

    sb, db = basis(s), basis(d)
    # Columns of each basis are x, y, z. rot = dest @ source.T.
    rot = tuple(
        tuple(sum(db[k][i] * sb[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )
    rotated = tuple(sum(rot[i][k] * sc[k] for k in range(3)) for i in range(3))
    trans = _sub(dc, rotated)
    return rot, trans


def _apply(rot, trans, p):
    return _add(
        (
            rot[0][0] * p[0] + rot[0][1] * p[1] + rot[0][2] * p[2],
            rot[1][0] * p[0] + rot[1][1] * p[1] + rot[1][2] * p[2],
            rot[2][0] * p[0] + rot[2][1] * p[1] + rot[2][2] * p[2],
        ),
        trans,
    )


def _dist(a, b):
    return math.sqrt(_dot(_sub(a, b), _sub(a, b)))


def _sphere(center, radius, n_az=16, n_pol=8):
    points = []
    for i in range(n_pol):
        pol = math.pi * (i + 0.5) / n_pol
        for j in range(n_az):
            az = 2 * math.pi * j / n_az
            points.append(
                _add(
                    center,
                    (
                        radius * math.sin(pol) * math.cos(az),
                        radius * math.sin(pol) * math.sin(az),
                        radius * math.cos(pol),
                    ),
                )
            )
    return points


def _angle(a, b, c):
    ba, bc = _sub(a, b), _sub(c, b)
    denom = math.sqrt(_dot(ba, ba) * _dot(bc, bc))
    cosv = max(-1.0, min(1.0, _dot(ba, bc) / denom))
    return math.degrees(math.acos(cosv))


def main():
    atoms = load_receptor()
    protein = [
        (resn, resi, name, element, xyz)
        for resn, resi, name, element, xyz in atoms
        if element != "H" and resn not in {"HOH", "WAT"} and name not in {"P0G", "BI167107"}
    ]
    # Drop the co-crystallized agonist if it is still in the file.
    protein = [row for row in protein if row[0] not in {"P0G", "LIG", "UNL"}]
    ogs = {resi: xyz for resn, resi, name, _, xyz in protein if resn == "SER" and name == "OG"}
    cbs = {resi: xyz for resn, resi, name, _, xyz in protein if resn == "SER" and name == "CB"}
    ods = [xyz for resn, resi, name, _, xyz in protein if resn == "ASP" and resi == "1113" and name in {"OD1", "OD2"}]
    cg = next(xyz for resn, resi, name, _, xyz in protein if resn == "ASP" and resi == "1113" and name == "CG")
    mid = _scale(_add(ods[0], ods[1]), 0.5)
    n_tgt = _add(mid, _scale(_unit(_sub(mid, cg)), 2.9))

    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261019
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=30, params=params))
    if not ids:
        raise SystemExit("no conformers")
    n_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N")
    phenolic = []
    phenol_carbon = {}
    phenol_h = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        carbons = [n for n in atom.GetNeighbors() if n.GetIsAromatic()]
        if carbons:
            phenolic.append(atom.GetIdx())
            phenol_carbon[atom.GetIdx()] = carbons[0].GetIdx()
            hydrogens = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "H"]
            phenol_h[atom.GetIdx()] = hydrogens[0]
    amine_h = [n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors() if n.GetSymbol() == "H"]
    amine_c = [n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors() if n.GetSymbol() == "C"]
    all_idx = [a.GetIdx() for a in mol.GetAtoms()]
    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H"]

    confs = []
    oo = []
    for cid in ids:
        conf = mol.GetConformer(int(cid))

        def xyz(idx, conf=conf):
            pos = conf.GetAtomPosition(int(idx))
            return (pos.x, pos.y, pos.z)

        pts = {i: xyz(i) for i in all_idx}
        confs.append(pts)
        oo.append(_dist(pts[phenolic[0]], pts[phenolic[1]]))

    grid = defaultdict(list)
    for row in protein:
        xyz = row[4]
        grid[(int(xyz[0] // 4), int(xyz[1] // 4), int(xyz[2] // 4))].append(row)

    def neighbors(point):
        ix, iy, iz = int(point[0] // 4), int(point[1] // 4), int(point[2] // 4)
        for a in (ix - 1, ix, ix + 1):
            for b in (iy - 1, iy, iy + 1):
                for c in (iz - 1, iz, iz + 1):
                    yield from grid.get((a, b, c), ())

    def clear_of_protein(point, ignore):
        best = 99.0
        for resn, resi, name, _, xyz in neighbors(point):
            if (resi, name) in ignore:
                continue
            best = min(best, _dist(point, xyz))
        return best

    target_oo = sorted(oo)[len(oo) // 2]
    best = None
    best_fit = None
    by_pair = {}
    trials = 0
    fitted = 0
    for radius in (2.7, 2.9, 3.1):
        clouds = {}
        for resi in ("1203", "1204", "1207"):
            clouds[resi] = [
                p
                for p in _sphere(ogs[resi], radius)
                if clear_of_protein(p, {(resi, "OG"), (resi, "CB")}) >= 2.3
            ]
        for pair in (("1203", "1204"), ("1203", "1207"), ("1204", "1207")):
            candidates = []
            for p1 in clouds[pair[0]]:
                for p2 in clouds[pair[1]]:
                    span = _dist(p1, p2)
                    if abs(span - target_oo) <= 0.45:
                        candidates.append((abs(span - target_oo), p1, p2))
            candidates.sort(key=lambda item: item[0])
            for _, p1, p2 in candidates[:12]:
                span = _dist(p1, p2)
                for conf_i, pts in enumerate(confs):
                    if abs(span - oo[conf_i]) > 0.45:
                        continue
                    for oa, ob in (phenolic, phenolic[::-1]):
                        trials += 1
                        src = (pts[n_idx], pts[oa], pts[ob])
                        dst = (n_tgt, p1, p2)
                        rot, trans = kabsch(src, dst)
                        mapped = tuple(_apply(rot, trans, p) for p in src)
                        rmsd = math.sqrt(sum(_dist(a, b) ** 2 for a, b in zip(mapped, dst)) / 3)
                        placed = {i: _apply(rot, trans, pts[i]) for i in all_idx}
                        clashes = []
                        closest = 99.0
                        n_clash = 0
                        for i in heavy_idx:
                            for resn, resi, name, _, xyz in neighbors(placed[i]):
                                d = _dist(placed[i], xyz)
                                if d < closest:
                                    closest = d
                                if d < 2.2:
                                    n_clash += 1
                                    if len(clashes) < 6:
                                        clashes.append([resn, resi, name, round(d, 2)])
                        n_od = min(_dist(placed[n_idx], ods[0]), _dist(placed[n_idx], ods[1]))
                        og_d = [
                            round(_dist(placed[oa], ogs[pair[0]]), 3),
                            round(_dist(placed[ob], ogs[pair[1]]), 3),
                        ]
                        phenol_angles = []
                        phenol_h_dist = []
                        for oxy, target in ((oa, ogs[pair[0]]), (ob, ogs[pair[1]])):
                            oxygen = placed[oxy]
                            hydrogen = placed[phenol_h[oxy]]
                            axis = _unit(_sub(oxygen, placed[phenol_carbon[oxy]]))
                            rel = _sub(hydrogen, oxygen)
                            proj = _scale(axis, _dot(rel, axis))
                            flipped = _add(oxygen, _sub(_scale(proj, 2), rel))
                            choices = (
                                (_angle(oxygen, hydrogen, target), _dist(hydrogen, target)),
                                (_angle(oxygen, flipped, target), _dist(flipped, target)),
                            )
                            chosen = max(choices, key=lambda item: (item[0] >= 120 and item[1] <= 2.5, item[0]))
                            phenol_angles.append(round(chosen[0], 1))
                            phenol_h_dist.append(round(chosen[1], 3))
                        amine_choices = []
                        for h in amine_h:
                            for od in ods:
                                amine_choices.append((_angle(placed[n_idx], placed[h], od), _dist(placed[h], od)))
                        amine_best = max(amine_choices, key=lambda item: (item[0] >= 120 and item[1] <= 2.5, item[0]))
                        nearest_od = min(ods, key=lambda od: _dist(placed[n_idx], od))
                        carbon_angles = [
                            round(_angle(placed[c], placed[n_idx], nearest_od), 1) for c in amine_c
                        ]
                        pointed = _add(placed[n_idx], _scale(_unit(_sub(nearest_od, placed[n_idx])), 1.01))
                        pointed_h = round(_dist(pointed, nearest_od), 3)
                        row = {
                            "clashes_under_2.2A": n_clash,
                            "closest_receptor_contact": round(closest, 3),
                            "three_point_rmsd": round(rmsd, 3),
                            "amine_N_to_nearest_OD": round(n_od, 3),
                            "phenol_O_to_assigned_OG": og_d,
                            "serines": list(pair),
                            "target_radius": radius,
                            "phenol_donor_angles": phenol_angles,
                            "phenol_H_to_OG": phenol_h_dist,
                            "amine_best_angle": round(amine_best[0], 1),
                            "amine_H_to_OD": round(amine_best[1], 3),
                            "amine_C_N_OD_angles": carbon_angles,
                            "amine_H_to_OD_if_pointed": pointed_h,
                            "conformer": conf_i,
                            "clash_examples": clashes,
                        }
                        key = (n_clash, rmsd, -closest)
                        if best is None or key < best[0]:
                            best = (key, row)
                        fit_key = (rmsd, n_clash, -closest)
                        if best_fit is None or fit_key < best_fit[0]:
                            best_fit = (fit_key, row)
                        if rmsd <= 0.5:
                            fitted += 1
                        pair_key = "-".join(pair)
                        pair_rank = (0 if rmsd <= 0.5 else 1, n_clash, rmsd)
                        if pair_key not in by_pair or pair_rank < by_pair[pair_key][0]:
                            by_pair[pair_key] = (pair_rank, row)
    payload = {
        "ligand": "CHEMBL776 orciprenaline",
        "smiles": SMILES,
        "receptor": "4LDE, side chains fixed, crystal ligand excluded",
        "phenolic_O_O_range": [round(min(oo), 3), round(max(oo), 3)],
        "ser_OG_distances": {
            "1203-1204": round(_dist(ogs["1203"], ogs["1204"]), 3),
            "1203-1207": round(_dist(ogs["1203"], ogs["1207"]), 3),
            "1204-1207": round(_dist(ogs["1204"], ogs["1207"]), 3),
        },
        "amine_site_to_OD": [round(_dist(n_tgt, ods[0]), 3), round(_dist(n_tgt, ods[1]), 3)],
        "trials_within_triangle": trials,
        "trials_with_rmsd_at_most_0.5": fitted,
        "least_clashes": None if best is None else best[1],
        "closest_triangle": None if best_fit is None else best_fit[1],
        "best_by_serine_pair": {key: value[1] for key, value in by_pair.items()},
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
