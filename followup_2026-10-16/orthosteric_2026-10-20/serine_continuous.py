#!/usr/bin/env python3
"""Continuous chi1 scan for Ser203/204/207. Criteria fixed before the run.

Accept a serine pair only when its OG-OG distance is inside 4.7-4.9 A.
Reject that angle pair when either OG is closer than 2.6 A to a backbone
heavy atom (N, CA, C, O) of any residue except its own. Backbone and CB
stay fixed. Chi1 is sampled every 1 degree.
"""

import json
import math
from pathlib import Path

ROOT = Path(
    r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs"
    r"\work\local_redock\receptors"
)
OUT = Path(__file__).resolve().parent / "serine_continuous.json"
WINDOW = (4.7, 4.9)
BACKBONE_CUTOFF = 2.6
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
RECEPTORS = (
    ("4LDE", "active_bi167107_4LDE.pdb"),
    ("4LDL", "active_hbi_4LDL.pdb"),
    ("3SN6", "active_gs_complex_3SN6.pdb"),
    ("5D5A", "inactive_carazolol_5D5A.pdb"),
    ("2RH1", "inactive_carazolol_2RH1.pdb"),
)
PAIRS = (("203", "204"), ("203", "207"), ("204", "207"))


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


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    return _scale(a, 1.0 / _norm(a))


def _dist(a, b):
    return _norm(_sub(a, b))


def _dihedral(p0, p1, p2, p3):
    b1, b2, b3 = _sub(p1, p0), _sub(p2, p1), _sub(p3, p2)
    n1, n2 = _cross(b1, b2), _cross(b2, b3)
    m1 = _cross(n1, _unit(b2))
    return math.degrees(math.atan2(_dot(m1, n2), _dot(n1, n2)))


def _rotate(point, origin, axis, degrees):
    k = _unit(axis)
    v = _sub(point, origin)
    theta = math.radians(degrees)
    c, s = math.cos(theta), math.sin(theta)
    rotated = _add(
        _add(_scale(v, c), _scale(_cross(k, v), s)),
        _scale(k, _dot(k, v) * (1.0 - c)),
    )
    return _add(origin, rotated)


def load(path):
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[16] not in (" ", "A"):
            continue
        atoms.append(
            (
                line[17:20].strip(),
                line[21].strip(),
                line[22:26].strip(),
                line[12:16].strip(),
                (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            )
        )
    return atoms


def serine_map(atoms):
    found = {}
    for resn, chain, resi, name, xyz in atoms:
        if resn != "SER" or name not in {"N", "CA", "CB", "OG"}:
            continue
        crystal = {"1203": "203", "1204": "204", "1207": "207"}.get(resi, resi)
        if crystal not in {"203", "204", "207"}:
            continue
        found.setdefault((chain, crystal, resi), {})[name] = xyz
    complete = [key for key, atoms_by_name in found.items() if set(atoms_by_name) >= {"N", "CA", "CB", "OG"}]
    chains = {key[0] for key in complete}
    if len(chains) != 1:
        raise SystemExit(f"expected one serine chain, got {sorted(chains)}")
    chain = chains.pop()
    side = {}
    for key, atoms_by_name in found.items():
        if key[0] != chain or key not in complete:
            continue
        crystal, file_resi = key[1], key[2]
        side[crystal] = atoms_by_name
        side[crystal]["file_resi"] = file_resi
        side[crystal]["chi1"] = _dihedral(
            atoms_by_name["N"], atoms_by_name["CA"], atoms_by_name["CB"], atoms_by_name["OG"]
        )
    return chain, side


def scan_receptor(pdb_id, filename):
    atoms = load(ROOT / filename)
    chain, side = serine_map(atoms)
    file_ids = {crystal: side[crystal]["file_resi"] for crystal in side}
    backbone = []
    for resn, atom_chain, resi, name, xyz in atoms:
        if atom_chain == chain and name in BACKBONE_NAMES:
            backbone.append((resn, resi, name, xyz))

    def og_circle(crystal):
        origin = side[crystal]["CA"]
        axis = _sub(side[crystal]["CB"], origin)
        base = side[crystal]["chi1"]
        return [_rotate(side[crystal]["OG"], origin, axis, chi - base) for chi in range(-180, 180)]

    circles = {crystal: og_circle(crystal) for crystal in ("203", "204", "207")}
    nearby = {}
    for crystal, circle in circles.items():
        own = file_ids[crystal]
        kept = []
        for resn, resi, name, xyz in backbone:
            if resi == own:
                continue
            if _dist(xyz, side[crystal]["CB"]) <= 4.5:
                kept.append((resn, resi, name, xyz))
        nearby[crystal] = kept

    def closest(og, crystal):
        best_d, best_atom = 99.0, None
        for resn, resi, name, xyz in nearby[crystal]:
            d = _dist(og, xyz)
            if d < best_d:
                best_d, best_atom = d, [resn, resi, name, round(d, 3)]
        return best_d, best_atom

    pairs = {}
    for a, b in PAIRS:
        n_window = 0
        n_clear = 0
        best = None
        min_d, max_d = 99.0, 0.0
        for i, og_a in enumerate(circles[a]):
            for j, og_b in enumerate(circles[b]):
                distance = _dist(og_a, og_b)
                min_d = min(min_d, distance)
                max_d = max(max_d, distance)
                if not (WINDOW[0] <= distance <= WINDOW[1]):
                    continue
                n_window += 1
                da, atom_a = closest(og_a, a)
                db, atom_b = closest(og_b, b)
                clearance = min(da, db)
                if clearance >= BACKBONE_CUTOFF:
                    n_clear += 1
                if best is None or clearance > best["backbone_clearance"]:
                    best = {
                        "chi1": {a: i - 180, b: j - 180},
                        "OG_OG": round(distance, 3),
                        "backbone_clearance": round(clearance, 3),
                        "closest_backbone": {a: atom_a, b: atom_b},
                    }
        pairs[f"{a}-{b}"] = {
            "OG_OG_min": round(min_d, 3),
            "OG_OG_max": round(max_d, 3),
            "window_angles": n_window,
            "window_angles_clear_of_backbone": n_clear,
            "best_in_window": best,
            "usable": n_clear > 0,
        }
    return {
        "file": filename,
        "chain": chain,
        "file_residue_numbers": file_ids,
        "chi1_crystal": {crystal: round(side[crystal]["chi1"], 1) for crystal in ("203", "204", "207")},
        "pairs": pairs,
        "any_usable_pair": any(row["usable"] for row in pairs.values()),
    }


def main():
    results = {pdb_id: scan_receptor(pdb_id, filename) for pdb_id, filename in RECEPTORS}
    payload = {
        "criteria_written_before_scan": {
            "accept_OG_OG": list(WINDOW),
            "reject_if_OG_to_nonself_backbone_heavy_atom_below": BACKBONE_CUTOFF,
            "backbone_atoms": sorted(BACKBONE_NAMES),
            "backbone_and_CB": "fixed",
            "chi1_step_degrees": 1,
        },
        "receptors": results,
        "structures_with_a_usable_pair": [pdb_id for pdb_id, row in results.items() if row["any_usable_pair"]],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    for pdb_id, row in results.items():
        print(pdb_id, "usable" if row["any_usable_pair"] else "none")
        for pair, item in row["pairs"].items():
            best = item["best_in_window"]
            clearance = None if best is None else best["backbone_clearance"]
            print(
                " ",
                pair,
                "min",
                item["OG_OG_min"],
                "window",
                item["window_angles"],
                "clear",
                item["window_angles_clear_of_backbone"],
                "best_clearance",
                clearance,
            )


if __name__ == "__main__":
    main()
