#!/usr/bin/env python3
"""Rotate only Ser203/204/207 (file 1203/1204/1207) on fixed 4LDE backbone.

Chi1 is the N-CA-CB-OG dihedral. Backbone and CB stay put. A pair is usable
when its OG-OG distance is 4.7-4.9 A and neither oxygen clashes with the pocket.
"""

import json
import math
from pathlib import Path

REC = Path(
    r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs"
    r"\work\local_redock\receptors\active_bi167107_4LDE.pdb"
)
OUT = Path(__file__).resolve().parent / "serine_rotamers_4LDE.json"
SERINES = ("1203", "1204", "1207")
PAIRS = (("1203", "1204"), ("1203", "1207"), ("1204", "1207"))
# Dunbrack-style serine wells: gauche+, trans, gauche-.
ROTAMERS = (("g+", 63.0), ("t", 178.0), ("g-", -65.0))
WINDOW = (4.7, 4.9)
CLASH = 2.4


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


def load():
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


def main():
    atoms = load()
    protein = [row for row in atoms if row[3] != "H"]
    side = {}
    for resi in SERINES:
        side[resi] = {
            name: xyz
            for resn, r, name, _, xyz in protein
            if resn == "SER" and r == resi and name in {"N", "CA", "CB", "OG"}
        }
        side[resi]["chi1"] = round(
            _dihedral(side[resi]["N"], side[resi]["CA"], side[resi]["CB"], side[resi]["OG"]),
            1,
        )

    def og_at(resi, chi):
        cur = side[resi]["chi1"]
        return _rotate(
            side[resi]["OG"],
            side[resi]["CA"],
            _sub(side[resi]["CB"], side[resi]["CA"]),
            chi - cur,
        )

    rotamer_ogs = {resi: {label: og_at(resi, chi) for label, chi in ROTAMERS} for resi in SERINES}

    def closest(og, resi, others):
        best = (99.0, None)
        for resn, r, name, _, xyz in protein:
            if r == resi:
                continue
            if r in SERINES and name == "OG":
                continue
            d = _dist(og, xyz)
            if d < best[0]:
                best = (d, [resn, r, name, round(d, 3)])
        for other_resi, other_og in others:
            d = _dist(og, other_og)
            if d < best[0] and d < WINDOW[0]:
                best = (d, ["SER", other_resi, "OG", round(d, 3)])
        return best

    def envelope(a, b):
        vals = []
        for i in range(36):
            oa = og_at(a, -180.0 + 10.0 * i)
            for j in range(36):
                ob = og_at(b, -180.0 + 10.0 * j)
                vals.append(_dist(oa, ob))
        return round(min(vals), 3), round(max(vals), 3)

    cb = {}
    crystal_og = {}
    reach = {}
    for a, b in PAIRS:
        cb[f"{a}-{b}"] = round(_dist(side[a]["CB"], side[b]["CB"]), 3)
        crystal_og[f"{a}-{b}"] = round(_dist(side[a]["OG"], side[b]["OG"]), 3)
        reach[f"{a}-{b}"] = envelope(a, b)

    hits = []
    nearest = []
    for la, oa in rotamer_ogs["1203"].items():
        for lb, ob in rotamer_ogs["1204"].items():
            for lc, oc in rotamer_ogs["1207"].items():
                ogs = {"1203": oa, "1204": ob, "1207": oc}
                labels = {"1203": la, "1204": lb, "1207": lc}
                clashes = {}
                for resi, og in ogs.items():
                    others = [(r, ogs[r]) for r in SERINES if r != resi]
                    clashes[resi] = closest(og, resi, others)
                distances = {f"{a}-{b}": round(_dist(ogs[a], ogs[b]), 3) for a, b in PAIRS}
                row = {
                    "rotamers": labels,
                    "OG_OG": distances,
                    "closest": {resi: clashes[resi][1] for resi in SERINES},
                }
                for pair, dist in distances.items():
                    a, b = pair.split("-")
                    clear = clashes[a][0] >= CLASH and clashes[b][0] >= CLASH
                    item = {
                        "pair": pair,
                        "OG_OG": dist,
                        "rotamers": {a: labels[a], b: labels[b]},
                        "clear_of_pocket": clear,
                        "closest": {a: clashes[a][1], b: clashes[b][1]},
                    }
                    if WINDOW[0] <= dist <= WINDOW[1] and clear:
                        hits.append(item)
                    nearest.append((abs(dist - 4.8), item))
    nearest.sort(key=lambda item: item[0])
    payload = {
        "receptor": "4LDE docking receptor, backbone and CB fixed",
        "chi1_crystal": {resi: side[resi]["chi1"] for resi in SERINES},
        "CB_CB": cb,
        "OG_OG_crystal": crystal_og,
        "OG_OG_if_chi1_spun_freely": {k: {"min": v[0], "max": v[1]} for k, v in reach.items()},
        "rotamer_chi1": {label: chi for label, chi in ROTAMERS},
        "window": list(WINDOW),
        "clash_cutoff": CLASH,
        "hits_in_window_and_clear": hits,
        "closest_rotamer_pairs": [item for _, item in nearest[:6]],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
