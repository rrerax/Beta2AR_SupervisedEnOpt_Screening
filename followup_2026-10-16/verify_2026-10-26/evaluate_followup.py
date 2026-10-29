#!/usr/bin/env python3
"""Paired scores and crystal-pattern distances. No new search.

The manual three-contact pose scored -4.376 kcal/mol. This script scores
the already docked flex modes with the same score_only call, after putting
each mode's Ser203 and Ser207 into the screen receptor. It also records
the phenolic oxygen-OG distances. A mode reaches the deposited pattern
only when two different phenolic oxygens are at least as close as the
longer 4LDO contact, 3.532 A, one to Ser203 and one to Ser207. That
comparison does not replace the 2.5 A rule.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import score_manual_pose as manual

HERE = Path(__file__).resolve().parent
POSE_ROOT = Path(r"D:\repositories\Enamine_GPCR_Library_plated\docking\flex_3sn6_ser203_207")
DEPOSITED_OG_MAX = 3.532
TYPE_ELEMENT = {"HD": "H", "HS": "H", "OA": "O", "OS": "O", "NA": "N", "N": "N", "A": "C", "C": "C", "SA": "S"}


def element(kind):
    return TYPE_ELEMENT.get(kind, "H" if kind.startswith("H") else kind[:1])


def models_of(path):
    current = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL"):
            yield current
            current = []
        elif current:
            current.append(line)


def atoms_of(lines):
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        kind = line.split()[-1]
        atoms.append(
            {
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "xyz": manual.xyz_of(line),
                "element": element(kind),
                "kind": kind,
                "charge": float(line.split()[-2]),
            }
        )
    return atoms


def affinity_of(lines):
    for line in lines:
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def analyze(lines):
    ligand_lines = []
    for line in lines:
        if line.startswith("BEGIN_RES"):
            break
        ligand_lines.append(line)
    ligand = atoms_of(ligand_lines)
    ser = {}
    current = None
    for line in lines:
        if line.startswith("BEGIN_RES"):
            current = int(line.split()[-1])
            ser[current] = {}
        elif current is not None and line.startswith(("ATOM", "HETATM")):
            atom = atoms_of([line])[0]
            ser[current][atom["name"]] = atom["xyz"]
    heavies = [atom for atom in ligand if atom["element"] != "H"]
    adj = {id(atom): [] for atom in heavies}
    for i, left in enumerate(heavies):
        for right in heavies[i + 1 :]:
            if manual.dist(left["xyz"], right["xyz"]) <= 1.7:
                adj[id(left)].append(right)
                adj[id(right)].append(left)
    phenols = []
    for atom in heavies:
        if atom["element"] != "O" or len(adj[id(atom)]) != 1:
            continue
        carbon = adj[id(atom)][0]
        if carbon["kind"] != "A" and carbon["element"] != "C":
            continue
        if carbon["kind"] == "C":
            continue
        hydrogens = [
            other
            for other in ligand
            if other["element"] == "H" and manual.dist(other["xyz"], atom["xyz"]) <= 1.25
        ]
        phenols.append((atom, hydrogens))
    nitrogens = [atom for atom in heavies if atom["element"] == "N"]
    amine_h = [
        other
        for nitrogen in nitrogens
        for other in ligand
        if other["element"] == "H" and manual.dist(other["xyz"], nitrogen["xyz"]) <= 1.25
    ]
    ods = None
    return {
        "ligand_lines": ligand_lines,
        "phenols": phenols,
        "nitrogens": nitrogens,
        "amine_h": amine_h,
        "ser": ser,
        "heavies": heavies,
        "ligand": ligand,
    }


def contact_table(parsed, ods):
    amine_hits = [
        manual.hbond(nitrogen["xyz"], hydrogen["xyz"], od)
        for nitrogen in parsed["nitrogens"]
        for hydrogen in parsed["amine_h"]
        for od in ods
    ]
    amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1])) if amine_hits else (False, None, None)
    phenols = []
    for atom, hydrogens in parsed["phenols"]:
        for resi in (203, 207):
            og = parsed["ser"][resi]["OG"]
            if hydrogens:
                best = max((manual.hbond(atom["xyz"], hydrogen["xyz"], og) for hydrogen in hydrogens), key=lambda hit: (hit[0], -hit[1]))
            else:
                best = (False, None, None)
            phenols.append(
                {
                    "oxygen": atom["name"],
                    "serine": resi,
                    "oxygen_og": round(manual.dist(atom["xyz"], og), 3),
                    "passes": best[0],
                    "h_distance": best[1],
                    "angle": best[2],
                }
            )
    pair = False
    best_pair = None
    n = len(parsed["phenols"])
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            left = manual.dist(parsed["phenols"][i][0]["xyz"], parsed["ser"][203]["OG"])
            right = manual.dist(parsed["phenols"][j][0]["xyz"], parsed["ser"][207]["OG"])
            def best_h(oxygen_atom, hydrogens, og):
                if not hydrogens:
                    return (False, None, None)
                return max(
                    (manual.hbond(oxygen_atom["xyz"], hydrogen["xyz"], og) for hydrogen in hydrogens),
                    key=lambda hit: (hit[0], -hit[1]),
                )

            if best_pair is None or left + right < best_pair[0]:
                h203 = best_h(parsed["phenols"][i][0], parsed["phenols"][i][1], parsed["ser"][203]["OG"])
                h207 = best_h(parsed["phenols"][j][0], parsed["phenols"][j][1], parsed["ser"][207]["OG"])
                best_pair = (
                    left + right,
                    round(left, 3),
                    round(right, 3),
                    parsed["phenols"][i][0]["name"],
                    parsed["phenols"][j][0]["name"],
                    {"passes": h203[0], "distance": h203[1], "angle": h203[2]},
                    {"passes": h207[0], "distance": h207[1], "angle": h207[2]},
                )
            if left <= DEPOSITED_OG_MAX and right <= DEPOSITED_OG_MAX:
                pair = True
    return {
        "amine_passes": bool(amine[0]),
        "amine_distance": amine[1],
        "amine_angle": amine[2],
        "phenols": phenols,
        "closest_pair_og": None
        if best_pair is None
        else {
            "sum": round(best_pair[0], 3),
            "ser203": best_pair[1],
            "ser207": best_pair[2],
            "atom203": best_pair[3],
            "atom207": best_pair[4],
            "h203": best_pair[5],
            "h207": best_pair[6],
        },
        "reaches_deposited_og": pair,
    }


def screen_ods():
    ods = []
    for line in manual.RECEPTOR.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[17:20] == "ASP" and line[21] == "R" and int(line[22:26]) == 113 and line[12:16].strip() in {"OD1", "OD2"}:
            ods.append(manual.xyz_of(line))
    if len(ods) != 2:
        raise SystemExit("Asp113 oxygens missing")
    return ods


def patch_receptor(ser, destination):
    lines = manual.RECEPTOR.read_text(encoding="utf-8", errors="ignore").splitlines()
    parsed = []
    for index, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            parsed.append(
                {
                    "index": index,
                    "name": line[12:16].strip(),
                    "resn": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "resi": int(line[22:26]),
                    "xyz": manual.xyz_of(line),
                    "element": element(line.split()[-1]),
                }
            )
    replacements = {}
    for resi, coords in ser.items():
        group = [atom for atom in parsed if atom["resn"] == "SER" and atom["chain"] == "R" and atom["resi"] == resi]
        og = next(atom for atom in group if atom["name"] == "OG")
        hydroxyl = min((atom for atom in group if atom["element"] == "H"), key=lambda atom: manual.dist(atom["xyz"], og["xyz"]))
        replacements[og["index"]] = coords["OG"]
        replacements[hydroxyl["index"]] = coords.get("H", hydroxyl["xyz"])
    destination.write_text(
        "\n".join(manual.with_xyz(line, replacements[index]) if index in replacements else line for index, line in enumerate(lines))
        + "\n",
        encoding="utf-8",
    )
    return destination


def main():
    log = (HERE / "manual_orciprenaline_3sn6_manual_chi1_score.log").read_text(encoding="utf-8")
    manual_affinity = None
    for line in log.splitlines():
        if "Estimated Free Energy of Binding" in line:
            manual_affinity = float(line.split(":")[1].split()[0])
    manual_json = json.loads((HERE / "manual_score.json").read_text(encoding="utf-8"))
    manual_json["vina"]["affinity"] = manual_affinity
    (HERE / "manual_score.json").write_text(json.dumps(manual_json, indent=2) + "\n", encoding="utf-8")

    ods = screen_ods()
    summary = {"manual_score_only": manual_affinity, "deposited_og_max_A": DEPOSITED_OG_MAX, "ligands": {}}
    for ligand_id in ("CHEMBL434", "CHEMBL776"):
        path = POSE_ROOT / f"{ligand_id}.pdbqt"
        rows = []
        chosen = {"energy": None, "geometry": None}
        for model_index, lines in enumerate(models_of(path), start=1):
            parsed = analyze(lines)
            geometry = contact_table(parsed, ods)
            energy = affinity_of(lines)
            row = {"model": model_index, "search_affinity": energy, **geometry}
            rows.append(row)
            if chosen["energy"] is None or energy < chosen["energy"]["search_affinity"]:
                chosen["energy"] = {**row, "ligand_lines": parsed["ligand_lines"], "ser": parsed["ser"]}
            key = (
                int(geometry["amine_passes"]) + sum(int(item["passes"]) for item in geometry["phenols"]) // 2,
                geometry["closest_pair_og"]["sum"] if geometry["closest_pair_og"] else 99,
            )
            # Geometry-best uses the October 26 key: hydrogen-bond count, then clearance is unavailable
            # here, so the stored flex decision remains the geometry record. Score the best-energy mode.
            _ = key
        energy_row = chosen["energy"]
        ligand_path = HERE / f"{ligand_id}_best_energy.pdbqt"
        kept = [line for line in energy_row["ligand_lines"] if not line.startswith(("MODEL", "ENDMDL"))]
        ligand_path.write_text("\n".join(kept) + "\n", encoding="utf-8")
        receptor_path = patch_receptor(energy_row["ser"], HERE / f"3sn6_{ligand_id}_best_energy.pdbqt")
        scored = manual.score_only(receptor_path, ligand_path)
        summary["ligands"][ligand_id] = {
            "best_energy_model": energy_row["model"],
            "search_affinity": energy_row["search_affinity"],
            "score_only": scored["affinity"],
            "reaches_deposited_og_in_any_mode": any(row["reaches_deposited_og"] for row in rows),
            "closest_pair_by_model": [
                {"model": row["model"], "search_affinity": row["search_affinity"], "pair": row["closest_pair_og"], "amine_passes": row["amine_passes"]}
                for row in rows
            ],
        }
        print(ligand_id, energy_row["search_affinity"], scored["affinity"], summary["ligands"][ligand_id]["reaches_deposited_og_in_any_mode"], flush=True)
    (HERE / "followup_eval.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("manual", manual_affinity)


if __name__ == "__main__":
    main()
