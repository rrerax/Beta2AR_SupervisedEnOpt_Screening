#!/usr/bin/env python3
"""Flexible Ser203/Ser207 on 3SN6. Rules fixed before the poses are read.

One receptor: 3SN6 chain R. Only Ser203 and Ser207 side chains move.
Asp113 stays fixed. Incomplete deposited residues are kept as they are.
Calipers are the October 22 exam values and are not widened.
A pose passes only when all three hydrogen bonds pass and every
ligand-receptor heavy-atom contact is at least 2.6 A.
Vina energy is recorded and is not the pass rule.
"""

from __future__ import annotations

import json
import math
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
VINA = REPO / "tools" / "vina.exe"
RIGID_SOURCE = REPO / "data" / "processed" / "receptors" / "active_gs_complex_3SN6.pdbqt"
SCORES = REPO / "followup_2026-10-16" / "gpr53_2026-10-26" / "amine_five_scores.csv"
MANIFEST = Path(r"D:\repositories\Enamine_GPCR_Library_plated\docking\ligand_manifest.csv")
POSE_ROOT = Path(r"D:\repositories\Enamine_GPCR_Library_plated\docking\flex_3sn6_ser203_207")

CONTROLS = ("CHEMBL434", "CHEMBL776")
FLEX = (("R", 203), ("R", 207))
ASP = ("R", 113)
BOX = {
    "center_x": 19.724,
    "center_y": 6.018,
    "center_z": -11.183,
    "size_x": 24.501,
    "size_y": 25.011,
    "size_z": 25.452,
}
DLOGP = 0.5
DMW = 50.0
DHBD = 1
DHBA = 1
HBOND_MAX_A = 2.5
HBOND_MIN_DEG = 120.0
CLEARANCE_MIN_A = 2.6
EXHAUSTIVENESS = 4
NUM_MODES = 9
SEED = 20260728
PARALLEL = 8

RULES = {
    "receptor": "3SN6 chain R, the only structure where a manual Ser203/Ser207 rotamer met the three-contact rule",
    "flexible_residues": ["Ser203", "Ser207"],
    "fixed_residues": ["Asp113 and every other residue, including truncated deposited side chains"],
    "calipers": {"abs_logp": DLOGP, "abs_mw": DMW, "abs_hbd": DHBD, "abs_hba": DHBA, "formal_charge": "equal"},
    "calipers_not_widened": True,
    "hbond": "H...acceptor <= 2.5 A and donor-H-acceptor >= 120 deg",
    "clearance": "every ligand-receptor heavy contact >= 2.6 A",
    "vina_energy_is_not_the_pass_rule": True,
    "pass": "both controls satisfy all three hydrogen bonds and the clearance, and both rank ahead of every matched negative on (hydrogen-bond count, clearance)",
    "exhaustiveness": EXHAUSTIVENESS,
    "num_modes": NUM_MODES,
    "seed": SEED,
    "box": BOX,
}


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def angle(donor, hydrogen, acceptor):
    ba = [donor[i] - hydrogen[i] for i in range(3)]
    bc = [acceptor[i] - hydrogen[i] for i in range(3)]
    denom = math.sqrt(sum(v * v for v in ba) * sum(v * v for v in bc))
    if denom == 0:
        return 0.0
    cosv = max(-1.0, min(1.0, sum(ba[i] * bc[i] for i in range(3)) / denom))
    return math.degrees(math.acos(cosv))


def element(atom):
    kind = atom["type"]
    if kind in {"HD", "HS", "H"}:
        return "H"
    if kind in {"OA", "OS", "O"}:
        return "O"
    if kind in {"NA", "NS", "N"}:
        return "N"
    if kind in {"A", "C", "CA"}:
        return "C"
    name = atom["name"]
    return "H" if name.startswith("H") else name[0]


def parse_atoms(text):
    atoms = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            {
                "line": line,
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip(),
                "resi": int(line[22:26]),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "type": line.split()[-1],
            }
        )
    return atoms


def neighbors(atoms):
    adj = {i: set() for i in range(len(atoms))}
    for i, atom in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            other = atoms[j]
            limit = 1.25 if element(atom) == "H" or element(other) == "H" else 1.70
            if dist(atom["xyz"], other["xyz"]) <= limit:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def select_negatives(scores, manifest):
    controls = scores[scores["ligand_id"].isin(CONTROLS)]
    if set(controls["ligand_id"]) != set(CONTROLS):
        raise SystemExit("both controls must be in the public score table")
    ref = controls.iloc[0]
    pool = scores[~scores["ligand_id"].isin(CONTROLS)].copy()
    keep = (
        (pool["crippen_logp"] - ref["crippen_logp"]).abs().le(DLOGP)
        & (pool["mw"] - ref["mw"]).abs().le(DMW)
        & (pool["hbd"] - ref["hbd"]).abs().le(DHBD)
        & (pool["hba"] - ref["hba"]).abs().le(DHBA)
        & (pool["formal_charge_pH74"] == ref["formal_charge_pH74"])
    )
    chosen = pool.loc[keep].copy()
    paths = manifest.set_index("ligand_id")["ligand_pdbqt"]
    chosen["ligand_pdbqt"] = chosen["ligand_id"].map(paths)
    if chosen["ligand_pdbqt"].isna().any():
        missing = chosen.loc[chosen["ligand_pdbqt"].isna(), "ligand_id"].tolist()
        raise SystemExit(f"missing pdbqt for {missing}")
    chosen = chosen.sort_values("ligand_id")
    return chosen, {
        "reference_logp": round(float(ref["crippen_logp"]), 3),
        "reference_mw": round(float(ref["mw"]), 3),
        "reference_hbd": int(ref["hbd"]),
        "reference_hba": int(ref["hba"]),
        "reference_charge": int(ref["formal_charge_pH74"]),
        "matched_negatives": int(len(chosen)),
        "calipers_widened": False,
    }


def top_profile(scores):
    top = scores.nsmallest(max(int(round(len(scores) * 0.01)), 1), "five_mean")
    controls = scores[scores["ligand_id"].isin(CONTROLS)]
    columns = ["five_mean", "crippen_logp", "mw", "hbd", "hba", "hac", "formal_charge_pH74"]

    def row(frame, label):
        med = frame[columns].median()
        return {
            "group": label,
            "n": int(len(frame)),
            **{name: round(float(med[name]), 3) for name in columns},
            "pains": int((frame["pains"].fillna("") != "").sum()),
            "brenk": int((frame["brenk"].fillna("") != "").sum()),
        }

    return [row(top, "top_1pct"), row(controls, "controls"), row(scores, "library")]


def write_receptor():
    atoms = parse_atoms(RIGID_SOURCE.read_text(errors="ignore"))
    flex_atoms = {key: [] for key in FLEX}
    rigid_lines = []
    for atom in atoms:
        key = (atom["chain"], atom["resi"])
        if key not in flex_atoms:
            rigid_lines.append(atom["line"])
            continue
        flex_atoms[key].append(atom)
    blocks = []
    for chain, resi in FLEX:
        group = flex_atoms[(chain, resi)]
        by_name = {atom["name"]: atom for atom in group}
        required = ["CA", "CB", "OG"]
        if any(name not in by_name for name in required):
            raise SystemExit(f"Ser {resi} is missing CA, CB, or OG")
        hydrogens = [atom for atom in group if element(atom) == "H"]
        hg = min(hydrogens, key=lambda atom: dist(atom["xyz"], by_name["OG"]["xyz"]))
        if dist(hg["xyz"], by_name["OG"]["xyz"]) > 1.25:
            raise SystemExit(f"Ser {resi} hydroxyl hydrogen was not found")
        side = [by_name["CA"], by_name["CB"], by_name["OG"], hg]
        keep_rigid = [atom for atom in group if atom not in side]
        rigid_lines.extend(atom["line"] for atom in keep_rigid)
        records = []
        for serial, atom in enumerate(side, start=1):
            records.append(
                "ATOM  {serial:5d} {name:4s} {resn:3s} {chain:1s}{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00    +0.000 {atype:>2s}".format(
                    serial=serial,
                    name=atom["line"][12:16],
                    resn="SER",
                    chain=chain,
                    resi=resi,
                    x=atom["xyz"][0],
                    y=atom["xyz"][1],
                    z=atom["xyz"][2],
                    atype=atom["type"],
                )
            )
        blocks.append(
            "\n".join(
                [
                    f"BEGIN_RES SER {chain} {resi}",
                    "REMARK  3 active torsions:",
                    "ROOT",
                    records[0],
                    "ENDROOT",
                    "BRANCH   1   2",
                    records[1],
                    "BRANCH   2   3",
                    records[2],
                    "BRANCH   3   4",
                    records[3],
                    "ENDBRANCH   3   4",
                    "ENDBRANCH   2   3",
                    "ENDBRANCH   1   2",
                    f"END_RES SER {chain} {resi}",
                ]
            )
        )
    rigid_path = HERE / "3sn6_rigid.pdbqt"
    flex_path = HERE / "3sn6_ser203_207.pdbqt"
    rigid_path.write_text("\n".join(rigid_lines) + "\n", encoding="utf-8")
    flex_path.write_text("\n".join(blocks) + "\n", encoding="utf-8")
    return rigid_path, flex_path


def dock_one(task):
    out_path = Path(task["out_path"])
    log_path = Path(task["log_path"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(VINA),
        "--receptor",
        task["rigid"],
        "--flex",
        task["flex"],
        "--ligand",
        task["ligand_pdbqt"],
        "--center_x",
        str(BOX["center_x"]),
        "--center_y",
        str(BOX["center_y"]),
        "--center_z",
        str(BOX["center_z"]),
        "--size_x",
        str(BOX["size_x"]),
        "--size_y",
        str(BOX["size_y"]),
        "--size_z",
        str(BOX["size_z"]),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--num_modes",
        str(NUM_MODES),
        "--cpu",
        "1",
        "--seed",
        str(SEED),
        "--out",
        str(out_path),
    ]
    with log_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    return {**task, "returncode": completed.returncode, "out_exists": out_path.exists()}


def models_of(path):
    text = Path(path).read_text(errors="ignore")
    chunks = []
    current = []
    score = None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = [line]
            score = None
            continue
        if line.startswith("ENDMDL"):
            current.append(line)
            chunks.append((score, "\n".join(current)))
            current = []
            continue
        if "VINA RESULT" in line:
            parts = line.split()
            try:
                score = float(parts[3])
            except (IndexError, ValueError):
                score = None
        if current:
            current.append(line)
    if not chunks and text.strip():
        chunks.append((None, text))
    return chunks


def pose_geometry(model_text, ligand_reference):
    model_atoms = parse_atoms(model_text)
    ligand = [atom for atom in model_atoms if atom["resn"] == "UNL"]
    flex = [atom for atom in model_atoms if atom["resn"] == "SER"]
    if not ligand:
        ligand = [atom for atom in model_atoms if atom["resn"] not in {"SER", "ASP"}]
    reference = parse_atoms(Path(ligand_reference).read_text(errors="ignore"))
    if len(ligand) != len(reference):
        typed = reference
        coords = ligand
        count = min(len(typed), len(coords))
        atoms = []
        for index in range(count):
            item = dict(typed[index])
            item["xyz"] = coords[index]["xyz"]
            atoms.append(item)
    else:
        atoms = []
        for ref, got in zip(reference, ligand):
            item = dict(ref)
            item["xyz"] = got["xyz"]
            atoms.append(item)
    link = neighbors(atoms)
    phenolic_h = []
    amine_h = []
    for index, atom in enumerate(atoms):
        if element(atom) != "H":
            continue
        for other_index in link[index]:
            parent = atoms[other_index]
            if element(parent) == "O" and any(atoms[n]["type"] == "A" for n in link[other_index]):
                phenolic_h.append((parent["xyz"], atom["xyz"]))
            if element(parent) == "N" and not any(atoms[n]["type"] == "A" for n in link[other_index]):
                amine_h.append((parent["xyz"], atom["xyz"]))
    serines = {}
    for resi in (203, 207):
        group = [atom for atom in flex if atom["resi"] == resi]
        og = next((atom["xyz"] for atom in group if atom["name"] == "OG"), None)
        serines[resi] = og
    asp = [atom["xyz"] for atom in parse_atoms(RIGID_SOURCE.read_text(errors="ignore")) if (atom["chain"], atom["resi"]) == ASP and atom["name"] in {"OD1", "OD2"}]
    def best_hbond(donors, acceptors):
        found = None
        for donor, hydrogen in donors:
            for acceptor in acceptors:
                if acceptor is None:
                    continue
                distance = dist(hydrogen, acceptor)
                value = angle(donor, hydrogen, acceptor)
                passes = distance <= HBOND_MAX_A and value >= HBOND_MIN_DEG
                item = (passes, round(distance, 3), round(value, 1))
                if found is None or (item[0], -item[1]) > (found[0], -found[1]):
                    found = item
        return found or (False, None, None)
    amine = best_hbond(amine_h, asp)
    phenol_count, phenol_one, phenol_two, map_one, map_two = 0, (False, None, None), (False, None, None), None, None
    if phenolic_h and all(value is not None for value in serines.values()):
        targets = (203, 207)
        choices = []
        if len(phenolic_h) == 1:
            choices = [((0, targets[0]),), ((0, targets[1]),)]
        else:
            choices = []
            for left in range(len(phenolic_h)):
                for right in range(len(phenolic_h)):
                    if left != right:
                        choices.append(((left, 203), (right, 207)))
        for choice in choices:
            scored = [best_hbond([phenolic_h[index]], [serines[resi]]) for index, resi in choice]
            count = sum(int(item[0]) for item in scored)
            one = scored[0]
            two = scored[1] if len(scored) > 1 else (False, None, None)
            if count > phenol_count:
                phenol_count = count
                phenol_one, phenol_two = one, two
                map_one = choice[0][1]
                map_two = choice[1][1] if len(choice) > 1 else None
    ligand_heavies = [atom["xyz"] for atom in atoms if element(atom) != "H"]
    receptor_heavies = []
    rigid = parse_atoms((HERE / "3sn6_rigid.pdbqt").read_text(errors="ignore"))
    receptor_heavies.extend(atom["xyz"] for atom in rigid if element(atom) != "H")
    receptor_heavies.extend(atom["xyz"] for atom in flex if element(atom) != "H")
    clearance = min(dist(left, right) for left in ligand_heavies for right in receptor_heavies)
    satisfied = int(amine[0]) + int(phenol_count)
    return {
        "hydrogen_bonds": satisfied,
        "amine_passes": bool(amine[0]),
        "amine_distance": amine[1],
        "amine_angle": amine[2],
        "phenol_bonds": int(phenol_count),
        "phenol_one": {"passes": bool(phenol_one[0]), "distance": phenol_one[1], "angle": phenol_one[2], "serine": map_one},
        "phenol_two": {"passes": bool(phenol_two[0]), "distance": phenol_two[1], "angle": phenol_two[2], "serine": map_two},
        "min_heavy_contact": round(float(clearance), 3),
        "full_pass": satisfied == 3 and clearance >= CLEARANCE_MIN_A,
    }


def main():
    scores = pd.read_csv(SCORES)
    manifest = pd.read_csv(MANIFEST)
    negatives, selection = select_negatives(scores, manifest)
    profile = top_profile(scores)
    pd.DataFrame(profile).to_csv(HERE / "top1_versus_controls.csv", index=False)
    negatives.to_csv(HERE / "matched_negatives.csv", index=False)
    rigid, flex = write_receptor()
    jobs = []
    control_paths = manifest.set_index("ligand_id")["ligand_pdbqt"]
    ordered = [(ligand_id, str(control_paths.loc[ligand_id]), "control") for ligand_id in CONTROLS]
    ordered.extend((row.ligand_id, row.ligand_pdbqt, "negative") for row in negatives.itertuples(index=False))
    for ligand_id, pdbqt, role in ordered:
        jobs.append(
            {
                "ligand_id": ligand_id,
                "role": role,
                "ligand_pdbqt": pdbqt,
                "rigid": str(rigid),
                "flex": str(flex),
                "out_path": str(POSE_ROOT / f"{ligand_id}.pdbqt"),
                "log_path": str(POSE_ROOT / f"{ligand_id}.log"),
            }
        )
    print(f"Docking {len(jobs)} ligands, flex Ser203/Ser207", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=PARALLEL) as executor:
        future_map = {executor.submit(dock_one, job): job for job in jobs}
        for future in as_completed(future_map):
            results.append(future.result())
            print(f"docked {len(results)}/{len(jobs)}", flush=True)
    rows = []
    for result in results:
        if result["returncode"] != 0 or not result["out_exists"]:
            rows.append({"ligand_id": result["ligand_id"], "role": result["role"], "status": "failed", "full_pass": False, "hydrogen_bonds": 0, "min_heavy_contact": None})
            continue
        best = None
        for score, model in models_of(result["out_path"]):
            geometry = pose_geometry(model, result["ligand_pdbqt"])
            geometry["vina_score"] = score
            key = (geometry["hydrogen_bonds"], geometry["min_heavy_contact"], -(score if score is not None else 0))
            if best is None or key > best[0]:
                best = (key, geometry)
        chosen = best[1]
        rows.append(
            {
                "ligand_id": result["ligand_id"],
                "role": result["role"],
                "status": "ready",
                "vina_score_of_geometry_pose": chosen["vina_score"],
                **{key: value for key, value in chosen.items() if key != "vina_score"},
            }
        )
    table = pd.DataFrame(rows)
    table["rank_key_hbonds"] = table["hydrogen_bonds"]
    table["rank_key_clearance"] = table["min_heavy_contact"]
    table = table.sort_values(["hydrogen_bonds", "min_heavy_contact", "ligand_id"], ascending=[False, False, True])
    table.to_csv(HERE / "flex_geometry.csv", index=False)
    controls = table[table["role"] == "control"]
    negatives_table = table[table["role"] == "negative"]
    both_pass = bool(len(controls) == 2 and controls["full_pass"].all())
    ahead = True
    if both_pass:
        for _, control in controls.iterrows():
            better = (
                (negatives_table["hydrogen_bonds"] > control["hydrogen_bonds"])
                | (
                    (negatives_table["hydrogen_bonds"] == control["hydrogen_bonds"])
                    & (negatives_table["min_heavy_contact"] >= control["min_heavy_contact"])
                )
            )
            if bool(better.any()):
                ahead = False
    passed = both_pass and ahead and int((table["status"] == "ready").sum()) == len(jobs)
    reason = "pass"
    if not passed:
        if int((table["status"] != "ready").sum()):
            reason = "docking_failed"
        elif not both_pass:
            reason = "three_contacts_not_simultaneous"
        else:
            reason = "controls_do_not_rank_ahead_of_matched_negatives"
    decision = {
        "rules_fixed_before_poses": True,
        "rules": RULES,
        "selection": selection,
        "top1_versus_controls": profile,
        "flex_pass": passed,
        "reason": reason,
        "model_trained": False,
        "library_redocked": False,
        "molecules": json.loads(table.to_json(orient="records")),
    }
    (HERE / "flex_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    print(json.dumps({"flex_pass": passed, "reason": reason, "selection": selection}, indent=2))


if __name__ == "__main__":
    main()
