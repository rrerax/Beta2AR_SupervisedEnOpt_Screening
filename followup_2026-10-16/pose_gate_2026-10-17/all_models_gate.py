#!/usr/bin/env python3
"""Same angle gate as geometry_gate.py, applied to every saved 4LDE mode."""

import json
from pathlib import Path

import geometry_gate as gate

ROOT = Path(__file__).resolve().parent
WS = ROOT.parents[1]
POSES = {
    "CHEMBL434": WS / "outputs" / "week2_2026-10-16" / "poses" / "active_bi167107_4LDE" / "CHEMBL434.pdbqt",
    "CHEMBL776": WS / "work" / "local_redock" / "poses_ex8" / "active_bi167107_4LDE" / "CHEMBL776.pdbqt",
}


def models(path):
    chunks = []
    current = []
    number = None
    score = None
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            current = []
            number = int(line.split()[1])
            score = None
            continue
        if line.startswith("REMARK VINA RESULT:") and number is not None:
            score = float(line.split()[3])
        if line.startswith(("ATOM", "HETATM")) and number is not None:
            current.append(line)
        if line.startswith("ENDMDL") and number is not None:
            chunks.append((number, score, current))
            number = None
    return chunks


def atoms_from(lines):
    text = ROOT / "_model_tmp.pdbqt"
    text.write_text("\n".join(lines) + "\nENDMDL\n", encoding="utf-8")
    return gate.parse_atoms(text)


def main():
    receptor = gate.parse_atoms(gate.REC)
    rows = []
    for ligand, path in POSES.items():
        for number, score, lines in models(path):
            report = gate.ligand_report(ligand, atoms_from(lines), receptor)
            rows.append(
                {
                    "ligand": ligand,
                    "model": number,
                    "vina_score": score,
                    "pass": report["pass"],
                    "amine_contact": report["amine_contact"],
                    "phenol_contacts": report["phenol_contacts"],
                }
            )
    passed = [row for row in rows if row["pass"]]
    by_ligand = {}
    for ligand in POSES:
        hits = [row["model"] for row in passed if row["ligand"] == ligand]
        by_ligand[ligand] = hits
    payload = {
        "receptor": "4LDE",
        "rule": "H...acceptor <= 2.5 A and donor-H-acceptor >= 120 deg; amine and both phenols",
        "n_models": {ligand: sum(row["ligand"] == ligand for row in rows) for ligand in POSES},
        "passing_models": by_ligand,
        "both_ligands_have_a_passing_pose": all(by_ligand[ligand] for ligand in POSES),
        "rows": rows,
    }
    out = ROOT / "all_models_gate.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (ROOT / "_model_tmp.pdbqt").unlink(missing_ok=True)
    brief = {
        "n_models": payload["n_models"],
        "passing_models": by_ligand,
        "both_ligands_have_a_passing_pose": payload["both_ligands_have_a_passing_pose"],
        "passed": [
            {
                "ligand": row["ligand"],
                "model": row["model"],
                "vina_score": row["vina_score"],
                "amine": row["amine_contact"],
                "phenols": row["phenol_contacts"],
            }
            for row in passed
        ],
    }
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
