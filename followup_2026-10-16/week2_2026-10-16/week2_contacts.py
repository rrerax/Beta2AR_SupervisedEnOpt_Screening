#!/usr/bin/env python3
"""Week 2: properties and pocket contacts for the buried actives.

CHEMBL776 poses are reused from the exhaustiveness-8 redock.
CHEMBL434 and the other low-ranked documented actives are docked
with the same box, seed, and Vina 1.2.5 settings.
"""
import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

WS = Path(r"D:\repositories\Codex\2026-08-09\d-desktop-project-ftsz-coumarin-htvs")
REPO = Path(
    r"D:\repositories\Codex\2026-08-04\rrerax-beta2ar-enopt-gpcr-screening-https"
    r"\work\Beta2AR_EnOpt_GPCR_Screening-main"
)
OUT = WS / "outputs" / "week2_2026-10-16"
POSE_OUT = OUT / "poses"
LOG_OUT = OUT / "logs"
LIG_DIR = OUT / "ligands"
for folder in (POSE_OUT, LOG_OUT, LIG_DIR):
    folder.mkdir(parents=True, exist_ok=True)

VINA = r"C:\Users\victor_rax\miniconda3\envs\beta2ar-enopt-screen\Library\bin\vina.exe"
OBABEL = r"C:\Users\victor_rax\miniconda3\envs\beta2ar-enopt-screen\Library\bin\obabel.exe"
REC_DIR = WS / "work" / "local_redock" / "receptors"
OLD_POSE = WS / "work" / "local_redock" / "poses_ex8"
SEED = 20260728

# Learned-rank extremes plus the worst oriented-mean documented actives.
TARGETS = {
    "CHEMBL776": "orciprenaline",
    "CHEMBL434": "catechol isostere",
    "CHEMBL24": "documented miss",
    "CHEMBL13": "worst oriented-mean active",
    "CHEMBL499": "low oriented-mean active",
    "CHEMBL29141": "low oriented-mean active",
}


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def parse_atoms(path, model=1):
    text = Path(path).read_text(errors="ignore").splitlines()
    atoms = []
    current = 1
    seen_model = False
    for line in text:
        if line.startswith("MODEL"):
            seen_model = True
            current = int(line.split()[1])
            continue
        if line.startswith("ENDMDL"):
            if current == model:
                break
            continue
        if seen_model and current != model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip()
        resi = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atype = line[77:79].strip() if len(line) >= 79 else ""
        atoms.append(
            {
                "name": name,
                "resn": resn,
                "chain": chain,
                "resi": resi,
                "xyz": (x, y, z),
                "atype": atype,
                "element": "".join(ch for ch in name if ch.isalpha())[:1],
            }
        )
    return atoms


def prepare_ligand(ligand_id, smiles):
    pdbqt = LIG_DIR / f"{ligand_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=1000) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=2000)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    sdf = LIG_DIR / f"{ligand_id}.sdf"
    Chem.MolToMolFile(mol, str(sdf))
    subprocess.run(
        [OBABEL, str(sdf), "-O", str(pdbqt), "--partialcharge", "gasteiger"],
        check=True,
        capture_output=True,
        text=True,
    )
    return pdbqt


def dock(ligand_id, ligand_pdbqt, row):
    conf = row["conformation_id"]
    out = POSE_OUT / conf / f"{ligand_id}.pdbqt"
    log = LOG_OUT / conf / f"{ligand_id}.log"
    out.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and log.exists() and "mode |" in log.read_text(errors="ignore"):
        return out
    cmd = [
        VINA,
        "--receptor",
        str(REC_DIR / f"{conf}.pdbqt"),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(row["center_x"]),
        "--center_y",
        str(row["center_y"]),
        "--center_z",
        str(row["center_z"]),
        "--size_x",
        str(row["size_x"]),
        "--size_y",
        str(row["size_y"]),
        "--size_z",
        str(row["size_z"]),
        "--exhaustiveness",
        "8",
        "--num_modes",
        "20",
        "--cpu",
        "2",
        "--seed",
        str(SEED),
        "--out",
        str(out),
    ]
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, check=True)
    return out


def pose_path(ligand_id, conf):
    new = POSE_OUT / conf / f"{ligand_id}.pdbqt"
    old = OLD_POSE / conf / f"{ligand_id}.pdbqt"
    if new.exists():
        return new
    if old.exists():
        return old
    return None


def best_score(path):
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def classify_ligand(atoms):
    heavy = [a for a in atoms if a["element"] != "H"]
    aromatic = [a for a in heavy if a["atype"] == "A"]
    roles = []
    for atom in heavy:
        bonded = [other for other in heavy if other is not atom and dist(atom["xyz"], other["xyz"]) < 1.7]
        if atom["element"] == "N":
            roles.append((atom, "amine_N"))
        elif atom["element"] == "O":
            if any(other["atype"] == "A" for other in bonded):
                roles.append((atom, "phenolic_O"))
            else:
                roles.append((atom, "aliphatic_O"))
    return roles, aromatic


def nearest(atom, candidates):
    if not candidates:
        return None, None
    best = min(candidates, key=lambda other: dist(atom["xyz"], other["xyz"]))
    return best, dist(atom["xyz"], best["xyz"])


def contacts(ligand_id, conf, receptor_atoms):
    path = pose_path(ligand_id, conf)
    if path is None:
        return None
    ligand_atoms = parse_atoms(path, model=1)
    roles, aromatic = classify_ligand(ligand_atoms)
    asp_o = [a for a in receptor_atoms if a["resn"] == "ASP" and a["name"] in {"OD1", "OD2"}]
    glu_o = [a for a in receptor_atoms if a["resn"] == "GLU" and a["name"] in {"OE1", "OE2"}]
    ser_og = [a for a in receptor_atoms if a["resn"] == "SER" and a["name"] == "OG"]
    asn_o = [a for a in receptor_atoms if a["resn"] in {"ASN", "GLN"} and a["name"] in {"OD1", "OE1"}]
    carbox = asp_o + glu_o
    ligand_heavy = [a for a in ligand_atoms if a["element"] != "H"]
    nearby = {}
    for lat in ligand_heavy:
        for rat in receptor_atoms:
            if rat["element"] == "H":
                continue
            if dist(lat["xyz"], rat["xyz"]) <= 4.0:
                key = f"{rat['resn']}{rat['resi']}"
                nearby[key] = nearby.get(key, 0) + 1
    row = {
        "ligand_id": ligand_id,
        "conformation_id": conf,
        "ex8_score": best_score(path),
        "nearby_residues": " ".join(sorted(nearby, key=nearby.get, reverse=True)[:12]),
    }
    for role_name in ("amine_N", "phenolic_O", "aliphatic_O"):
        chosen = [atom for atom, role in roles if role == role_name]
        row[f"n_{role_name}"] = len(chosen)
        if not chosen:
            row[f"{role_name}_nearest_ASP_GLU"] = None
            row[f"{role_name}_nearest_SER"] = None
            continue
        asp_hits = []
        ser_hits = []
        for atom in chosen:
            asp, asp_d = nearest(atom, carbox)
            ser, ser_d = nearest(atom, ser_og)
            if asp is not None:
                asp_hits.append((asp_d, f"{asp['resn']}{asp['resi']}"))
            if ser is not None:
                ser_hits.append((ser_d, f"{ser['resn']}{ser['resi']}"))
        if asp_hits:
            d, name = min(asp_hits)
            row[f"{role_name}_nearest_ASP_GLU"] = f"{name}:{d:.2f}"
        else:
            row[f"{role_name}_nearest_ASP_GLU"] = None
        if ser_hits:
            d, name = min(ser_hits)
            row[f"{role_name}_nearest_SER"] = f"{name}:{d:.2f}"
        else:
            row[f"{role_name}_nearest_SER"] = None
    # Aromatic: ligand aromatic carbons within 4.5 A of a PHE/TYR/TRP ring atom.
    ring_res = {}
    ring_atoms = [
        a
        for a in receptor_atoms
        if a["resn"] in {"PHE", "TYR", "TRP"}
        and a["name"] in {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CE3", "CZ2", "CZ3", "CH2"}
    ]
    for lat in aromatic:
        for rat in ring_atoms:
            if dist(lat["xyz"], rat["xyz"]) <= 4.5:
                ring_res[f"{rat['resn']}{rat['resi']}"] = True
    row["aromatic_residues"] = " ".join(sorted(ring_res))
    _ = asn_o
    return row


def main():
    manifest = pd.read_csv(REPO / "results/tables/ligand_manifest.csv").set_index("ligand_id")
    receptors = pd.read_csv(REPO / "results/tables/receptor_manifest.csv")
    receptors = receptors[receptors["status"] == "ready"]
    scores = pd.read_csv(REPO / "results/tables/docking_scores.csv")
    wide = scores.pivot_table(
        index="ligand_id", columns="conformation_id", values="best_score_kcal_mol", aggfunc="first"
    )
    rank = pd.read_csv(WS / "outputs/decoy_validation/library_ranking_with_decoy_model_v2.csv").set_index(
        "ligand_id"
    )
    wide["rank_oriented_mean"] = (-wide.mean(axis=1)).rank(ascending=False, method="min").astype(int)

    prop_rows = []
    for ligand_id, note in TARGETS.items():
        smiles = manifest.loc[ligand_id, "smiles"]
        mol = Chem.MolFromSmiles(smiles)
        prop_rows.append(
            {
                "ligand_id": ligand_id,
                "note": note,
                "smiles": smiles,
                "mw": round(Descriptors.MolWt(mol), 2),
                "logp": round(Crippen.MolLogP(mol), 2),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "hac": rdMolDescriptors.CalcNumHeavyAtoms(mol),
                "rank_oriented_mean": int(wide.loc[ligand_id, "rank_oriented_mean"]),
                "rank_learned": int(rank.loc[ligand_id, "enopt_decoy_rank"]),
                "ex4_mean": round(float(wide.loc[ligand_id].drop(labels=["rank_oriented_mean"]).mean()), 3),
            }
        )
        if ligand_id != "CHEMBL776":
            print("preparing", ligand_id, flush=True)
            lig = prepare_ligand(ligand_id, smiles)
            for _, row in receptors.iterrows():
                print(" docking", ligand_id, row["conformation_id"], flush=True)
                dock(ligand_id, lig, row)

    contact_rows = []
    for ligand_id in TARGETS:
        for _, row in receptors.iterrows():
            conf = row["conformation_id"]
            receptor_atoms = parse_atoms(REC_DIR / f"{conf}.pdbqt")
            item = contacts(ligand_id, conf, receptor_atoms)
            if item is None:
                continue
            item["ex4_score"] = round(float(wide.loc[ligand_id, conf]), 3)
            contact_rows.append(item)

    props = pd.DataFrame(prop_rows).sort_values("rank_learned")
    contact_df = pd.DataFrame(contact_rows)
    props.to_csv(OUT / "week2_properties.csv", index=False)
    contact_df.to_csv(OUT / "week2_contacts.csv", index=False)
    print(props.to_string(index=False))
    print(contact_df.to_string(index=False))


if __name__ == "__main__":
    main()
