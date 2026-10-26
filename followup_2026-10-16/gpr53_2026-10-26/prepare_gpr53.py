#!/usr/bin/env python3
"""Prepare the plated GPR-53 library for a new beta2-AR docking screen.

The protonation rule is fixed at pH 7.4 before the library is read.
Aliphatic amines, amidines, and guanidines are protonated. Carboxylic acids
and sulfonic acids are deprotonated. Anilines, amides, sulfonamides,
pyridines, imidazoles, and phenols are left neutral. PAINS and Brenk are
recorded on the desalted neutral structure and do not remove compounds.
CHEMBL434 and CHEMBL776 are added as controls when the library does not
already contain them. No model is trained and no docking is started.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

ROOT = Path(r"D:\repositories\Enamine_GPCR_Library_plated")
SOURCE = ROOT / "Enamine_GPCR_Library_plated_53440cmpds_20261004.csv"
OUT = ROOT / "prepared"

CONTROLS = (
    ("CHEMBL434", "CC(C)NCC(O)c1ccc(O)c(O)c1"),
    ("CHEMBL776", "CC(C)NCC(O)c1cc(O)cc(O)c1"),
)
CHECKS = (
    ("CC(C)NCC(O)c1ccc(O)c(O)c1", 1, True, "catechol_A(92)"),
    ("CC(C)NCC(O)c1cc(O)cc(O)c1", 1, True, ""),
    ("Nc1ccccc1", 0, False, ""),
    ("CC(=O)O", -1, False, ""),
    ("c1ccncc1", 0, False, ""),
    ("CC(=O)N", 0, False, ""),
    ("CCN(C)C", 1, True, ""),
    ("Oc1ccccc1", 0, False, ""),
    ("CC(=N)N", 1, False, ""),
)


def catalog(kind: FilterCatalogParams.FilterCatalogs) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(kind)
    return FilterCatalog(params)


PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
BRENK = catalog(FilterCatalogParams.FilterCatalogs.BRENK)


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())


def is_imine_carbon(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 7 and bond.GetBondType() == Chem.BondType.DOUBLE:
            return True
    return False


def is_carbonyl(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() in (8, 16) and bond.GetBondType() == Chem.BondType.DOUBLE:
            return True
    return False


def is_sulfonyl(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 16:
        return False
    double_oxygen = 0
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
            double_oxygen += 1
    return double_oxygen >= 2


def acidic_oxygens(mol: Chem.Mol) -> list[int]:
    chosen = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() < 1:
            continue
        if atom.GetDegree() != 1:
            continue
        parent = atom.GetNeighbors()[0]
        if is_carbonyl(parent):
            chosen.append(atom.GetIdx())
        elif is_sulfonyl(parent):
            chosen.append(atom.GetIdx())
    return chosen


def basic_amines(mol: Chem.Mol) -> list[int]:
    chosen = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        if any(bond.GetBondType() != Chem.BondType.SINGLE for bond in atom.GetBonds()):
            continue
        if not 1 <= atom.GetDegree() <= 3:
            continue
        blocked = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                blocked = True
            elif is_carbonyl(neighbor) or is_sulfonyl(neighbor) or is_imine_carbon(neighbor):
                blocked = True
        if not blocked:
            chosen.append(atom.GetIdx())
    return chosen


def basic_imines(mol: Chem.Mol) -> list[int]:
    pattern = Chem.MolFromSmarts("[CX3;!$(C-a)](=[NX2;!$(N-a);+0])[#7]")
    chosen = []
    for match in mol.GetSubstructMatches(pattern):
        imine = match[1]
        atom = mol.GetAtomWithIdx(imine)
        if atom.GetFormalCharge() == 0:
            chosen.append(imine)
    return chosen


def apply_ph74(mol: Chem.Mol) -> tuple[str, int, int]:
    writable = Chem.RWMol(mol)
    for index in acidic_oxygens(writable):
        writable.GetAtomWithIdx(index).SetFormalCharge(-1)
    amines = basic_amines(writable)
    for index in amines:
        writable.GetAtomWithIdx(index).SetFormalCharge(1)
    for index in basic_imines(writable):
        atom = writable.GetAtomWithIdx(index)
        if atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(1)
    prepared = writable.GetMol()
    Chem.SanitizeMol(prepared)
    charge = Chem.GetFormalCharge(prepared)
    smiles = Chem.MolToSmiles(prepared, canonical=True)
    return smiles, charge, len(amines)


def flags(mol: Chem.Mol) -> tuple[str, str]:
    pains = "|".join(match.GetDescription() for match in PAINS.GetMatches(mol))
    brenk = "|".join(match.GetDescription() for match in BRENK.GetMatches(mol))
    return pains, brenk


def describe(mol: Chem.Mol) -> dict[str, float | int]:
    return {
        "mw": round(float(Descriptors.MolWt(mol)), 3),
        "crippen_logp": round(float(Crippen.MolLogP(mol)), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hac": int(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
    }


def prepare_smiles(smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    parent = largest_fragment(mol)
    Chem.SanitizeMol(parent)
    pains, brenk = flags(parent)
    properties = describe(parent)
    protonated, charge, n_basic = apply_ph74(parent)
    return {
        "smiles_desalted": Chem.MolToSmiles(parent, canonical=True),
        "smiles_pH74": protonated,
        "formal_charge_pH74": charge,
        "n_basic_amine": n_basic,
        "pains": pains,
        "brenk": brenk,
        **properties,
    }


def check_rule() -> None:
    for smiles, charge, basic, pains_name in CHECKS:
        row = prepare_smiles(smiles)
        if row is None:
            raise SystemExit(f"protonation check did not parse {smiles}")
        if row["formal_charge_pH74"] != charge or bool(row["n_basic_amine"]) != basic:
            raise SystemExit(f"protonation check failed for {smiles}: {row}")
        if pains_name and pains_name not in row["pains"]:
            raise SystemExit(f"PAINS check failed for {smiles}: {row['pains']}")
        if not pains_name and "catechol_A(92)" in row["pains"]:
            raise SystemExit(f"unexpected catechol flag for {smiles}")


def read_library() -> list[dict[str, str]]:
    rows = []
    with SOURCE.open(newline="", encoding="utf-8-sig") as handle:
        first = handle.readline()
        if not first.lower().startswith("sep="):
            handle.seek(0)
        for record in csv.DictReader(handle):
            rows.append({"ligand_id": record["Catalog ID"].strip(), "smiles": record["SMILES"].strip()})
    return rows


def main() -> None:
    check_rule()
    OUT.mkdir(parents=True, exist_ok=True)
    library = read_library()
    seen = set()
    prepared_rows = []
    failed = []
    by_parent = {}
    for record in library:
        ligand_id = record["ligand_id"]
        if ligand_id in seen:
            continue
        seen.add(ligand_id)
        prepared = prepare_smiles(record["smiles"])
        if prepared is None:
            failed.append(ligand_id)
            continue
        row = {"ligand_id": ligand_id, "role": "library", "smiles_as_supplied": record["smiles"], **prepared}
        prepared_rows.append(row)
        by_parent.setdefault(prepared["smiles_desalted"], ligand_id)

    control_status = []
    for ligand_id, smiles in CONTROLS:
        prepared = prepare_smiles(smiles)
        if prepared is None:
            raise SystemExit(f"control {ligand_id} did not parse")
        existing = by_parent.get(prepared["smiles_desalted"])
        if existing:
            for row in prepared_rows:
                if row["ligand_id"] == existing:
                    row["role"] = "library_control"
                    row["control_id"] = ligand_id
            control_status.append({"control_id": ligand_id, "disposition": "already_in_library", "ligand_id": existing})
            continue
        prepared_rows.append(
            {"ligand_id": ligand_id, "role": "control", "smiles_as_supplied": smiles, "control_id": ligand_id, **prepared}
        )
        control_status.append({"control_id": ligand_id, "disposition": "appended", "ligand_id": ligand_id})

    fieldnames = [
        "ligand_id", "role", "control_id", "smiles_as_supplied", "smiles_desalted", "smiles_pH74",
        "formal_charge_pH74", "n_basic_amine", "mw", "crippen_logp", "hbd", "hba", "hac", "pains", "brenk",
    ]
    for row in prepared_rows:
        row.setdefault("control_id", "")
    table_path = OUT / "gpr53_pH74_table.csv"
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(prepared_rows)
    docking_path = OUT / "gpr53_docking_input.csv"
    with docking_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand_id", "smiles", "role", "n_basic_amine", "formal_charge_pH74"])
        writer.writeheader()
        for row in prepared_rows:
            writer.writerow({
                "ligand_id": row["ligand_id"],
                "smiles": row["smiles_pH74"],
                "role": row["role"],
                "n_basic_amine": row["n_basic_amine"],
                "formal_charge_pH74": row["formal_charge_pH74"],
            })

    charges = Counter(row["formal_charge_pH74"] for row in prepared_rows)
    summary = {
        "source": SOURCE.name,
        "source_rows": len(library),
        "unique_catalog_ids": len(seen),
        "prepared_rows": len(prepared_rows),
        "parse_failures": failed,
        "with_basic_amine": sum(row["n_basic_amine"] > 0 for row in prepared_rows),
        "without_basic_amine": sum(row["n_basic_amine"] == 0 for row in prepared_rows),
        "pains_flagged": sum(bool(row["pains"]) for row in prepared_rows),
        "brenk_flagged": sum(bool(row["brenk"]) for row in prepared_rows),
        "formal_charge_counts": {str(key): charges[key] for key in sorted(charges)},
        "controls": control_status,
        "pains_brenk_remove_compounds": False,
        "model_trained": False,
        "docking_started": False,
        "protonation": "pH 7.4 fixed rule: aliphatic amine, amidine, and guanidine protonated; carboxylic and sulfonic acids deprotonated",
    }
    (OUT / "gpr53_prep_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "prepared_rows", "with_basic_amine", "without_basic_amine", "pains_flagged",
        "brenk_flagged", "formal_charge_counts", "controls", "parse_failures",
    )}, indent=2))


if __name__ == "__main__":
    main()
