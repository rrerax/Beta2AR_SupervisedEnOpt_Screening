from __future__ import annotations

import argparse
import csv
import signal
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

from common import PROJECT_ROOT, ensure_dir, load_config, project_path


RDLogger.DisableLog("rdApp.*")


class LigandPreparationTimeout(RuntimeError):
    pass


def _timeout_handler(signum, frame):
    raise LigandPreparationTimeout("ligand_preparation_timeout")


def shard_name(index: int, shard_size: int) -> str:
    shard_start = (index // shard_size) * shard_size
    shard_end = shard_start + shard_size - 1
    return f"{shard_start:05d}_{shard_end:05d}"


def make_3d_sdf(smiles: str, ligand_id: str, output_sdf: Path, seed: int) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "rdkit_parse_failed"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=int(seed), useRandomCoords=True)
    if status != 0:
        return "embed_failed"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=300)
    mol.SetProp("_Name", ligand_id)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()
    return "sdf_ready"


def convert_to_pdbqt(input_sdf: Path, output_pdbqt: Path, log_path: Path) -> str:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = ["obabel", str(input_sdf), "-O", str(output_pdbqt), "--partialcharge", "gasteiger"]
    with log_path.open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    if completed.returncode != 0 or not output_pdbqt.exists() or output_pdbqt.stat().st_size == 0:
        return "pdbqt_failed"
    return "ready"


def prepare_one(task: dict[str, object]) -> dict[str, str]:
    ligand_id = str(task["ligand_id"])
    smiles = str(task["smiles"])
    index = int(task["index"])
    shard_size = int(task["shard_size"])
    seed = int(task["seed"])
    force = bool(task["force"])
    timeout_seconds = int(task["timeout_seconds"])
    shard = shard_name(index, shard_size)
    processed_root = Path(str(task["processed_root"]))
    log_root = Path(str(task["log_root"]))
    sdf_path = processed_root / "sdf" / shard / f"{ligand_id}.sdf"
    pdbqt_path = processed_root / "pdbqt" / shard / f"{ligand_id}.pdbqt"
    log_path = log_root / shard / f"{ligand_id}.log"

    use_alarm = hasattr(signal, "SIGALRM")
    previous_handler = None
    if use_alarm:
        previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(timeout_seconds)
    try:
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0 and not force:
            status = "ready"
        else:
            sdf_status = make_3d_sdf(smiles, ligand_id, sdf_path, seed + index)
            status = convert_to_pdbqt(sdf_path, pdbqt_path, log_path) if sdf_status == "sdf_ready" else sdf_status
    except LigandPreparationTimeout:
        status = "prep_timeout"
    finally:
        if use_alarm:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, previous_handler)

    return {
        "ligand_id": ligand_id,
        "smiles": smiles,
        "sdf_path": str(sdf_path) if sdf_path.exists() else "",
        "ligand_pdbqt": str(pdbqt_path) if pdbqt_path.exists() else "",
        "status": status,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/beta2ar_screen.yml")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--parallel-jobs", type=int, default=None)
    parser.add_argument("--per-ligand-timeout", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--manifest", default="results/tables/ligand_manifest.csv")
    parser.add_argument("--processed-root", default="data/processed/ligands")
    parser.add_argument("--log-root", default="logs/ligands")
    args = parser.parse_args()

    config = load_config(args.config)
    ligands = pd.read_csv(project_path(config["ligands"]["input_csv"]))
    full_limit = int(config["ligands"].get("full_limit", len(ligands)))
    limit = args.limit or full_limit
    ligands = ligands.head(limit).copy()
    shard_size = int(config["ligands"].get("shard_size", 1000))
    seed = int(config["docking"].get("seed", 20260728))
    id_column = config["ligands"].get("id_column", "ligand_id")
    smiles_column = config["ligands"].get("smiles_column", "smiles")
    parallel_jobs = args.parallel_jobs or min(int(config["docking"].get("parallel_jobs", 32)), 128)

    processed_root = project_path(args.processed_root)
    log_root = project_path(args.log_root)
    ensure_dir(processed_root / "sdf")
    ensure_dir(processed_root / "pdbqt")
    ensure_dir(log_root)
    manifest_path = project_path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["ligand_id", "smiles", "sdf_path", "ligand_pdbqt", "status"]

    tasks = [
        {
            "index": int(index),
            "ligand_id": str(row[id_column]).replace("/", "_").replace(" ", "_"),
            "smiles": str(row[smiles_column]),
            "shard_size": shard_size,
            "seed": seed,
            "force": args.force,
            "timeout_seconds": args.per_ligand_timeout,
            "processed_root": str(processed_root),
            "log_root": str(log_root),
        }
        for index, row in ligands.iterrows()
    ]
    print(f"Preparing {len(tasks)} ligands with {parallel_jobs} workers", flush=True)
    rows: list[dict[str, str]] = []
    if parallel_jobs <= 1:
        for task in tasks:
            rows.append(prepare_one(task))
            if (len(rows) % 500) == 0:
                print(f"prepared {len(rows)} / {len(tasks)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
            future_map = {executor.submit(prepare_one, task): task for task in tasks}
            for future in as_completed(future_map):
                rows.append(future.result())
                if (len(rows) % 500) == 0:
                    print(f"prepared {len(rows)} / {len(tasks)}", flush=True)

    rows.sort(key=lambda row: row["ligand_id"])

    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    ready_count = sum(1 for row in rows if row["status"] == "ready")
    print(f"Wrote {manifest_path} with {ready_count}/{len(rows)} ready ligands")


if __name__ == "__main__":
    main()
