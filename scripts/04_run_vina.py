from __future__ import annotations

import argparse
import csv
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from common import ensure_dir, load_config, path_for_table, project_path


SCORE_PATTERN = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")


def parse_vina_log(log_path: Path) -> float | None:
    if not log_path.exists():
        return None
    for line in log_path.read_text(errors="ignore").splitlines():
        match = SCORE_PATTERN.match(line)
        if match:
            return float(match.group(1))
    return None


def run_one(task: dict[str, object]) -> dict[str, object]:
    out_path = Path(str(task["out_path"]))
    log_path = Path(str(task["log_path"]))
    if out_path.exists() and log_path.exists():
        score = parse_vina_log(log_path)
        if score is not None:
            return {**task, "best_score_kcal_mol": score, "status": "ready"}

    out_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "vina",
        "--receptor",
        str(task["receptor_pdbqt"]),
        "--ligand",
        str(task["ligand_pdbqt"]),
        "--center_x",
        str(task["center_x"]),
        "--center_y",
        str(task["center_y"]),
        "--center_z",
        str(task["center_z"]),
        "--size_x",
        str(task["size_x"]),
        "--size_y",
        str(task["size_y"]),
        "--size_z",
        str(task["size_z"]),
        "--exhaustiveness",
        str(task["exhaustiveness"]),
        "--num_modes",
        str(task["num_modes"]),
        "--cpu",
        str(task["cpu_per_job"]),
        "--seed",
        str(task["seed"]),
        "--out",
        str(out_path),
    ]
    with log_path.open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    score = parse_vina_log(log_path)
    status = "ready" if completed.returncode == 0 and score is not None else "failed"
    return {**task, "best_score_kcal_mol": score if score is not None else "", "status": status}


def existing_keys(path: Path) -> set[tuple[str, str]]:
    if not path.exists():
        return set()
    keys: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("status") == "ready":
                keys.add((row["ligand_id"], row["conformation_id"]))
    return keys


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/beta2ar_screen.yml")
    parser.add_argument("--limit-ligands", type=int, default=None)
    parser.add_argument("--parallel-jobs", type=int, default=None)
    parser.add_argument("--exhaustiveness", type=int, default=None)
    parser.add_argument("--score-file", default="results/tables/docking_scores.csv")
    parser.add_argument("--work-dir", default="work/docking")
    parser.add_argument("--log-dir", default="logs/vina")
    parser.add_argument("--ligand-manifest", default="results/tables/ligand_manifest.csv")
    parser.add_argument("--receptor-manifest", default="results/tables/receptor_manifest.csv")
    args = parser.parse_args()

    config = load_config(args.config)
    docking_cfg = config["docking"]
    parallel_jobs = args.parallel_jobs or int(docking_cfg.get("parallel_jobs", 32))
    exhaustiveness = args.exhaustiveness or int(docking_cfg.get("exhaustiveness", 8))
    num_modes = int(docking_cfg.get("num_modes", 5))
    cpu_per_job = int(docking_cfg.get("cpu_per_job", 1))
    seed = int(docking_cfg.get("seed", 20260728))

    receptor_manifest = pd.read_csv(project_path(args.receptor_manifest))
    ligand_manifest = pd.read_csv(project_path(args.ligand_manifest))
    receptor_manifest = receptor_manifest[receptor_manifest["status"] == "ready"].copy()
    ligand_manifest = ligand_manifest[ligand_manifest["status"] == "ready"].copy()
    if args.limit_ligands:
        ligand_manifest = ligand_manifest.head(args.limit_ligands).copy()

    if receptor_manifest.empty:
        raise ValueError("No ready receptors found. Run scripts/02_fetch_receptors.py first.")
    if ligand_manifest.empty:
        raise ValueError("No ready ligands found. Run scripts/03_prepare_ligands.py first.")

    work_dir = ensure_dir(args.work_dir)
    log_dir = ensure_dir(args.log_dir)
    score_path = project_path(args.score_file)
    score_path.parent.mkdir(parents=True, exist_ok=True)
    done = existing_keys(score_path)
    fieldnames = [
        "ligand_id",
        "conformation_id",
        "pdb_id",
        "state_label",
        "best_score_kcal_mol",
        "out_path",
        "log_path",
        "status",
    ]
    write_header = not score_path.exists()

    tasks: list[dict[str, object]] = []
    for _, ligand in ligand_manifest.iterrows():
        ligand_id = str(ligand["ligand_id"])
        ligand_pdbqt = project_path(str(ligand["ligand_pdbqt"]))
        for _, receptor in receptor_manifest.iterrows():
            conformation_id = str(receptor["conformation_id"])
            if (ligand_id, conformation_id) in done:
                continue
            out_path = work_dir / conformation_id / f"{ligand_id}.pdbqt"
            log_path = log_dir / conformation_id / f"{ligand_id}.log"
            tasks.append(
                {
                    "ligand_id": ligand_id,
                    "conformation_id": conformation_id,
                    "pdb_id": receptor["pdb_id"],
                    "state_label": receptor["state_label"],
                    "receptor_pdbqt": project_path(str(receptor["receptor_pdbqt"])),
                    "ligand_pdbqt": ligand_pdbqt,
                    "center_x": receptor["center_x"],
                    "center_y": receptor["center_y"],
                    "center_z": receptor["center_z"],
                    "size_x": receptor["size_x"],
                    "size_y": receptor["size_y"],
                    "size_z": receptor["size_z"],
                    "exhaustiveness": exhaustiveness,
                    "num_modes": num_modes,
                    "cpu_per_job": cpu_per_job,
                    "seed": seed,
                    "out_path": out_path,
                    "log_path": log_path,
                }
            )

    print(f"Docking tasks to run: {len(tasks)}")
    if not tasks:
        return

    completed_count = 0
    with score_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
            future_map = {executor.submit(run_one, task): task for task in tasks}
            for future in as_completed(future_map):
                result = future.result()
                writer.writerow(
                    {
                        "ligand_id": result["ligand_id"],
                        "conformation_id": result["conformation_id"],
                        "pdb_id": result["pdb_id"],
                        "state_label": result["state_label"],
                        "best_score_kcal_mol": result["best_score_kcal_mol"],
                        "out_path": path_for_table(result["out_path"]),
                        "log_path": path_for_table(result["log_path"]),
                        "status": result["status"],
                    }
                )
                completed_count += 1
                if completed_count % 100 == 0:
                    handle.flush()
                    print(f"completed {completed_count}/{len(tasks)}", flush=True)

    print(f"Wrote scores to {score_path}")


if __name__ == "__main__":
    main()
