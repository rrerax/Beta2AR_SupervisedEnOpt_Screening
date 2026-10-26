from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def ensure_tool_path() -> None:
    """Put local and conda binaries on PATH for obabel and vina."""
    candidates = [
        PROJECT_ROOT / "tools",
        Path(sys.executable).resolve().parent / "Library" / "bin",
    ]
    current = os.environ.get("PATH", "")
    present = {part.lower() for part in current.split(os.pathsep)}
    prefix = []
    for bin_dir in candidates:
        if not bin_dir.is_dir():
            continue
        entry = str(bin_dir)
        if entry.lower() not in present:
            prefix.append(entry)
            present.add(entry.lower())
    if prefix:
        os.environ["PATH"] = os.pathsep.join(prefix) + os.pathsep + current


ensure_tool_path()


def path_for_table(path: str | Path) -> str:
    """Store a repo-relative path when possible, otherwise an absolute path."""
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(resolved)


def project_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def load_config(config_path: str | Path = "configs/beta2ar_screen.yml") -> dict[str, Any]:
    with project_path(config_path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ensure_dir(path: str | Path) -> Path:
    path = project_path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def receptor_entries(config: dict[str, Any]) -> list[dict[str, Any]]:
    receptors = config.get("receptors", {})
    if "entries" in receptors:
        return list(receptors["entries"])
    return [
        {
            "pdb_id": pdb_id,
            "conformation_id": str(pdb_id).upper(),
            "receptor_chain": None,
            "native_ligand_resname": None,
            "state_label": "unspecified",
        }
        for pdb_id in receptors.get("pdb_ids", [])
    ]


def run_command(command: list[str], log_path: str | Path | None = None) -> int:
    if log_path is None:
        completed = subprocess.run(command, text=True)
        return completed.returncode
    log_path = project_path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return completed.returncode


def write_csv(path: str | Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> Path:
    path = project_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path

