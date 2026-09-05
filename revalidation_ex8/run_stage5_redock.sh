#!/usr/bin/env bash
set -euo pipefail
REPO="${1:?usage: run_stage5_redock.sh <repo_dir> [cores]}"
CORES="${2:-$(nproc)}"
cd "$REPO"

# copy helpers next to the repo scripts (run from repo root)
HERE="$(cd "$(dirname "$0")" && pwd)"
cp -n "$HERE/make_validation_inputs.py" "$HERE/eval_ex8.py" . || true

# 0) micromamba env under repo/.tools + repo/.mamba_vina (from environment-vina.yml)
bash scripts/00_setup_env.sh
MICRO="$PWD/.tools/bin/micromamba"

echo "== install autodock-vina (linux binary from conda-forge) =="
"$MICRO" install -y -p "$PWD/.mamba_vina" -c conda-forge autodock-vina

echo "== receptors: 5 pockets (2RH1 5D5A 3SN6 4LDE 4LDL) =="
"$MICRO" run -p "$PWD/.mamba_vina" python scripts/02_fetch_receptors.py --config configs/beta2ar_screen.yml

echo "== merged validation ligand csv (206 actives + 3000 decoys) + config =="
"$MICRO" run -p "$PWD/.mamba_vina" python make_validation_inputs.py --repo .

echo "== prepare ligands on $CORES cores =="
"$MICRO" run -p "$PWD/.mamba_vina" python scripts/03_prepare_ligands.py --config configs/validation_redock.yml --parallel-jobs "$CORES"

echo "== ex8 docking on $CORES cores (resume-safe; ctrl-c and rerun to continue) =="
"$MICRO" run -p "$PWD/.mamba_vina" python scripts/04_run_vina.py \
  --config configs/validation_redock.yml \
  --exhaustiveness 8 \
  --parallel-jobs "$CORES" \
  --score-file results/tables/docking_scores_ex8_validation.csv \
  --work-dir work/docking_ex8 \
  --log-dir logs/vina_ex8

echo "== evaluation: ex8 OOF vs ex4 references (raw 0.696 / raw+ha 0.751) =="
"$MICRO" run -p "$PWD/.mamba_vina" python eval_ex8.py \
  --repo . --score results/tables/docking_scores_ex8_validation.csv

echo "DONE -> results/tables/docking_scores_ex8_validation.csv"
echo "resume: just rerun this script; finished docks are skipped."
