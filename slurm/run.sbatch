#!/usr/bin/env bash
# Run one script on the cluster:  sbatch slurm/run.sbatch <script.py|script.R> [args]
# Adjust the partition and resources to your site; MFX_PYTHON / MFX_RSCRIPT select the interpreters.
#SBATCH --job-name=mfx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mfx_%j.log
set -euo pipefail
export MFX_ROOT="${MFX_ROOT:-$PWD}"
export PYTHONPATH="$MFX_ROOT/src:${PYTHONPATH:-}" PYTHONNOUSERSITE=1 HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
script=$1; shift
case "$script" in
  *.R)  "${MFX_RSCRIPT:-Rscript}" "$script" "$@" ;;
  *.sh) bash "$script" "$@" ;;
  *)   "${MFX_PYTHON:-python}" -u "$script" "$@" ;;
esac
