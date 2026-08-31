#!/usr/bin/env bash
#SBATCH --job-name=04_copy_number
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.log
#
# Stage 4: copy-number & conservation ranking, --top-n 100.
# Measured cost (10-candidate probe): fixed makeblastdb+blast+index dominates;
# per-candidate mafft is sub-second with the fast preset. Full/top-100 runs land
# well under an hour, so `short` (<=11h) with a 2h ceiling is comfortable.
#
# Submit from the repo root:   sbatch scripts/04_copy_number_ranking.sbatch
set -euo pipefail

# Slurm runs from the submit dir; make logs/ exist for --output before anything.
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

# Activate the env exactly as you do interactively.
eval "$(conda shell.bash hook)"
conda activate tropi-primers

# Threads come from the allocation: 04_copy_number_ranking.py --threads defaults
# to $SLURM_CPUS_PER_TASK, so -c 8 above is the single source of truth.
./scripts/04_copy_number_ranking.py \
    data/raw/tropi_assembly.fasta \
    results/candidates/unique_candidates.fasta \
    --top-n 100
