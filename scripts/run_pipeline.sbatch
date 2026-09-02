#!/usr/bin/env bash
#SBATCH --job-name=tropi_pipeline
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.log
#
# Runs the rest of the pipeline as ONE batch job (never on the login node).
# run_pipeline.sh skips any stage whose output already exists, so with stage 4
# done it proceeds 5 -> 6 -> 5L -> 6L -> 7 (report). Stage 6L (BLAST of every
# LAMP oligo component vs all off-target genomes) is the heaviest; the -c 8 /
# 24G / 4h sizing is set for it, comfortably inside the `short` partition (<=11h).
#
# Submit from the repo root:  sbatch scripts/run_pipeline.sbatch
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

eval "$(conda shell.bash hook)"
conda activate tropi-primers

# Pin every stage's thread count to the allocation. Stage 4 reads --threads /
# $SLURM_CPUS_PER_TASK on its own, but stages 3/6L still read $THREADS, and the
# runner would otherwise default to nproc (all cores on the node). Exporting it
# keeps the whole pipeline honest to the 8-core reservation.
export THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Off-target genomes are flat .fna in data/reference/ (Apis x2, Varroa x2); the
# congener markers (T_clareae_markers.fasta) are intentionally NOT in the *.fna
# glob — they are mito fragments, not a genome, so isPcr/BLAST validation skips
# them by design. Expect exactly 4 off-target lines in the stage 6 and 6L logs.
./scripts/run_pipeline.sh
