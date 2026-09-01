#!/usr/bin/env bash
#SBATCH --job-name=06_ispcr
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.log
#
# Stage 6 only: in-silico PCR validation (prefilter BLAST + isPcr vs target and
# all off-target genomes). Runs on a compute node, never the login node.
#
# Submit from the repo root:  sbatch scripts/06_ispcr_validation.sbatch
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

eval "$(conda shell.bash hook)"
conda activate tropi-primers

# Prefilter BLAST reads $THREADS; pin it to the allocation.
export THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Off-targets: flat .fna in data/reference (Apis x2, Varroa x2). Congener markers
# (T_clareae_markers.fasta) are intentionally excluded by the *.fna glob.
# Override the site ceiling if needed, e.g.: MAX_PRIMER_SITES=100 sbatch ...
./scripts/06_ispcr_validation.sh \
    results/candidates/primers.tsv \
    data/raw/tropi_assembly.fasta \
    data/reference/*.fna
