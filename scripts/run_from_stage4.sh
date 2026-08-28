#!/usr/bin/env bash
# Resume the pipeline from stage 4. Re-running is safe: each stage is skipped if
# its output already exists, so a killed run can be restarted without redoing work.
set -euo pipefail

# Run from the repo root regardless of where this is invoked from.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

# Activate the env only if it isn't already active and we can find a conda.
# On a cluster you have usually already done `conda activate` (or loaded a
# module), so this must not assume a particular install location.
if [ -z "${CONDA_PREFIX:-}" ] || [ "$(basename "${CONDA_PREFIX:-}")" != "tropi-primers" ]; then
  for _c in "${CONDA_ROOT:-}" "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" \
            "$(conda info --base 2>/dev/null)"; do
    if [ -n "$_c" ] && [ -f "$_c/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      . "$_c/etc/profile.d/conda.sh"
      conda activate tropi-primers && break
    fi
  done
fi
command -v blastn >/dev/null 2>&1 || {
  echo "tropi-primers env not active and could not be activated." >&2
  echo "Run 'conda activate tropi-primers' (or load your site's modules) first." >&2
  exit 1
}

# Default to all available cores; override with THREADS=n.
if [ -z "${THREADS:-}" ]; then
  THREADS="$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) )"
fi
export THREADS

LOG="$REPO/logs"; mkdir -p "$LOG"
R="results/candidates"
mkdir -p "$R" data/interim
ts() { date '+%H:%M:%S'; }
have() { [ -s "$1" ]; }

if have "$R/conserved_cores.fasta"; then
  echo "[$(ts)] STAGE 4: skipped (conserved_cores.fasta exists)"
else
  echo "[$(ts)] STAGE 4: copy-number & conservation ranking"
  ./scripts/04_copy_number_ranking.py data/raw/tropi_assembly.fasta "$R/unique_candidates.fasta" 2>&1 | tee "$LOG/stage04.log"
fi

if have "$R/primers.tsv"; then
  echo "[$(ts)] STAGE 5: skipped (primers.tsv exists)"
else
  echo "[$(ts)] STAGE 5: PCR primer design"
  ./scripts/05_primer_design.py "$R/conserved_cores.fasta" 2>&1 | tee "$LOG/stage05.log"
fi

if have "$R/validated_primers.tsv"; then
  echo "[$(ts)] STAGE 6: skipped (validated_primers.tsv exists)"
else
  echo "[$(ts)] STAGE 6: in-silico PCR validation"
  ./scripts/06_ispcr_validation.sh "$R/primers.tsv" data/raw/tropi_assembly.fasta data/reference/*.fna 2>&1 | tee "$LOG/stage06.log"
fi

if have "$R/lamp_primers.tsv"; then
  echo "[$(ts)] STAGE 5L: skipped (lamp_primers.tsv exists)"
else
  echo "[$(ts)] STAGE 5L: LAMP primer design"
  ./scripts/05L_lamp_primer_design.py "$R/conserved_cores.fasta" 2>&1 | tee "$LOG/stage05L.log"
fi

if have "$R/lamp_validated.tsv"; then
  echo "[$(ts)] STAGE 6L: skipped (lamp_validated.tsv exists)"
else
  echo "[$(ts)] STAGE 6L: LAMP specificity validation"
  ./scripts/06L_lamp_validation.py "$R/lamp_primers.tsv" data/raw/tropi_assembly.fasta data/reference/*.fna 2>&1 | tee "$LOG/stage06L.log"
fi

echo "[$(ts)] DONE. Results in $REPO/results/candidates/"
echo "  PCR  : primers.tsv -> validated_primers.tsv"
echo "  LAMP : lamp_primers.tsv -> lamp_validated.tsv (+ lamp_rejected.tsv)"
