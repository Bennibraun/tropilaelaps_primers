#!/usr/bin/env bash
# Stage 3c: re-cut the stage-3 screen at different thresholds WITHOUT re-BLASTing,
# and apply a minimum candidate length.
#
# Stage 3 records all off-target evidence in offtarget_summary.tsv; the
# disqualification threshold is a parameter, not a decision baked into the search.
# This re-applies it, so exploring thresholds costs seconds instead of hours.
#
# MIN_LEN exists because a candidate shorter than an assay footprint cannot host
# one no matter how unique it is: a qPCR pair needs ~110bp (70bp product + two
# ~20bp primers), a LAMP set needs ~180bp. The stage-2 TRF track emits tens of
# thousands of short monomers (median 33bp) that are unusable by construction and
# make stage 4's self-BLAST quadratic. Filtering them here is removing
# never-viable candidates, not loosening specificity.
#
# Usage:
#   DQ_COV_FRAC=0.25 MIN_LEN=110 scripts/03c_recut.sh
set -euo pipefail

OUT="results/candidates"
SUMMARY="$OUT/offtarget_summary.tsv"
CANDIDATES="${CANDIDATES:-data/interim/repeats/candidates_raw.fasta}"

[ -s "$SUMMARY" ] || { echo "no $SUMMARY — run scripts/03_specificity_screen.sh first" >&2; exit 1; }

DQ_COV_FRAC="${DQ_COV_FRAC:-0.25}"
MIN_LEN="${MIN_LEN:-110}"

echo ">> re-cutting at DQ_COV_FRAC=$DQ_COV_FRAC, MIN_LEN=${MIN_LEN}bp"

# candidates with off-target coverage at or above the threshold are disqualified.
# NOTE: exact ID match — an earlier prefix-matching bug inflated survivor counts.
awk -v F="$DQ_COV_FRAC" -F'\t' 'NR>1 && $6>0 && $7>=F { print $1 }' "$SUMMARY" \
  | sort -u > "$OUT/disqualified_by_assemblies.txt"
cp "$OUT/disqualified_by_assemblies.txt" "$OUT/disqualified.txt"

# survivors = everything not disqualified, then length-filtered
seqkit grep -v -n -f "$OUT/disqualified.txt" "$CANDIDATES" 2>/dev/null \
  | seqkit seq -m "$MIN_LEN" > "$OUT/unique_candidates.fasta"

n_in=$(grep -c '^>' "$CANDIDATES")
n_dq=$(wc -l < "$OUT/disqualified.txt")
n_out=$(grep -c '^>' "$OUT/unique_candidates.fasta" || echo 0)
n_rm=$(grep -c '^>RM_' "$OUT/unique_candidates.fasta" || echo 0)
n_trf=$(grep -c '^>TRF' "$OUT/unique_candidates.fasta" || echo 0)

echo
echo "  input candidates      : $n_in"
echo "  disqualified (hits)   : $n_dq"
echo "  surviving >= ${MIN_LEN}bp : $n_out   (RepeatModeler: $n_rm, TRF: $n_trf)"
echo
seqkit stats "$OUT/unique_candidates.fasta"
echo
echo "Surviving candidates -> $OUT/unique_candidates.fasta"
