#!/usr/bin/env bash
# Stage 3: screen candidate sequences against off-target *assembled genomes*.
# Keeps only candidates with NO meaningful hit in Apis / Varroa / other off-targets.
# For off-targets that have only raw WGS reads (e.g. a Tropilaelaps congener),
# use 03b_specificity_screen_reads.sh instead — it writes to the same
# results/candidates/disqualified.txt so the two screens compose.
#
# Usage: scripts/03_specificity_screen.sh candidates.fasta data/reference/*/*.fna
set -euo pipefail

CANDIDATES="${1:?need candidates FASTA}"
shift
OFFTARGETS=("$@")
OUT="results/candidates"
mkdir -p "$OUT" data/interim/blastdb

# Thresholds — tune & document. Conservative = err toward calling things NON-unique.
MIN_IDENT=75      # % identity to count as an off-target hit
MIN_LEN=30        # bp alignment length to count as a hit (< typical primer footprint)

hits="$OUT/offtarget_hits.tsv"
: > "$hits"

for ref in "${OFFTARGETS[@]}"; do
  db="data/interim/blastdb/$(basename "${ref%.*}")"
  [ -f "${db}.nsq" ] || makeblastdb -in "$ref" -dbtype nucl -out "$db" >/dev/null
  echo ">> screening vs $(basename "$ref")"
  # permissive blastn: catch even weak similarity
  blastn -query "$CANDIDATES" -db "$db" \
    -task blastn -word_size 7 -evalue 1 \
    -perc_identity "$MIN_IDENT" \
    -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore' \
    | awk -v L="$MIN_LEN" '$4>=L' \
    | sed "s|^|$(basename "$ref")\t|" >> "$hits"
done

# candidates with any recorded hit are disqualified
cut -f2 "$hits" | sort -u > "$OUT/disqualified.txt"
seqkit grep -v -f "$OUT/disqualified.txt" "$CANDIDATES" > "$OUT/unique_candidates.fasta"

echo "Off-target hits logged to $hits"
echo "Surviving unique candidates -> $OUT/unique_candidates.fasta"
seqkit stats "$OUT/unique_candidates.fasta"
