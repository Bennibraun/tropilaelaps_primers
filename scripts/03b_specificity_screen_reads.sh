#!/usr/bin/env bash
# Stage 3b: screen candidate sequences against off-target *raw WGS reads*.
#
# Use this when an off-target has sequencing reads but no assembly — e.g. a
# Tropilaelaps congener (T. clareae) if/when WGS is obtained (from SRA or your own
# sequencer). No assembly needed: we map the reads TO the candidates and reject any
# candidate that gets covered. Unassembled reads are actually preferable here
# because they retain repeat copies that assemblies collapse.
#
# Logic: a candidate is disqualified if off-target reads map across enough of it to
# support a primer/probe footprint. We require both a minimum breadth of coverage
# and a minimum depth, so a few spurious reads don't kill a good candidate.
#
# Usage:
#   Paired:  scripts/03b_specificity_screen_reads.sh candidates.fasta OFFNAME reads_R1.fq.gz reads_R2.fq.gz
#   Single:  scripts/03b_specificity_screen_reads.sh candidates.fasta OFFNAME reads.fq.gz
#
# Requires: minimap2, samtools, bedtools, seqkit
set -euo pipefail

CANDIDATES="${1:?need candidates FASTA}"
OFFNAME="${2:?need a label for this off-target read set, e.g. T_clareae}"
shift 2
READS=("$@")
[ ${#READS[@]} -ge 1 ] || { echo "need at least one reads file" >&2; exit 1; }

OUT="results/candidates"
WORK="data/interim/readscreen/$OFFNAME"
mkdir -p "$OUT" "$WORK"

# --- disqualification thresholds (tune & document) ---
MIN_BREADTH=0.60   # fraction of a candidate covered by >=1 off-target read
MIN_DEPTH=2        # min reads at a base for it to count as "covered"
# Rationale: a real shared region gets broad, repeated coverage. A stray
# low-complexity read gives a thin, narrow smear that stays under these.

# minimap2 preset for short reads; use map-ont for Nanopore off-target reads.
PRESET="${MM2_PRESET:-sr}"

echo ">> mapping $OFFNAME reads (${#READS[@]} file(s), preset=$PRESET) onto candidates"
minimap2 -ax "$PRESET" "$CANDIDATES" "${READS[@]}" 2>"$WORK/minimap2.log" \
  | samtools sort -o "$WORK/aln.bam" -
samtools index "$WORK/aln.bam"

# per-base depth, then per-candidate breadth at >= MIN_DEPTH
samtools depth -a "$WORK/aln.bam" > "$WORK/depth.tsv"

# candidate lengths
seqkit fx2tab -nl "$CANDIDATES" | awk 'BEGIN{OFS="\t"}{print $1,$2}' > "$WORK/lengths.tsv"

# compute breadth = covered_bases(>=MIN_DEPTH) / length, flag if >= MIN_BREADTH
awk -v D="$MIN_DEPTH" -v B="$MIN_BREADTH" '
  FNR==NR { len[$1]=$2; next }                 # lengths.tsv
  $3>=D   { cov[$1]++ }                          # depth.tsv, covered bases
  END {
    for (c in len) {
      br = (len[c]>0) ? cov[c]/len[c] : 0
      status = (br>=B) ? "DISQUALIFIED" : "ok"
      printf "%s\t%s\t%d\t%.3f\t%s\n", c, "'"$OFFNAME"'", len[c], br, status
      if (status=="DISQUALIFIED") print c > "/dev/stderr"
    }
  }
' "$WORK/lengths.tsv" "$WORK/depth.tsv" \
  > "$WORK/breadth.tsv" 2> "$WORK/disqualified_by_reads.txt"

# append to the shared disqualified list so downstream stages see it
sort -u "$WORK/disqualified_by_reads.txt" >> "$OUT/disqualified.txt"
sort -u -o "$OUT/disqualified.txt" "$OUT/disqualified.txt"

echo "Per-candidate breadth -> $WORK/breadth.tsv"
echo "Newly disqualified by $OFFNAME reads:"
cat "$WORK/disqualified_by_reads.txt" || true

# refresh the surviving set (respects hits from the assembly screen too)
if [ -s "$OUT/disqualified.txt" ]; then
  seqkit grep -v -f "$OUT/disqualified.txt" "$CANDIDATES" > "$OUT/unique_candidates.fasta"
else
  cp "$CANDIDATES" "$OUT/unique_candidates.fasta"
fi
echo "Surviving unique candidates -> $OUT/unique_candidates.fasta"
seqkit stats "$OUT/unique_candidates.fasta"
