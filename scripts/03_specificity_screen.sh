#!/usr/bin/env bash
# Stage 3: screen candidate sequences against off-target *assembled genomes*.
#
# SCORE, DON'T FILTER. Earlier versions collapsed all BLAST evidence into a binary
# kill list. That threw away the information needed to tell a candidate with one
# marginal 30bp hit apart from one sharing 80% identity across its whole length —
# wildly different risks. We now:
#   - screen PERMISSIVELY (SCREEN_IDENT/SCREEN_LEN) so the evidence record is rich
#   - record every hit, plus a per-candidate summary (worst identity, longest hit,
#     aggregate fraction of the candidate covered by off-target homology)
#   - apply a SEPARATE, tunable disqualification threshold (DQ_*) to produce the
#     surviving set — rerun with different DQ_* values without re-BLASTing.
# Evidence lives in offtarget_summary.tsv; the cut is a parameter, not a decision
# baked into the search.
#
# NOTE: PCR specificity ultimately lives in the ~40bp of primer footprint, not in
# the candidate region as a whole — a shared fragment is only fatal if a primer
# sits on it. Stage 6 (isPcr) is the authoritative test. This stage ranks risk.
#
# For off-targets that have only raw WGS reads, use 03b_specificity_screen_reads.sh
# — it writes to the same results/candidates/disqualified.txt so the two compose.
#
# Usage: scripts/03_specificity_screen.sh candidates.fasta data/reference/*.fna
set -euo pipefail

CANDIDATES="${1:?need candidates FASTA}"
shift
OFFTARGETS=("$@")
OUT="results/candidates"
mkdir -p "$OUT" data/interim/blastdb

# --- screening sensitivity: what gets RECORDED (permissive on purpose) ---
SCREEN_IDENT="${SCREEN_IDENT:-70}"   # % identity to record a hit
SCREEN_LEN="${SCREEN_LEN:-25}"       # bp alignment length to record a hit

# --- disqualification: what gets CUT (tune after inspecting the summary) ---
DQ_IDENT="${DQ_IDENT:-75}"           # a hit must reach this identity to count toward DQ
DQ_LEN="${DQ_LEN:-30}"               # ...and this length
DQ_COV_FRAC="${DQ_COV_FRAC:-0.0}"    # ...and candidate must have >= this fraction of its
                                     # length covered by such hits. 0.0 = any single
                                     # qualifying hit disqualifies (original behavior).

THREADS="${THREADS:-4}"

hits="$OUT/offtarget_hits.tsv"
: > "$hits"

for ref in "${OFFTARGETS[@]}"; do
  [ -e "$ref" ] || continue
  db="data/interim/blastdb/$(basename "${ref%.*}")"
  [ -f "${db}.nsq" ] || { echo ">> makeblastdb $(basename "$ref")"; makeblastdb -in "$ref" -dbtype nucl -out "$db" >/dev/null; }
  echo ">> screening vs $(basename "$ref")"
  # permissive blastn: catch even weak similarity. qstart/qend retained so we can
  # compute per-candidate covered fraction (and, later, feed primer3 excluded regions).
  blastn -query "$CANDIDATES" -db "$db" \
    -task blastn -word_size 7 -evalue 1 \
    -perc_identity "$SCREEN_IDENT" \
    -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore qstart qend qlen' \
    | awk -v L="$SCREEN_LEN" '$4>=L' \
    | sed "s|^|$(basename "$ref")\t|" >> "$hits"
done

echo ">> summarizing per-candidate off-target risk"
# columns in $hits: ref qseqid sseqid pident length mismatch evalue bitscore qstart qend qlen
summary="$OUT/offtarget_summary.tsv"
sort -k2,2 "$hits" | awk -v DI="$DQ_IDENT" -v DL="$DQ_LEN" -F'\t' '
  {
    c=$2; ident=$4; len=$5; qs=$9; qe=$10; qlen=$11
    n[c]++; qlength[c]=qlen
    if (ident>maxid[c]) maxid[c]=ident
    if (len>maxlen[c]) maxlen[c]=len
    refs[c]=refs[c] (index(refs[c],$1)?"":(refs[c]==""?"":",") $1)
    if (ident>=DI && len>=DL) {
      dq[c]++
      s=(qs<qe?qs:qe); e=(qs<qe?qe:qs)
      # accumulate covered intervals per candidate for a union length
      iv[c]=iv[c] s","e";"
    }
  }
  END {
    for (c in n) {
      # union of DQ-qualifying intervals
      cov=0
      if (c in iv) {
        m=split(iv[c], arr, ";")
        # simple sweep: collect, sort by start
        k=0; delete S; delete E
        for (i=1;i<=m;i++) { if (arr[i]=="") continue; split(arr[i], p, ","); k++; S[k]=p[1]+0; E[k]=p[2]+0 }
        for (i=1;i<k;i++) for (j=1;j<=k-i;j++) if (S[j]>S[j+1]) { t=S[j];S[j]=S[j+1];S[j+1]=t; t=E[j];E[j]=E[j+1];E[j+1]=t }
        cs=S[1]; ce=E[1]
        for (i=2;i<=k;i++) {
          if (S[i]<=ce) { if (E[i]>ce) ce=E[i] }
          else { cov+=ce-cs+1; cs=S[i]; ce=E[i] }
        }
        if (k>0) cov+=ce-cs+1
      }
      frac=(qlength[c]>0)?cov/qlength[c]:0
      printf "%s\t%d\t%d\t%.1f\t%d\t%d\t%.4f\t%s\n", c, qlength[c], n[c], maxid[c], maxlen[c], (c in dq?dq[c]:0), frac, refs[c]
    }
  }' | sort -k7,7gr -k4,4gr > "$summary.body"

printf "candidate_id\tcand_len\tn_hits\tmax_pident\tmax_hit_len\tn_dq_hits\tdq_covered_frac\toff_targets\n" > "$summary"
cat "$summary.body" >> "$summary"; rm -f "$summary.body"

echo ">> applying disqualification threshold (DQ_IDENT=$DQ_IDENT DQ_LEN=$DQ_LEN DQ_COV_FRAC=$DQ_COV_FRAC)"
awk -v F="$DQ_COV_FRAC" -F'\t' 'NR>1 && $6>0 && $7>=F { print $1 }' "$summary" \
  | sort -u > "$OUT/disqualified_by_assemblies.txt"

# shared disqualified list (03b appends to this too)
sort -u "$OUT/disqualified_by_assemblies.txt" > "$OUT/disqualified.txt"

if [ -s "$OUT/disqualified.txt" ]; then
  seqkit grep -v -f "$OUT/disqualified.txt" "$CANDIDATES" > "$OUT/unique_candidates.fasta"
else
  cp "$CANDIDATES" "$OUT/unique_candidates.fasta"
fi

# everything, annotated — candidates that were cut are still here with their evidence
seqkit grep -f <(awk -F'\t' 'NR>1{print $1}' "$summary") "$CANDIDATES" > "$OUT/candidates_with_hits.fasta" 2>/dev/null || true

n_all=$(grep -c '^>' "$CANDIDATES" || echo 0)
n_dq=$(wc -l < "$OUT/disqualified.txt")
echo
echo "Per-hit evidence      -> $hits"
echo "Per-candidate summary -> $summary   (sorted by dq_covered_frac, then identity)"
echo "Disqualified          -> $OUT/disqualified.txt   ($n_dq of $n_all)"
echo "Surviving candidates  -> $OUT/unique_candidates.fasta"
seqkit stats "$OUT/unique_candidates.fasta"
echo
echo "To re-cut WITHOUT re-BLASTing, edit thresholds and rerun the awk on $summary,"
echo "e.g. DQ_COV_FRAC=0.25 to only kill candidates with >=25% of their length shared."
