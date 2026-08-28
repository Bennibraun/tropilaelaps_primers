#!/usr/bin/env bash
# Stage 2: build the raw candidate consensus set from two sources:
#   (a) a RepeatModeler families FASTA you already generated (e.g. on a cluster,
#       or wherever RepeatModeler + RepeatClassifier were run) — this is the
#       de novo repeat family library, headers typically "familyName#Class/Subclass".
#   (b) TRF, run here directly on the assembly, to catch tandem repeats/satellites
#       whose periodic structure RepeatModeler's classifier can miss or lump in
#       with "Unspecified".
# RepeatModeler itself is NOT run by this repo — it's one of the heaviest/most
# fragile bioconda installs and a multi-hour-to-multi-day job, so it's assumed to
# already be done by the time you have its output. TRF is fast and safe to run
# locally on a several-hundred-Mb genome.
#
# Usage: scripts/02_repeat_discovery.sh data/raw/tropi_assembly.fasta path/to/repeatmodeler-families.fa
set -euo pipefail

ASSEMBLY="${1:?need path to the T. mercedesae assembly FASTA}"
RM_FAMILIES="${2:?need path to the RepeatModeler families FASTA (<db>-families.fa)}"

WORK="data/interim/repeats"
OUT="data/interim/repeats"
mkdir -p "$WORK"

# --- (a) RepeatModeler families: just re-tag headers so origin is traceable ---
# RepeatModeler/RepeatClassifier headers look like:
#   >rnd-1_family-41#LINE/R1 ( RepeatScout Family Size = 34, ... )
# i.e. "familyName#Class/Subclass" followed by free-text provenance we don't
# want in the ID. We keep family name + class (useful later for ranking) but:
#  - take only the first whitespace token after '#' (drops the "( ... )" provenance)
#  - replace '/' in the class (e.g. "LINE/R1") with '_' — candidate IDs get used
#    as filesystem path components downstream (scripts/04), and RepeatModeler
#    classes routinely contain '/'.
echo ">> tagging RepeatModeler families ($RM_FAMILIES)"
awk '
  /^>/ {
    n = split(substr($0,2), a, "#")
    if (n > 1) {
      split(a[2], b, /[ \t]/)
      cls = b[1]
    } else {
      cls = "unclassified"
    }
    gsub(/[[:space:]\/]/, "_", a[1])
    gsub(/[[:space:]\/]/, "_", cls)
    print ">RM_" a[1] "__" cls
    next
  }
  { print }
' "$RM_FAMILIES" > "$WORK/rm_candidates.fasta"

# --- (b) TRF: whole-genome tandem repeat / satellite scan ---
# Params: match mismatch delta PM PI minscore maxperiod (standard TRF defaults,
# maxperiod widened to 2000 to allow larger satellite monomers).
TRF_PARAMS="2 7 7 80 10 50 2000"
MIN_PERIOD=20   # bp; drop microsatellite-like short-period repeats (low specificity, prone to homoplasy)
MIN_COPIES=5    # tandem copies at the locus; we want genuinely repetitive loci, not incidental duplications

echo ">> TRF: whole-genome tandem repeat scan (parallelized over ${THREADS:-4} workers)"
# TRF is single-threaded, so we split the assembly into chunks and run one TRF
# per chunk in parallel, then concatenate. -ngs output is a flat per-sequence
# stream (@seqname header followed by that sequence's repeats), so concatenating
# chunk outputs is safe as long as no sequence is split across chunks — seqkit
# split2 splits on sequence boundaries, never mid-sequence.
# NOTE: TRF's exit code is NOT a success/failure code (historically it returns a
# repeat-count-derived value), so `set -e` must not treat a nonzero exit as
# failure — we check the output file instead, not $?.
NCHUNK="${THREADS:-4}"
CHUNKDIR="$WORK/trf_chunks"
rm -rf "$CHUNKDIR"; mkdir -p "$CHUNKDIR"

echo "   splitting assembly into $NCHUNK parts"
seqkit split2 -p "$NCHUNK" -O "$CHUNKDIR" -f "$ASSEMBLY" >/dev/null 2>&1

pids=()
for chunk in "$CHUNKDIR"/*.fa "$CHUNKDIR"/*.fasta; do
  [ -e "$chunk" ] || continue
  ( trf "$chunk" $TRF_PARAMS -h -ngs > "${chunk}.trf" 2>/dev/null || true ) &
  pids+=($!)
done
echo "   ${#pids[@]} TRF workers running"
for pid in "${pids[@]}"; do wait "$pid"; done

cat "$CHUNKDIR"/*.trf > "$WORK/trf_ngs.txt" 2>/dev/null || true

if [ ! -s "$WORK/trf_ngs.txt" ]; then
  echo "TRF produced no output — check trf is installed and the assembly path is correct." >&2
  exit 1
fi

# -ngs columns (per TRF README): start end period copies consensus_size
# pct_matches pct_indels score A C G T entropy consensus_seq repeat_seq
# Verify against `trf --help` / TRF README if your installed version differs.
awk '
  /^@/ { seq = substr($1,2); next }
  NF>=14 { print seq"\t"$1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$14 }
' "$WORK/trf_ngs.txt" > "$WORK/trf_hits.tsv"   # seq start end period copies score consensus

awk -v p="$MIN_PERIOD" -v c="$MIN_COPIES" -F'\t' '$4>=p && $5>=c' \
  "$WORK/trf_hits.tsv" > "$WORK/trf_hits.filtered.tsv"

awk -F'\t' '{ printf(">TRF_%s_%s_%s_p%s\n%s\n", $1,$2,$3,$4,$7) }' \
  "$WORK/trf_hits.filtered.tsv" > "$WORK/trf_candidates_raw.fasta"

echo ">> deduping identical TRF consensus sequences (same monomer found at multiple loci)"
seqkit rmdup -s "$WORK/trf_candidates_raw.fasta" -o "$WORK/trf_candidates.fasta" 2>"$WORK/trf_rmdup.log"

# --- merge ---
cat "$WORK/rm_candidates.fasta" "$WORK/trf_candidates.fasta" > "$OUT/candidates_raw.fasta"

echo "RepeatModeler families : $(grep -c '^>' "$WORK/rm_candidates.fasta") families"
echo "TRF satellite/tandem candidates (deduped) : $(grep -c '^>' "$WORK/trf_candidates.fasta")"
echo "Combined raw candidate set -> $OUT/candidates_raw.fasta"
seqkit stats "$OUT/candidates_raw.fasta"
echo
echo "Next: scripts/03_specificity_screen.sh $OUT/candidates_raw.fasta data/reference/*/*.fna"
