#!/usr/bin/env bash
# Stage 6: in-silico PCR validation.
# Every primer pair from stage 5 must:
#   - amplify the T. mercedesae assembly (confirms the pair actually works on-target)
#   - amplify NOTHING in any off-target genome (Apis, Varroa, ...)
# Any off-target product kills that pair. Survivors -> validated_primers.tsv,
# the final shortlist handed to wet-lab (Stage 7 in docs/plan.md).
#
# Usage: scripts/06_ispcr_validation.sh primers.tsv data/raw/tropi_assembly.fasta data/reference/*/*.fna
set -euo pipefail

PRIMERS="${1:?need primers.tsv (from scripts/05_primer_design.py)}"
TARGET="${2:?need path to the T. mercedesae assembly FASTA}"
shift 2
OFFTARGETS=("$@")

OUT="results/candidates"
WORK="data/interim/ispcr"
mkdir -p "$OUT" "$WORK/2bit"

MAX_SIZE=300  # upper bound on amplicon size isPcr will report; well above our 70-150bp design window

# Fail loudly if the primer table is missing or has no data rows. Otherwise the
# awk below writes an EMPTY isPcr primer file, isPcr "succeeds" searching with
# nothing, and every pair looks like it survived-then-failed validation — a fake
# scientific result caused by a plumbing error. (This actually happened: a re-run
# with an empty primers.tsv clobbered a good primer file and produced 0 hits.)
if [ ! -s "$PRIMERS" ] || [ "$(awk 'NR>1{c++} END{print c+0}' "$PRIMERS")" -eq 0 ]; then
  echo "ERROR: $PRIMERS is missing or has no primer rows — refusing to run isPcr" >&2
  echo "       (re-run stage 5 first; this stage will not overwrite its inputs with empty data)" >&2
  exit 1
fi

# isPcr primer file: name<TAB>fwd_primer<TAB>rev_primer (extra columns ignored)
PRIMER_FILE_ALL="$WORK/primers_ispcr_all.tsv"
awk -F'\t' 'NR>1 { print $1"_pair"$2"\t"$3"\t"$6 }' "$PRIMERS" > "$PRIMER_FILE_ALL"

# --- prefilter: drop pairs that prime at too many genomic sites ---------------
# A pair whose primers occur at thousands of tandem sites (satellite/repeat
# families dominate the top of stage 4's copy-number ranking) is both a bad assay
# target (smear, not a clean amplicon) AND crashes isPcr: enumerating amplicons
# across thousands of interleaved fwd/rev sites overflows its coordinate bin
# ("start 0, end 0 out of range in findBin"). So we count each primer's perfect
# full-length occurrences via blastn (fast, threaded, reuses stage-4 DB) and drop
# any pair whose fwd OR rev exceeds MAX_PRIMER_SITES. Rejects are logged with
# counts — nothing is silently discarded.
# NOTE on the count vs. what isPcr actually sees: blastn here counts PERFECT
# FULL-LENGTH primer occurrences, but isPcr triggers an alignment on an 11-mer
# tile (default -tileSize=11) and only requires a 15bp perfect 3' match
# (-minPerfect). So isPcr enumerates from MANY more seed sites than this count,
# and pairs near the threshold can still explode into a huge fwd/rev product
# search. This prefilter is a coarse guard, NOT the primary crash fix — the .ooc
# over-used-tile mask below is. Default 150 keeps single- through moderate-copy
# targets while dropping the dense-repeat pairs. Override via MAX_PRIMER_SITES.
MAX_PRIMER_SITES="${MAX_PRIMER_SITES:-150}"
DB="data/interim/copy_number/assembly_db"
if [ ! -f "${DB}.nsq" ]; then
  DB="$WORK/target_blastdb"
  [ -f "${DB}.nsq" ] || makeblastdb -in "$TARGET" -dbtype nucl -out "$DB" >/dev/null
fi

# one FASTA of every primer sequence: <pairname>_F / <pairname>_R
awk -F'\t' '{ print ">"$1"_F\n"$2"\n>"$1"_R\n"$3 }' "$PRIMER_FILE_ALL" > "$WORK/all_primers.fa"
blastn -query "$WORK/all_primers.fa" -db "$DB" \
  -task blastn-short -word_size 7 -perc_identity 100 -qcov_hsp_perc 100 \
  -num_threads "${THREADS:-4}" -outfmt '6 qseqid' \
  | sort | uniq -c | awk '{ print $2"\t"$1 }' > "$WORK/per_primer_hits.tsv"

# per-pair max(fwd,rev); split into keep vs reject at MAX_PRIMER_SITES
REJECT="$OUT/high_copy_rejected.tsv"
PRIMER_FILE="$WORK/primers_ispcr.tsv"
{ printf 'pair_name\tfwd_sites\trev_sites\tmax_sites\n'; } > "$REJECT"
awk -F'\t' -v MAX="$MAX_PRIMER_SITES" -v KEEP="$PRIMER_FILE" -v REJ="$REJECT" '
  FNR==NR { h[$1]=$2; next }                                  # per_primer_hits
  { f=h[$1"_F"]+0; r=h[$1"_R"]+0; m=(f>r?f:r)
    if (m>MAX) { printf "%s\t%d\t%d\t%d\n", $1,f,r,m >> REJ }
    else       { print $0 > KEEP }                            # keep original 3-col isPcr row
  }
' "$WORK/per_primer_hits.tsv" "$PRIMER_FILE_ALL"
[ -f "$PRIMER_FILE" ] || : > "$PRIMER_FILE"

n_all=$(wc -l < "$PRIMER_FILE_ALL")
n_rej=$(( $(wc -l < "$REJECT") - 1 ))
n_pairs=$(wc -l < "$PRIMER_FILE")
echo ">> prefilter: $n_rej/$n_all pairs prime at >${MAX_PRIMER_SITES} genomic sites (logged -> $REJECT, excluded from isPcr)"
if [ "$n_pairs" -eq 0 ]; then
  echo "ERROR: all pairs were filtered out (none <= ${MAX_PRIMER_SITES} sites) — raise MAX_PRIMER_SITES or revisit candidates" >&2
  exit 1
fi
echo ">> $n_pairs primer pairs to validate"

# --- isPcr with an over-used-tile (.ooc) mask --------------------------------
# The crash this stage kept hitting —
#   "start 0, end 0 out of range in findBin (max is 512M)"
# — is NOT a genome-size ceiling. findBin's 512Mb limit is on a SINGLE feature's
# end coordinate; the largest scaffold here is ~10Mb and a 449Mb sub-chunk still
# crashed, so total size was never the cause. The real trigger is amplicon
# enumeration blowing up over a repeat-dense genome: isPcr seeds an alignment at
# every 11bp tile (-tileSize=11) needing only a 15bp perfect 3' match, so
# satellite tiles spawn a combinatorial fwd/rev product search whose coordinate
# arithmetic collapses to the 0,0 interval that trips findBin.
#
# The intended isPcr mechanism for exactly this is an .ooc "over-used tile" file:
# tiles occurring more than -repMatch times are ignored as seeds, which is how
# UCSC runs isPcr/BLAT against whole genomes. We build one per genome from the
# genome itself (makeOoc needs the complete genome, which is what we pass) and
# feed it to every isPcr call. This removes the repeat-driven enumeration at the
# source, so no genome-size chunking is needed.
REP_MATCH="${REP_MATCH:-1024}"   # isPcr default for tileSize 11; lower = mask more repeats

# Run isPcr on one whole FASTA (single 2bit, .ooc-masked) -> data-row PSL at $2.
ispcr_genome() {
  local fasta="$1" out_psl="$2"
  local name; name="$(basename "${fasta%.*}")"
  local total_mb; total_mb=$(seqkit stats -T "$fasta" | awk 'NR==2{printf "%d", $5/1000000}')

  local tb="$WORK/2bit/${name}.2bit"
  [ -f "$tb" ] || faToTwoBit "$fasta" "$tb" >/dev/null

  # Build the over-used-tile mask once per genome. makeOoc requires three
  # positional args (database query output) but exits after reading only the
  # database, so the /dev/null query/output are never opened. The .ooc bakes in
  # tileSize (default 11), which must match the search below — both use the
  # default, so they agree.
  local ooc="$WORK/2bit/${name}.ooc"
  if [ ! -f "$ooc" ]; then
    isPcr "$tb" /dev/null /dev/null -makeOoc="$ooc" -repMatch="$REP_MATCH" >/dev/null 2>&1 \
      || echo "   ($name: warning — makeOoc failed, running without .ooc mask)" >&2
  fi

  # (empty-array expansion under `set -u` throws on the cluster's older bash,
  # so pass the .ooc flag by conditional, not an array splat.)
  if [ -f "$ooc" ]; then
    isPcr "$tb" "$PRIMER_FILE" "$WORK/${name}.psl" \
      -out=psl -maxSize="$MAX_SIZE" -ooc="$ooc"
  else
    isPcr "$tb" "$PRIMER_FILE" "$WORK/${name}.psl" \
      -out=psl -maxSize="$MAX_SIZE"
  fi
  awk '$1 ~ /^[0-9]+$/' "$WORK/${name}.psl" > "$out_psl"
  echo "   ($name: ${total_mb}Mb, ooc-masked)" >&2
}

echo ">> isPcr vs target assembly (must amplify)"
ispcr_genome "$TARGET" "$WORK/target_hits.psl"
# PSL col 10 = qName = the primer-pair name (one row per amplicon it produces).
# Data rows are already filtered by ispcr_genome, so no header-skipping needed.
awk '{ print $10 }' "$WORK/target_hits.psl" | sort -u > "$WORK/amplifies_target.txt"

# --- on-target product COUNT per pair ---
# A pair sitting on a high-copy repeat can prime at many loci in the right
# orientation and spacing, giving multiple products: smeared bands and unreliable
# qPCR quantification. High copy number helps sensitivity but hurts assay
# cleanliness, so count products rather than just checking the pair amplifies.
#   1 product = clean/ideal | 2-5 = tolerable for presence/absence | >5 = smear risk
# target_hits.psl is already filtered to data rows by ispcr_genome; col 10 = pair.
awk '{ print $10 }' "$WORK/target_hits.psl" \
  | sort | uniq -c | awk '{ print $2"\t"$1 }' | sort -k2,2nr > "$WORK/target_product_counts.tsv"
awk -F'\t' '$2==1{a++} $2>=2&&$2<=5{b++} $2>5{c++} END{
  printf "   on-target products: %d pairs=1 (clean) | %d pairs=2-5 | %d pairs>5 (smear risk)\n", a+0,b+0,c+0 }' \
  "$WORK/target_product_counts.tsv"

: > "$WORK/offtarget_hits.tsv"
for ref in "${OFFTARGETS[@]}"; do
  name="$(basename "${ref%.*}")"
  echo ">> isPcr vs off-target: $name (must amplify nothing)"
  ispcr_genome "$ref" "$WORK/${name}_hits.psl"   # chunk-aware; data-row PSL
  awk -v n="$name" '{ print n"\t"$10 }' "$WORK/${name}_hits.psl" >> "$WORK/offtarget_hits.tsv"
done
cut -f2 "$WORK/offtarget_hits.tsv" | sort -u > "$WORK/hits_any_offtarget.txt"

# survivors: amplify target AND absent from every off-target hit list
comm -23 "$WORK/amplifies_target.txt" "$WORK/hits_any_offtarget.txt" > "$WORK/validated_pair_names.txt"

VALIDATED="$OUT/validated_primers.tsv"
# append n_target_products so the shortlist shows assay cleanliness alongside each pair
{ head -n1 "$PRIMERS" | tr -d '\n'; printf '\tn_target_products\n'; } > "$VALIDATED"
awk -F'\t' 'NR>1 { print $1"_pair"$2"\t"$0 }' "$PRIMERS" \
  | grep -Ff "$WORK/validated_pair_names.txt" -w \
  | awk -F'\t' -v CF="$WORK/target_product_counts.tsv" '
      BEGIN { while ((getline line < CF) > 0) { split(line, pp, "\t"); cnt[pp[1]]=pp[2] } }
      { key=$1; row=$0; sub(/^[^\t]*\t/, "", row); print row"\t"(key in cnt?cnt[key]:0) }
    ' >> "$VALIDATED" || true

echo
echo "Off-target hit log       -> $WORK/offtarget_hits.tsv"
echo "Validated primer pairs   -> $VALIDATED"
echo "$(($(wc -l < "$VALIDATED") - 1)) / $n_pairs pairs survived (amplify target, silent on all off-targets)"
