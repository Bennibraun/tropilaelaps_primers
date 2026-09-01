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
# Default 500: from the observed hit distribution this keeps single- through
# high-copy targets (~560 pairs) while dropping the smear/crash zone (>500 sites,
# incl. all the >2000-site pairs that overflow isPcr's coordinate bin). Override
# via the MAX_PRIMER_SITES env var (e.g. 100 for cleaner single-product assays).
MAX_PRIMER_SITES="${MAX_PRIMER_SITES:-500}"
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

to_2bit() {
  local fasta="$1"
  local name; name="$(basename "${fasta%.*}")"
  local twobit="$WORK/2bit/${name}.2bit"
  [ -f "$twobit" ] || faToTwoBit "$fasta" "$twobit" >/dev/null
  echo "$twobit"
}

echo ">> isPcr vs target assembly (must amplify)"
target_2bit="$(to_2bit "$TARGET")"
isPcr "$target_2bit" "$PRIMER_FILE" "$WORK/target_hits.psl" -out=psl -maxSize="$MAX_SIZE"
# PSL output has a 5-line header block; filter to data rows only (col 1 = match
# count, always numeric on a real row) rather than assuming a fixed line count.
awk '$1 ~ /^[0-9]+$/ { print $10 }' "$WORK/target_hits.psl" | sort -u > "$WORK/amplifies_target.txt"

# --- on-target product COUNT per pair ---
# A pair sitting on a high-copy repeat can prime at many loci in the right
# orientation and spacing, giving multiple products: smeared bands and unreliable
# qPCR quantification. High copy number helps sensitivity but hurts assay
# cleanliness, so count products rather than just checking the pair amplifies.
#   1 product = clean/ideal | 2-5 = tolerable for presence/absence | >5 = smear risk
awk '$1 ~ /^[0-9]+$/ { print $10 }' "$WORK/target_hits.psl" \
  | sort | uniq -c | awk '{ print $2"\t"$1 }' | sort -k2,2nr > "$WORK/target_product_counts.tsv"
awk -F'\t' '$2==1{a++} $2>=2&&$2<=5{b++} $2>5{c++} END{
  printf "   on-target products: %d pairs=1 (clean) | %d pairs=2-5 | %d pairs>5 (smear risk)\n", a+0,b+0,c+0 }' \
  "$WORK/target_product_counts.tsv"

: > "$WORK/offtarget_hits.tsv"
for ref in "${OFFTARGETS[@]}"; do
  name="$(basename "${ref%.*}")"
  echo ">> isPcr vs off-target: $name (must amplify nothing)"
  twobit="$(to_2bit "$ref")"
  isPcr "$twobit" "$PRIMER_FILE" "$WORK/${name}_hits.psl" -out=psl -maxSize="$MAX_SIZE"
  awk -v n="$name" '$1 ~ /^[0-9]+$/ { print n"\t"$10 }' "$WORK/${name}_hits.psl" >> "$WORK/offtarget_hits.tsv"
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
