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

# --- chunk-aware isPcr -------------------------------------------------------
# This isPcr build has a CUMULATIVE genome-size ceiling: a single .2bit above
# ~552Mb overflows its coordinate binning and dies with
#   "start 0, end 0 out of range in findBin (max is 512M)".
# Verified empirically: 552Mb and 128Mb subsets run clean, the full 680Mb
# assembly crashes, and no single scaffold exceeds 10Mb (so it is genuinely the
# TOTAL, not any one sequence). glibc 2.17 (CentOS 7) rules out a newer 64-bit
# build, so we partition any oversized genome into sub-ceiling bins of WHOLE
# scaffolds, run isPcr per bin, and concatenate the PSLs. Whole-scaffold bins
# mean no amplicon can straddle a bin boundary (our products are <=300bp; the
# smallest scaffolds are 1kb), and PSL coordinates are per-scaffold so merging
# needs no offset arithmetic.
MAX_2BIT_MB="${MAX_2BIT_MB:-450}"   # safe cap under the observed ~552Mb ceiling

# Run isPcr on one FASTA (chunking if over the cap) -> merged data-row PSL at $2.
ispcr_genome() {
  local fasta="$1" out_psl="$2"
  local name; name="$(basename "${fasta%.*}")"
  local total_mb; total_mb=$(seqkit stats -T "$fasta" | awk 'NR==2{printf "%d", $5/1000000}')

  if [ "$total_mb" -le "$MAX_2BIT_MB" ]; then
    local tb="$WORK/2bit/${name}.2bit"
    [ -f "$tb" ] || faToTwoBit "$fasta" "$tb" >/dev/null
    isPcr "$tb" "$PRIMER_FILE" "$WORK/${name}.psl" -out=psl -maxSize="$MAX_SIZE"
    awk '$1 ~ /^[0-9]+$/' "$WORK/${name}.psl" > "$out_psl"
    echo "   ($name: ${total_mb}Mb, 1 chunk)" >&2
    return
  fi

  # Oversized: split scaffolds into bins <= MAX_2BIT_MB by greedy packing.
  local chunkdir="$WORK/chunks/${name}"
  mkdir -p "$chunkdir"
  # names + lengths, assign each scaffold to a bin whose running total stays
  # under the cap; write one scaffold-name list per bin.
  seqkit fx2tab -nl "$fasta" | awk -v CAP="$((MAX_2BIT_MB*1000000))" -v D="$chunkdir" '
    BEGIN { bin=0; sum=0 }
    { if (sum+$2 > CAP && sum>0) { bin++; sum=0 }
      print $1 > (D"/bin"bin".ids"); sum+=$2 }
    END { print bin+1 }' > "$chunkdir/_nbins"
  local nbins; nbins=$(cat "$chunkdir/_nbins")
  echo "   ($name: ${total_mb}Mb -> $nbins chunks of <= ${MAX_2BIT_MB}Mb)" >&2

  : > "$out_psl"
  local b
  for idlist in "$chunkdir"/bin*.ids; do
    b="$(basename "${idlist%.ids}")"
    seqkit grep -f "$idlist" "$fasta" > "$chunkdir/$b.fa"
    faToTwoBit "$chunkdir/$b.fa" "$chunkdir/$b.2bit" >/dev/null
    isPcr "$chunkdir/$b.2bit" "$PRIMER_FILE" "$chunkdir/$b.psl" -out=psl -maxSize="$MAX_SIZE"
    awk '$1 ~ /^[0-9]+$/' "$chunkdir/$b.psl" >> "$out_psl"
  done
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
