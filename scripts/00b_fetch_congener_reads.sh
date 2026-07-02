#!/usr/bin/env bash
# Fetch Tropilaelaps congener WGS reads from SRA, if/when any exist.
#
# As of 2026-07-02 there are NO public WGS runs for T. clareae (or any congener) —
# only T. mercedesae genomic reads + RNA-seq + 16S amplicon. Re-check periodically:
#   esearch -db sra -query 'Tropilaelaps[Organism] NOT mercedesae' | efetch -format runinfo
# and drop any real WGS run accessions (SRR/ERR/DRR) into RUNS below.
#
# Requires: sra-tools (prefetch, fasterq-dump)  (conda install -c bioconda sra-tools)
set -euo pipefail

OUTDIR="data/reference/congener_reads"
mkdir -p "$OUTDIR"

# Add congener WGS run accessions here once they exist (or from your own sequencer,
# in which case skip this script and point 03b straight at your FASTQs).
RUNS=(
  # "SRRXXXXXXX"
)

if [ ${#RUNS[@]} -eq 0 ]; then
  cat >&2 <<'MSG'
No congener WGS runs configured.
Reason: none are public yet. Options:
  1) Re-run the esearch check above and add any new SRR/ERR/DRR WGS runs.
  2) Sequence a congener (T. clareae) yourself — even low-coverage Illumina is
     enough to SCREEN against; no assembly required. Put the FASTQs anywhere and
     run: scripts/03b_specificity_screen_reads.sh candidates.fasta T_clareae R1.fq.gz R2.fq.gz
MSG
  exit 0
fi

for run in "${RUNS[@]}"; do
  echo ">> $run"
  prefetch "$run" -O "$OUTDIR"
  fasterq-dump "$run" -O "$OUTDIR" --split-files
  gzip -f "$OUTDIR/$run"_*.fastq
done
echo "Congener reads -> $OUTDIR"
