#!/usr/bin/env bash
# Stage 1: QC the incoming T. mercedesae assembly before anything else runs on it.
# Repeats are often collapsed in assemblies (see docs/plan.md) — a mediocre N50
# doesn't block the pipeline, this is a sanity check / record, not a gate.
#
# Usage: scripts/01_assembly_qc.sh data/raw/tropi_assembly.fasta
set -euo pipefail

ASSEMBLY="${1:?need path to the T. mercedesae assembly FASTA}"
OUT="data/interim/qc"
mkdir -p "$OUT"

echo ">> seqkit stats"
seqkit stats -a "$ASSEMBLY" | tee "$OUT/seqkit_stats.txt"

echo ">> per-contig length/GC"
seqkit fx2tab -nlg "$ASSEMBLY" > "$OUT/contig_len_gc.tsv"

# Optional: BUSCO completeness (arachnid lineage). Off by default — it downloads a
# lineage dataset and runs augustus/metaeuk, which is heavy for a laptop, so it's
# deliberately kept out of env/environment.yml. Install separately if you want it:
#   conda install -c bioconda busco
# Then: RUN_BUSCO=1 scripts/01_assembly_qc.sh data/raw/tropi_assembly.fasta
if [ "${RUN_BUSCO:-0}" = "1" ]; then
  if ! command -v busco >/dev/null 2>&1; then
    echo "RUN_BUSCO=1 but busco is not installed (see comment above)." >&2
    exit 1
  fi
  echo ">> BUSCO (arachnida_odb10)"
  busco -i "$ASSEMBLY" -o busco_tropi --out_path "$OUT" \
    -l arachnida_odb10 -m genome -c "${THREADS:-4}"
fi

echo "QC summary   -> $OUT/seqkit_stats.txt"
echo "Contig table -> $OUT/contig_len_gc.tsv"
