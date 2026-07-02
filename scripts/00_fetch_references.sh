#!/usr/bin/env bash
# Fetch off-target reference genomes into data/reference/.
# Fill in the accessions (see docs/references.md) before running.
# Requires: NCBI 'datasets' CLI  (conda install -c conda-forge ncbi-datasets-cli)
set -euo pipefail

REFDIR="data/reference"
mkdir -p "$REFDIR"

# Off-target genome accessions — confirmed vs NCBI Datasets 2026-07-02.
# See docs/references.md for the full table & rationale.
declare -A GENOMES=(
  ["Apis_mellifera"]="GCF_003254395.2"     # host honey bee (RefSeq, chromosome)
  ["Apis_cerana"]="GCF_029169275.1"        # Asian honey bee, natural tropi host (RefSeq)
  ["Varroa_destructor"]="GCF_002443255.2"  # co-occurring mite (RefSeq)
  ["Varroa_jacobsoni"]="GCF_002532875.2"   # Asian-range Varroa (RefSeq)
  # Optional sanity cross-check — public T. mercedesae assembly (NOT the off-target):
  # ["Tmercedesae_public"]="GCA_002081605.1"
)

for name in "${!GENOMES[@]}"; do
  acc="${GENOMES[$name]}"
  echo ">> $name  ($acc)"
  datasets download genome accession "$acc" --include genome \
    --filename "$REFDIR/${name}.zip"
  unzip -o "$REFDIR/${name}.zip" -d "$REFDIR/${name}"
done

# --- Tropilaelaps congener markers (NO genome exists for T. clareae) ---
# Species-level specificity was locked, but there is no congener assembly to
# subtract against — only ~31 GenBank records, mostly mitochondrial markers.
# Pull them so candidates can at least be checked for divergence at these loci.
# Requires: entrez-direct  (conda install -c bioconda entrez-direct)
echo ">> Tropilaelaps clareae markers (Entrez)"
if command -v esearch >/dev/null 2>&1; then
  esearch -db nuccore -query "Tropilaelaps clareae" \
    | efetch -format fasta > "$REFDIR/T_clareae_markers.fasta" || \
    echo "   (Entrez fetch failed — pull manually; see docs/references.md)" >&2
else
  echo "   entrez-direct not installed; skipping congener markers." >&2
fi

echo "Done. Record MD5s:"
echo "  find $REFDIR \\( -name '*.fna' -o -name '*.fasta' \\) -exec md5sum {} \\; >> $REFDIR/manifest.tsv"
