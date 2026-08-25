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
  if [ -f "$REFDIR/${name}.fna" ]; then
    echo ">> $name  ($acc) — already fetched, skipping"
    continue
  fi
  echo ">> $name  ($acc)"
  datasets download genome accession "$acc" --include genome \
    --filename "$REFDIR/${name}.zip"
  # NCBI's zip nests the genomic FASTA several directories deep
  # (ncbi_dataset/data/<accession>/<accession>_..._genomic.fna) — flatten it
  # to $REFDIR/<name>.fna so downstream scripts can just glob data/reference/*.fna.
  unzip -o "$REFDIR/${name}.zip" -d "$REFDIR/${name}_raw" >/dev/null
  fna="$(find "$REFDIR/${name}_raw" -iname '*_genomic.fna' | head -n1)"
  [ -n "$fna" ] || { echo "   no genomic FASTA found in ${name}.zip" >&2; exit 1; }
  mv "$fna" "$REFDIR/${name}.fna"
  rm -rf "$REFDIR/${name}_raw" "$REFDIR/${name}.zip"
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

echo ">> recording MD5s -> $REFDIR/manifest.tsv"
find "$REFDIR" -maxdepth 1 \( -name '*.fna' -o -name '*.fasta' \) -exec md5sum {} \; > "$REFDIR/manifest.tsv"

echo "Done."
