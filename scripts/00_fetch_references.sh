#!/usr/bin/env bash
# Fetch off-target reference genomes into data/reference/.
# Fill in the accessions (see docs/references.md) before running.
# Requires: NCBI 'datasets' CLI  (conda install -c conda-forge ncbi-datasets-cli)
set -euo pipefail

REFDIR="data/reference"
mkdir -p "$REFDIR"

# accession list: edit these once confirmed in docs/references.md
declare -A GENOMES=(
  # ["Apis_mellifera"]="GCF_003254395.2"
  # ["Apis_cerana"]="GCF_XXXXXXXXX.X"
  # ["Varroa_destructor"]="GCF_002443255.1"
  # ["Varroa_jacobsoni"]="GCA_XXXXXXXXX.X"
)

if [ ${#GENOMES[@]} -eq 0 ]; then
  echo "No accessions set yet. Edit GENOMES[] in this script (see docs/references.md)." >&2
  exit 1
fi

for name in "${!GENOMES[@]}"; do
  acc="${GENOMES[$name]}"
  echo ">> $name  ($acc)"
  datasets download genome accession "$acc" --include genome \
    --filename "$REFDIR/${name}.zip"
  unzip -o "$REFDIR/${name}.zip" -d "$REFDIR/${name}"
done

echo "Done. Record MD5s: find $REFDIR -name '*.fna' -exec md5sum {} \; >> data/reference/manifest.tsv"
