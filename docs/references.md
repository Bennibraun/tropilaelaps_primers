# Off-target reference genomes

Record every downloaded genome here with its exact accession and source so the
specificity screen is reproducible. Fill in accessions when fetched.

Accessions confirmed against NCBI Datasets (checked 2026-07-02).

| Role | Species | Assembly accession | Level | Size | Notes |
|------|---------|--------------------|-------|------|-------|
| Target | _Tropilaelaps mercedesae_ | in-house (incoming) | — | — | primary subject. Public alt: **GCA_002081605.1** (scaffold, 352 Mb, XJTLU "Wuxi") — useful as a sanity cross-check |
| Off-target (host) | _Apis mellifera_ | **GCF_003254395.2** | Chromosome | 225 Mb | main honey bee host — RefSeq |
| Off-target (host) | _Apis cerana_ | **GCF_029169275.1** | Chromosome | 223 Mb | Asian honey bee, natural tropi host — RefSeq |
| Off-target (mite) | _Varroa destructor_ | **GCF_002443255.2** | Scaffold | 369 Mb | co-occurring parasitic mite — RefSeq |
| Off-target (mite) | _Varroa jacobsoni_ | **GCF_002532875.2** | Scaffold | 364 Mb | Asian-range Varroa — RefSeq |
| **Specificity — congener** | _Tropilaelaps clareae_ | **NO GENOME** — 31 nuccore records (25 mito/marker) | — | — | **See gap note below.** Pull marker set via Entrez, not a genome |
| Debris flora (opt.) | _Vairimorpha (Nosema)_, _Galleria_ (wax moth), hive microbes | TBD | — | — | realistic sample background |

## ⚠ Species-specificity gap (important)

We locked **species-level specificity** (unique to _T. mercedesae_ vs. other
_Tropilaelaps_), but **there is no public _T. clareae_ (or any congener) genome
assembly** — only ~31 GenBank nucleotide records, mostly mitochondrial markers
(COI, 16S, cytb, ITS). Consequences:

- We **cannot** do a whole-genome subtraction against a congener. True genome-wide
  species specificity is not verifiable in-silico with current public data.
- **Workable fallback:** avoid designing on any repeat/region whose consensus
  matches the available _T. clareae_ marker sequences, and **prioritize candidate
  families that fall in or near loci where congener sequence exists** (so we can at
  least confirm divergence there). Then push species-vs-congener discrimination to
  **wet-lab validation** with _T. clareae_ gDNA if obtainable.
- Pull the congener marker set with:
  `datasets download` won't work (no genome) — use Entrez instead:
  `esearch -db nuccore -query "Tropilaelaps clareae" | efetch -format fasta > data/reference/T_clareae_markers.fasta`
  (requires `entrez-direct`, on bioconda).

## Fetch method
Genomes: NCBI `datasets` CLI (`scripts/00_fetch_references.sh`).
Congener markers: `entrez-direct` (esearch/efetch).
Store everything in `data/reference/` (gitignored). Log MD5s in `manifest.tsv`.
