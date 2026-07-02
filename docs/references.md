# Off-target reference genomes

Record every downloaded genome here with its exact accession and source so the
specificity screen is reproducible. Fill in accessions when fetched.

| Role | Species | Assembly accession | Source | Notes |
|------|---------|--------------------|--------|-------|
| Target | _Tropilaelaps mercedesae_ | (incoming) | in-house assembly | primary subject |
| Off-target (host) | _Apis mellifera_ | TBD (RefSeq) | NCBI | main honey bee host |
| Off-target (host) | _Apis cerana_ | TBD | NCBI | Asian honey bee — natural tropi host |
| Off-target (mite) | _Varroa destructor_ | TBD (RefSeq) | NCBI | co-occurring parasitic mite |
| Off-target (mite) | _Varroa jacobsoni_ | TBD | NCBI/GenBank | Asian-range Varroa |
| Specificity | _Tropilaelaps clareae_ / spp. | TBD if available | NCBI | species-level discrimination |
| Debris flora (opt.) | _Vairimorpha (Nosema)_, wax moth, hive microbes | TBD | NCBI | realistic sample background |

## Fetch method
Use NCBI `datasets` CLI or direct FTP. Store genomes in `data/reference/`,
never committed (see `.gitignore`). Log MD5 checksums in `manifest.tsv`.
