# tropi_primers

Discovery of PCR primer/probe targets **unique to _Tropilaelaps mercedesae_** for
environmental surveillance of its geographic spread. Samples are "hive debris"
(gunk from the base of beehives), so the assay must detect tropi-mite DNA against
a heavy background of honey bee and _Varroa_ DNA.

## Design goal

Find a sequence that is:

1. **Present** in the _T. mercedesae_ genome — ideally **high-copy / repeated**,
   so a small amount of mite DNA in a debris sample gives strong signal.
2. **Absent** in the two dominant off-targets: **_Apis_** (honey bee) and
   **_Varroa destructor / V. jacobsoni_**.
3. **Conserved** across _T. mercedesae_ populations (so we don't miss variants in
   the field) but **specific** to the species vs. related _Tropilaelaps_ (a nice-to-have,
   secondary).
4. Amenable to a robust PCR/qPCR assay (amplicon ~70–150 bp for qPCR, primer Tm ~60 °C).

The ideal target is a **species-specific tandem repeat / satellite DNA** — high
copy number maximizes sensitivity, and satellites evolve fast enough to be
lineage-specific. Fallback targets: multi-copy gene families, mitochondrial
regions, or any unique single-copy region if no clean repeat emerges.

## Pipeline overview

```
[T. mercedesae assembly]          [off-target genomes: Apis, Varroa, other Tropilaelaps]
        |                                          |
        v                                          |
(1) Repeat discovery                               |
    RepeatModeler / TRF / Red  --> candidate repeat & satellite families
        |                                          |
        v                                          v
(2) Specificity screen  <---- blastn / nucmer / minimap2 ---->  REMOVE anything
    keep only candidates with NO significant hit in off-targets
        |
        v
(3) Copy-number & conservation ranking
    count genomic occurrences; check they cluster in consensus (conserved core)
        |
        v
(4) Primer/probe design (primer3) on the conserved core of surviving candidates
        |
        v
(5) In-silico PCR validation (isPcr) against tropi (should amplify)
    AND against every off-target genome (must NOT amplify)
        |
        v
(6) Shortlist -> wet-lab validation
```

Each numbered stage is a script in `scripts/` and a rule in the Snakemake
workflow (`workflow/Snakefile`, added once inputs land).

## Repo layout

```
data/
  raw/          # the T. mercedesae assembly (gitignored; tracked via manifest)
  reference/    # off-target genomes: Apis mellifera, Varroa spp., other Tropilaelaps
  interim/      # derived files: repeat libs, blast DBs, alignments
results/
  candidates/   # surviving unique sequences + designed primers (small files tracked)
scripts/        # numbered pipeline steps
notebooks/      # exploratory analysis / QC plots
docs/           # methods notes, off-target accession list, design decisions
env/            # conda environment spec
```

## Off-target reference set (to download into data/reference/)

Minimum viable set — expand as needed. See `docs/references.md` for accessions.

- _Apis mellifera_ (host honey bee) — RefSeq assembly
- _Apis cerana_ (Asian honey bee, the natural tropi host) — RefSeq
- _Varroa destructor_ — RefSeq
- _Varroa jacobsoni_ — GenBank
- _Tropilaelaps clareae_ / other _Tropilaelaps_ spp. if available (for species-level specificity)
- Optional negative controls likely in debris: _Nosema/Vairimorpha_, common hive fungi/bacteria, wax moth (_Galleria_)

## Status

- [x] Repo + structure + environment
- [ ] Acquire off-target reference genomes (`scripts/00_fetch_references.sh`)
- [ ] Receive T. mercedesae assembly -> QC (`scripts/01_assembly_qc.sh`)
- [ ] Repeat discovery
- [ ] Specificity screen — assembled off-targets (`03`) + read-based off-targets (`03b`)
- [ ] Primer design + in-silico PCR
- [ ] Wet-lab shortlist

See `docs/plan.md` for the detailed methods rationale.
