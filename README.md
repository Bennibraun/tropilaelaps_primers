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
[+ RepeatModeler families, run externally]
        |                                          |
        v                                          |
(1) Assembly QC        (2) Repeat discovery         |
    seqkit stats            RepeatModeler families (external input) + TRF (local)
        |                       --> candidate repeat & satellite consensus set
        |                                          |
        v                                          v
(3) Specificity screen  <---- blastn / minimap2 ---->  REMOVE anything
    keep only candidates with NO significant hit in off-targets
        |
        v
(4) Copy-number & conservation ranking
    blastn self-mapping + mafft; count genomic occurrences, find conserved core
        |
        v
(5) Primer/probe design (primer3) on the conserved core of surviving candidates
        |
        v
(6) In-silico PCR validation (isPcr) against tropi (should amplify)
    AND against every off-target genome (must NOT amplify)
        |
        v
(7) Shortlist -> wet-lab validation
```

Each numbered stage is a script in `scripts/` and a rule in the Snakemake
workflow (`workflow/Snakefile`, added once inputs land). RepeatModeler itself
is not run by this repo — see `docs/plan.md` Stage 2 and "Compute footprint".

## Repo layout

```
data/
  raw/          # the T. mercedesae assembly + RepeatModeler families FASTA (gitignored; tracked via manifest)
  reference/    # off-target genomes: Apis mellifera, Varroa spp., other Tropilaelaps
  interim/      # derived files: repeat libs, blast DBs, alignments, 2bit files
results/
  candidates/   # surviving unique sequences + designed primers (run output, gitignored)
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
- [x] Pipeline scripts written end-to-end (`scripts/00`–`06`)
- [ ] Acquire off-target reference genomes (`scripts/00_fetch_references.sh`)
- [ ] Receive T. mercedesae assembly + RepeatModeler families -> QC (`scripts/01_assembly_qc.sh`)
- [ ] Repeat discovery (`scripts/02_repeat_discovery.sh`)
- [ ] Specificity screen — assembled off-targets (`03`) + read-based off-targets (`03b`)
- [ ] Copy-number & conservation ranking (`scripts/04_copy_number_ranking.py`)
- [ ] Primer design (`scripts/05_primer_design.py`) + in-silico PCR (`scripts/06_ispcr_validation.sh`)
- [ ] Wet-lab shortlist

None of the scripts have been run yet (no assembly in hand) — they're written
against documented tool behavior but unverified end-to-end. Run stage 0 now
(doesn't need the tropi genome) so off-targets are ready when the assembly lands.

See `docs/plan.md` for the detailed methods rationale.
