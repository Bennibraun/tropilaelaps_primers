# Detailed methods rationale

The problem in one line: **from a mixed environmental sample (hive debris), detect
the presence of _T. mercedesae_ via PCR without cross-reacting to _Apis_ or _Varroa_.**

This drives every design choice below.

## Why target a repeat/satellite

- **Sensitivity.** Debris is dominated by bee and _Varroa_ DNA; tropi DNA may be a
  tiny fraction. A target present in thousands of copies per genome gives many more
  template molecules per mite than a single-copy gene — the difference between a
  Ct of 25 and "no amplification."
- **Specificity.** Satellite DNA and repeat families turn over rapidly between
  lineages, so they are frequently species- or genus-specific. This is exactly the
  property that lets us exclude _Apis_/_Varroa_.
- **Precedent.** Species-specific satellite qPCR is an established approach for
  insect/mite detection and eDNA-style surveillance.

We still keep single-copy unique regions as a fallback, because a clean, verified
unique region beats a repeat that turns out to be shared.

## Stage-by-stage

### 0. References
Download off-target genomes first, so the specificity filter is ready before we
even look at candidates. Script: `00_fetch_references.sh` (uses NCBI datasets /
direct FTP). Record exact accessions + MD5s in `docs/references.md` and a
`manifest.tsv` so the analysis is reproducible.

### 1. Assembly QC
When the tropi assembly arrives: `seqkit stats`, N50, total length, GC, and a
quick contamination sniff (a lot of the "gunk" DNA can end up in a raw assembly).
BUSCO (arachnid lineage) optional to gauge completeness. We don't need a perfect
assembly — repeats are often collapsed in assemblies, which actually understates
copy number, so real-genome copy number is a floor, not a ceiling.

### 2. Repeat discovery
Run in parallel and merge:
- **RepeatModeler** (or **Red** for speed) → de novo repeat family consensus library.
- **TRF** → tandem repeats / satellite monomers directly (report period, copy number).
- Optionally **k-mer counting** (e.g. very high-frequency k-mers) as an
  orthogonal way to find high-copy motifs without assembling repeat families.

Output: a FASTA of candidate repeat **consensus** sequences + a table of
(family, monomer length, estimated copy number, genomic span).

### 3. Specificity screen (the make-or-break step)
For every candidate consensus:
- `blastn` against each off-target genome (Apis, Varroa, other Tropilaelaps) with a
  deliberately **permissive** setting (short word size, low identity threshold) —
  we want to catch even weak similarity and be conservative about what we call "absent."
- Cross-check with `nucmer` (MUMmer) whole-genome alignment to catch diverged
  homology that blast might miss.
- **Keep only candidates with no meaningful off-target hit.** Define "meaningful"
  explicitly (e.g. no alignment ≥ X bp at ≥ Y% identity, especially none spanning
  the region a primer pair would sit on). Document the threshold.

Also screen against nt/NCBI later for peace of mind, but local genomes are the
authority for the two species that will actually be in the sample.

### 4. Copy-number & conservation ranking
Map candidate back to the tropi assembly to (a) confirm high copy number and (b)
extract all copies and align them. We want the **conserved core** of the repeat —
the stretch that is near-identical across all copies — because primers must sit on
invariant bases to hybridize to every copy in every field population. Rank by:
copy number ↑, core conservation ↑, off-target distance ↑, GC/complexity sane.

### 5. Primer / probe design
Run `primer3` on the conserved core:
- qPCR-friendly: amplicon 70–150 bp, primer Tm ≈ 60 °C, GC 40–60%, avoid runs/hairpins.
- If going probe-based (TaqMan), design a probe on an internal conserved segment.
- Generate several pairs per candidate.

### 6. In-silico PCR
`isPcr` (or primer-blast style):
- Against tropi assembly → must produce the expected product (ideally many, if repeat).
- Against every off-target genome → must produce **nothing**. Any product here kills
  the pair.
Survivors go to `results/candidates/` with their sequences, coordinates, predicted
amplicon, and the full off-target-clearance evidence.

### 7. Wet-lab handoff
Shortlist ~5–10 pairs spanning different candidate families (don't put all eggs in
one repeat). Provide: sequences, expected amplicon size/Tm, positive control
(tropi gDNA), negative controls (Apis, Varroa gDNA, and a no-mite debris sample),
and predicted cross-reactivity notes.

## Open questions to resolve with the user
- Endpoint PCR + gel, or qPCR/TaqMan? (affects amplicon length & probe design)
- Is a single _T. mercedesae_ assembly enough, or do we have/plan multiple
  populations to confirm conservation across the geographic range we're surveilling?
- Do we care about distinguishing _T. mercedesae_ from _T. clareae_, or is
  genus-level "Tropilaelaps present" acceptable for the surveillance question?
- Availability of a _Varroa jacobsoni_ and _Apis cerana_ assembly (Asian-range
  off-targets matter most since that's tropi's home range).
