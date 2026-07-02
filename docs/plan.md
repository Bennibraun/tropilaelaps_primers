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

## Decisions (locked)

- **Assay format:** not decided yet → design **flexibly**. Keep amplicons in the
  70–150 bp qPCR window (works for SYBR/TaqMan/endpoint), and additionally flag a
  conserved internal segment on each survivor as a **candidate TaqMan probe site**,
  so we don't have to redesign if we go probe-based. No format is foreclosed.
- **Specificity level:** **species-specific** — target must be unique to
  _T. mercedesae_ vs. other _Tropilaelaps_ (esp. _T. clareae_). This makes the
  off-target set critical: **other _Tropilaelaps_ genomes are required off-targets**,
  not just Apis/Varroa. If no _T. clareae_ assembly exists publicly, this becomes a
  known gap to flag (may need in-silico or wet-lab congener check on the shortlist).
- **Tropi genomes for conservation:** unsure how many → proceed on the **single
  in-house assembly**, but treat cross-population conservation as an **explicit
  wet-lab validation risk**. Design primers on the most conserved core of each
  repeat family (lowest per-copy variance) to maximize the odds the target holds
  across the surveilled range. Revisit computationally if more assemblies arrive.

## Consequence for the off-target set
Because we committed to species-level specificity, acquiring **any _Tropilaelaps_
congener sequence** (whole genome, or even marker/satellite reads) is now a
priority in Stage 0, alongside Apis and Varroa. Track availability in
`docs/references.md`.

### Reality check (accessions confirmed 2026-07-02)
Off-target genomes exist and are locked in (Apis mellifera GCF_003254395.2,
Apis cerana GCF_029169275.1, Varroa destructor GCF_002443255.2, Varroa jacobsoni
GCF_002532875.2 — all RefSeq). **But there is NO public _T. clareae_ or any
_Tropilaelaps_ congener genome** — only ~31 GenBank records (25 mitochondrial
markers). So:

- Full genome-wide species-vs-congener subtraction is **not possible in-silico**
  with current public data. Honest framing: our pipeline delivers **verified
  exclusion of Apis and Varroa** (which is what actually matters for the debris
  sample) plus **best-effort congener divergence at known marker loci**, with true
  species-level specificity confirmed in the wet lab against _T. clareae_ gDNA.
- Practically this makes the assay **"detects _Tropilaelaps_, validated as
  _T. mercedesae_ at tested loci"** until a congener genome or wet-lab cross-test
  closes the gap. Worth deciding with the user whether that's acceptable for the
  surveillance question, or whether sourcing _T. clareae_ material is warranted.
