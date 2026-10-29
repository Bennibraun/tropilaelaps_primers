# tropi_primers

Find PCR/qPCR (and LAMP) primers that detect *Tropilaelaps mercedesae* DNA in
hive debris (a mix including *Apis* and potentially *Varroa* material) without cross-reacting
with either. The ideal target is a high-copy repeat/satellite unique to tropi,
since more template copies per mite genome means more sensitivity from a trace
amount of debris DNA.

## Genomes

| Role | Species | Accession |
|---|---|---|
| Target | *Tropilaelaps mercedesae* | in-house assembly |
| Off-target | *Apis mellifera* | GCF_003254395.2 |
| Off-target | *Apis cerana* | GCF_029169275.1 |
| Off-target | *Varroa destructor* | GCF_002443255.2 |
| Off-target | *Varroa jacobsoni* | GCF_002532875.2 |

## Scripts

- `00_fetch_references.sh` - download the genomes above.
- `01_assembly_qc.sh` - seqkit stats / sanity check on the tropi assembly.
- `02_repeat_discovery.sh` - RepeatModeler families (external input) + local TRF, merged into one candidate set.
- `02b_kmer_discovery.py` - **WIP**, alternate track: canonical k-mers common in tropi and absent from every off-target, independent of any repeat annotation. Catches repetitive/specific sequence RepeatModeler or TRF miss or misclassify. Not yet run through primer design/validation.
- `03_specificity_screen.sh` / `03b_specificity_screen_reads.sh` / `03c_recut.sh` - drop any candidate with a hit in an off-target genome (or off-target reads, for a congener with only WGS reads); `03c` re-applies the threshold without redoing the blast/minimap runs.
- `04_copy_number_ranking.py` - self-map each candidate back onto the assembly to get real copy number and a conserved core (MSA-based) for primer design.
- `05_primer_design.py` - primer3 on each conserved core.
- `05L_lamp_primer_design.py` - LAMP track: F3/B3/FIP/BIP/LF/LB oligo sets on the same conserved cores.
- `06_pcr_validation_v2.py` - in-silico PCR (BLASTN-short + explicit product enumeration in both primer orientations) against tropi and every off-target.
- `06L_lamp_validation.py` - specificity check for LAMP sets (co-location of binding components; a coarser proxy than PCR product enumeration).
- `06Lb_lamp_independent_verify.py` - independent exact raw-FASTA verification of LAMP component binding and co-location for the clean shortlist.
- `06b_independent_verify.py` - non-BLAST re-check of the shortlist: literal and fuzzy substring search directly against genome text.
- `07_report.py` - builds `results/report.md` from everything above.
- `run_pipeline.sh` - runs stages 4 through 7 in order, skipping any stage whose output already exists.

## Reusable full Snakemake workflow

`Snakefile` runs the complete analysis from a target genome FASTA:

1. Assembly QC and reference acquisition
2. TRF repeat discovery, with optional RepeatModeler families
3. Candidate off-target screening
4. Copy-number ranking and conserved-core extraction
5. PCR primer design and validation, including independent verification
6. LAMP primer design and validation, including independent verification
7. Combined report generation

The default off-target panel is the two *Apis* and two *Varroa* assemblies
listed in `scripts/00_fetch_references.sh`. Configure the target assembly and
optional name in `config/full.yaml`, then run:

```bash
snakemake --cores 8 --configfile config/full.yaml \
  --config target_assembly=/path/to/target.fasta target_name=my_species
```

For another project, copy `config/full.yaml`, replace the target path and name,
and use `--configfile` to select it. Use Snakemake's `--directory` option to
keep each analysis in its own output directory. If a RepeatModeler families FASTA is
available, provide `repeatmodeler_families`; otherwise the workflow uses its
TRF-only candidate discovery path:

```bash
snakemake --directory analyses/my_species --cores 8 \
  --snakefile /path/to/tropilaelaps_primers/Snakefile \
  --configfile /path/to/tropilaelaps_primers/config/full.yaml \
  --config target_assembly=/path/to/my_species.fasta target_name=my_species
```

Each `--directory` is an independent analysis; outputs are written below its
`results/` and `data/interim/` directories.

## Status

- PCR and LAMP tracks (repeat-based, stages 2–7): done, validated, shortlists in `results/`.
- K-mer track (`02b`): candidates generated, not yet run through primer design/validation.
- Species-level specificity vs. *T. clareae*: unverifiable in-silico (no genome exists) - needs wet-lab cross-test.
- Nothing here has been run on a bench yet. Everything is a computational prediction.
