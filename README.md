# Genome primer discovery pipeline

This repository contains a reusable Snakemake workflow for discovering
repeat-derived PCR and LAMP primer candidates from any genome assembly.

## Inputs

Configure these values in `config/full.yaml` or copy
`config/full.example.yaml` and edit it:

- `target_assembly`: target genome FASTA.
- `off_targets`: optional list of FASTA files that must not amplify.
- `target_name`: label used in the report.
- `output_dir`: directory for all generated data and results.
- `repeatmodeler_families`: optional RepeatModeler family FASTA.

The workflow does not download genomes or assume a species-specific reference
panel.

## Run

```bash
snakemake --cores 8 \
  --snakefile /path/to/genome-primer-pipeline/Snakefile \
  --configfile /path/to/config.yaml
```

Alternatively, override values on the command line:

```bash
snakemake --cores 8 --snakefile Snakefile --configfile config/full.yaml \
  --config target_assembly=/path/to/target.fasta \
  target_name=my_genome \
  off_targets='["/path/to/negative1.fasta", "/path/to/negative2.fasta"]' \
  output_dir=/path/to/analysis
```

All generated files are written below `output_dir`:

- `data/interim/` for intermediate files.
- `results/candidates/` for candidate and primer tables.
- `results/report.md` for the final report.

Use Snakemake's `--directory` instead of `output_dir` when preferred. Do not
commit input genomes, intermediate files, results, logs, or indexes.

## Stages

1. Assembly quality control.
2. Repeat discovery using TRF and optional RepeatModeler families.
3. Candidate screening against configured off-target genomes.
4. Copy-number ranking and conserved-core extraction.
5. PCR primer design and validation.
6. LAMP primer design and validation.
7. Independent verification and report generation.

The supported entrypoint is `Snakefile`; individual scripts are retained for
advanced or partial-stage use.
