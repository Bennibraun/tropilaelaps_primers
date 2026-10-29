"""Full reusable primer-discovery workflow.

The only required biological input is a target genome FASTA:

    snakemake --cores 8 --config target_assembly=/path/to/target.fasta

Optional:

    --config target_name=my_species

The workflow automatically downloads the default Apis/Varroa off-target panel,
uses TRF-only repeat discovery unless repeatmodeler_families is supplied, and
runs the PCR and LAMP tracks through their independent checks and report.
"""
from pathlib import Path

configfile: "config/full.yaml"

ROOT = Path(workflow.basedir)
if not config.get("target_assembly"):
    raise ValueError(
        "Set the target genome FASTA with "
        "--config target_assembly=/path/to/genome.fasta"
    )
TARGET = str(Path(config["target_assembly"]).resolve())
NAME = config.get("target_name") or Path(TARGET).stem
THREADS = int(config.get("threads", 8))
RM = config.get("repeatmodeler_families", "")
REFS = [
    "data/reference/Apis_cerana.fna",
    "data/reference/Apis_mellifera.fna",
    "data/reference/Varroa_destructor.fna",
    "data/reference/Varroa_jacobsoni.fna",
]
REF_MARKER = "data/reference/.references_ready"
R = "results/candidates"
W = "data/interim"


rule all:
    input:
        "data/interim/qc/seqkit_stats.txt",
        f"{R}/unique_candidates.fasta",
        f"{R}/ranked_candidates.tsv",
        f"{R}/primers.tsv",
        f"{R}/validated_primers.tsv",
        f"{W}/pcr_validation/independent_verification.tsv",
        f"{R}/lamp_primers.tsv",
        f"{R}/lamp_validated.tsv",
        f"{W}/lamp_validation/independent_verification.tsv",
        "results/report.md",


rule fetch_references:
    output:
        refs=REFS,
        marker=REF_MARKER,
    conda: "env/environment.yml"
    shell:
        "mkdir -p data/reference && {ROOT}/scripts/00_fetch_references.sh && touch {output}"


rule assembly_qc:
    input: assembly=TARGET
    output:
        stats="data/interim/qc/seqkit_stats.txt",
        contigs="data/interim/qc/contig_len_gc.tsv",
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        "THREADS={threads} {ROOT}/scripts/01_assembly_qc.sh {input.assembly}"


rule discover_repeats:
    input:
        assembly=TARGET,
        qc="data/interim/qc/seqkit_stats.txt",
    output: candidates="data/interim/repeats/candidates_raw.fasta"
    params: rm=RM
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        r"""
        THREADS={threads} {ROOT}/scripts/02_repeat_discovery.sh {input.assembly} {params.rm}
        test -s {output.candidates}
        """


rule screen_candidate_specificity:
    input:
        candidates="data/interim/repeats/candidates_raw.fasta",
        refs=REFS,
    output:
        unique=f"{R}/unique_candidates.fasta",
        summary=f"{R}/offtarget_summary.tsv",
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        r"""
        THREADS={threads} {ROOT}/scripts/03_specificity_screen.sh {input.candidates} {input.refs}
        """


rule rank_copy_number:
    input:
        assembly=TARGET,
        candidates=f"{R}/unique_candidates.fasta",
    output:
        ranked=f"{R}/ranked_candidates.tsv",
        cores=f"{R}/conserved_cores.fasta",
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        "python {ROOT}/scripts/04_copy_number_ranking.py {input.assembly} {input.candidates} "
        "--outdir {R} --threads {threads}"


rule design_pcr:
    input: cores=f"{R}/conserved_cores.fasta"
    output: primers=f"{R}/primers.tsv"
    conda: "env/environment.yml"
    shell: "python {ROOT}/scripts/05_primer_design.py {input.cores} -o {output.primers}"


rule validate_pcr:
    input:
        primers=f"{R}/primers.tsv",
        target=TARGET,
        refs=REFS,
    output:
        validated=f"{R}/validated_primers.tsv",
        summary=f"{W}/pcr_validation/rejection_summary.tsv",
        products=f"{W}/pcr_validation/offtarget_products.tsv",
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        "THREADS={threads} python {ROOT}/scripts/06_pcr_validation_v2.py {input.primers} "
        "{input.target} {input.refs} -o {output.validated} --work {W}/pcr_validation"


rule independently_verify_pcr:
    input:
        primers=f"{R}/validated_primers.tsv",
        target=TARGET,
        refs=REFS,
    output: verification=f"{W}/pcr_validation/independent_verification.tsv"
    conda: "env/environment.yml"
    shell:
        "python {ROOT}/scripts/06b_independent_verify.py {input.primers} {input.target} "
        "{input.refs} -o {output.verification}"


rule design_lamp:
    input: cores=f"{R}/conserved_cores.fasta"
    output:
        primers=f"{R}/lamp_primers.tsv",
        rejected=f"{R}/lamp_rejected.tsv",
    conda: "env/environment.yml"
    shell:
        "python {ROOT}/scripts/05L_lamp_primer_design.py {input.cores} -o {output.primers} "
        "--rejects {output.rejected}"


rule validate_lamp:
    input:
        primers=f"{R}/lamp_primers.tsv",
        target=TARGET,
        refs=REFS,
    output: validated=f"{R}/lamp_validated.tsv"
    threads: THREADS
    conda: "env/environment.yml"
    shell:
        "THREADS={threads} python {ROOT}/scripts/06L_lamp_validation.py {input.primers} "
        "{input.target} {input.refs} -o {output.validated}"


rule independently_verify_lamp:
    input:
        validated=f"{R}/lamp_validated.tsv",
        designed=f"{R}/lamp_primers.tsv",
        target=TARGET,
        refs=REFS,
    output: verification=f"{W}/lamp_validation/independent_verification.tsv"
    conda: "env/environment.yml"
    shell:
        "python {ROOT}/scripts/06Lb_lamp_independent_verify.py {input.validated} "
        "{input.target} {input.refs} --designed-tsv {input.designed} "
        "-o {output.verification}"


rule report:
    input:
        ranked=f"{R}/ranked_candidates.tsv",
        primers=f"{R}/primers.tsv",
        validated=f"{R}/validated_primers.tsv",
        pcr_summary=f"{W}/pcr_validation/rejection_summary.tsv",
        pcr_products=f"{W}/pcr_validation/offtarget_products.tsv",
        pcr_independent=f"{W}/pcr_validation/independent_verification.tsv",
        lamp=f"{R}/lamp_primers.tsv",
        lamp_validated=f"{R}/lamp_validated.tsv",
    output: report="results/report.md"
    conda: "env/environment.yml"
    shell:
        "python {ROOT}/scripts/07_report.py --candidates-dir {R} "
        "--pcr-work-dir {W}/pcr_validation --title {NAME} -o {output.report}"
