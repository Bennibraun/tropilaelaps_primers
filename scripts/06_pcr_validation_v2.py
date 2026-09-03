#!/usr/bin/env python3
"""
Stage 6: explicit in-silico PCR validation.

For every primer pair from Stage 5:

  TARGET GENOME
    - must produce >=1 PCR product
    - product must be 70-150 bp
    - multiple products are accepted (PASS_HIGH_COPY) only when they are all
      the same size within MULTI_COPY_MAX_SIZE_RANGE bp -- the signature of a
      real high-copy repeat amplifying cleanly at every copy. Multiple
      products of visibly different sizes are flagged as
      TARGET_MULTIPLE_PRODUCTS and rejected (nonspecific/multi-locus).

  OFF-TARGET GENOMES
    - must produce ZERO PCR-compatible products
    - any plausible product rejects the pair

Instead of UCSC isPcr, this script:

  1. Uses BLASTN-short to map every primer independently.
  2. Keeps full-length, sufficiently good primer alignments.
  3. Explicitly pairs forward/reverse binding sites on each contig.
  4. Applies PCR orientation and product-size constraints itself.

This avoids isPcr's seed enumeration / coordinate-bin behavior on
repeat-rich genomes.

Primer binding criteria
-----------------------
By default:

  - alignment covers the complete primer
  - >=90% sequence identity
  - <=2 mismatches total
  - <=1 mismatch in the final 5 bases at the primer's 3' end
  - no gaps

These are configurable from the command line.

Usage
-----
scripts/06_pcr_validation.py \
    results/candidates/primers.tsv \
    data/raw/tropi_assembly.fasta \
    data/reference/*/*.fna

Requirements
------------
  - BLAST+ (blastn, makeblastdb)
  - Python 3.9+

Outputs
-------
results/candidates/validated_primers.tsv

data/interim/pcr_validation/
    target_binding_sites.tsv
    target_products.tsv
    offtarget_binding_sites.tsv
    offtarget_products.tsv
    rejection_summary.tsv
    primer_binding_summary.tsv
"""


import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

MIN_PRODUCT = 70
MAX_PRODUCT = 150

# Primer alignment criteria.
MIN_IDENTITY = 90.0
MAX_MISMATCHES = 2
MAX_3P_MISMATCHES = 1
THREE_PRIME_WINDOW = 5

# A primer occurring at enormous numbers of sites is not a useful assay
# primer, regardless of whether some of those sites happen to pair.
MAX_PRIMER_SITES = 150

# Multiple target products are only accepted (PASS_HIGH_COPY) when they are
# all within this many bp of each other -- see the classification note below.
# 5bp tolerates a little indel variation between diverged repeat copies while
# still rejecting the genuinely scattered-size multi-locus failure mode.
MULTI_COPY_MAX_SIZE_RANGE = 5

# BLAST word size. 7 is appropriate for short primers.
BLAST_WORD_SIZE = 7


# ---------------------------------------------------------------------------
# FASTA / utilities
# ---------------------------------------------------------------------------

def fasta_lengths(path):
    """Return {sequence_name: length} for a FASTA file."""
    lengths = {}
    name = None
    length = 0

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())

    if name is not None:
        lengths[name] = length

    return lengths


def check_program(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required program '{name}' was not found in PATH"
        )
    return path


def run_command(cmd, *, stdout=None, stderr=None):
    """Run command and raise a useful error if it fails."""
    print("  $ " + " ".join(str(x) for x in cmd), file=sys.stderr)

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=stdout,
            stderr=stderr,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Command failed with exit status {e.returncode}: "
            f"{' '.join(str(x) for x in cmd)}"
        ) from e


# ---------------------------------------------------------------------------
# Primer input
# ---------------------------------------------------------------------------

def read_primers(path):
    """
    Read Stage 5 primers.tsv.

    Expected columns:
      candidate_id
      pair_rank
      fwd_seq
      ...
      rev_seq
      ...
      product_size
      ...
    """
    primers = []

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {
            "candidate_id",
            "pair_rank",
            "fwd_seq",
            "rev_seq",
            "product_size",
        }

        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path} is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            pair_name = f"{row['candidate_id']}_pair{row['pair_rank']}"

            primers.append({
                "pair_name": pair_name,
                "candidate_id": row["candidate_id"],
                "pair_rank": row["pair_rank"],
                "fwd_seq": row["fwd_seq"].upper(),
                "rev_seq": row["rev_seq"].upper(),
                "designed_product_size": row["product_size"],
                "row": row,
            })

    if not primers:
        raise RuntimeError(
            f"{path} contains no primer rows"
        )

    return primers


def write_primer_fasta(primers, path):
    """
    Write each primer as a BLAST query.

    The reverse primer is deliberately queried in its normal primer
    5'->3' orientation. BLAST reports whether the genomic hit is on
    the forward or reverse strand.
    """
    with open(path, "w") as fh:
        for p in primers:
            fh.write(f">{p['pair_name']}__F\n{p['fwd_seq']}\n")
            fh.write(f">{p['pair_name']}__R\n{p['rev_seq']}\n")


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def make_blast_db(genome, db_prefix):
    """Create a nucleotide BLAST database."""
    nsq = Path(str(db_prefix) + ".nsq")

    if nsq.exists():
        return

    print(
        f"  building BLAST database for {genome}",
        file=sys.stderr,
    )

    run_command([
        "makeblastdb",
        "-in", str(genome),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ])


def run_blast(
    primer_fasta,
    db_prefix,
    output,
):
    """
    Map primers against a genome.

    We request:
      qseqid
      sseqid
      pident
      length
      mismatch
      gapopen
      qstart
      qend
      sstart
      send
      evalue
      bitscore
      qlen
      qseq
      sseq

    qseq/sseq are retained so that the 3' mismatch criterion can be
    evaluated explicitly.
    """

    outfmt = (
        "6 "
        "qseqid "
        "sseqid "
        "pident "
        "length "
        "mismatch "
        "gapopen "
        "qstart "
        "qend "
        "sstart "
        "send "
        "evalue "
        "bitscore "
        "qlen "
        "qseq "
        "sseq"
    )

    with open(output, "w") as out:
        run_command(
            [
                "blastn",
                "-query", str(primer_fasta),
                "-db", str(db_prefix),
                "-task", "blastn-short",
                "-word_size", str(BLAST_WORD_SIZE),
                "-ungapped",
                "-dust", "no",
                "-outfmt", outfmt,
            ],
            stdout=out,
        )


# ---------------------------------------------------------------------------
# BLAST parsing / primer-site filtering
# ---------------------------------------------------------------------------

def count_3prime_mismatches(qseq, sseq, window):
    qseq = qseq.upper()
    sseq = sseq.upper()

    if len(qseq) != len(sseq):
        return None

    window = min(window, len(qseq))

    return sum(
        a != b
        for a, b in zip(qseq[-window:], sseq[-window:])
    )

def parse_blast_hits(
    blast_file,
    primer_lookup,
    min_identity,
    max_mismatches,
    max_3p_mismatches,
    three_prime_window,
):
    """
    Parse BLAST output and return:

        sites[primer_name] = [
            site dictionaries...
        ]

    Only full-length, ungapped, sufficiently good alignments are retained.
    """

    sites = defaultdict(list)

    with open(blast_file) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) != 15:
                raise RuntimeError(
                    f"Malformed BLAST line {line_number} in {blast_file}: "
                    f"expected 15 columns, got {len(fields)}"
                )

            (
                qseqid,
                sseqid,
                pident,
                length,
                mismatch,
                gapopen,
                qstart,
                qend,
                sstart,
                send,
                evalue,
                bitscore,
                qlen,
                qseq,
                sseq,
            ) = fields

            pident = float(pident)
            length = int(length)
            mismatches = int(mismatch)
            gapopen = int(gapopen)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            qlen = int(qlen)

            if qseqid not in primer_lookup:
                continue

            # Full-length query coverage is mandatory.
            if qstart != 1 or qend != qlen:
                continue

            if length != qlen:
                continue

            if gapopen != 0:
                continue

            if pident < min_identity:
                continue

            if mismatches > max_mismatches:
                continue

            p3_mismatches = count_3prime_mismatches(
                qseq,
                sseq,
                three_prime_window,
            )

            if p3_mismatches is None:
                continue

            if p3_mismatches > max_3p_mismatches:
                continue

            if sstart <= send:
                strand = "+"
                start = sstart
                end = send
            else:
                strand = "-"
                start = send
                end = sstart

            sites[qseqid].append({
                "primer_name": qseqid,
                "contig": sseqid,
                "start": start,
                "end": end,
                "strand": strand,
                "identity": pident,
                "mismatches": mismatches,
                "mismatches_3p": p3_mismatches,
                "length": length,
                "evalue": evalue,
                "bitscore": bitscore,
            })

    return sites


# ---------------------------------------------------------------------------
# PCR product enumeration
# ---------------------------------------------------------------------------

def enumerate_products(
    f_sites,
    r_sites,
    min_product,
    max_product,
):
    """
    Find all PCR-compatible F/R site combinations.

    Required geometry (either orientation is a valid amplicon, since a
    repeat copy can be inserted on either genomic strand):

        F --->                <--- R      (copy on the "+" reference strand)
        R --->                <--- F      (copy inverted relative to reference)

    Both sites must be on the same contig.

    Product coordinates include both primer sequences.

    Returns a list of product dictionaries.
    """

    by_contig_f = defaultdict(list)
    by_contig_r = defaultdict(list)

    for site in f_sites:
        by_contig_f[site["contig"]].append(site)

    for site in r_sites:
        by_contig_r[site["contig"]].append(site)

    products = []

    for contig in set(by_contig_f) & set(by_contig_r):
        forwards = by_contig_f[contig]
        reverses = by_contig_r[contig]

        # Sort to make the search deterministic.
        forwards.sort(key=lambda x: x["start"])
        reverses.sort(key=lambda x: x["start"])

        for f in forwards:
            for r in reverses:

                # Two valid geometries produce a real amplicon:
                #  1. F on "+", R on "-", R downstream of F (copy in the
                #     same orientation as the reference strand).
                #  2. R on "+", F on "-", F downstream of R (copy inserted
                #     in the opposite orientation).
                # Checking only (1) silently drops every product at a
                # repeat copy inserted in reverse orientation.
                if f["strand"] == "+" and r["strand"] == "-":
                    upstream, downstream = f, r
                elif r["strand"] == "+" and f["strand"] == "-":
                    upstream, downstream = r, f
                else:
                    continue

                if downstream["start"] <= upstream["start"]:
                    continue

                product_size = downstream["end"] - upstream["start"] + 1

                if product_size < min_product:
                    continue

                if product_size > max_product:
                    continue

                products.append({
                    "contig": contig,
                    "start": upstream["start"],
                    "end": downstream["end"],
                    "product_size": product_size,

                    "f_identity": f["identity"],
                    "f_mismatches": f["mismatches"],
                    "f_3p_mismatches": f["mismatches_3p"],

                    "r_identity": r["identity"],
                    "r_mismatches": r["mismatches"],
                    "r_3p_mismatches": r["mismatches_3p"],
                })

    return products


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

BINDING_FIELDS = [
    "pair_name",
    "primer",
    "primer_type",
    "contig",
    "start",
    "end",
    "strand",
    "identity",
    "mismatches",
    "mismatches_3p",
    "length",
    "evalue",
    "bitscore",
]


def write_binding_sites(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=BINDING_FIELDS,
            delimiter="\t",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


PRODUCT_FIELDS = [
    "pair_name",
    "genome",
    "contig",
    "start",
    "end",
    "product_size",
    "f_identity",
    "f_mismatches",
    "f_3p_mismatches",
    "r_identity",
    "r_mismatches",
    "r_3p_mismatches",
]


def write_products(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=PRODUCT_FIELDS,
            delimiter="\t",
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------

def validate_genome(
    genome_name,
    genome,
    primers,
    primer_lookup,
    work,
    min_identity,
    max_mismatches,
    max_3p_mismatches,
    min_product,
    max_product,
    max_primer_sites,
    three_prime_window,
):
    """
    BLAST all primers against one genome and enumerate PCR products.

    Returns:
        sites
        products_by_pair
        overused_primers
    """

    db_prefix = work / f"db_{genome_name}"
    primer_fasta = work / "all_primers.fa"
    blast_file = work / f"{genome_name}.blast.tsv"

    make_blast_db(genome, db_prefix)

    print(
        f"  mapping primers against {genome_name}",
        file=sys.stderr,
    )

    if blast_file.exists():
        print(f"  using existing BLAST output: {blast_file}")
    else:
        run_blast(
            primer_fasta,
            db_prefix,
            blast_file,
        )

    sites = parse_blast_hits(
        blast_file,
        primer_lookup,
        min_identity,
        max_mismatches,
        max_3p_mismatches,
        three_prime_window,
    )

    # Map query names to their actual pair and primer type.
    pair_sites = defaultdict(lambda: {"F": [], "R": []})

    for query_name, hit_sites in sites.items():
        pair_name, primer_type = query_name.rsplit("__", 1)

        for site in hit_sites:
            pair_sites[pair_name][primer_type].append(site)

    # Reject individual primers that have an excessive number of sites.
    # Importantly, this happens before product enumeration.
    overused = set()

    for pair_name, typed_sites in pair_sites.items():
        for primer_type in ("F", "R"):
            n = len(typed_sites[primer_type])

            if n > max_primer_sites:
                overused.add((pair_name, primer_type))

    products_by_pair = defaultdict(list)

    for primer in primers:
        pair_name = primer["pair_name"]

        f_overused = (pair_name, "F") in overused
        r_overused = (pair_name, "R") in overused

        if f_overused or r_overused:
            continue

        typed = pair_sites[pair_name]

        products = enumerate_products(
            typed["F"],
            typed["R"],
            min_product,
            max_product,
        )

        for product in products:
            product["pair_name"] = pair_name
            product["genome"] = genome_name

        products_by_pair[pair_name].extend(products)

    return sites, products_by_pair, overused


def binding_rows_for_genome(
    genome_name,
    primers,
    sites,
):
    rows = []

    for primer in primers:
        pair_name = primer["pair_name"]

        for primer_type, seq in (
            ("F", primer["fwd_seq"]),
            ("R", primer["rev_seq"]),
        ):
            query_name = f"{pair_name}__{primer_type}"

            for site in sites.get(query_name, []):
                rows.append({
                    "pair_name": pair_name,
                    "primer": query_name,
                    "primer_type": primer_type,
                    "contig": site["contig"],
                    "start": site["start"],
                    "end": site["end"],
                    "strand": site["strand"],
                    "identity": site["identity"],
                    "mismatches": site["mismatches"],
                    "mismatches_3p": site["mismatches_3p"],
                    "length": site["length"],
                    "evalue": site["evalue"],
                    "bitscore": site["bitscore"],
                })

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument(
        "primers",
        help="Stage 5 primers.tsv",
    )

    ap.add_argument(
        "target",
        help="target genome FASTA",
    )

    ap.add_argument(
        "offtargets",
        nargs="*",
        help="off-target genome FASTAs",
    )

    ap.add_argument(
        "-o",
        "--out",
        default="results/candidates/validated_primers.tsv",
        help="validated primer output",
    )

    ap.add_argument(
        "--work",
        default="data/interim/pcr_validation",
        help="working directory",
    )

    ap.add_argument(
        "--min-product",
        type=int,
        default=MIN_PRODUCT,
        help=f"minimum PCR product size (default: {MIN_PRODUCT})",
    )

    ap.add_argument(
        "--max-product",
        type=int,
        default=MAX_PRODUCT,
        help=f"maximum PCR product size (default: {MAX_PRODUCT})",
    )

    ap.add_argument(
        "--min-identity",
        type=float,
        default=MIN_IDENTITY,
        help=f"minimum primer identity (default: {MIN_IDENTITY})",
    )

    ap.add_argument(
        "--max-mismatches",
        type=int,
        default=MAX_MISMATCHES,
        help=f"maximum total primer mismatches (default: {MAX_MISMATCHES})",
    )

    ap.add_argument(
        "--max-3p-mismatches",
        type=int,
        default=MAX_3P_MISMATCHES,
        help=f"maximum mismatches in 3' window (default: {MAX_3P_MISMATCHES})",
    )

    ap.add_argument(
        "--three-prime-window",
        type=int,
        default=THREE_PRIME_WINDOW,
        help=f"3' mismatch window size (default: {THREE_PRIME_WINDOW})",
    )

    ap.add_argument(
        "--max-primer-sites",
        type=int,
        default=MAX_PRIMER_SITES,
        help=f"maximum genomic sites per primer (default: {MAX_PRIMER_SITES})",
    )

    args = ap.parse_args()


    check_program("blastn")
    check_program("makeblastdb")

    three_prime_window=args.three_prime_window
    primers_path = Path(args.primers)
    target = Path(args.target)
    offtargets = [Path(x) for x in args.offtargets]
    work = Path(args.work)
    out_path = Path(args.out)

    work.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    primers = read_primers(primers_path)

    print(
        f">> loaded {len(primers)} primer pairs",
        file=sys.stderr,
    )

    # ------------------------------------------------------------------
    # Write all primers once. Same query file is used for every genome.
    # ------------------------------------------------------------------

    primer_fasta = work / "all_primers.fa"
    write_primer_fasta(primers, primer_fasta)

    primer_lookup = {}

    for p in primers:
        primer_lookup[f"{p['pair_name']}__F"] = p
        primer_lookup[f"{p['pair_name']}__R"] = p

    # ------------------------------------------------------------------
    # TARGET
    # ------------------------------------------------------------------

    print(
        ">> validating against target genome",
        file=sys.stderr,
    )

    target_name = target.stem

    (
        target_sites,
        target_products,
        target_overused,
    ) = validate_genome(
        target_name,
        target,
        primers,
        primer_lookup,
        work,
        args.min_identity,
        args.max_mismatches,
        args.max_3p_mismatches,
        args.min_product,
        args.max_product,
        args.max_primer_sites,
        args.three_prime_window,
    )

    target_binding_rows = binding_rows_for_genome(
        target_name,
        primers,
        target_sites,
    )

    write_binding_sites(
        work / "target_binding_sites.tsv",
        target_binding_rows,
    )

    target_product_rows = []

    for pair_name, products in target_products.items():
        target_product_rows.extend(products)

    write_products(
        work / "target_products.tsv",
        target_product_rows,
    )

    # ------------------------------------------------------------------
    # OFF-TARGETS
    # ------------------------------------------------------------------

    all_offtarget_products = []
    all_offtarget_binding_rows = []

    # pair -> total number of off-target products
    off_product_counts = defaultdict(int)

    # pair -> genomes with products
    off_product_genomes = defaultdict(set)

    for offtarget in offtargets:

        name = offtarget.stem

        print(
            f">> validating against off-target: {name}",
            file=sys.stderr,
        )

        (
            sites,
            products,
            overused,
        ) = validate_genome(
            name,
            offtarget,
            primers,
            primer_lookup,
            work,
            args.min_identity,
            args.max_mismatches,
            args.max_3p_mismatches,
            args.min_product,
            args.max_product,
            args.max_primer_sites,
            args.three_prime_window,
        )

        all_offtarget_binding_rows.extend(
            binding_rows_for_genome(
                name,
                primers,
                sites,
            )
        )

        for pair_name, product_list in products.items():
            if not product_list:
                # products_by_pair always has an entry per pair (even an
                # empty one, from the unconditional .extend() in
                # validate_genome), so skipping empty lists here is required
                # -- otherwise every off-target genome gets recorded as a
                # cross-reactivity hit for every pair, regardless of whether
                # any product was actually found.
                continue
            off_product_counts[pair_name] += len(product_list)
            off_product_genomes[pair_name].add(name)
            all_offtarget_products.extend(product_list)

    write_binding_sites(
        work / "offtarget_binding_sites.tsv",
        all_offtarget_binding_rows,
    )

    write_products(
        work / "offtarget_products.tsv",
        all_offtarget_products,
    )

    # ------------------------------------------------------------------
    # Classification
    # ------------------------------------------------------------------

    classifications = {}

    for p in primers:
        pair = p["pair_name"]

        n_target_sites_f = sum(
            1
            for s in target_sites.get(f"{pair}__F", [])
        )

        n_target_sites_r = sum(
            1
            for s in target_sites.get(f"{pair}__R", [])
        )

        n_target_products = len(
            target_products.get(pair, [])
        )

        n_off_products = off_product_counts.get(pair, 0)

        f_overused = (pair, "F") in target_overused
        r_overused = (pair, "R") in target_overused

        pair_products = target_products.get(pair, [])
        product_sizes = [pr["product_size"] for pr in pair_products]
        size_range = (max(product_sizes) - min(product_sizes)) if product_sizes else 0

        if f_overused or r_overused:
            status = "REJECT_HIGH_COPY_TARGET"

        elif n_target_products == 0:
            status = "REJECT_NO_TARGET_PRODUCT"

        elif n_off_products > 0:
            status = "REJECT_OFFTARGET"

        elif n_target_products == 1:
            status = "PASS"

        elif size_range <= MULTI_COPY_MAX_SIZE_RANGE:
            # Many products, but all (near-)identical in size: this is the
            # signature of a good high-copy repeat assay (design goal, see
            # docs/plan.md) -- one clean band/Ct from many template copies,
            # not a smear. A small range is tolerated because divergent copies
            # of a real repeat can carry a base or two of indel variation
            # between the primer sites; a pair whose products vary widely in
            # size is the actual multi-locus/nonspecific failure mode this
            # rejects instead.
            status = "PASS_HIGH_COPY"

        else:
            status = "TARGET_MULTIPLE_PRODUCTS"

        classifications[pair] = {
            "status": status,
            "target_f_sites": n_target_sites_f,
            "target_r_sites": n_target_sites_r,
            "target_products": n_target_products,
            "target_product_size_range": size_range,
            "offtarget_products": n_off_products,
            "offtarget_genomes": ",".join(
                sorted(off_product_genomes.get(pair, set()))
            ),
        }

    # ------------------------------------------------------------------
    # Write rejection summary
    # ------------------------------------------------------------------

    summary_path = work / "rejection_summary.tsv"

    summary_fields = [
        "pair_name",
        "candidate_id",
        "pair_rank",
        "designed_product_size",
        "target_f_sites",
        "target_r_sites",
        "target_products",
        "target_product_size_range",
        "offtarget_products",
        "offtarget_genomes",
        "status",
    ]

    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=summary_fields,
            delimiter="\t",
        )
        writer.writeheader()

        for p in primers:
            pair = p["pair_name"]
            c = classifications[pair]

            writer.writerow({
                "pair_name": pair,
                "candidate_id": p["candidate_id"],
                "pair_rank": p["pair_rank"],
                "designed_product_size": p["designed_product_size"],
                **c,
            })

    # ------------------------------------------------------------------
    # Write primer binding summary
    # ------------------------------------------------------------------

    binding_summary_path = work / "primer_binding_summary.tsv"

    with open(binding_summary_path, "w", newline="") as fh:
        fields = [
            "pair_name",
            "primer_type",
            "target_sites",
            "target_overused",
        ]

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()

        for p in primers:
            pair = p["pair_name"]

            for primer_type in ("F", "R"):
                query = f"{pair}__{primer_type}"

                n = len(target_sites.get(query, []))
                overused = (pair, primer_type) in target_overused

                writer.writerow({
                    "pair_name": pair,
                    "primer_type": primer_type,
                    "target_sites": n,
                    "target_overused": int(overused),
                })

    # ------------------------------------------------------------------
    # Write final validated primer table
    # ------------------------------------------------------------------

    validated_fields = list(primers[0]["row"].keys()) + [
        "validation_status",
        "target_f_sites",
        "target_r_sites",
        "target_products",
        "target_product_size_range",
        "offtarget_products",
        "offtarget_genomes",
    ]

    PASS_STATUSES = {"PASS", "PASS_HIGH_COPY"}

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=validated_fields,
            delimiter="\t",
        )
        writer.writeheader()

        for p in primers:
            pair = p["pair_name"]
            c = classifications[pair]

            # Only actual PASS pairs (single-product, or same-size multi-copy)
            # enter the final validated table.
            if c["status"] not in PASS_STATUSES:
                continue

            row = dict(p["row"])

            row.update({
                "validation_status": c["status"],
                "target_f_sites": c["target_f_sites"],
                "target_r_sites": c["target_r_sites"],
                "target_products": c["target_products"],
                "target_product_size_range": c["target_product_size_range"],
                "offtarget_products": c["offtarget_products"],
                "offtarget_genomes": c["offtarget_genomes"],
            })

            writer.writerow(row)

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------

    counts = defaultdict(int)

    for c in classifications.values():
        counts[c["status"]] += 1

    print(file=sys.stderr)
    print(">> PCR validation complete", file=sys.stderr)

    for status in (
        "PASS",
        "PASS_HIGH_COPY",
        "TARGET_MULTIPLE_PRODUCTS",
        "REJECT_NO_TARGET_PRODUCT",
        "REJECT_OFFTARGET",
        "REJECT_HIGH_COPY_TARGET",
    ):
        print(
            f"   {status}: {counts[status]}",
            file=sys.stderr,
        )

    print(
        f"   validated pairs -> {out_path}",
        file=sys.stderr,
    )

    print(
        f"   binding sites   -> {work / 'target_binding_sites.tsv'}",
        file=sys.stderr,
    )

    print(
        f"   target products -> {work / 'target_products.tsv'}",
        file=sys.stderr,
    )

    print(
        f"   off-target products -> {work / 'offtarget_products.tsv'}",
        file=sys.stderr,
    )

    print(
        f"   rejection summary  -> {summary_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()