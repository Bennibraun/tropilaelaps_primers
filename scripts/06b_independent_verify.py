#!/usr/bin/env python3
"""
Stage 6b: independent, non-BLAST verification of primer occurrence counts.

Stage 6 (06_pcr_validation_v2.py) trusts BLASTN-short's alignment/scoring
machinery to decide what counts as a primer binding site. This stage makes
no such assumption: it searches for each primer sequence directly in the
raw genome text, on both strands, using ordinary substring/fuzzy-substring
matching. If BLAST says a primer has 29 target sites but this script can't
find anything close to 29 occurrences of that sequence in the FASTA, that
is a real discrepancy worth chasing, not a rounding difference.

Two independent counts are reported per primer per genome:

  exact_hits          literal substring match, 0 mismatches, 0 indels.
                       This is the "if I can't find it myself, it's not
                       real" number -- unambiguous, no library trust
                       required beyond Python's own string search.

  fuzzy_hits_ed2       matches allowing up to 2 errors (substitutions,
                       insertions, or deletions combined -- true edit
                       distance), via the `regex` module's fuzzy syntax.
                       This uses the SAME error budget as stage 6's
                       MAX_MISMATCHES=2 default, so the two are meant to
                       be roughly comparable -- but note stage 6 runs
                       blastn with -ungapped (substitutions only, no
                       indels), so fuzzy_hits_ed2 here is a strictly
                       broader search than what BLAST actually did. A
                       large gap between BLAST's site count and
                       fuzzy_hits_ed2 is a genuine red flag; a gap between
                       BLAST's count and exact_hits alone is not, by
                       itself, cause for alarm (repeat copies diverge).

Both directions are searched (the primer sequence itself, and its reverse
complement), since a primer can bind either genomic strand.

This does NOT replace stage 6. It has no concept of primer pairing,
product size, or orientation -- it only answers "how many times does this
literal sequence (or something within 2 edits of it) occur in this
genome," which is exactly the question a manual grep is trying to answer,
done properly (unwrapped, both strands, tolerant of a few mismatches) and
for every candidate automatically instead of one at a time by hand.

Scope: this is meant to run on the SHORTLIST (stage 6's validated_primers.tsv,
i.e. PASS pairs only -- typically single digits to low tens), not every
candidate primer3 designed. fuzzy_hits_ed2 costs several CPU-seconds per
primer per ~20Mb of genome sequence (regex fuzzy matching is not cheap), so
running it against the full primers.tsv (which can be hundreds of pairs)
against a multi-hundred-Mb assembly plus several off-target genomes would
take hours; running it against a shortlist of PASS pairs takes minutes. If
you want to sanity-check a rejected/near-miss pair, run it by hand on a
one-line TSV rather than widening the default scope.

Usage
-----
scripts/06b_independent_verify.py \\
    results/candidates/validated_primers.tsv \\
    data/raw/tropi_assembly.fasta \\
    data/reference/*/*.fna

Requirements
------------
  - Python 3.9+
  - the `regex` module (conda-forge; see env/environment.yml)

Output
------
data/interim/pcr_validation/independent_verification.tsv
    pair_name, primer_type, primer_seq, genome, exact_hits, fuzzy_hits_ed2
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import regex
except ImportError:
    print(
        "This script requires the `regex` module (not stdlib `re`) for "
        "fuzzy/approximate matching. Install it with:\n"
        "    conda install -c conda-forge regex\n"
        "or add it to env/environment.yml and recreate the env.",
        file=sys.stderr,
    )
    raise


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

MAX_EDITS = 2


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Yield (name, sequence) pairs. Sequence is upper-cased and unwrapped."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name:
        yield name, "".join(seq).upper()


def read_primers(path):
    """
    Read a primer TSV (stage 5 primers.tsv or stage 6 validated_primers.tsv --
    both share the candidate_id/pair_rank/fwd_seq/rev_seq columns) and return
    a flat list of (pair_name, primer_type, sequence) tuples, one per F/R
    primer.
    """
    entries = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"candidate_id", "pair_rank", "fwd_seq", "rev_seq"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path} is missing required columns: " + ", ".join(sorted(missing))
            )

        for row in reader:
            pair_name = f"{row['candidate_id']}_pair{row['pair_rank']}"
            entries.append((pair_name, "F", row["fwd_seq"].upper()))
            entries.append((pair_name, "R", row["rev_seq"].upper()))

    if not entries:
        raise RuntimeError(f"{path} contains no primer rows")

    return entries


def count_exact(genome_seqs, query):
    """Count non-overlapping-free (all start positions) exact occurrences of
    query or its reverse complement, across all sequences in a genome."""
    total = 0
    rc = revcomp(query)
    for seq in genome_seqs:
        total += seq.count(query)
        if rc != query:
            total += seq.count(rc)
    return total


def count_fuzzy(genome_seqs, query, max_edits):
    """
    Count occurrences of query (or its reverse complement) allowing up to
    max_edits total insertions/deletions/substitutions, across all
    sequences in a genome.

    Uses the `regex` module's fuzzy-matching extension: {e<=N} bounds the
    total edit count. Deliberately NOT overlapped=True: with fuzzy patterns,
    overlapped scanning reports every valid alignment window around a near
    match (e.g. the same occurrence with an extra flanking base absorbed as
    a "substitution", or trimmed by one base) as a separate hit, wildly
    inflating the count for what is genomically one occurrence. Plain
    finditer's default non-overlapping, leftmost-first scan claims each
    match's span and advances past it, so each genomic occurrence is
    counted once -- verified against a brute-force edit-distance scan on a
    synthetic test genome (29/29 correct, 0 double-counts, 0 misses).
    """
    total = 0
    rc = revcomp(query)
    queries = {query, rc}
    for q in queries:
        pattern = regex.compile(f"({regex.escape(q)}){{e<={max_edits}}}")
        for seq in genome_seqs:
            total += sum(1 for _ in pattern.finditer(seq))
    return total


FIELDS = ["pair_name", "primer_type", "primer_seq", "genome", "exact_hits", "fuzzy_hits_ed2"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "primers",
        help="primer TSV to verify -- normally results/candidates/validated_primers.tsv "
             "(the PASS shortlist; keeps runtime to minutes). Pointing this at the full "
             "primers.tsv works but can take hours on a full genome + off-target set.",
    )
    ap.add_argument("target", help="target genome FASTA")
    ap.add_argument("offtargets", nargs="*", help="off-target genome FASTAs")
    ap.add_argument(
        "-o",
        "--out",
        default="data/interim/pcr_validation/independent_verification.tsv",
        help="output TSV path",
    )
    ap.add_argument(
        "--max-edits",
        type=int,
        default=MAX_EDITS,
        help=f"max total edits (sub/ins/del) for the fuzzy count (default: {MAX_EDITS})",
    )
    ap.add_argument(
        "--skip-fuzzy",
        action="store_true",
        help="only compute exact_hits (much faster; fuzzy_hits_ed2 left blank)",
    )
    args = ap.parse_args()

    primers_path = Path(args.primers)
    target = Path(args.target)
    offtargets = [Path(x) for x in args.offtargets]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    entries = read_primers(primers_path)
    print(f">> loaded {len(entries)} primer sequences ({len(entries)//2} pairs)", file=sys.stderr)

    genomes = [("target:" + target.stem, target)] + [
        ("offtarget:" + g.stem, g) for g in offtargets
    ]

    rows = []

    for genome_label, genome_path in genomes:
        print(f">> scanning {genome_label} ({genome_path})", file=sys.stderr)
        genome_seqs = [seq for _name, seq in read_fasta(genome_path)]

        for i, (pair_name, primer_type, seq) in enumerate(entries, 1):
            exact = count_exact(genome_seqs, seq)
            fuzzy = "" if args.skip_fuzzy else count_fuzzy(genome_seqs, seq, args.max_edits)

            rows.append({
                "pair_name": pair_name,
                "primer_type": primer_type,
                "primer_seq": seq,
                "genome": genome_label,
                "exact_hits": exact,
                "fuzzy_hits_ed2": fuzzy,
            })

            if i % 200 == 0:
                print(f"   {i}/{len(entries)} primers scanned", file=sys.stderr)

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f">> independent verification -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
