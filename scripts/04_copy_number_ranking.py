#!/usr/bin/env python3
"""Stage 4: copy-number & conservation ranking.

For each candidate that survived the specificity screen (stage 3/3b), map it
back onto the T. mercedesae assembly to (a) confirm/count real genomic copies
and (b) extract those copies and find their conserved core — the near-invariant
stretch primers must sit on to hybridize to every copy in every field
population (see docs/plan.md, Stage 4).

Written in Python (not bash, unlike the other stages) because it needs to parse
a multiple-sequence alignment and compute per-column conservation — awk isn't a
good fit for that part.

Usage:
    scripts/04_copy_number_ranking.py <assembly.fasta> <unique_candidates.fasta>

Requires on PATH: makeblastdb, blastn, seqkit, mafft
Requires Python package: biopython
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import AlignIO

MIN_COPY_IDENT = 80.0   # % identity to count a blast hit as a real copy (looser than the
                         # stage-3 off-target screen, since we WANT to find diverged same-family copies)
MIN_COPY_LEN_FRAC = 0.5 # hit must cover at least this fraction of the candidate length
MAX_COPIES_FOR_MSA = 50 # cap copies fed to mafft per candidate; keeps this laptop-feasible
                         # on high-copy satellites (true copy number is still reported in full)
MAX_HITS_PER_CANDIDATE = 500  # blastn -max_target_seqs; see blast_copies(). Candidates
                              # at this ceiling have n_copies reported as a floor.
CORE_MIN_IDENTITY = 0.90  # per-column agreement fraction required for a base to count as "core"
CORE_MIN_LEN = 20         # minimum core window length to bother reporting

REQUIRED_TOOLS = ["makeblastdb", "blastn", "seqkit", "mafft"]


def default_threads():
    """Thread count from the Slurm allocation, else the machine's CPUs, else 4.
    Under sbatch, $SLURM_CPUS_PER_TASK reflects `-c/--cpus-per-task`, so the job
    uses exactly what it reserved without anyone exporting an env var."""
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit():
        return int(slurm)
    return os.cpu_count() or 4


def check_tools():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        sys.exit(f"Missing required tools on PATH: {', '.join(missing)}")


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, **kw)


def seq_lengths(fasta_path):
    out = subprocess.run(["seqkit", "fx2tab", "-nl", str(fasta_path)],
                          check=True, text=True, capture_output=True).stdout
    lengths = {}
    for line in out.splitlines():
        name, length = line.rstrip("\n").split("\t")
        lengths[name] = int(length)
    return lengths


def build_self_blastdb(assembly, workdir):
    db = workdir / "assembly_db"
    if not (workdir / "assembly_db.nsq").exists():
        run(["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db)],
            stdout=subprocess.DEVNULL)
    return db


def blast_copies(candidates, db, workdir, threads):
    out = workdir / "self_hits.tsv"
    # Reuse a completed hit table if one is already present and non-empty. BLAST
    # is the fixed, expensive part of this stage (the loop is what we optimize
    # with --top-n), so a restart after a killed/timed-out run should not redo it.
    # To force a fresh BLAST, delete self_hits.tsv first.
    if out.exists() and out.stat().st_size > 0:
        print(f">> reusing existing BLAST hits: {out} "
              f"({out.stat().st_size} bytes) — delete it to force a fresh run",
              file=sys.stderr)
        return out
    # -max_target_seqs caps hits per candidate. Without it, self-BLASTing a
    # repeat library against a repeat-rich genome is effectively quadratic: a
    # high-copy family matches thousands of loci, and an unbounded run produced a
    # 4GB hit table that had not finished. We only need enough copies to (a)
    # estimate copy number and (b) fill an MSA capped at MAX_COPIES_FOR_MSA, so
    # a generous cap loses nothing downstream. Reported n_copies becomes a floor
    # for families that hit the cap — noted in the output.
    run([
        "blastn", "-query", str(candidates), "-db", str(db),
        "-task", "blastn", "-word_size", "11",
        "-perc_identity", str(MIN_COPY_IDENT),
        "-num_threads", str(threads),
        "-max_target_seqs", str(MAX_HITS_PER_CANDIDATE),
        "-dust", "yes",
        "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send",
    ], stdout=open(out, "w"))
    return out


# --- locus clustering -------------------------------------------------------
#
# A single HSP covering >=MIN_COPY_LEN_FRAC of the query used to be required
# for a hit to count as a "copy" (see git history / docs/plan.md notes). That
# silently dropped every genuinely diverged copy of a candidate: RepeatModeler
# consensus sequences are built by averaging an alignment of many real copies,
# so an individual copy that diverges from the consensus in different places
# often only aligns to it in fragments -- each fragment covers well under 50%
# of the query on its own, even though the copy is real and the fragments
# together cover most of it. The fix: cluster HSPs into genomic loci first,
# then judge each *locus* (union of its fragments' query coverage) against the
# length-fraction threshold, instead of judging each HSP in isolation.
#
# The tricky part is not merging together truly distinct copies. A tandem
# satellite has many real copies sitting right next to each other in the
# genome with ~zero gap between them -- genomic adjacency alone can't tell
# "one copy split by an indel" from "two adjacent copies". The discriminator
# used here is the QUERY coordinates: fragments of one diverged copy cover
# *different, largely non-overlapping* parts of the query (the parts that
# still match, on either side of wherever it diverged); adjacent hits from
# genuinely separate tandem copies instead each re-cover the *same* region of
# the query (the same repeat unit), so their query ranges overlap heavily.
# Only the first pattern gets merged.
LOCUS_MAX_GENOMIC_GAP = 1000  # bp; max genomic gap between HSPs to treat as one insertion
LOCUS_MAX_QUERY_OVERLAP_FRAC = 0.2  # max reciprocal query-range overlap to still merge


def _query_overlap_frac(a, b):
    """Fraction of the shorter [qstart,qend] interval that a and b share."""
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    overlap = max(0, hi - lo + 1)
    shorter = min(a[1] - a[0] + 1, b[1] - b[0] + 1)
    return overlap / shorter if shorter else 0.0


def cluster_loci(raw_hits, qlen):
    """
    Group per-candidate HSPs into genomic loci and keep the ones whose
    combined (unioned) query coverage reaches MIN_COPY_LEN_FRAC of the query.

    raw_hits: list of (chrom, qstart, qend, start0, end, strand), 1-based
              qstart/qend, 0-based half-open start0/end.

    Returns: list of (chrom, start0, end, strand) loci -- same shape
    write_bed_and_extract() already expects.
    """
    loci = []

    by_group = defaultdict(list)
    for chrom, qstart, qend, start0, end, strand in raw_hits:
        by_group[(chrom, strand)].append((start0, end, qstart, qend))

    for (chrom, strand), group_hits in by_group.items():
        group_hits.sort(key=lambda h: h[0])

        clusters = []  # each: {"start0", "end", "qspans": [(qstart,qend), ...]}
        for start0, end, qstart, qend in group_hits:
            placed = False
            for cl in clusters:
                gap = start0 - cl["end"]
                if gap > LOCUS_MAX_GENOMIC_GAP:
                    continue
                overlaps_existing = any(
                    _query_overlap_frac((qstart, qend), qs) > LOCUS_MAX_QUERY_OVERLAP_FRAC
                    for qs in cl["qspans"]
                )
                if overlaps_existing:
                    continue
                cl["start0"] = min(cl["start0"], start0)
                cl["end"] = max(cl["end"], end)
                cl["qspans"].append((qstart, qend))
                placed = True
                break
            if not placed:
                clusters.append({"start0": start0, "end": end, "qspans": [(qstart, qend)]})

        for cl in clusters:
            qspans = sorted(cl["qspans"])
            covered = 0
            cur_lo, cur_hi = qspans[0]
            for lo, hi in qspans[1:]:
                if lo <= cur_hi + 1:
                    cur_hi = max(cur_hi, hi)
                else:
                    covered += cur_hi - cur_lo + 1
                    cur_lo, cur_hi = lo, hi
            covered += cur_hi - cur_lo + 1

            if covered >= MIN_COPY_LEN_FRAC * qlen:
                loci.append((chrom, cl["start0"], cl["end"], strand))

    return loci


def parse_hits(hits_path):
    """candidate -> list of raw (chrom, qstart, qend, start0, end, strand) HSPs."""
    hits = defaultdict(list)
    with open(hits_path) as fh:
        for line in fh:
            qseqid, sseqid, pident, length, qlen, qstart, qend, sstart, send = (
                line.rstrip("\n").split("\t")
            )
            qstart, qend = int(qstart), int(qend)
            s, e = int(sstart), int(send)
            strand = "+" if s <= e else "-"
            start0, end = (s - 1, e) if s <= e else (e - 1, s)
            hits[qseqid].append((sseqid, qstart, qend, start0, end, strand))
    return hits


def write_bed_and_extract(cand_id, copies, assembly, workdir):
    bed = workdir / f"{cand_id}.copies.bed"
    with open(bed, "w") as fh:
        for i, (chrom, start0, end, strand) in enumerate(copies[:MAX_COPIES_FOR_MSA]):
            fh.write(f"{chrom}\t{start0}\t{end}\t{cand_id}__copy{i}\t0\t{strand}\n")
    fasta = workdir / f"{cand_id}.copies.fasta"
    run(["seqkit", "subseq", "--bed", str(bed), str(assembly)], stdout=open(fasta, "w"))
    return fasta


def align_and_find_core(copies_fasta, single_seq, workdir, cand_id, threads):
    n = sum(1 for _ in open(copies_fasta) if _.startswith(">"))
    if n <= 1:
        return single_seq, 1.0
    aln_path = workdir / f"{cand_id}.aln.fasta"
    # Fast, threaded alignment. `--auto` chooses an accurate but O(N^2 L^2)
    # strategy (L-INS-i) for these repeat monomers, which is what made the
    # per-candidate loop overrun the 11h wall silently. `--retree 1` pins the
    # fast progressive method (FFT-NS-2); copies of one repeat family are near-
    # identical, so the accurate refinement buys nothing here.
    run(["mafft", "--retree", "1", "--thread", str(threads), "--quiet", str(copies_fasta)],
        stdout=open(aln_path, "w"))
    aln = AlignIO.read(aln_path, "fasta")
    ncols = aln.get_alignment_length()

    # Per column: agreement fraction AND the majority (consensus) base. We need
    # the consensus base per column, not a single copy's base, so the reported
    # core is genuinely conserved — see the core-extraction note below.
    agreement = []
    consensus_bases = []
    for col in range(ncols):
        bases = [rec.seq[col].upper() for rec in aln if rec.seq[col] != "-"]
        if not bases:
            agreement.append(0.0)
            consensus_bases.append("-")
            continue
        counts = defaultdict(int)
        for b in bases:
            counts[b] += 1
        top_base, top_count = max(counts.items(), key=lambda kv: kv[1])
        agreement.append(top_count / len(bases))
        consensus_bases.append(top_base)

    # longest run of columns meeting CORE_MIN_IDENTITY
    best_start, best_len, cur_start, cur_len = 0, 0, 0, 0
    for i, a in enumerate(agreement):
        if a >= CORE_MIN_IDENTITY:
            if cur_len == 0:
                cur_start = i
            cur_len += 1
            if cur_len > best_len:
                best_start, best_len = cur_start, cur_len
        else:
            cur_len = 0
    if best_len < CORE_MIN_LEN:
        return None, max(agreement) if agreement else 0.0

    # Core = column-wise CONSENSUS across the conserved window, NOT aln[0].
    # Using aln[0] (one arbitrary copy) sliced by alignment columns could return
    # bases where that copy disagrees with the family, letting primers sit on
    # non-conserved positions. The consensus is the invariant target we want.
    core_seq = "".join(consensus_bases[best_start:best_start + best_len]).replace("-", "")
    mean_ident = sum(agreement[best_start:best_start + best_len]) / best_len
    return core_seq, mean_ident


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("assembly")
    ap.add_argument("candidates")
    ap.add_argument("--outdir", default="results/candidates")
    ap.add_argument("--top-n", type=int, default=None,
                    help="Run the expensive MSA/conserved-core step only on the "
                         "top N candidates by copy number (BLAST-derived, cheap). "
                         "The full copy-number table still covers all candidates; "
                         "candidates outside the top N are marked 'core not computed "
                         "(outside --top-n)'. Omit to process every candidate.")
    ap.add_argument("--threads", type=int, default=default_threads(),
                    help="Threads for blastn and mafft. Defaults to the Slurm "
                         "allocation ($SLURM_CPUS_PER_TASK) if set, else the CPU "
                         "count, else 4.")
    args = ap.parse_args()
    threads = str(args.threads)

    check_tools()
    assembly = Path(args.assembly)
    candidates = Path(args.candidates)
    workdir = Path("data/interim/copy_number")
    workdir.mkdir(parents=True, exist_ok=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lengths = seq_lengths(candidates)
    seqs = {}
    name = None
    for line in open(candidates):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line

    print(">> building self-blastdb from assembly", file=sys.stderr)
    db = build_self_blastdb(assembly, workdir)
    print(">> mapping candidates back onto the assembly", file=sys.stderr)
    hits_path = blast_copies(candidates, db, workdir, threads)
    raw_hits = parse_hits(hits_path)

    # Copy number is BLAST-derived and cheap, so compute it for ALL candidates
    # first. The expensive part (per-candidate MAFFT + core) is what blew the
    # wall clock, so when --top-n is set we run it only on the highest-copy
    # families — which is exactly what a high-copy PCR target needs anyway.
    #
    # Cluster fragmented HSPs into loci before counting -- see cluster_loci()
    # docstring. This is the expensive-ish step here (pure Python, no
    # subprocess), but still cheap relative to BLAST/MAFFT.
    print(">> clustering fragmented BLAST hits into genomic loci", file=sys.stderr)
    hits = {cid: cluster_loci(raw_hits.get(cid, []), lengths[cid]) for cid in lengths}
    copy_counts = {cid: len(hits.get(cid, [])) for cid in lengths}
    if args.top_n is not None:
        # candidates with >=1 copy, most copies first; ties broken by longer candidate
        ranked_ids = sorted(
            (cid for cid in lengths if copy_counts[cid] > 0),
            key=lambda cid: (-copy_counts[cid], -lengths[cid]),
        )
        selected = set(ranked_ids[:args.top_n])
        print(f">> --top-n {args.top_n}: running MSA/core on {len(selected)} of "
              f"{len(lengths)} candidates (highest copy number); "
              f"copy-number table still covers all.", file=sys.stderr)
    else:
        selected = None  # process everything

    rows = []
    core_records = []
    total = len(lengths)
    for idx, (cand_id, length) in enumerate(lengths.items(), 1):
        copies = hits.get(cand_id, [])
        n_copies = len(copies)
        # Progress so a long MSA loop is never a silent black box (the 11h
        # timeout printed nothing after the blast line). Flush so tee/Slurm logs
        # show it live.
        print(f"[{idx}/{total}] {cand_id}: {n_copies} copies",
              file=sys.stderr, flush=True)
        if n_copies == 0:
            rows.append({"candidate_id": cand_id, "n_copies": 0, "core_len": 0,
                         "core_identity": "", "note": "no locus reached the coverage threshold (unique single-copy candidate)"})
            core_records.append((cand_id, seqs[cand_id]))
            continue

        # Skip the expensive core step for candidates outside the top-N, but keep
        # their copy count in the table so nothing is silently dropped.
        if selected is not None and cand_id not in selected:
            rows.append({"candidate_id": cand_id, "n_copies": n_copies, "core_len": 0,
                         "core_identity": "", "note": "core not computed (outside --top-n)"})
            continue

        copies_fasta = write_bed_and_extract(cand_id, copies, assembly, workdir)
        core_seq, core_ident = align_and_find_core(copies_fasta, seqs[cand_id], workdir, cand_id, threads)
        if core_seq is None:
            rows.append({"candidate_id": cand_id, "n_copies": n_copies, "core_len": 0,
                         "core_identity": f"{core_ident:.3f}",
                         "note": "no conserved core >= threshold — copies too divergent"})
            continue
        capped = " (n_copies is a floor: hit cap reached)" if n_copies >= MAX_HITS_PER_CANDIDATE else ""
        rows.append({"candidate_id": cand_id, "n_copies": n_copies, "core_len": len(core_seq),
                     "core_identity": f"{core_ident:.3f}", "note": capped.strip()})
        core_records.append((cand_id, core_seq))

    rows.sort(key=lambda r: (-r["n_copies"], -(float(r["core_identity"]) if r["core_identity"] else 0)))

    ranked_path = outdir / "ranked_candidates.tsv"
    with open(ranked_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["candidate_id", "n_copies", "core_len", "core_identity", "note"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    cores_path = outdir / "conserved_cores.fasta"
    with open(cores_path, "w") as fh:
        for cand_id, seq in core_records:
            if seq:
                fh.write(f">{cand_id}\n{seq}\n")

    print(f"Ranked candidates -> {ranked_path}", file=sys.stderr)
    print(f"Conserved cores   -> {cores_path}", file=sys.stderr)
    print("Note: candidates whose copies are too divergent to yield a conserved core "
          "are ranked but excluded from conserved_cores.fasta (nothing to design primers on).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
