#!/usr/bin/env python3
"""
Stage 2b: repeat-agnostic k-mer based candidate discovery.

PARALLEL TRACK to scripts/02_repeat_discovery.sh (RepeatModeler + TRF), feeding
the same downstream stages (05_primer_design.py, 06_pcr_validation_v2.py) with
a different, complementary discovery method. Neither replaces the other.

Why: RepeatModeler/TRF only find sequence that got *classified* as a repeat
family or tandem array. That classification step can miss real, useful signal
-- diverged satellite variants RepeatModeler's classifier lumps into
"Unspecified" or drops, tandem arrays whose true monomer is shorter than the
window TRF happened to report (see docs/plan.md's stage-4 note), or repetitive
sequence that just never got annotated at all. This script skips annotation
entirely and asks the direct question instead: which exact k-mers are common
in the T. mercedesae assembly and completely absent from every off-target
genome? That is precisely the property a specific, high-copy primer site
needs, independent of whether anything ever called the region a "repeat."

Method
------
  1. Count canonical k-mers (both strands) in the tropi assembly (jellyfish).
  2. Keep the ones repeated >= --min-tropi-copies times (still cheap: a few
     million k-mers out of ~680M genomic positions, not the whole genome).
  3. Count canonical k-mers in each off-target genome, and query the
     high-copy tropi set against each -- keep only k-mers with
     <= --max-offtarget-count hits in EVERY off-target. This queries the small
     high-copy set, not the full genome, against each off-target: cheap.
  4. Find every genomic position of the surviving k-mers with a *single*
     full-genome query pass (one linear scan, not one per off-target) against
     a tiny jellyfish DB built from just the survivors.
  5. Merge nearby surviving positions on each contig into windows (a real
     repeat unit shows up as a run of consecutive/overlapping surviving
     k-mers, not an isolated one). Windows below --min-window-len are dropped
     -- too short to be a useful primer3 template regardless of specificity.
  6. Group windows into families by their highest-copy constituent k-mer (the
     "anchor"): windows sharing an anchor are the same repeat instance seen at
     different genomic loci. One representative window per family becomes the
     candidate sequence -- no MSA/consensus step is needed the way stage 4
     needs one, because family membership here is already defined by exact
     k-mer identity, which is stronger evidence of conservation than an
     alignment consensus.

Output
------
results/candidates/kmer_candidates.fasta   -- one representative sequence per family
results/candidates/kmer_ranked.tsv         -- family_id, n_loci, window_len, anchor_kmer, anchor_tropi_copies

Feed kmer_candidates.fasta into scripts/05_primer_design.py and
scripts/06_pcr_validation_v2.py exactly like conserved_cores.fasta.

Requires: jellyfish (bioconda) on PATH.

Usage
-----
scripts/02b_kmer_discovery.py data/raw/tropi_assembly.fasta data/reference/*.fna
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_K = 20
DEFAULT_MIN_TROPI_COPIES = 5
DEFAULT_MAX_OFFTARGET_COUNT = 0
DEFAULT_MERGE_GAP = 30       # bp; merge surviving k-mer windows into one locus if closer than this
DEFAULT_MIN_WINDOW_LEN = 120  # bp; drop loci too short to be a useful primer3 template
DEFAULT_MAX_WINDOW_LEN = 400  # bp; cap window growth -- see note in the merge loop below


def check_tools():
    for t in ("jellyfish",):
        if shutil.which(t) is None:
            sys.exit(f"Required program '{t}' was not found in PATH (conda install -c bioconda jellyfish)")


def run(cmd, **kw):
    print("  $ " + " ".join(str(x) for x in cmd), file=sys.stderr)
    return subprocess.run(cmd, check=True, text=True, **kw)


def fasta_lengths(path):
    """{record_name: length}, in file order. Also flags whether any N was seen."""
    lengths = {}
    order = []
    name, length = None, 0
    has_n = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                    order.append(name)
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
                if "N" in line or "n" in line:
                    has_n = True
        if name is not None:
            lengths[name] = length
            order.append(name)
    return lengths, order, has_n


def load_fasta(path):
    """{record_name: sequence}. Used once, to extract candidate windows."""
    seqs = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def jellyfish_count(fasta, out_jf, k, hash_size, threads):
    run(["jellyfish", "count", "-m", str(k), "-s", str(hash_size), "-C",
         "-t", str(threads), "-o", str(out_jf), str(fasta)])


def jellyfish_dump_highcopy(jf_path, min_copies):
    """Yield (kmer, count) for k-mers with count >= min_copies."""
    out = subprocess.run(
        ["jellyfish", "dump", "-c", "-t", "-L", str(min_copies), str(jf_path)],
        check=True, text=True, capture_output=True,
    ).stdout
    for line in out.splitlines():
        kmer, count = line.split("\t")
        yield kmer, int(count)


def write_fasta(seqs, path):
    """seqs: iterable of (name, seq)."""
    with open(path, "w") as fh:
        for name, seq in seqs:
            fh.write(f">{name}\n{seq}\n")


def jellyfish_query_stream(query_fasta, jf_path):
    """
    Run `jellyfish query -s query_fasta jf_path`, streaming stdout.

    Yields (kmer, count) in the same order as k-mer windows appear in
    query_fasta (jellyfish's own ordering guarantee, verified empirically --
    one line per valid k-mer window, skipping any window that would contain
    an ambiguous base, no output at all for a record shorter than k).
    """
    proc = subprocess.Popen(
        ["jellyfish", "query", "-s", str(query_fasta), str(jf_path)],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    for line in proc.stdout:
        kmer, count = line.rstrip("\n").split(" ")
        yield kmer, int(count)
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"jellyfish query failed (exit {ret}) on {jf_path}")


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(kmer):
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assembly", help="T. mercedesae assembly FASTA")
    ap.add_argument("offtargets", nargs="+", help="off-target genome FASTAs")
    ap.add_argument("-k", "--kmer-len", type=int, default=DEFAULT_K)
    ap.add_argument("--min-tropi-copies", type=int, default=DEFAULT_MIN_TROPI_COPIES,
                     help="minimum genome-wide count in tropi to keep a k-mer as "
                          "'repetitive' before off-target screening (default: %(default)s)")
    ap.add_argument("--max-offtarget-count", type=int, default=DEFAULT_MAX_OFFTARGET_COUNT,
                     help="max allowed count in any single off-target genome (default: %(default)s, "
                          "i.e. must be completely absent)")
    ap.add_argument("--merge-gap", type=int, default=DEFAULT_MERGE_GAP,
                     help="max bp gap between surviving k-mer positions to merge into one "
                          "locus/window (default: %(default)s)")
    ap.add_argument("--min-window-len", type=int, default=DEFAULT_MIN_WINDOW_LEN,
                     help="drop merged windows shorter than this -- too short to be a "
                          "useful primer3 template regardless of specificity (default: %(default)s)")
    ap.add_argument("--max-window-len", type=int, default=DEFAULT_MAX_WINDOW_LEN,
                     help="cap merged window growth at this length (default: %(default)s). "
                          "Without a cap, a densely-packed tandem satellite (copies spaced well "
                          "under --merge-gap apart) collapses into one multi-hundred-kb blob "
                          "instead of being counted as many loci -- the same undercounting "
                          "failure mode stage 4's MIN_COPY_LEN_FRAC bug had, reintroduced here. "
                          "Capping splits a dense array into many capped windows instead, which "
                          "the anchor-based family grouping (step 6) then correctly tallies as "
                          "many loci of one family.")
    ap.add_argument("--threads", type=int, default=min(32, os.cpu_count() or 4))
    ap.add_argument("--workdir", default="data/interim/kmer_discovery")
    ap.add_argument("--outdir", default="results/candidates")
    args = ap.parse_args()

    check_tools()
    k = args.kmer_len
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    assembly = Path(args.assembly)
    offtargets = [Path(x) for x in args.offtargets]

    # ------------------------------------------------------------------
    # 1. Count k-mers in the tropi assembly.
    # ------------------------------------------------------------------
    tropi_jf = work / "tropi.jf"
    if tropi_jf.exists():
        print(f">> reusing existing {tropi_jf}", file=sys.stderr)
    else:
        print(">> counting k-mers in the tropi assembly", file=sys.stderr)
        # Hash size: genome is ~680Mb, so up to ~680M distinct k-mer slots in
        # the worst case (no repeats at all). Oversize a bit to avoid a costly
        # in-flight resize; RAM is not the constraint here (a few GB either way).
        jellyfish_count(assembly, tropi_jf, k, hash_size=1_000_000_000, threads=args.threads)

    # ------------------------------------------------------------------
    # 2. Keep only the high-copy tropi k-mers -- this is the expensive-ish
    #    genome-wide step done exactly once, and even so it's just a hash
    #    dump, not a genome rescan.
    # ------------------------------------------------------------------
    print(f">> extracting tropi k-mers with count >= {args.min_tropi_copies}", file=sys.stderr)
    tropi_highcopy = dict(jellyfish_dump_highcopy(tropi_jf, args.min_tropi_copies))
    print(f"   {len(tropi_highcopy):,} distinct high-copy k-mers", file=sys.stderr)
    if not tropi_highcopy:
        sys.exit("No k-mers passed --min-tropi-copies -- nothing to screen. "
                  "Lower --min-tropi-copies and retry.")

    highcopy_fa = work / "tropi_highcopy.fa"
    write_fasta(((f"k{i}", kmer) for i, kmer in enumerate(tropi_highcopy)), highcopy_fa)

    # ------------------------------------------------------------------
    # 3. Screen the (small) high-copy set against every off-target genome.
    #    Off-target counting is genome-sized but the genomes here are small
    #    (200-360Mb); querying is proportional to the high-copy set size, not
    #    genome size, so this stays cheap regardless of off-target count.
    # ------------------------------------------------------------------
    surviving = dict(tropi_highcopy)  # kmer -> tropi count; shrinks as off-targets are checked
    for off in offtargets:
        if not surviving:
            break
        name = off.stem
        off_jf = work / f"off_{name}.jf"
        if off_jf.exists():
            print(f">> reusing existing {off_jf}", file=sys.stderr)
        else:
            print(f">> counting k-mers in off-target: {name}", file=sys.stderr)
            jellyfish_count(off, off_jf, k, hash_size=500_000_000, threads=args.threads)

        print(f">> screening {len(surviving):,} candidate k-mers against {name}", file=sys.stderr)
        # Query only the k-mers still alive -- write a fresh fasta each round
        # so a k-mer eliminated by an earlier off-target isn't re-queried.
        alive_fa = work / f"alive_before_{name}.fa"
        write_fasta(((f"k{i}", kmer) for i, kmer in enumerate(surviving)), alive_fa)
        eliminated = 0
        next_surviving = {}
        kmers_in_order = list(surviving.keys())
        for (kmer_seq, count), orig_kmer in zip(
                jellyfish_query_stream(alive_fa, off_jf), kmers_in_order):
            if count <= args.max_offtarget_count:
                next_surviving[orig_kmer] = surviving[orig_kmer]
            else:
                eliminated += 1
        surviving = next_surviving
        alive_fa.unlink(missing_ok=True)
        print(f"   {eliminated:,} eliminated, {len(surviving):,} remain", file=sys.stderr)

    print(f">> {len(surviving):,} k-mers are high-copy in tropi and absent from "
          f"all {len(offtargets)} off-target genome(s)", file=sys.stderr)
    if not surviving:
        sys.exit("No k-mers survived the off-target screen. Nothing to report.")

    # ------------------------------------------------------------------
    # 4. Locate every genomic occurrence of the survivors with ONE full
    #    tropi-genome pass (not one per off-target). Build a tiny jellyfish DB
    #    from just the survivors so the query only matches what we care about.
    # ------------------------------------------------------------------
    survivors_fa = work / "survivors.fa"
    write_fasta(((f"s{i}", kmer) for i, kmer in enumerate(surviving)), survivors_fa)
    survivors_jf = work / "survivors.jf"
    jellyfish_count(survivors_fa, survivors_jf, k,
                     hash_size=max(1000, len(surviving) * 2), threads=1)

    print(">> locating survivor k-mers across the tropi assembly (single full-genome pass)",
          file=sys.stderr)
    lengths, record_order, has_n = fasta_lengths(assembly)
    if has_n:
        print("   WARNING: assembly contains ambiguous (N) bases. jellyfish silently skips "
              "any k-mer window spanning one, which breaks the simple position-by-window-index "
              "bookkeeping this script uses. Coordinates below may be wrong near N runs.",
              file=sys.stderr)

    # positions[contig] = sorted list of 0-based start offsets with a surviving k-mer
    positions = {c: [] for c in record_order}
    rec_iter = iter(record_order)
    cur_rec = next(rec_iter, None)
    cur_windows_left = max(0, lengths[cur_rec] - k + 1) if cur_rec else 0
    cur_offset = 0

    for kmer, count in jellyfish_query_stream(assembly, survivors_jf):
        while cur_windows_left == 0:
            cur_rec = next(rec_iter, None)
            if cur_rec is None:
                raise RuntimeError("jellyfish query produced more windows than the "
                                    "assembly's record lengths account for")
            cur_windows_left = max(0, lengths[cur_rec] - k + 1)
            cur_offset = 0
        if count > 0:
            positions[cur_rec].append(cur_offset)
        cur_offset += 1
        cur_windows_left -= 1

    n_loci_raw = sum(len(v) for v in positions.values())
    print(f"   {n_loci_raw:,} raw survivor positions across {len(record_order):,} contigs",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # 5. Merge nearby positions per contig into windows.
    # ------------------------------------------------------------------
    windows = []  # (contig, start0, end)  end exclusive, covers [start, last_kmer_start+k)
    for contig, starts in positions.items():
        if not starts:
            continue
        starts.sort()
        win_start = starts[0]
        win_last = starts[0]
        for s in starts[1:]:
            if s - win_last <= args.merge_gap and (s + k - win_start) <= args.max_window_len:
                win_last = s
            else:
                windows.append((contig, win_start, win_last + k))
                win_start = s
                win_last = s
        windows.append((contig, win_start, win_last + k))

    windows = [(c, s, e) for c, s, e in windows if (e - s) >= args.min_window_len]
    print(f"   {len(windows):,} merged windows >= {args.min_window_len}bp", file=sys.stderr)
    if not windows:
        sys.exit(f"No windows survived --min-window-len {args.min_window_len}. "
                  "Try lowering it or --merge-gap.")

    # ------------------------------------------------------------------
    # 6. Group windows into families by their highest-copy constituent k-mer
    #    (the "anchor"), extract one representative sequence per family.
    # ------------------------------------------------------------------
    print(">> extracting window sequences and grouping into families", file=sys.stderr)
    genome = load_fasta(assembly)

    family_loci = {}       # anchor_kmer -> list of (contig, start0, end)
    family_best_seq = {}   # anchor_kmer -> (len, seq) of the longest representative seen

    for contig, start0, end in windows:
        seq = genome[contig][start0:end]
        # anchor = the survivor k-mer within this window with the highest tropi copy number
        best_anchor, best_count = None, -1
        for i in range(0, len(seq) - k + 1):
            kmer = canonical(seq[i:i + k])
            c = surviving.get(kmer)
            if c is not None and c > best_count:
                best_anchor, best_count = kmer, c
        if best_anchor is None:
            continue  # shouldn't happen, but don't crash a big run over one edge case
        family_loci.setdefault(best_anchor, []).append((contig, start0, end))
        cur_best = family_best_seq.get(best_anchor)
        if cur_best is None or len(seq) > cur_best[0]:
            family_best_seq[best_anchor] = (len(seq), seq)

    families = sorted(family_loci.items(), key=lambda kv: -len(kv[1]))
    print(f"   {len(families):,} distinct families", file=sys.stderr)

    # ------------------------------------------------------------------
    # Write outputs.
    # ------------------------------------------------------------------
    cand_fasta = outdir / "kmer_candidates.fasta"
    ranked_tsv = outdir / "kmer_ranked.tsv"

    with open(cand_fasta, "w") as ffa, open(ranked_tsv, "w") as ftsv:
        ftsv.write("candidate_id\tn_loci\twindow_len\tanchor_kmer\tanchor_tropi_copies\n")
        for anchor, loci in families:
            n_loci = len(loci)
            seq_len, seq = family_best_seq[anchor]
            cand_id = f"KMER_{anchor}_n{n_loci}"
            ffa.write(f">{cand_id}\n{seq}\n")
            ftsv.write(f"{cand_id}\t{n_loci}\t{seq_len}\t{anchor}\t{surviving[anchor]}\n")

    print(f"\nCandidates -> {cand_fasta}", file=sys.stderr)
    print(f"Ranked table -> {ranked_tsv}", file=sys.stderr)
    print(f"Next: scripts/05_primer_design.py {cand_fasta} -o {outdir}/kmer_primers.tsv",
          file=sys.stderr)
    print(f"      scripts/06_pcr_validation_v2.py {outdir}/kmer_primers.tsv "
          f"{assembly} {' '.join(str(o) for o in offtargets)} "
          f"-o {outdir}/kmer_validated.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
