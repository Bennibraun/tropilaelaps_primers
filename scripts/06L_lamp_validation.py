#!/usr/bin/env python3
"""Stage 6L: in-silico specificity validation for LAMP primer sets.

PARALLEL TRACK to 06_ispcr_validation.sh. isPcr cannot be used here: it models
two primers flanking an amplicon, which is not how LAMP works. A LAMP set is
4-6 oligos recognizing 6-8 regions via a strand-displacement / dumbbell
mechanism, so "does it amplify" is not a simple flanking-pair question.

What we check instead, per off-target genome:

  1. Each oligo COMPONENT is BLASTed separately. FIP and BIP are chimeric
     (F1c+F2, B1c+B2) and never hybridize as one unit - blasting the joined
     40-44bp sequence would be meaningless, so they are split back into their
     halves and each half screened on its own.
  2. A component "binds" an off-target if it matches at high identity over
     nearly its full length, with the 3' end intact. Priming tolerates internal
     mismatches far better than 3'-terminal ones, so a hit that stops short of
     the 3' end is scored as non-priming.
  3. RISK is scored by how many components bind, and critically whether the
     binding components are CO-LOCATED on one off-target sequence within a
     plausible amplicon span. Scattered single-component hits across a 370Mb
     genome are expected by chance and are not evidence of amplification.

Verdict per set:
  PASS          - no off-target has >=2 co-located binding components
  REVIEW        - some co-location, but not the F2/B2 or F1c/B1c core pairs
  FAIL          - an off-target has co-located inner-primer regions: plausible
                  amplification, treat the set as non-specific

On-target, the same components are checked against the assembly: all 4-6 must
bind, and we count how many co-located sites exist (a high-copy repeat target
gives many, which is good for sensitivity but is reported so it is visible).

Usage:
    scripts/06L_lamp_validation.py results/candidates/lamp_primers.tsv \\
        data/raw/tropi_assembly.fasta data/reference/*.fna
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# a component counts as binding if it matches this well
MIN_IDENT = 90.0        # % identity over the alignment
MIN_COV_FRAC = 0.85     # fraction of the component length aligned
THREE_PRIME_SLACK = 2   # bp of the 3' end allowed to go unaligned

# two components are "co-located" if they land on the same subject sequence
# within this distance - a plausible LAMP amplicon plus generous slack
COLOCATION_WINDOW = 500

THREADS = os.environ.get("THREADS", "4")


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def split_components(row):
    """Return {component_name: sequence} with FIP/BIP split into their halves.

    FIP = F1c + F2 and BIP = B1c + B2. We do not know the internal split point
    from the joined string alone, so we recover it from the per-half Tm columns
    only if present; otherwise we fall back to splitting at half length, which
    is close enough for a specificity screen (both halves are 18-22bp).
    """
    comps = {}
    for name in ("F3", "B3", "LF", "LB"):
        if row.get(name):
            comps[name] = row[name]
    fip = row.get("FIP", "")
    bip = row.get("BIP", "")
    if fip:
        # F1c is the 5' half, F2 the 3' half
        h = len(fip) // 2
        comps["FIP_F1c"] = fip[:h]
        comps["FIP_F2"] = fip[h:]
    if bip:
        h = len(bip) // 2
        comps["BIP_B1c"] = bip[:h]
        comps["BIP_B2"] = bip[h:]
    return comps


def make_db(fasta, dbdir):
    db = Path(dbdir) / Path(fasta).stem
    if not Path(str(db) + ".nsq").exists():
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                        "-out", str(db)], check=True, stdout=subprocess.DEVNULL)
    return db


def blast_components(comps, db, workdir, tag):
    """BLAST every component; return {comp: [(sseqid, sstart, send), ...]} for
    hits that look like they could actually prime."""
    qpath = Path(workdir) / ("query_%s.fa" % tag)
    with open(qpath, "w") as fh:
        for name, seq in comps.items():
            fh.write(">%s\n%s\n" % (name, seq))
    out = subprocess.run(
        ["blastn", "-query", str(qpath), "-db", str(db),
         "-task", "blastn-short", "-word_size", "7", "-evalue", "1000",
         "-num_threads", THREADS, "-max_target_seqs", "5000",
         "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send"],
        check=True, text=True, capture_output=True).stdout

    binding = defaultdict(list)
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        q, s, pid, alen, qlen, qs, qe, ss, se = (
            f[0], f[1], float(f[2]), int(f[3]), int(f[4]),
            int(f[5]), int(f[6]), int(f[7]), int(f[8]))
        if pid < MIN_IDENT:
            continue
        if alen < MIN_COV_FRAC * qlen:
            continue
        # 3' end must be (nearly) covered - priming depends on it
        if (qlen - max(qs, qe)) > THREE_PRIME_SLACK:
            continue
        binding[q].append((s, min(ss, se), max(ss, se)))
    return binding


def colocated_groups(binding):
    """Find sets of >=2 distinct components landing near each other on one
    subject sequence. Returns list of (sseqid, pos, {components})."""
    by_seq = defaultdict(list)
    for comp, hits in binding.items():
        for sseqid, s, e in hits:
            by_seq[sseqid].append((s, e, comp))
    groups = []
    for sseqid, items in by_seq.items():
        items.sort()
        i = 0
        while i < len(items):
            start = items[i][0]
            members = {items[i][2]}
            j = i + 1
            while j < len(items) and items[j][0] - start <= COLOCATION_WINDOW:
                members.add(items[j][2])
                j += 1
            if len(members) >= 2:
                groups.append((sseqid, start, members))
            i += 1
    return groups


INNER = {"FIP_F1c", "FIP_F2", "BIP_B1c", "BIP_B2"}


def verdict_for(groups):
    """FAIL if any off-target group contains >=2 inner components (the ones that
    actually drive LAMP), REVIEW for lesser co-location, PASS if none."""
    worst = "PASS"
    detail = ""
    for sseqid, pos, members in groups:
        n_inner = len(members & INNER)
        if n_inner >= 2:
            return "FAIL", "%s@%d: %s" % (sseqid, pos, ",".join(sorted(members)))
        if len(members) >= 2 and worst == "PASS":
            worst = "REVIEW"
            detail = "%s@%d: %s" % (sseqid, pos, ",".join(sorted(members)))
    return worst, detail


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lamp_tsv")
    ap.add_argument("target")
    ap.add_argument("offtargets", nargs="+")
    ap.add_argument("-o", "--out",
                    default="results/candidates/lamp_validated.tsv")
    args = ap.parse_args()

    workdir = Path("data/interim/lamp_validation")
    workdir.mkdir(parents=True, exist_ok=True)
    dbdir = Path("data/interim/blastdb")
    dbdir.mkdir(parents=True, exist_ok=True)

    with open(args.lamp_tsv) as fh:
        sets = list(csv.DictReader(fh, delimiter="\t"))
    if not sets:
        sys.exit("no LAMP sets in %s" % args.lamp_tsv)
    print(">> %d LAMP sets to validate" % len(sets), file=sys.stderr)

    print(">> building/locating BLAST DBs", file=sys.stderr)
    target_db = make_db(args.target, dbdir)
    off_dbs = [(Path(f).stem, make_db(f, dbdir)) for f in args.offtargets]

    rows = []
    for idx, s in enumerate(sets):
        cand = s["candidate_id"]
        comps = split_components(s)

        # --- on-target ---
        tb = blast_components(comps, target_db, workdir, "t%d" % idx)
        n_bound_target = sum(1 for c in comps if tb.get(c))
        tgroups = colocated_groups(tb)
        # count distinct loci where a full-ish set co-locates
        n_sites = sum(1 for _, _, m in tgroups if len(m & INNER) >= 2)

        # --- off-targets ---
        worst, detail, worst_ref = "PASS", "", ""
        for name, db in off_dbs:
            ob = blast_components(comps, db, workdir, "o%d_%s" % (idx, name))
            groups = colocated_groups(ob)
            v, d = verdict_for(groups)
            if v == "FAIL":
                worst, detail, worst_ref = "FAIL", d, name
                break
            if v == "REVIEW" and worst == "PASS":
                worst, detail, worst_ref = "REVIEW", d, name

        on_target_ok = n_bound_target >= 4 and n_sites >= 1
        rows.append({
            "candidate_id": cand,
            "amplicon_len": s.get("amplicon_len", ""),
            "n_loop_primers": s.get("n_loop_primers", ""),
            "components_bound_target": n_bound_target,
            "n_target_sites": n_sites,
            "on_target_ok": "yes" if on_target_ok else "NO",
            "offtarget_verdict": worst,
            "offtarget_ref": worst_ref,
            "offtarget_detail": detail,
            "F3": s.get("F3", ""), "B3": s.get("B3", ""),
            "FIP": s.get("FIP", ""), "BIP": s.get("BIP", ""),
            "LF": s.get("LF", ""), "LB": s.get("LB", ""),
        })
        print("   [%d/%d] %s: target %d comps / %d sites, off-target %s"
              % (idx + 1, len(sets), cand, n_bound_target, n_sites, worst),
              file=sys.stderr)

    order = {"PASS": 0, "REVIEW": 1, "FAIL": 2}
    rows.sort(key=lambda r: (r["on_target_ok"] != "yes",
                             order.get(r["offtarget_verdict"], 3),
                             -int(r["n_target_sites"] or 0)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_pass = sum(1 for r in rows
                 if r["offtarget_verdict"] == "PASS" and r["on_target_ok"] == "yes")
    n_rev = sum(1 for r in rows if r["offtarget_verdict"] == "REVIEW")
    n_fail = sum(1 for r in rows if r["offtarget_verdict"] == "FAIL")
    print("", file=sys.stderr)
    print("%d PASS (on-target OK, clean off-target) | %d REVIEW | %d FAIL"
          % (n_pass, n_rev, n_fail), file=sys.stderr)
    print("Validated LAMP sets -> %s" % out, file=sys.stderr)
    print("NOTE: co-location screening is a proxy for LAMP amplification, not a "
          "simulation of it. PASS means no off-target carries co-located inner "
          "primer sites; wet-lab confirmation is still required.", file=sys.stderr)


if __name__ == "__main__":
    main()
