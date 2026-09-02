#!/usr/bin/env python3
"""Stage 5L: LAMP primer set design on each candidate's conserved core.

PARALLEL TRACK to 05_primer_design.py (PCR/qPCR) - both consume the same
conserved_cores.fasta from stage 4. Neither replaces the other; the locked
"format-flexible" decision means we carry both until an assay format is chosen.

LAMP geometry (Notomi et al. 2000). Six regions on the plus strand, in order:

    5'--[F3]--[F2]---[F1]--------[B1c]---[B2c]--[B3c]--3'
             |<------ FIP ----->|       |<--- BIP --->|

  F3   : outer forward primer
  FIP  : F1c + F2   (F1c = reverse complement of F1, joined 5'->3' to F2)
  BIP  : B1c + B2   (B1c is plus-strand B1c region; B2 = revcomp of B2c)
  B3   : outer backward primer (revcomp of B3c)
  LF/LB: optional loop primers (accelerate the reaction ~2x); designed when the
         F2-F1 / B1c-B2c gaps admit a valid oligo.

Because a LAMP set spans ~200-280bp of contiguous conserved sequence (vs ~40bp
of primer footprint for PCR), core length - not copy number - is the binding
constraint. Cores too short for a valid set are reported, not silently dropped.

No standalone LAMP designer exists in bioconda (the 'lamps' package is an
unrelated Qt lab-data GUI; PrimerExplorer V5 is web-only and unscriptable), so
this implements the geometry directly and takes thermodynamics from primer3-py
(SantaLucia nearest-neighbour), the same engine behind stage 5.

Usage:
    scripts/05L_lamp_primer_design.py results/candidates/conserved_cores.fasta
"""
import argparse
import csv
import sys
from pathlib import Path

import primer3

# --- design constraints (Eiken/PrimerExplorer conventions) ---
# Tm targets differ by role: inner primers (F1c/B1c) must fold back first, so
# they run hotter than the outer displacement primers (F3/B3).
TM_INNER = (64.0, 66.0, 68.0)   # min, opt, max  -- F1c/B1c, F2/B2
TM_OUTER = (57.0, 59.0, 61.0)   # min, opt, max  -- F3/B3
TM_LOOP = (63.0, 65.0, 67.0)    # min, opt, max  -- LF/LB

LEN_INNER = (18, 22)   # F1c/B1c, F2/B2 component length
LEN_OUTER = (17, 21)   # F3/B3
LEN_LOOP = (18, 22)

GC_RANGE = (40.0, 65.0)
MAX_POLY_X = 4          # reject runs like AAAAA (LAMP is very prone to these misfiring)
MAX_HAIRPIN_TM = 50.0   # secondary structure ceiling for any single oligo
MAX_DIMER_TM = 45.0     # heterodimer ceiling between primers in a set

# --- spacing constraints (bp between region ends) ---
F2_F1_GAP = (20, 60)     # distance F2 end -> F1 start (loop region, hosts LF)
B1_B2_GAP = (20, 60)     # distance B1c end -> B2c start (hosts LB)
F3_F2_GAP = (0, 20)      # F3 sits just outside F2
B2_B3_GAP = (0, 20)
AMPLICON_MAX = 280       # F3 start -> B3c end; classic LAMP works best <=280bp
AMPLICON_MIN = 180
MAX_COMBINATIONS = 40000  # per-candidate cap on region placements evaluated

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def gc(s):
    if not s:
        return 0.0
    return 100.0 * (s.count("G") + s.count("C") + s.count("g") + s.count("c")) / len(s)


def has_poly_run(s, n=MAX_POLY_X):
    run, prev = 1, ""
    for ch in s.upper():
        if ch == prev:
            run += 1
            if run > n:
                return True
        else:
            run, prev = 1, ch
    return False


_OLIGO_CACHE = {}


def ok_oligo(seq, tm_range, len_range):
    """Basic per-oligo QC: length, GC, poly-run, Tm window, hairpin.

    Memoized: the nested region search evaluates heavily-overlapping windows, so
    the same substring is otherwise re-tested thousands of times (profiling
    showed ~1.5M calls / 1.1M calc_tm for a single 400bp core).
    """
    ck = (seq, tm_range[0], tm_range[2], len_range[0], len_range[1])
    if ck in _OLIGO_CACHE:
        return _OLIGO_CACHE[ck]
    r = _ok_oligo_uncached(seq, tm_range, len_range)
    _OLIGO_CACHE[ck] = r
    return r


def _ok_oligo_uncached(seq, tm_range, len_range):
    if not (len_range[0] <= len(seq) <= len_range[1]):
        return None
    if "N" in seq:
        return None
    if not (GC_RANGE[0] <= gc(seq) <= GC_RANGE[1]):
        return None
    if has_poly_run(seq):
        return None
    tm = primer3.calc_tm(seq)
    if not (tm_range[0] <= tm <= tm_range[2]):
        return None
    try:
        if primer3.calc_hairpin(seq).tm > MAX_HAIRPIN_TM:
            return None
    except Exception:
        pass
    return tm


def best_oligo(template, lo, hi, tm_range, len_range):
    """Scan window template[lo:hi] for the oligo closest to optimal Tm.

    Returns (seq, start, end, tm) or None.
    """
    best = None
    lo = max(0, lo)
    hi = min(len(template), hi)
    for L in range(len_range[0], len_range[1] + 1):
        for s in range(lo, hi - L + 1):
            sub = template[s:s + L]
            tm = ok_oligo(sub, tm_range, len_range)
            if tm is None:
                continue
            score = abs(tm - tm_range[1])
            if best is None or score < best[0]:
                best = (score, sub, s, s + L, tm)
    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def all_oligos(template, lo, hi, tm_range, len_range, limit=40):
    """All viable oligos in template[lo:hi], best-Tm first.

    best_oligo() returns only the single closest-to-optimal oligo, which makes
    the nested placement search brittle: one locally-optimal pick upstream can
    leave no legal placement downstream. LAMP has six interdependent regions, so
    we need to keep alternatives and let the search back off to them.
    """
    lo = max(0, lo)
    hi = min(len(template), hi)
    out = []
    for L_ in range(len_range[0], len_range[1] + 1):
        for s_ in range(lo, hi - L_ + 1):
            sub = template[s_:s_ + L_]
            tm = ok_oligo(sub, tm_range, len_range)
            if tm is None:
                continue
            out.append((abs(tm - tm_range[1]), sub, s_, s_ + L_, tm))
    out.sort(key=lambda x: x[0])
    return [(seq, st, en, tm) for _, seq, st, en, tm in out[:limit]]


_DIMER_CACHE = {}


def _pair_tm(a, b):
    """Cached duplex Tm. The nested search re-tests the same oligo pairs many
    times over; each primer3 duplex call is comparatively expensive."""
    k = (a, b) if a <= b else (b, a)
    if k in _DIMER_CACHE:
        return _DIMER_CACHE[k]
    try:
        tm = (primer3.calc_homodimer(a).tm if a == b
              else primer3.calc_heterodimer(a, b).tm)
    except Exception:
        tm = -100.0
    _DIMER_CACHE[k] = tm
    return tm


def dimer_ok(oligos):
    """Reject sets whose members cross-hybridize. LAMP runs 4-6 primers at once,
    so pairwise heterodimers matter far more than in a 2-primer PCR."""
    names = list(oligos)
    for i in range(len(names)):
        for j in range(i, len(names)):
            a, b = oligos[names[i]], oligos[names[j]]
            if not a or not b:
                continue
            tm = _pair_tm(a, b)
            if tm > MAX_DIMER_TM:
                return False, "%s/%s dimer Tm %.1f" % (names[i], names[j], tm)
    return True, ""


def design_lamp(cand_id, seq):
    """Place F3/F2/F1 ... B1c/B2c/B3c left-to-right and return the best set.

    Search strategy: regions are placed in order, but at each step we keep a
    list of viable oligos (all_oligos) rather than only the best-Tm one, and
    iterate over them. A LAMP set has six interdependent regions packed into
    ~200-280bp, so a greedy single-choice search almost always paints itself
    into a corner - an oligo that is optimal locally can leave no legal
    placement for everything downstream.
    """
    seq = seq.upper()
    _OLIGO_CACHE.clear()
    _DIMER_CACHE.clear()
    n = len(seq)
    if n < AMPLICON_MIN:
        return None, "core %dbp < %dbp minimum for a LAMP set" % (n, AMPLICON_MIN)

    best_set = None
    # windows are generous: each region only needs to START within its window,
    # and we let the Tm filter inside all_oligos do the real selection.
    W = 45
    # Bound the search. The region space is combinatorial and a core with no
    # valid set will otherwise burn minutes proving it. Candidates that admit a
    # set almost always find one early (the oligo lists are Tm-sorted), so a cap
    # costs little sensitivity and turns worst-case minutes into seconds.
    tried = 0

    for f3_seq, f3_s, f3_e, f3_tm in all_oligos(seq, 0, min(60, n), TM_OUTER, LEN_OUTER, limit=12):
        for f2_seq, f2_s, f2_e, f2_tm in all_oligos(
                seq, f3_e + F3_F2_GAP[0], min(f3_e + F3_F2_GAP[1] + W, n),
                TM_INNER, LEN_INNER, limit=12):
            for f1_seq, f1_s, f1_e, f1_tm in all_oligos(
                    seq, f2_e + F2_F1_GAP[0], min(f2_e + F2_F1_GAP[1] + W, n),
                    TM_INNER, LEN_INNER, limit=10):
                for b1c_seq, b1_s, b1_e, b1_tm in all_oligos(
                        seq, f1_e, min(f1_e + 25 + W, n),
                        TM_INNER, LEN_INNER, limit=10):
                    for b2c_seq, b2_s, b2_e, b2_tm in all_oligos(
                            seq, b1_e + B1_B2_GAP[0], min(b1_e + B1_B2_GAP[1] + W, n),
                            TM_INNER, LEN_INNER, limit=10):
                        for b3c_seq, b3_s, b3_e, b3_tm in all_oligos(
                                seq, b2_e + B2_B3_GAP[0], min(b2_e + B2_B3_GAP[1] + W, n),
                                TM_OUTER, LEN_OUTER, limit=8):
                            tried += 1
                            if tried > MAX_COMBINATIONS:
                                if best_set:
                                    return best_set[1], ""
                                return None, ("no valid LAMP set in %dbp core "
                                              "(search cap %d reached)"
                                              % (n, MAX_COMBINATIONS))
                            amp = b3_e - f3_s
                            if not (AMPLICON_MIN <= amp <= AMPLICON_MAX):
                                continue

                            # assemble: FIP = F1c + F2 ; BIP = B1c + B2
                            fip = revcomp(f1_seq) + f2_seq
                            bip = b1c_seq + revcomp(b2c_seq)
                            b3_primer = revcomp(b3c_seq)

                            # loop primers sit in the F2-F1 and B1c-B2c gaps
                            lf_seq = ""
                            if f1_s - f2_e >= LEN_LOOP[0]:
                                lf = all_oligos(seq, f2_e, f1_s, TM_LOOP, LEN_LOOP, limit=1)
                                if lf:
                                    lf_seq = revcomp(lf[0][0])
                            lb_seq = ""
                            if b2_s - b1_e >= LEN_LOOP[0]:
                                lb = all_oligos(seq, b1_e, b2_s, TM_LOOP, LEN_LOOP, limit=1)
                                if lb:
                                    lb_seq = lb[0][0]

                            oligos = {"F3": f3_seq, "B3": b3_primer,
                                      "FIP": fip, "BIP": bip}
                            if lf_seq:
                                oligos["LF"] = lf_seq
                            if lb_seq:
                                oligos["LB"] = lb_seq
                            ok, _why = dimer_ok(oligos)
                            if not ok:
                                continue

                            score = (-(bool(lf_seq) + bool(lb_seq)), amp)
                            row = {
                                "candidate_id": cand_id, "core_len": n,
                                "F3": f3_seq, "F3_tm": round(f3_tm, 1),
                                "B3": b3_primer, "B3_tm": round(b3_tm, 1),
                                "FIP": fip, "FIP_F1c_tm": round(f1_tm, 1),
                                "FIP_F2_tm": round(f2_tm, 1),
                                "BIP": bip, "BIP_B1c_tm": round(b1_tm, 1),
                                "BIP_B2_tm": round(b2_tm, 1),
                                "LF": lf_seq, "LB": lb_seq,
                                "n_loop_primers": bool(lf_seq) + bool(lb_seq),
                                "amplicon_len": amp,
                                "F3_start": f3_s, "B3c_end": b3_e,
                            }
                            if best_set is None or score < best_set[0]:
                                best_set = (score, row)
                            if lf_seq and lb_seq:
                                return best_set[1], ""
    if best_set:
        return best_set[1], ""
    return None, "no valid LAMP set in %dbp core (Tm/GC/spacing/dimer constraints)" % n


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cores_fasta")
    ap.add_argument("-o", "--out", default="results/candidates/lamp_primers.tsv")
    ap.add_argument("--rejects", default="results/candidates/lamp_rejected.tsv")
    args = ap.parse_args()

    fields = ["candidate_id", "core_len", "F3", "F3_tm", "B3", "B3_tm",
              "FIP", "FIP_F1c_tm", "FIP_F2_tm", "BIP", "BIP_B1c_tm", "BIP_B2_tm",
              "LF", "LB", "n_loop_primers", "amplicon_len", "F3_start", "B3c_end"]
    rows, rejects = [], []
    n_cand = 0
    for cand_id, seq in read_fasta(args.cores_fasta):
        n_cand += 1
        row, why = design_lamp(cand_id, seq)
        if row:
            rows.append(row)
        else:
            rejects.append({"candidate_id": cand_id, "core_len": len(seq), "reason": why})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (-r["n_loop_primers"], r["amplicon_len"]))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    rej = Path(args.rejects)
    with open(rej, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["candidate_id", "core_len", "reason"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rejects)

    n_full = sum(1 for r in rows if r["n_loop_primers"] == 2)
    print("%d/%d candidates yielded a valid LAMP set (%d with both loop primers)"
          % (len(rows), n_cand, n_full), file=sys.stderr)
    print("LAMP sets -> %s" % out, file=sys.stderr)
    print("Rejected  -> %s  (mostly cores too short - LAMP needs %d-%dbp vs ~40bp "
          "for a PCR pair)" % (rej, AMPLICON_MIN, AMPLICON_MAX), file=sys.stderr)
    print("Next: scripts/06L_lamp_validation.py %s <assembly> data/reference/*.fna"
          % out, file=sys.stderr)


if __name__ == "__main__":
    main()
