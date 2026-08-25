#!/usr/bin/env python3
"""Stage 5: primer3 design on each candidate's conserved core.

qPCR-friendly settings per docs/plan.md: amplicon 70-150 bp, primer Tm ~60C,
GC 40-60%, avoid homopolymer runs/hairpins. Also designs an internal oligo
(candidate TaqMan probe site) on every core, per the locked "format-flexible"
decision — no assay format is foreclosed, so we flag a probe candidate even
though the format isn't chosen yet.

If a core is too short for a 70-150bp product with flanking primers, primer3
simply returns zero pairs for it and it's skipped — no special-casing needed.

Usage:
    scripts/05_primer_design.py results/candidates/conserved_cores.fasta

Requires Python package: primer3-py (>=2.0; uses the design_primers(seq_args,
global_args) two-dict API — if your installed version only has the older
camelCase designPrimers(...), update primer3-py rather than patching this).
"""
import argparse
import csv
import sys
from pathlib import Path

import primer3

GLOBAL_ARGS = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_SELF_ANY": 8,
    "PRIMER_MAX_SELF_END": 3,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]],
    "PRIMER_NUM_RETURN": 3,
    "PRIMER_PICK_INTERNAL_OLIGO": 1,
    "PRIMER_INTERNAL_OPT_TM": 68.0,
    "PRIMER_INTERNAL_MIN_TM": 65.0,
    "PRIMER_INTERNAL_MAX_TM": 72.0,
    "PRIMER_INTERNAL_MAX_POLY_X": 4,
}


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


def design_for_candidate(cand_id, seq):
    seq_args = {"SEQUENCE_ID": cand_id, "SEQUENCE_TEMPLATE": seq}
    result = primer3.bindings.design_primers(seq_args, GLOBAL_ARGS)
    n = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    rows = []
    for i in range(n):
        row = {
            "candidate_id": cand_id,
            "pair_rank": i,
            "fwd_seq": result[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "fwd_tm": round(result[f"PRIMER_LEFT_{i}_TM"], 1),
            "fwd_gc": round(result[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1),
            "rev_seq": result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "rev_tm": round(result[f"PRIMER_RIGHT_{i}_TM"], 1),
            "rev_gc": round(result[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
            "product_size": result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "probe_seq": result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", ""),
            "probe_tm": round(result[f"PRIMER_INTERNAL_{i}_TM"], 1) if f"PRIMER_INTERNAL_{i}_TM" in result else "",
        }
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cores_fasta")
    ap.add_argument("-o", "--out", default="results/candidates/primers.tsv")
    args = ap.parse_args()

    fieldnames = ["candidate_id", "pair_rank", "fwd_seq", "fwd_tm", "fwd_gc",
                  "rev_seq", "rev_tm", "rev_gc", "product_size", "probe_seq", "probe_tm"]
    all_rows = []
    n_candidates = 0
    n_with_pairs = 0
    for cand_id, seq in read_fasta(args.cores_fasta):
        n_candidates += 1
        rows = design_for_candidate(cand_id, seq)
        if rows:
            n_with_pairs += 1
        all_rows.extend(rows)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    print(f"{n_with_pairs}/{n_candidates} candidates yielded at least one primer pair "
          f"(short/uninformative cores are skipped, not an error)", file=sys.stderr)
    print(f"Primer pairs -> {out_path}", file=sys.stderr)
    print(f"Next: scripts/06_ispcr_validation.sh {out_path} <tropi_assembly> data/reference/*/*.fna", file=sys.stderr)


if __name__ == "__main__":
    main()
