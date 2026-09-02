#!/usr/bin/env python3
"""
Stage 7: build a human-readable summary report of the full pipeline run.

Pulls together the per-stage tables that stages 4-6 already write to
results/candidates/ and produces one Markdown report with:

  - a pipeline funnel (candidates in/out at each stage)
  - a ranked shortlist of PASS primer pairs, joined back to their
    copy-number / conservation evidence from stage 4
  - flagged pairs that amplify the target but were rejected (off-target
    hit, multi-product, or high-copy-on-target) so you can see what was
    close but didn't make it
  - per-candidate off-target hit counts, broken out by genome, so you can
    see *which* off-target genome (if any) came close

This script does not re-run BLAST or primer3. It only reads the TSVs
stages 4-6 already wrote. If a file is missing, that section is skipped
with a note rather than the whole report failing -- useful if you're
generating a report mid-pipeline (e.g. LAMP stages not run).

Usage
-----
scripts/07_report.py [--candidates-dir results/candidates] [--out results/report.md]

Inputs (all optional; missing ones are noted, not fatal)
----------------------------------------------------------
results/candidates/ranked_candidates.tsv       (stage 4)
results/candidates/primers.tsv                 (stage 5)
results/candidates/validated_primers.tsv       (stage 6)
data/interim/pcr_validation/rejection_summary.tsv   (stage 6, all pairs incl. rejects)
data/interim/pcr_validation/offtarget_products.tsv  (stage 6, per-genome detail)

Output
------
results/report.md
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    """Return list[dict] for a TSV, or None if the file doesn't exist."""
    if not path.exists():
        return None
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_int(x, default="?"):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return default


def fmt_pct(x, default="?"):
    """Format a 0-1 fraction (as written by stage 4's core_identity column) as a percent."""
    try:
        return f"{float(x) * 100:.1f}%"
    except (TypeError, ValueError):
        return default


def md_table(rows, columns, headers=None):
    """Render a list of dicts as a Markdown table over the given columns."""
    headers = headers or columns
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(str(row.get(c, "")) for c in columns) + " |"
        )
    return "\n".join(lines)


def build_funnel(ranked, primers, validated_all):
    """Candidate/pair counts at each stage, for a quick top-of-report summary."""
    lines = []

    if ranked is not None:
        n_cand = len(ranked)
        n_repeat = sum(1 for r in ranked if fmt_int(r.get("n_copies"), "0") != "0"
                       and r.get("n_copies", "0") != "0")
        n_core = sum(1 for r in ranked if r.get("core_len") not in (None, "", "0"))
        lines.append(f"- **Stage 4 (copy-number ranking):** {n_cand:,} candidates scored; "
                      f"{n_core:,} had a conserved core computed.")
    else:
        lines.append("- **Stage 4:** ranked_candidates.tsv not found — skipped.")

    if primers is not None:
        pairs_per_cand = defaultdict(int)
        for p in primers:
            pairs_per_cand[p["candidate_id"]] += 1
        lines.append(f"- **Stage 5 (primer design):** {len(primers):,} primer pairs designed "
                      f"across {len(pairs_per_cand):,} candidates.")
    else:
        lines.append("- **Stage 5:** primers.tsv not found — skipped.")

    if validated_all is not None:
        counts = defaultdict(int)
        for r in validated_all:
            counts[r.get("status", "UNKNOWN")] += 1
        total = len(validated_all)
        lines.append(f"- **Stage 6 (PCR validation):** {total:,} pairs tested.")
        for status in (
            "PASS",
            "TARGET_MULTIPLE_PRODUCTS",
            "REJECT_NO_TARGET_PRODUCT",
            "REJECT_OFFTARGET",
            "REJECT_HIGH_COPY_TARGET",
        ):
            if counts.get(status):
                lines.append(f"    - {status}: {counts[status]:,}")
        other = total - sum(counts.get(s, 0) for s in (
            "PASS", "TARGET_MULTIPLE_PRODUCTS", "REJECT_NO_TARGET_PRODUCT",
            "REJECT_OFFTARGET", "REJECT_HIGH_COPY_TARGET",
        ))
        if other:
            lines.append(f"    - other/unrecognized status: {other:,}")
    else:
        lines.append("- **Stage 6:** rejection_summary.tsv not found — skipped.")

    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates-dir", default="results/candidates",
                     help="directory with stage 4-6 outputs (default: results/candidates)")
    ap.add_argument("--pcr-work-dir", default="data/interim/pcr_validation",
                     help="stage 6 working directory (default: data/interim/pcr_validation)")
    ap.add_argument("-o", "--out", default="results/report.md",
                     help="output Markdown report path")
    ap.add_argument("--top", type=int, default=20,
                     help="max number of PASS pairs to list in the shortlist (default: 20)")
    args = ap.parse_args()

    cdir = Path(args.candidates_dir)
    wdir = Path(args.pcr_work_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ranked = read_tsv(cdir / "ranked_candidates.tsv")
    primers = read_tsv(cdir / "primers.tsv")
    validated_pass = read_tsv(cdir / "validated_primers.tsv")
    validated_all = read_tsv(wdir / "rejection_summary.tsv")
    offtarget_products = read_tsv(wdir / "offtarget_products.tsv")

    ranked_by_id = {r["candidate_id"]: r for r in ranked} if ranked else {}

    # ------------------------------------------------------------------
    # Off-target hit counts, per candidate x genome (for the "how close
    # did this get" table -- useful even for pairs that passed, since a
    # near-miss on a related species is still worth knowing about).
    # ------------------------------------------------------------------
    off_by_pair_genome = defaultdict(lambda: defaultdict(int))
    if offtarget_products:
        for row in offtarget_products:
            off_by_pair_genome[row["pair_name"]][row["genome"]] += 1

    # ------------------------------------------------------------------
    # Shortlist: PASS pairs, joined with stage-4 evidence, sorted by
    # copy number (more copies = more sensitive assay) then core identity.
    # ------------------------------------------------------------------
    shortlist_rows = []
    if validated_pass:
        for p in validated_pass:
            cand = ranked_by_id.get(p["candidate_id"], {})
            shortlist_rows.append({
                "pair_name": f"{p['candidate_id']}_pair{p['pair_rank']}",
                "candidate_id": p["candidate_id"],
                "n_copies": cand.get("n_copies", "?"),
                "core_identity": fmt_pct(cand.get("core_identity")) if cand.get("core_identity") not in (None, "") else "?",
                "product_size": p.get("product_size", "?"),
                "fwd_seq": p.get("fwd_seq", ""),
                "fwd_tm": p.get("fwd_tm", ""),
                "rev_seq": p.get("rev_seq", ""),
                "rev_tm": p.get("rev_tm", ""),
                "target_products": p.get("target_products", ""),
                "_sort_copies": int(cand["n_copies"]) if str(cand.get("n_copies", "")).isdigit() else -1,
                "_sort_identity": float(cand["core_identity"]) if cand.get("core_identity") not in (None, "") else -1,
            })
        shortlist_rows.sort(key=lambda r: (-r["_sort_copies"], -r["_sort_identity"]))

    # ------------------------------------------------------------------
    # Near-misses: pairs that DID amplify the target (target_products >= 1)
    # but were rejected for another reason. These are candidates worth a
    # second look if the shortlist ends up thin.
    # ------------------------------------------------------------------
    near_miss_rows = []
    if validated_all:
        for r in validated_all:
            if r.get("status") == "PASS":
                continue
            if fmt_int(r.get("target_products"), "0") == "0" or r.get("target_products", "0") == "0":
                continue
            cand = ranked_by_id.get(r.get("candidate_id"), {})
            genomes = r.get("offtarget_genomes", "")
            near_miss_rows.append({
                "pair_name": r["pair_name"],
                "status": r["status"],
                "n_copies": cand.get("n_copies", "?"),
                "target_products": r.get("target_products", ""),
                "offtarget_products": r.get("offtarget_products", ""),
                "offtarget_genomes": genomes or "-",
            })

    # ------------------------------------------------------------------
    # Write report
    # ------------------------------------------------------------------
    lines = []
    lines.append("# Tropilaelaps mercedesae primer pipeline — results report")
    lines.append("")
    lines.append("Auto-generated by `scripts/07_report.py`. Source tables: "
                  f"`{cdir}/`, `{wdir}/`.")
    lines.append("")

    lines.append("## Pipeline funnel")
    lines.append("")
    lines.append(build_funnel(ranked, primers, validated_all))
    lines.append("")

    lines.append("## Shortlist: PASS primer pairs")
    lines.append("")
    if shortlist_rows:
        lines.append(
            f"{len(shortlist_rows)} pair(s) passed validation "
            f"(exactly one target product, zero off-target products, "
            f"not high-copy-on-target-flagged). "
            f"Showing top {min(args.top, len(shortlist_rows))}, ranked by target "
            f"copy number then conserved-core identity."
        )
        lines.append("")
        lines.append(md_table(
            shortlist_rows[:args.top],
            columns=["pair_name", "n_copies", "core_identity", "product_size",
                     "fwd_seq", "fwd_tm", "rev_seq", "rev_tm"],
            headers=["pair", "target copies", "core identity", "amplicon (bp)",
                     "fwd primer", "fwd Tm", "rev primer", "rev Tm"],
        ))
    else:
        lines.append("No PASS pairs found (or validated_primers.tsv is missing/empty). "
                      "See the near-miss table below for what came closest.")
    lines.append("")

    lines.append("## Near misses (amplified target, rejected for another reason)")
    lines.append("")
    if near_miss_rows:
        lines.append(md_table(
            near_miss_rows,
            columns=["pair_name", "status", "n_copies", "target_products",
                     "offtarget_products", "offtarget_genomes"],
            headers=["pair", "rejection reason", "target copies", "target products",
                     "off-target products", "off-target genome(s)"],
        ))
    else:
        lines.append("None, or rejection_summary.tsv not found.")
    lines.append("")

    lines.append("## Notes")
    lines.append("")
    lines.append(
        "- `n_copies` / `core_identity` come from stage 4 (genome self-mapping); "
        "they describe the *candidate repeat region*, not the primer pair itself."
    )
    lines.append(
        "- `target_products` / off-target counts come from stage 6 in-silico PCR "
        "(BLASTN-short + explicit product enumeration). This is a computational "
        "prediction, not a wet-lab result — see the validation guidance below "
        "before ordering primers."
    )
    lines.append("")

    out_path.write_text("\n".join(lines) + "\n")

    print(f">> report written -> {out_path}", file=sys.stderr)
    print(f"   shortlist: {len(shortlist_rows)} PASS pairs", file=sys.stderr)
    print(f"   near misses: {len(near_miss_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
