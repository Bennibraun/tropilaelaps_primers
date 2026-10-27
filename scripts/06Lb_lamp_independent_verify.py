#!/usr/bin/env python3
"""Stage 6Lb: independent raw-FASTA verification of LAMP component binding.

This is deliberately independent of BLAST. It searches every LAMP component
and its reverse complement directly in unwrapped FASTA sequences, then
reconstructs the stage 6L co-location rule from those coordinates. It verifies
sequence binding and genomic co-location; it is not a biochemical simulation.

The default scope is the stage 6L shortlist: rows with on_target_ok=yes and
offtarget_verdict=PASS. Use --all to include REVIEW, FAIL, and on-target
failures.
"""
import argparse
import csv
import sys
from collections import deque
from collections import defaultdict
from pathlib import Path

COLOCATION_WINDOW = 500
INNER = {"FIP_F1c", "FIP_F2", "BIP_B1c", "BIP_B2"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def occurrences(sequence, query):
    """Return all overlapping exact occurrences of query in sequence."""
    starts = []
    start = sequence.find(query)
    while start >= 0:
        starts.append(start)
        start = sequence.find(query, start + 1)
    return starts


def build_automaton(rows):
    """Build an Aho-Corasick matcher for every component in every set."""
    transitions = [{}]
    failures = [0]
    outputs = [[]]
    for row_index, row in enumerate(rows):
        for name, query in components(row).items():
            for strand, oriented in (("+", query), ("-", revcomp(query))):
                node = 0
                for base in oriented:
                    node = transitions[node].setdefault(base, len(transitions))
                    if node == len(failures):
                        transitions.append({})
                        failures.append(0)
                        outputs.append([])
                outputs[node].append((row_index, name, strand, len(query)))

    queue = deque()
    for node in transitions[0].values():
        queue.append(node)
    while queue:
        node = queue.popleft()
        for base, child in transitions[node].items():
            queue.append(child)
            fallback = failures[node]
            while fallback and base not in transitions[fallback]:
                fallback = failures[fallback]
            failures[child] = transitions[fallback].get(base, 0)
            outputs[child].extend(outputs[failures[child]])
    return transitions, failures, outputs


def scan_genome(genome, automaton, row_count):
    transitions, failures, outputs = automaton
    hits_by_row = [defaultdict(list) for _ in range(row_count)]
    for subject, sequence in read_fasta(genome):
        node = 0
        for end, base in enumerate(sequence):
            while node and base not in transitions[node]:
                node = failures[node]
            node = transitions[node].get(base, 0)
            for row_index, name, strand, length in outputs[node]:
                start = end - length + 1
                hits_by_row[row_index][name].append(
                    (subject, start, end + 1, strand)
                )
    return hits_by_row


def colocated_groups(hits):
    by_subject = defaultdict(list)
    for component, locations in hits.items():
        for subject, start, end, strand in locations:
            by_subject[subject].append((start, end, component, strand))

    groups = []
    for subject, locations in by_subject.items():
        locations.sort()
        for i, (start, _end, component, strand) in enumerate(locations):
            members = {(component, strand)}
            j = i + 1
            while j < len(locations) and locations[j][0] - start <= COLOCATION_WINDOW:
                members.add((locations[j][2], locations[j][3]))
                j += 1
            components = {name for name, _strand in members}
            if len(components) >= 2:
                groups.append((subject, start, components))
    return groups


def components(row):
    result = {name: row[name] for name in ("F3", "B3", "LF", "LB") if row.get(name)}
    required = ("FIP_F1c", "FIP_F2", "BIP_B1c", "BIP_B2")
    missing = [name for name in required if not row.get(name)]
    if missing:
        raise RuntimeError(
            f"{row.get('candidate_id', '?')}: missing LAMP component columns: "
            + ", ".join(missing)
        )
    result.update({name: row[name] for name in required})
    return result


def classify(hits):
    groups = colocated_groups(hits)
    inner_groups = [group for group in groups if len(group[2] & INNER) >= 2]
    return groups, inner_groups


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("validated_tsv")
    ap.add_argument("target")
    ap.add_argument("offtargets", nargs="+")
    ap.add_argument("-o", "--out", default="data/interim/lamp_validation/independent_verification.tsv")
    ap.add_argument(
        "--designed-tsv",
        default="",
        help="stage 5L output containing split FIP/BIP components; "
        "defaults to lamp_primers.tsv beside validated_tsv",
    )
    ap.add_argument("--all", action="store_true", help="verify every row, not only clean PASS rows")
    args = ap.parse_args()

    with open(args.validated_tsv, newline="") as fh:
        source_rows = list(csv.DictReader(fh, delimiter="\t"))
    if not source_rows:
        raise RuntimeError(f"no rows in {args.validated_tsv}")
    designed_path = Path(args.designed_tsv) if args.designed_tsv else (
        Path(args.validated_tsv).with_name("lamp_primers.tsv")
    )
    with open(designed_path, newline="") as fh:
        designed_rows = {
            row["candidate_id"]: row for row in csv.DictReader(fh, delimiter="\t")
        }
    rows = []
    for row in source_rows:
        design = designed_rows.get(row["candidate_id"])
        if design is None:
            raise RuntimeError(
                f"{row['candidate_id']}: not found in {designed_path}"
            )
        merged = dict(design)
        merged.update(row)
        if args.all or (
            row.get("on_target_ok") == "yes"
            and row.get("offtarget_verdict") == "PASS"
        ):
            rows.append(merged)
    if not rows:
        raise RuntimeError("no clean LAMP rows selected; use --all to verify every row")

    genomes = [("target", Path(args.target))] + [
        (f"offtarget:{Path(path).stem}", Path(path)) for path in args.offtargets
    ]
    output = []

    automaton = build_automaton(rows)

    hits_by_genome = {}
    for label, genome in genomes:
        print(f">> scanning {label}", file=sys.stderr)
        hits_by_genome[label] = scan_genome(genome, automaton, len(rows))

    for index, row in enumerate(rows, 1):
        comps = components(row)
        target_hits = hits_by_genome["target"][index - 1]
        target_groups, target_inner = classify(target_hits)
        target_bound = sum(bool(target_hits.get(name)) for name in comps)
        for label, genome in genomes:
            hits = hits_by_genome[label][index - 1]
            groups, inner_groups = (
                (target_groups, target_inner)
                if label == "target" else classify(hits)
            )
            output.append({
                "candidate_id": row["candidate_id"],
                "genome": label,
                "components_total": len(comps),
                "components_bound": sum(bool(hits.get(name)) for name in comps),
                "inner_components_bound": sum(bool(hits.get(name)) for name in INNER),
                "co_located_groups": len(groups),
                "inner_co_located_groups": len(inner_groups),
                "target_on_target_ok": "yes" if target_bound >= 4 and target_inner else "NO",
                "stage6L_offtarget_verdict": row.get("offtarget_verdict", ""),
                "independent_verdict": (
                    "TARGET_OK" if label == "target" and target_bound >= 4 and target_inner
                    else "TARGET_NOT_OK" if label == "target"
                    else "FAIL" if inner_groups
                    else "REVIEW" if groups
                    else "PASS"
                ),
            })
        if index % 10 == 0:
            print(f">> verified {index}/{len(rows)} LAMP sets", file=sys.stderr)

    fields = list(output[0])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(f">> independent LAMP verification -> {out}", file=sys.stderr)
    print(f">> verified {len(rows)} LAMP sets across {len(genomes)} genomes", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
