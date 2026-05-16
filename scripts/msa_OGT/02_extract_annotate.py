#!/usr/bin/env python3
"""
Step 02: extract UniRef sequences and annotate every hit with OGTFinder.

Outputs:
  - hits_all.fa                 all extracted returned hits
  - hits_metadata.tsv           all extracted hits with search + OGT metadata
  - hits_filtered.fa            post-search candidate FASTA for class filtering
  - hits_filtered_metadata.tsv  metadata rows passing default candidate filters
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


SEARCH_FIELDS = [
    "id", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
    "tstart", "tend", "evalue", "bits", "qcov", "tcov", "theader",
]

BASE_FIELDS = [
    "id", "accession", "taxid", "organism", "genus", "length",
    "pident", "qcov", "tcov", "evalue", "bits", "passes_default_filters",
    "filter_reasons",
]
METADATA_FIELDS = BASE_FIELDS + cfg.OGTFinderLookup.fields


def load_search_hits(path: Path):
    hits = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(SEARCH_FIELDS):
                continue
            row = dict(zip(SEARCH_FIELDS, parts))
            row["pident"] = float(row["pident"])
            row["qcov"] = float(row["qcov"])
            row["tcov"] = float(row["tcov"])
            row["bits"] = float(row["bits"])
            row["evalue"] = float(row["evalue"])
            sid = row["id"]
            if sid not in hits or row["bits"] > hits[sid]["bits"]:
                hits[sid] = row
    return hits


def extract_from_uniref(hit_ids, output_fa: Path):
    source_fa = cfg.INTER_DIR / f"{cfg.ACTIVE_DB_NAME}_photolyase_subset.fa"
    if not source_fa.exists():
        source_fa = cfg.ACTIVE_DB_FASTA
    hit_set = set(hit_ids)
    found = 0
    print(f"Extracting {len(hit_set)} sequences from {source_fa}")
    with cfg.open_fasta(source_fa) as fin, open(output_fa, "w") as fout:
        writing = False
        for line in fin:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                writing = sid in hit_set
                if writing:
                    found += 1
                    fout.write(line)
                    if found % 1000 == 0:
                        print(f"  found {found}/{len(hit_set)}")
            elif writing:
                fout.write(line)
    print(f"Extracted {found}/{len(hit_set)} sequences")


def filter_reasons(row, seq):
    reasons = []
    length = len(seq.replace("-", "").replace("*", ""))
    if length < cfg.MIN_LENGTH:
        reasons.append("short")
    if length > cfg.MAX_LENGTH:
        reasons.append("long")
    if row["pident"] < cfg.POST_MIN_IDENTITY:
        reasons.append("low_identity")
    if row["pident"] > cfg.POST_MAX_IDENTITY:
        reasons.append("too_similar")
    if row["qcov"] < cfg.POST_MIN_QCOV:
        reasons.append("low_qcov")
    return reasons


def main():
    search_tsv = cfg.INTER_DIR / "mmseqs_search_results.tsv"
    all_fa = cfg.INTER_DIR / "hits_all.fa"
    all_meta = cfg.INTER_DIR / "hits_metadata.tsv"
    filtered_fa = cfg.INTER_DIR / "hits_filtered.fa"
    filtered_meta = cfg.INTER_DIR / "hits_filtered_metadata.tsv"

    if not search_tsv.exists():
        sys.exit(f"ERROR: missing {search_tsv}; run 01_mmseqs_search.py first")

    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    hits = load_search_hits(search_tsv)
    print(f"Unique search targets: {len(hits)}")

    if not all_fa.exists():
        extract_from_uniref(hits.keys(), all_fa)
    else:
        n = sum(1 for line in open(all_fa) if line.startswith(">"))
        print(f"Using existing {all_fa.name} ({n} sequences)")

    lookup = cfg.OGTFinderLookup()
    entries = cfg.read_fasta(all_fa)
    annotated_rows = []
    filtered_entries = []
    filtered_rows = []
    reason_counts = {}

    for sid, hdr, seq in entries:
        h = hits.get(sid)
        if not h:
            continue
        clean_seq = seq.replace("-", "").replace("*", "")
        taxid = cfg.extract_taxid(hdr)
        organism = cfg.extract_organism(hdr)
        reasons = filter_reasons(h, clean_seq)
        passes = not reasons
        for reason in reasons:
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
        ogt = lookup.lookup(taxid, organism)
        row = {
            "id": sid,
            "accession": cfg.extract_acc(sid),
            "taxid": taxid,
            "organism": organism,
            "genus": cfg.extract_genus(hdr),
            "length": str(len(clean_seq)),
            "pident": f"{h['pident']:.3f}",
            "qcov": f"{h['qcov']:.3f}",
            "tcov": f"{h['tcov']:.3f}",
            "evalue": f"{h['evalue']:.3g}",
            "bits": f"{h['bits']:.1f}",
            "passes_default_filters": "yes" if passes else "no",
            "filter_reasons": ";".join(reasons),
            **ogt,
        }
        annotated_rows.append(row)
        if passes:
            filtered_entries.append((sid, hdr, clean_seq))
            filtered_rows.append(row)

    with open(all_meta, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=METADATA_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated_rows)
    with open(filtered_meta, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=METADATA_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(filtered_rows)
    cfg.write_fasta(filtered_entries, filtered_fa)

    labelled = [r for r in annotated_rows if r["ogt"]]
    ranges = [r for r in labelled if r["ogt_has_range"] == "yes"]
    print(f"Saved: {all_fa} ({len(entries)} sequences)")
    print(f"Saved: {all_meta} ({len(annotated_rows)} rows)")
    print(f"Saved: {filtered_fa} ({len(filtered_entries)} candidates)")
    print(f"Saved: {filtered_meta}")
    print(f"OGTFinder-labelled hits: {len(labelled)}; with range-derived OGT: {len(ranges)}")
    if reason_counts:
        print("Filter removals:")
        for reason, count in sorted(reason_counts.items()):
            print(f"  {reason}: {count}")


if __name__ == "__main__":
    main()

