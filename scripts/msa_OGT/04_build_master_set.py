#!/usr/bin/env python3
"""
Step 04: build the OGT-annotated master homolog set.

This joins class-I confirmed sequences to the taxid/OGTFinder metadata from
step 02 and writes one authoritative FASTA/TSV pair for sampling and alignment.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def load_tsv(path: Path, key: str):
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {row[key]: row for row in csv.DictReader(f, delimiter="\t")}


def main():
    in_fa = cfg.INTER_DIR / "classI_confirmed.fa"
    hit_meta = cfg.INTER_DIR / "hits_metadata.tsv"
    class_report = cfg.INTER_DIR / "classification_report.tsv"
    out_fa = cfg.INTER_DIR / "master_homologs.fa"
    out_tsv = cfg.INTER_DIR / "master_homologs.tsv"

    if not in_fa.exists():
        sys.exit("ERROR: run 03_classify_classI.py first")
    if not hit_meta.exists():
        sys.exit("ERROR: run 02_extract_annotate.py first")

    metadata = load_tsv(hit_meta, "id")
    class_meta = load_tsv(class_report, "id")
    entries = cfg.read_fasta(in_fa)
    cfg.write_fasta(entries, out_fa)

    fields = [
        "id", "accession", "source", "taxid", "organism", "genus", "length",
        "pident", "qcov", "tcov", "evalue", "bits",
        "assigned_class", "best_ref", "class_pident", "class_bits",
    ] + cfg.OGTFinderLookup.fields

    rows = []
    for sid, hdr, seq in entries:
        base = metadata.get(sid, {})
        cls = class_meta.get(sid, {})
        if not base:
            organism = cfg.extract_organism(hdr)
            base = {
                "id": sid,
                "accession": cfg.extract_acc(sid),
                "taxid": cfg.extract_taxid(hdr),
                "organism": organism,
                "genus": cfg.extract_genus(hdr),
                "length": str(len(seq.replace("-", "").replace("*", ""))),
                **{field: "" for field in cfg.OGTFinderLookup.fields},
            }
        row = {field: "" for field in fields}
        row.update({k: base.get(k, "") for k in row})
        row["id"] = sid
        row["source"] = "uniref_homolog_classI" if sid.startswith("UniRef") else "query"
        row["assigned_class"] = cls.get("assigned_class", "query" if row["source"] == "query" else "")
        row["best_ref"] = cls.get("best_ref", "")
        row["class_pident"] = cls.get("class_pident", "")
        row["class_bits"] = cls.get("class_bits", "")
        rows.append(row)

    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    labelled = [row for row in rows if row["ogt"]]
    range_labelled = [row for row in labelled if row["ogt_has_range"] == "yes"]
    regimes = {"psychro": 0, "meso": 0, "thermo": 0}
    for row in labelled:
        if row["regime"] in regimes:
            regimes[row["regime"]] += 1

    print(f"Saved: {out_fa} ({len(entries)} sequences)")
    print(f"Saved: {out_tsv}")
    print(
        "OGTFinder labels: "
        f"{len(labelled)} total; ranges={len(range_labelled)}; "
        f"psychro={regimes['psychro']}, meso={regimes['meso']}, thermo={regimes['thermo']}"
    )


if __name__ == "__main__":
    main()

