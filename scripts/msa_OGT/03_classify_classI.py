#!/usr/bin/env python3
"""
Step 03: classify filtered candidates as class-I CPD photolyase.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

MMSEQS = cfg.resolve_bin("mmseqs")


def build_reference_db(work_dir: Path):
    ref_fa = work_dir / "classification_refs.fa"
    entries = []
    labels = {}
    for path, label, prefix in [
        (cfg.UNIPROT_CLASSI_CPD, "class_I", "CLASS_I"),
        (cfg.CLASSII_KNOWN, "class_II", "CLASS_II"),
        (cfg.OTHER_KNOWN, "other", "OTHER"),
    ]:
        if not path.exists():
            continue
        for sid, hdr, seq in cfg.read_fasta(path):
            ref_hdr = f"{prefix}|{hdr}"
            entries.append((sid, ref_hdr, seq.replace("-", "")))
            labels[ref_hdr.split()[0]] = label
    cfg.write_fasta(entries, ref_fa)
    return ref_fa, labels


def classify(candidates_fa: Path, work_dir: Path):
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_fa, labels = build_reference_db(work_dir)
    if not labels:
        return {}

    query_db = work_dir / "query"
    ref_db = work_dir / "ref"
    result_db = work_dir / "result"
    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    subprocess.run([MMSEQS, "createdb", str(candidates_fa), str(query_db)],
                   check=True, capture_output=True)
    subprocess.run([MMSEQS, "createdb", str(ref_fa), str(ref_db)],
                   check=True, capture_output=True)
    subprocess.run([
        MMSEQS, "search", str(query_db), str(ref_db), str(result_db), str(tmp_dir),
        "-s", "7.5", "-e", "1e-5", "--max-seqs", "50",
        "--threads", str(cfg.N_THREADS),
    ], check=True, capture_output=True)

    out_tsv = work_dir / "classification_hits.tsv"
    subprocess.run([
        MMSEQS, "convertalis", str(query_db), str(ref_db), str(result_db), str(out_tsv),
        "--format-output", "query,target,pident,evalue,bits",
    ], check=True, capture_output=True)

    calls = {}
    with open(out_tsv) as f:
        for line in f:
            qid, tid, pident, _evalue, bits = line.rstrip("\n").split("\t")[:5]
            bits = float(bits)
            if qid in calls and bits <= calls[qid]["class_bits"]:
                continue
            calls[qid] = {
                "assigned_class": labels.get(tid, "unknown"),
                "best_ref": tid,
                "class_pident": float(pident),
                "class_bits": bits,
            }
    return calls


def main():
    in_fa = cfg.INTER_DIR / "hits_filtered.fa"
    in_meta = cfg.INTER_DIR / "hits_filtered_metadata.tsv"
    out_fa = cfg.INTER_DIR / "classI_confirmed.fa"
    out_report = cfg.INTER_DIR / "classification_report.tsv"

    if not in_fa.exists() or not in_meta.exists():
        sys.exit("ERROR: run 02_extract_annotate.py first")
    if out_fa.exists():
        n = sum(1 for line in open(out_fa) if line.startswith(">"))
        print(f"Classification already exists: {out_fa} ({n} class-I sequences)")
        return

    entries = cfg.read_fasta(in_fa)
    print(f"Candidates to classify: {len(entries)}")
    calls = classify(in_fa, cfg.INTER_DIR / "classify_work")
    if not calls:
        print("WARNING: no reference labels available; keeping all candidates")

    kept = []
    counts = defaultdict(int)
    rows = []
    for sid, hdr, seq in entries:
        call = calls.get(sid, {
            "assigned_class": "no_hit",
            "best_ref": "-",
            "class_pident": 0.0,
            "class_bits": 0.0,
        })
        counts[call["assigned_class"]] += 1
        rows.append({"id": sid, **call})
        if not calls or call["assigned_class"] == "class_I":
            kept.append((sid, hdr, seq))

    query_entry = cfg.read_fasta(cfg.QUERY_FASTA)[0]
    if not any(sid == query_entry[0] for sid, _, _ in kept):
        kept.insert(0, query_entry)

    cfg.write_fasta(kept, out_fa)
    with open(out_report, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["id", "assigned_class", "best_ref", "class_pident", "class_bits"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            row = dict(row)
            row["class_pident"] = f"{row['class_pident']:.1f}"
            row["class_bits"] = f"{row['class_bits']:.1f}"
            writer.writerow(row)

    print("Classification counts:")
    for cls, count in sorted(counts.items()):
        print(f"  {cls}: {count}")
    print(f"Saved: {out_fa} ({len(kept)} sequences, including query)")
    print(f"Saved: {out_report}")


if __name__ == "__main__":
    main()
