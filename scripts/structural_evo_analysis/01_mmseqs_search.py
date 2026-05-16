#!/usr/bin/env python3
"""Step 01: search a query protein against a sequence database with MMseqs2."""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


FIELDS = [
    "target", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
    "tstart", "tend", "evalue", "bits", "qcov", "tcov", "theader",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", type=Path, required=True, help="Single-sequence query FASTA.")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR)
    parser.add_argument("--db-name", default=cfg.DEFAULT_DB_NAME)
    parser.add_argument("--db-fasta", type=Path, default=cfg.DEFAULT_DB_FASTA)
    parser.add_argument("--db-mmseqs", type=Path, default=cfg.DEFAULT_DB_MMSEQS)
    parser.add_argument("--sensitivity", type=float, default=float(os.environ.get("SEA_MMSEQS_S", "7.5")))
    parser.add_argument("--evalue", default=os.environ.get("SEA_MMSEQS_E", "1e-5"))
    parser.add_argument("--min-seq-id", type=float, default=float(os.environ.get("SEA_MIN_SEQ_ID", "0.05")))
    parser.add_argument("--coverage", type=float, default=float(os.environ.get("SEA_COVERAGE", "0.50")))
    parser.add_argument("--max-seqs", type=int, default=int(os.environ.get("SEA_MAX_SEQS", "50000")))
    parser.add_argument("--min-length", type=int, default=int(os.environ.get("SEA_MIN_LENGTH", "0")))
    parser.add_argument("--max-length", type=int, default=int(os.environ.get("SEA_MAX_LENGTH", "1000000")))
    parser.add_argument("--post-min-qcov", type=float, default=float(os.environ.get("SEA_POST_MIN_QCOV", "0.50")))
    parser.add_argument("--post-min-identity", type=float, default=float(os.environ.get("SEA_POST_MIN_IDENTITY", "5.0")))
    parser.add_argument("--post-max-identity", type=float, default=float(os.environ.get("SEA_POST_MAX_IDENTITY", "99.9")))
    parser.add_argument("--threads", type=int, default=cfg.N_THREADS)
    parser.add_argument("--force", action="store_true", help="Overwrite existing search TSV and derived FASTAs.")
    return parser.parse_args()


def require_mmseqs() -> str:
    mmseqs = cfg.resolve_bin("mmseqs")
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit("ERROR: mmseqs not found. Activate ./envs/ankros or add mmseqs to PATH.")
    return mmseqs


def create_db(mmseqs: str, db_fasta: Path, db_mmseqs: Path) -> None:
    if Path(f"{db_mmseqs}.dbtype").exists():
        return
    if not db_fasta.exists():
        raise SystemExit(f"ERROR: database FASTA not found: {db_fasta}")
    db_mmseqs.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([mmseqs, "createdb", str(db_fasta), str(db_mmseqs)], check=True)


def run_search(args: argparse.Namespace, mmseqs: str, out_tsv: Path) -> None:
    work = args.out_dir / "mmseqs_work"
    tmp = work / "tmp"
    query_db = work / "query_db"
    result_db = work / "search_result"
    work.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)
    subprocess.run([mmseqs, "createdb", str(args.query), str(query_db)], check=True)
    subprocess.run([
        mmseqs, "search",
        str(query_db), str(args.db_mmseqs), str(result_db), str(tmp),
        "-s", str(args.sensitivity),
        "-e", str(args.evalue),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.coverage),
        "--cov-mode", "1",
        "--max-seqs", str(args.max_seqs),
        "--threads", str(args.threads),
    ], check=True)
    subprocess.run([
        mmseqs, "convertalis",
        str(query_db), str(args.db_mmseqs), str(result_db), str(out_tsv),
        "--format-output", ",".join(FIELDS),
    ], check=True)
    shutil.rmtree(tmp, ignore_errors=True)


def load_best_hits(path: Path) -> dict[str, dict[str, str]]:
    hits: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(FIELDS):
                continue
            row = dict(zip(FIELDS, parts))
            sid = row["target"]
            bits = float(row["bits"])
            if sid not in hits or bits > float(hits[sid]["bits"]):
                hits[sid] = row
    return hits


def extract_hits(hit_ids: set[str], db_fasta: Path, out_fa: Path) -> int:
    found = 0
    writing = False
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with cfg.open_text(db_fasta) as source, out_fa.open("w") as out:
        for line in source:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                writing = sid in hit_ids
                if writing:
                    found += 1
                    out.write(line)
            elif writing:
                out.write(line)
    return found


def filter_reasons(row: dict[str, str], seq: str, args: argparse.Namespace) -> list[str]:
    reasons = []
    length = len(cfg.clean_sequence(seq))
    pident = float(row["pident"])
    qcov = float(row["qcov"])
    if length < args.min_length:
        reasons.append("short")
    if length > args.max_length:
        reasons.append("long")
    if pident < args.post_min_identity:
        reasons.append("low_identity")
    if pident > args.post_max_identity:
        reasons.append("too_similar")
    if qcov < args.post_min_qcov:
        reasons.append("low_qcov")
    return reasons


def write_outputs(args: argparse.Namespace, hits_tsv: Path, all_fa: Path) -> None:
    query_sid, query_header, query_seq = cfg.query_entry(args.query)
    hits = load_best_hits(hits_tsv)
    if args.force or not all_fa.exists():
        found = extract_hits(set(hits), args.db_fasta, all_fa)
        print(f"Extracted {found}/{len(hits)} database hits into {all_fa}")
    entries = cfg.read_fasta(all_fa)
    rows = []
    filtered = [(query_sid, query_header, query_seq)]
    for sid, header, seq in entries:
        row = hits.get(sid)
        if row is None:
            continue
        reasons = filter_reasons(row, seq, args)
        accession = cfg.normalize_accession(sid)
        clean = cfg.clean_sequence(seq)
        meta = {
            "id": sid,
            "accession": accession,
            "header": header,
            "length": len(clean),
            "pident": f"{float(row['pident']):.3f}",
            "qcov": f"{float(row['qcov']):.3f}",
            "tcov": f"{float(row['tcov']):.3f}",
            "evalue": f"{float(row['evalue']):.3g}",
            "bits": f"{float(row['bits']):.1f}",
            "passes_default_filters": "yes" if not reasons else "no",
            "filter_reasons": ";".join(reasons),
        }
        rows.append(meta)
        if not reasons:
            filtered.append((sid, header, clean))
    query_meta = {
        "id": query_sid,
        "accession": cfg.normalize_accession(query_sid),
        "header": query_header,
        "length": len(query_seq),
        "pident": "100.000",
        "qcov": "1.000",
        "tcov": "1.000",
        "evalue": "0",
        "bits": "",
        "passes_default_filters": "query",
        "filter_reasons": "",
    }
    fields = [
        "id", "accession", "header", "length", "pident", "qcov", "tcov",
        "evalue", "bits", "passes_default_filters", "filter_reasons",
    ]
    cfg.write_tsv(args.out_dir / "hits_metadata.tsv", rows, fields)
    cfg.write_tsv(args.out_dir / "repset_metadata.tsv", [query_meta] + rows, fields)
    cfg.write_fasta(filtered, args.out_dir / "repset.fa")
    print(f"Saved: {args.out_dir / 'hits_metadata.tsv'} ({len(rows)} rows)")
    print(f"Saved: {args.out_dir / 'repset.fa'} ({len(filtered)} sequences including query)")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.out_dir / "mmseqs_search_results.tsv"
    all_fa = args.out_dir / "hits_all.fa"
    if out_tsv.exists() and not args.force:
        print(f"Using existing search results: {out_tsv}")
    else:
        mmseqs = require_mmseqs()
        create_db(mmseqs, args.db_fasta, args.db_mmseqs)
        run_search(args, mmseqs, out_tsv)
    write_outputs(args, out_tsv, all_fa)


if __name__ == "__main__":
    main()
