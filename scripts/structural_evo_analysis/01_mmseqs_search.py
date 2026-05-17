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
    parser.add_argument(
        "--db-name",
        default=None,
        help="Database name or UniRef shorthand. Examples: uniref50, uniref90, uniref100, 50, 90, 100.",
    )
    parser.add_argument(
        "--db-fasta",
        type=Path,
        default=None,
        help="Protein FASTA to search. Defaults to $UNIREF_DIR/<db-name>/<db-name>.fasta.gz.",
    )
    parser.add_argument(
        "--db-mmseqs",
        type=Path,
        default=None,
        help="MMseqs database prefix. Defaults beside the selected FASTA as <fasta-stem>_db.",
    )
    parser.add_argument(
        "--ogt-metadata",
        type=Path,
        default=cfg.DEFAULT_OGT_TSV,
        help="OGT metadata TSV used only with --join-ogt. Supports data/ogt_taxid_summary.tsv or an external raw OGTFinder table.",
    )
    parser.add_argument("--join-ogt", action="store_true", help="Join OGT/regime metadata into hit tables.")
    parser.add_argument("--no-ogt", action="store_true", help="Deprecated no-op; OGT joins are disabled unless --join-ogt is set.")
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
    parser.add_argument(
        "--max-repset-seqs",
        type=int,
        default=int(os.environ.get("SEA_MAX_REPSET_SEQS", "500")),
        help="Maximum representative-set size including query. Use 0 to keep all filtered hits.",
    )
    parser.add_argument(
        "--identity-bin-width",
        type=float,
        default=float(os.environ.get("SEA_IDENTITY_BIN_WIDTH", "5.0")),
        help="Query-identity bin width for diverse default hit selection.",
    )
    parser.add_argument("--threads", type=int, default=cfg.N_THREADS)
    parser.add_argument("--force", action="store_true", help="Overwrite existing search TSV and derived FASTAs.")
    args = parser.parse_args()
    db_name_from_cli = args.db_name is not None
    db_fasta_from_cli = args.db_fasta is not None
    args.db_name = cfg.normalize_db_name(args.db_name or cfg.DEFAULT_DB_NAME)
    if args.db_fasta is None:
        if db_name_from_cli:
            args.db_fasta = cfg.default_db_fasta(args.db_name)
        else:
            args.db_fasta = Path(os.environ.get("SEA_DB_FASTA", cfg.default_db_fasta(args.db_name)))
    if args.db_mmseqs is None:
        if db_name_from_cli or db_fasta_from_cli:
            args.db_mmseqs = cfg.default_db_mmseqs(args.db_name, args.db_fasta)
        else:
            args.db_mmseqs = Path(os.environ.get("SEA_DB_MMSEQS", cfg.default_db_mmseqs(args.db_name, args.db_fasta)))
    return args


def require_mmseqs() -> str:
    mmseqs = cfg.resolve_bin("mmseqs")
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit("ERROR: mmseqs not found. Activate ./envs/structural_evo or add mmseqs to PATH.")
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


def select_diverse_hits(rows: list[dict[str, object]], args: argparse.Namespace) -> list[dict[str, object]]:
    max_hits = args.max_repset_seqs - 1 if args.max_repset_seqs > 0 else 0
    if max_hits <= 0 or len(rows) <= max_hits:
        return rows
    width = max(args.identity_bin_width, 0.1)
    bins: dict[int, list[dict[str, object]]] = {}
    for row in rows:
        pident = float(row["pident"])
        bin_id = int(pident // width)
        bins.setdefault(bin_id, []).append(row)
    for bin_rows in bins.values():
        bin_rows.sort(key=lambda row: (-float(row["bits"]), -float(row["qcov"]), str(row["id"])))

    selected = []
    ordered_bins = sorted(bins)
    while len(selected) < max_hits and ordered_bins:
        next_bins = []
        for bin_id in ordered_bins:
            bin_rows = bins[bin_id]
            if bin_rows and len(selected) < max_hits:
                selected.append(bin_rows.pop(0))
            if bin_rows:
                next_bins.append(bin_id)
        ordered_bins = next_bins
    return selected


def write_outputs(args: argparse.Namespace, hits_tsv: Path, all_fa: Path) -> None:
    query_sid, query_header, query_seq = cfg.query_entry(args.query)
    hits = load_best_hits(hits_tsv)
    if args.force or not all_fa.exists():
        found = extract_hits(set(hits), args.db_fasta, all_fa)
        print(f"Extracted {found}/{len(hits)} database hits into {all_fa}")
    entries = cfg.read_fasta(all_fa)
    lookup = cfg.OGTFinderLookup(args.ogt_metadata) if args.join_ogt and not args.no_ogt else None
    rows = []
    passed_entries = []
    for sid, header, seq in entries:
        row = hits.get(sid)
        if row is None:
            continue
        reasons = filter_reasons(row, seq, args)
        accession = cfg.normalize_accession(sid)
        taxid = cfg.extract_taxid(header)
        organism = cfg.extract_organism(header)
        clean = cfg.clean_sequence(seq)
        meta = {
            "id": sid,
            "accession": accession,
            "taxid": taxid,
            "organism": organism,
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
        if lookup is not None:
            meta.update(lookup.lookup(taxid, organism))
        rows.append(meta)
        if not reasons:
            passed_entries.append({
                **meta,
                "sid": sid,
                "fasta_header": header,
                "sequence": clean,
            })
    selected_entries = select_diverse_hits(passed_entries, args)
    selected_ids = {entry["sid"] for entry in selected_entries}
    filtered = [(query_sid, query_header, query_seq)] + [
        (entry["sid"], entry["fasta_header"], entry["sequence"]) for entry in selected_entries
    ]
    for row in rows:
        if row["passes_default_filters"] == "yes":
            row["selected_for_msa"] = "yes" if row["id"] in selected_ids else "no"
            row["selection_reason"] = "diverse_identity_subset" if row["id"] in selected_ids else "not_selected_diverse_subset"
        else:
            row["selected_for_msa"] = "no"
            row["selection_reason"] = "filtered"
    query_meta = {
        "id": query_sid,
        "accession": cfg.normalize_accession(query_sid),
        "taxid": cfg.extract_taxid(query_header),
        "organism": cfg.extract_organism(query_header),
        "header": query_header,
        "length": len(query_seq),
        "pident": "100.000",
        "qcov": "1.000",
        "tcov": "1.000",
        "evalue": "0",
        "bits": "",
        "passes_default_filters": "query",
        "filter_reasons": "",
        "selected_for_msa": "yes",
        "selection_reason": "query",
    }
    if lookup is not None:
        query_meta.update(lookup.lookup(query_meta["taxid"], query_meta["organism"]))
    fields = [
        "id", "accession", "taxid", "organism", "header", "length", "pident", "qcov", "tcov",
        "evalue", "bits", "passes_default_filters", "filter_reasons", "selected_for_msa", "selection_reason",
    ] + ([] if lookup is None else cfg.OGTFinderLookup.fields)
    cfg.write_tsv(args.out_dir / "hits_metadata.tsv", rows, fields)
    cfg.write_tsv(args.out_dir / "repset_metadata.tsv", [query_meta] + rows, fields)
    cfg.write_fasta(filtered, args.out_dir / "repset.fa")
    print(f"Saved: {args.out_dir / 'hits_metadata.tsv'} ({len(rows)} rows)")
    print(f"Selected {len(filtered) - 1}/{len(passed_entries)} filtered hits for diverse MSA subset")
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
