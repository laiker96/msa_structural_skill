#!/usr/bin/env python3
"""Step 06: download AlphaFold DB models for query/homolog accessions."""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", nargs="*", help="Explicit UniProt accessions.")
    parser.add_argument("--metadata", type=Path, default=cfg.OUTPUT_DIR / "repset_metadata.tsv")
    parser.add_argument("--column", default="accession")
    parser.add_argument("--list", type=Path, help="One accession per line.")
    parser.add_argument("--query-accession", default="", help="Optional query UniProt accession to include.")
    parser.add_argument("--dest", type=Path, default=cfg.STRUCTURE_DIR / "afdb")
    parser.add_argument("--manifest", type=Path, default=cfg.OUTPUT_DIR / "afdb_downloads" / "download_manifest.tsv")
    parser.add_argument("--afdb-version", default="6")
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--strict-missing", action="store_true")
    return parser.parse_args()


def afdb_url(acc: str, version: str) -> str:
    return f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{version}.pdb"


def out_path(acc: str, dest: Path, version: str) -> Path:
    return dest / f"AF-{acc}-F1-model_v{version}.pdb"


def read_list(path: Path) -> list[str]:
    return [
        line.split("#", 1)[0].strip()
        for line in path.read_text().splitlines()
        if line.split("#", 1)[0].strip()
    ]


def read_metadata(path: Path, column: str) -> list[str]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if column not in (reader.fieldnames or []):
            raise SystemExit(f"ERROR: {path} has no column {column}; available: {reader.fieldnames}")
        return [row.get(column, "").strip() for row in reader]


def collect_accessions(args: argparse.Namespace) -> list[str]:
    values = list(args.accessions)
    if args.query_accession:
        values.append(args.query_accession)
    if args.list:
        values.extend(read_list(args.list))
    values.extend(read_metadata(args.metadata, args.column))
    return cfg.unique_preserve_order([cfg.normalize_accession(value) for value in values])


def fetch_one(acc: str, args: argparse.Namespace) -> dict[str, str]:
    url = afdb_url(acc, args.afdb_version)
    out = out_path(acc, args.dest, args.afdb_version)
    row = {"accession": acc, "status": "", "path": str(out), "url": url, "error": ""}
    if out.exists() and out.stat().st_size > 0 and not args.force:
        row["status"] = "skip"
        return row
    args.dest.mkdir(parents=True, exist_ok=True)
    tmp_name = ""
    try:
        with tempfile.NamedTemporaryFile(dir=args.dest, prefix=f".{out.name}.", delete=False) as tmp:
            tmp_name = tmp.name
            with urllib.request.urlopen(url, timeout=args.timeout) as response:
                shutil.copyfileobj(response, tmp)
        Path(tmp_name).replace(out)
        row["status"] = "ok"
    except urllib.error.HTTPError as exc:
        if tmp_name:
            Path(tmp_name).unlink(missing_ok=True)
        row["status"] = "not-in-afdb" if exc.code == 404 else "error"
        row["error"] = "" if exc.code == 404 else f"HTTP {exc.code}"
    except Exception as exc:
        if tmp_name:
            Path(tmp_name).unlink(missing_ok=True)
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    cfg.write_tsv(path, rows, ["accession", "status", "path", "url", "error"])


def main() -> None:
    args = parse_args()
    accessions = collect_accessions(args)
    if not accessions:
        raise SystemExit("ERROR: no accessions resolved.")
    if args.dry_run:
        rows = [
            {
                "accession": acc,
                "status": "dry-run",
                "path": str(out_path(acc, args.dest, args.afdb_version)),
                "url": afdb_url(acc, args.afdb_version),
                "error": "",
            }
            for acc in accessions
        ]
        write_manifest(args.manifest, rows)
        print(f"Dry run only. Wrote: {args.manifest}")
        return
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_acc = {executor.submit(fetch_one, acc, args): acc for acc in accessions}
        for future in as_completed(future_to_acc):
            row = future.result()
            rows.append(row)
            print(f"[{row['status']}] {row['accession']}", flush=True)
    order = {acc: idx for idx, acc in enumerate(accessions)}
    rows.sort(key=lambda row: order[row["accession"]])
    write_manifest(args.manifest, rows)
    counts = {status: sum(row["status"] == status for row in rows) for status in ["ok", "skip", "not-in-afdb", "error"]}
    print(f"Saved: {args.manifest}")
    print(
        "Summary: "
        f"downloaded={counts['ok']} already-present={counts['skip']} "
        f"not-in-AFDB={counts['not-in-afdb']} errors={counts['error']}"
    )
    if args.strict_missing and counts["not-in-afdb"]:
        raise SystemExit("ERROR: at least one accession is absent from AlphaFold DB.")
    if counts["error"]:
        raise SystemExit("ERROR: at least one AFDB download failed.")


if __name__ == "__main__":
    main()
