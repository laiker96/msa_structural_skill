#!/usr/bin/env python3
"""Download AlphaFold DB v6 models for the MSA representative set.

By default this utility reads the QC-filtered representative metadata from
``results/msa_OGT/repset_metadata_qc.tsv`` and downloads every non-query
UniProt accession into ``structures/afdb``. Missing AFDB entries are recorded
as ``not-in-afdb`` and do not make the command fail unless
``--strict-missing`` is used.

Examples:
  python scripts/msa_OGT/13_download_afdb.py
  python scripts/msa_OGT/13_download_afdb.py --pre-qc-repset
  python scripts/msa_OGT/13_download_afdb.py --dry-run
  python scripts/msa_OGT/13_download_afdb.py Q6XKP4 R4YNG8
  python scripts/msa_OGT/13_download_afdb.py --tsv results/msa_OGT/repset_metadata.tsv --column accession
"""
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


ROOT = Path(__file__).resolve().parents[2]
MSA_OUT = ROOT / "results" / "msa_OGT"
DEFAULT_REPSET_METADATA = MSA_OUT / "repset_metadata_qc.tsv"
PRE_QC_REPSET_METADATA = MSA_OUT / "repset_metadata.tsv"
AFDB_DIR = ROOT / "structures" / "afdb"
DEFAULT_MANIFEST = MSA_OUT / "afdb_downloads" / "download_manifest.tsv"
AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"


def _out_path(acc: str, dest: Path) -> Path:
    return dest / f"AF-{acc}-F1-model_v6.pdb"


def normalize_accession(value: str) -> str:
    """Return a UniProt accession-like value from metadata/header identifiers."""
    text = (value or "").strip()
    if not text or text == "photoHymenobact":
        return ""
    if text.startswith("UniRef") and "_" in text:
        text = text.split("_", 1)[1]
    if "|" in text:
        parts = [part for part in text.split("|") if part]
        if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
            text = parts[1]
        else:
            text = parts[-1]
    if text == "photoHymenobact":
        return ""
    return text


def unique_accessions(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        acc = normalize_accession(value)
        if not acc or acc in seen:
            continue
        seen.add(acc)
        out.append(acc)
    return out


def fetch_one(acc: str, dest: Path = AFDB_DIR, timeout: float = 60.0, force: bool = False) -> dict[str, str]:
    """Download one accession if missing and return a manifest row."""
    url = AFDB_URL.format(acc=acc)
    out = _out_path(acc, dest)
    row = {"accession": acc, "status": "", "path": str(out), "url": url, "error": ""}
    if out.exists() and out.stat().st_size > 0 and not force:
        row["status"] = "skip"
        return row

    dest.mkdir(parents=True, exist_ok=True)
    tmp_name = ""
    try:
        with tempfile.NamedTemporaryFile(dir=dest, prefix=f".{out.name}.", delete=False) as tmp:
            tmp_name = tmp.name
            with urllib.request.urlopen(url, timeout=timeout) as response:
                shutil.copyfileobj(response, tmp)
        Path(tmp_name).replace(out)
        row["status"] = "ok"
    except urllib.error.HTTPError as exc:
        if tmp_name:
            Path(tmp_name).unlink(missing_ok=True)
        if exc.code == 404:
            row["status"] = "not-in-afdb"
        else:
            row["status"] = "error"
            row["error"] = f"HTTP {exc.code}"
    except Exception as exc:
        if tmp_name:
            Path(tmp_name).unlink(missing_ok=True)
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def fetch_records(
    accessions: list[str],
    n_workers: int = 20,
    dest: Path = AFDB_DIR,
    progress: bool = True,
    timeout: float = 60.0,
    force: bool = False,
) -> list[dict[str, str]]:
    """Download many AFDB models in parallel and return per-accession rows."""
    accessions = unique_accessions(accessions)
    dest.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futures = {ex.submit(fetch_one, acc, dest, timeout, force): acc for acc in accessions}
        for fut in as_completed(futures):
            row = fut.result()
            rows.append(row)
            if progress:
                tag = {
                    "ok": "[get ]",
                    "skip": "[skip]",
                    "not-in-afdb": "[miss]",
                    "error": "[err ]",
                }[row["status"]]
                print(f"{tag} {row['accession']}", flush=True)
    return sorted(rows, key=lambda item: accessions.index(item["accession"]))


def fetch_many(
    accessions: list[str],
    n_workers: int = 20,
    dest: Path = AFDB_DIR,
    progress: bool = True,
    timeout: float = 60.0,
    force: bool = False,
) -> dict[str, int]:
    """Download many AFDB models and return status counts."""
    return count_statuses(fetch_records(accessions, n_workers, dest, progress, timeout, force))


def count_statuses(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = {"ok": 0, "skip": 0, "not-in-afdb": 0, "error": 0}
    for row in rows:
        counts[row["status"]] += 1
    return counts


def read_list(path: Path) -> list[str]:
    out: list[str] = []
    for line in path.read_text().splitlines():
        text = line.split("#", 1)[0].strip()
        if text:
            out.append(text)
    return out


def read_tsv_column(path: Path, column: str) -> list[str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if column not in (reader.fieldnames or []):
            raise SystemExit(f"{path}: no column '{column}'. Have: {reader.fieldnames}")
        return [row[column].strip() for row in reader if row[column].strip()]


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "status", "path", "url", "error"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def collect_accessions(args: argparse.Namespace) -> tuple[list[str], str]:
    accs: list[str] = list(args.accessions)
    source = "explicit"
    if args.list:
        accs += read_list(args.list)
        source = str(args.list)
    if args.tsv:
        accs += read_tsv_column(args.tsv, args.column)
        source = f"{args.tsv}:{args.column}"
    if not accs:
        repset = PRE_QC_REPSET_METADATA if args.pre_qc_repset else args.repset_metadata
        accs = read_tsv_column(repset, args.column)
        source = f"{repset}:{args.column}"
    return unique_accessions(accs), source


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("accessions", nargs="*", help="Explicit UniProt accessions")
    parser.add_argument("--list", type=Path, help="File with one accession per line")
    parser.add_argument("--tsv", type=Path, help="TSV to read accessions from")
    parser.add_argument("--column", default="accession", help="Column in --tsv/repset metadata")
    parser.add_argument(
        "--repset-metadata",
        type=Path,
        default=DEFAULT_REPSET_METADATA,
        help="Default repset metadata when no explicit accessions/list/TSV are given",
    )
    parser.add_argument(
        "--pre-qc-repset",
        action="store_true",
        help="Use results/msa_OGT/repset_metadata.tsv instead of repset_metadata_qc.tsv",
    )
    parser.add_argument("--dest", type=Path, default=AFDB_DIR, help="Destination directory")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Per-accession status TSV")
    parser.add_argument("--workers", type=int, default=20, help="Parallel workers")
    parser.add_argument("--timeout", type=float, default=60.0, help="HTTP timeout per accession in seconds")
    parser.add_argument("--force", action="store_true", help="Re-download existing AFDB files instead of skipping them")
    parser.add_argument("--dry-run", action="store_true", help="Resolve accessions without downloading")
    parser.add_argument(
        "--strict-missing",
        action="store_true",
        help="Exit nonzero if any accession is absent from AlphaFold DB",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    accessions, source = collect_accessions(args)
    if not accessions:
        raise SystemExit("No AFDB accessions resolved.")

    print(f"AFDB accession source: {source}")
    print(f"Resolved accessions: {len(accessions)}")
    print(f"Destination: {args.dest}")
    print(f"Manifest: {args.manifest}")

    if args.dry_run:
        rows = [
            {
                "accession": acc,
                "status": "dry-run",
                "path": str(_out_path(acc, args.dest)),
                "url": AFDB_URL.format(acc=acc),
                "error": "",
            }
            for acc in accessions
        ]
        write_manifest(args.manifest, rows)
        print(f"Dry run only. Wrote: {args.manifest}")
        return

    rows = fetch_records(accessions, args.workers, args.dest, True, args.timeout, args.force)
    write_manifest(args.manifest, rows)
    counts = count_statuses(rows)
    print(
        "\nSummary: "
        f"downloaded={counts['ok']}  already-present={counts['skip']}  "
        f"not-in-AFDB={counts['not-in-afdb']}  error={counts['error']}  "
        f"of {len(accessions)} accessions"
    )
    print(f"Wrote: {args.manifest}")
    fail_for_missing = args.strict_missing and counts["not-in-afdb"]
    sys.exit(1 if (counts["error"] or fail_for_missing) else 0)


if __name__ == "__main__":
    main()
