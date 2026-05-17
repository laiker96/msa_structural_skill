#!/usr/bin/env python3
"""Prepare or validate a local UniRef FASTA/MMseqs database for step 01."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


UNIPROT_UNIREF_BASE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a local UniRef database location, optionally download the "
            "FASTA from UniProt, and optionally build the MMseqs database prefix."
        )
    )
    parser.add_argument(
        "--db",
        "--db-name",
        dest="db_name",
        default=cfg.DEFAULT_DB_NAME,
        help="UniRef database: 50, 90, 100, uniref50, uniref90, or uniref100.",
    )
    parser.add_argument(
        "--uniref-dir",
        type=Path,
        default=cfg.UNIREF_DIR,
        help="Directory containing UniRef subdirectories. Defaults to $UNIREF_DIR or ~/databases.",
    )
    parser.add_argument(
        "--db-fasta",
        type=Path,
        default=None,
        help="Exact FASTA path. Defaults to <uniref-dir>/<db>/<db>.fasta.gz.",
    )
    parser.add_argument(
        "--db-mmseqs",
        type=Path,
        default=None,
        help="Exact MMseqs DB prefix. Defaults beside the FASTA as <fasta-stem>_db.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the UniRef FASTA if it is missing. Use only after user approval.",
    )
    parser.add_argument(
        "--create-mmseqs",
        action="store_true",
        help="Create the MMseqs database prefix if it is missing.",
    )
    parser.add_argument("--force-download", action="store_true", help="Redownload an existing FASTA.")
    parser.add_argument("--dry-run", action="store_true", help="Print the plan without downloading or creating DBs.")
    return parser.parse_args()


def default_db_fasta(uniref_dir: Path, db_name: str) -> Path:
    return uniref_dir / db_name / f"{db_name}.fasta.gz"


def url_for(db_name: str) -> str:
    return f"{UNIPROT_UNIREF_BASE_URL}/{db_name}/{db_name}.fasta.gz"


def mmseqs_exists(prefix: Path) -> bool:
    return Path(f"{prefix}.dbtype").exists()


def require_mmseqs() -> str:
    mmseqs = cfg.resolve_bin("mmseqs")
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit("ERROR: mmseqs not found. Run setup_envs.sh or add mmseqs to PATH.")
    return mmseqs


def download_fasta(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        print(f"FASTA already exists: {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"Downloading {url}")
    print(f"Destination: {dest}")
    try:
        with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    except urllib.error.URLError as exc:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: UniRef download failed: {exc}") from exc
    tmp.replace(dest)


def create_mmseqs_db(db_fasta: Path, db_mmseqs: Path) -> None:
    if mmseqs_exists(db_mmseqs):
        print(f"MMseqs database already exists: {db_mmseqs}")
        return
    if not db_fasta.exists():
        raise SystemExit(f"ERROR: FASTA not found, cannot create MMseqs DB: {db_fasta}")
    db_mmseqs.parent.mkdir(parents=True, exist_ok=True)
    mmseqs = require_mmseqs()
    subprocess.run([mmseqs, "createdb", str(db_fasta), str(db_mmseqs)], check=True)


def shell_quote(path: Path) -> str:
    text = str(path)
    return "'" + text.replace("'", "'\"'\"'") + "'"


def print_summary(db_name: str, db_fasta: Path, db_mmseqs: Path, url: str) -> None:
    print(f"UniRef database: {db_name}")
    print(f"FASTA: {db_fasta}")
    print(f"MMseqs prefix: {db_mmseqs}")
    print(f"Source URL: {url}")
    print("")
    print("Use these exports for the pipeline:")
    print(f"export SEA_DB={db_name}")
    print(f"export SEA_DB_FASTA={shell_quote(db_fasta)}")
    print(f"export SEA_DB_MMSEQS={shell_quote(db_mmseqs)}")


def main() -> None:
    args = parse_args()
    db_name = cfg.normalize_db_name(args.db_name)
    if db_name not in {"uniref50", "uniref90", "uniref100"}:
        raise SystemExit(f"ERROR: unsupported UniRef database name for download helper: {db_name}")
    db_fasta = args.db_fasta or default_db_fasta(args.uniref_dir, db_name)
    db_mmseqs = args.db_mmseqs or cfg.default_db_mmseqs(db_name, db_fasta)
    url = url_for(db_name)

    print_summary(db_name, db_fasta, db_mmseqs, url)

    if args.dry_run:
        return

    if not db_fasta.exists():
        if not args.download:
            print("")
            print("Database FASTA is missing.")
            print("Provide an existing path with --db-fasta, or rerun with --download after user approval.")
            raise SystemExit(2)
        download_fasta(url, db_fasta, force=args.force_download)
    elif args.download and args.force_download:
        download_fasta(url, db_fasta, force=True)

    if args.create_mmseqs:
        create_mmseqs_db(db_fasta, db_mmseqs)
    elif not mmseqs_exists(db_mmseqs):
        print("")
        print("MMseqs database prefix is missing.")
        print("Run again with --create-mmseqs before step 01, or let step 01 create it during the search.")


if __name__ == "__main__":
    main()
