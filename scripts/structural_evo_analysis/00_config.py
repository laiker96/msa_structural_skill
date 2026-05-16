"""Shared helpers for the general structural-evolution MSA pipeline."""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(os.environ.get("SEA_OUT_DIR", PROJECT_ROOT / "results" / "structural_evo_analysis"))
STRUCTURE_DIR = Path(os.environ.get("SEA_STRUCTURE_DIR", PROJECT_ROOT / "structures" / "structural_evo_analysis"))
UNIREF_DIR = Path(os.environ.get("UNIREF_DIR", Path.home() / "databases"))
DEFAULT_DB_NAME = os.environ.get("SEA_DB", os.environ.get("UNIREF_DB", "uniref90"))
DEFAULT_DB_FASTA = Path(
    os.environ.get("SEA_DB_FASTA", UNIREF_DIR / DEFAULT_DB_NAME / f"{DEFAULT_DB_NAME}.fasta.gz")
)
DEFAULT_DB_MMSEQS = Path(
    os.environ.get("SEA_DB_MMSEQS", UNIREF_DIR / DEFAULT_DB_NAME / f"{DEFAULT_DB_NAME}_db")
)
N_THREADS = int(os.environ.get("SEA_THREADS", "8"))

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


def resolve_bin(name: str) -> str:
    candidate = Path(sys.executable).parent / name
    if candidate.exists():
        return str(candidate)
    project_env = PROJECT_ROOT / "envs" / "ankros" / "bin" / name
    if project_env.exists():
        return str(project_env)
    found = shutil.which(name)
    return found or name


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    entries: list[tuple[str, str, str]] = []
    header = None
    parts: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header.split()[0], header, "".join(parts)))
                header = line[1:].strip()
                parts = []
            elif header is not None:
                parts.append(line)
    if header is not None:
        entries.append((header.split()[0], header, "".join(parts)))
    return entries


def write_fasta(entries: list[tuple[str, str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for _sid, header, seq in entries:
            handle.write(f">{header}\n")
            clean = seq.replace("*", "")
            for start in range(0, len(clean), 80):
                handle.write(clean[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_sequence(seq: str) -> str:
    return "".join(aa for aa in seq.upper().replace("*", "") if aa in AA_SET)


def normalize_accession(value: str) -> str:
    text = (value or "").strip()
    if not text:
        return ""
    if text.startswith("UniRef") and "_" in text:
        text = text.split("_", 1)[1]
    if "|" in text:
        parts = [part for part in text.split("|") if part]
        if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
            text = parts[1]
        else:
            text = parts[-1]
    match = re.search(r"\b([A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z0-9]{10})\b", text)
    return match.group(1) if match else text


def unique_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def query_entry(path: Path) -> tuple[str, str, str]:
    entries = read_fasta(path)
    if len(entries) != 1:
        raise SystemExit(f"ERROR: query FASTA must contain exactly one sequence: {path}")
    sid, header, seq = entries[0]
    clean = clean_sequence(seq)
    if not clean:
        raise SystemExit(f"ERROR: query sequence has no standard amino acids: {path}")
    return sid, header, clean

