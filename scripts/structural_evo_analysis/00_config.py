"""Shared helpers for the general structural-evolution MSA pipeline."""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
import statistics
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = Path(os.environ.get("SEA_OUT_DIR", PROJECT_ROOT / "results" / "structural_evo_analysis"))
STRUCTURE_DIR = Path(os.environ.get("SEA_STRUCTURE_DIR", PROJECT_ROOT / "structures" / "structural_evo_analysis"))
DEFAULT_QUERY_FASTA = PROJECT_ROOT / "sequences" / "photoHymenobact.fa"
DEFAULT_OGT_TSV = PROJECT_ROOT / "data" / "growth_temp_dataset_OGTFinder.tsv"
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
    project_env = PROJECT_ROOT / "envs" / "structural_evo" / "bin" / name
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


def extract_taxid(header: str) -> str:
    match = re.search(r"\bTaxID=(\d+)\b", header)
    return match.group(1) if match else ""


def extract_organism(header: str) -> str:
    if "Tax=" in header:
        text = header.split("Tax=", 1)[1]
        if " TaxID=" in text:
            text = text.split(" TaxID=", 1)[0]
        return " ".join(text.split())
    if "[" in header:
        return " ".join(header.split("[")[-1].rstrip("]").split())
    return ""


def parse_temp_value(raw: str) -> tuple[float | None, str, bool]:
    text = (raw or "").strip().replace("−", "-")
    if not text:
        return None, "missing", False
    range_match = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*-\s*(-?\d+(?:\.\d+)?)\s*", text)
    if range_match:
        lo = float(range_match.group(1))
        hi = float(range_match.group(2))
        return (lo + hi) / 2.0, "range_midpoint", True
    point_match = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*", text)
    if point_match:
        return float(point_match.group(1)), "point", False
    nums = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", text)]
    if nums:
        return statistics.median(nums), "numeric_median_fallback", len(nums) > 1
    return None, "unparsed", False


def regime_of(ogt: float | None, psychro_max: float = 20.0, thermo_min: float = 45.0) -> str:
    if ogt is None:
        return ""
    if ogt < psychro_max:
        return "psychro"
    if ogt >= thermo_min:
        return "thermo"
    return "meso"


class OGTFinderLookup:
    """Taxid-first lookup for the committed growth-temperature table."""

    fields = [
        "ogt", "regime", "ogt_match_type", "ogt_taxid", "ogt_species_id",
        "ogt_species", "ogt_raw_temps", "ogt_sources", "ogt_types",
        "ogt_source_ids", "ogt_parse_modes", "ogt_has_range", "ogt_range_count",
        "ogt_row_count",
    ]

    def __init__(self, path: Path = DEFAULT_OGT_TSV):
        self.path = path
        self.taxid_to_species_id: dict[str, str] = {}
        self.optimum_by_taxid: dict[str, list[dict[str, object]]] = {}
        self.optimum_by_species_id: dict[str, list[dict[str, object]]] = {}
        self.optimum_by_species_name: dict[str, list[dict[str, object]]] = {}
        if path.exists():
            self._load()

    def _load(self) -> None:
        with self.path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                taxid = (row.get("ncbiTaxID_new") or "").strip()
                species_id = (row.get("species_id") or "").strip()
                species = " ".join((row.get("species") or "").split())
                if taxid and species_id:
                    self.taxid_to_species_id.setdefault(taxid, species_id)
                if (row.get("Type") or "").strip().lower() != "optimum":
                    continue
                temp, parse_mode, is_range = parse_temp_value(row.get("Temp", ""))
                if temp is None:
                    continue
                record = {
                    "taxid": taxid,
                    "species_id": species_id,
                    "species": species,
                    "temp": temp,
                    "raw_temp": row.get("Temp", ""),
                    "source": row.get("Source", ""),
                    "type": row.get("Type", ""),
                    "source_id": row.get("Source_ID", ""),
                    "parse_mode": parse_mode,
                    "is_range": is_range,
                }
                if taxid:
                    self.optimum_by_taxid.setdefault(taxid, []).append(record)
                if species_id:
                    self.optimum_by_species_id.setdefault(species_id, []).append(record)
                if species:
                    self.optimum_by_species_name.setdefault(species, []).append(record)

    def lookup(self, taxid: str, organism: str = "") -> dict[str, str]:
        taxid = (taxid or "").strip()
        organism = " ".join((organism or "").split())
        records: list[dict[str, object]] = []
        match_type = "none"
        match_taxid = ""

        if taxid and taxid in self.optimum_by_taxid:
            records = self.optimum_by_taxid[taxid]
            match_type = "exact_taxid"
            match_taxid = taxid
        elif taxid:
            species_id = self.taxid_to_species_id.get(taxid, taxid)
            if species_id in self.optimum_by_species_id:
                records = self.optimum_by_species_id[species_id]
                match_type = "species_taxid"
                match_taxid = taxid

        if not records and organism:
            records = self.optimum_by_species_name.get(organism, [])
            match_type = "name_exact" if records else "none"
            if not records:
                parts = organism.split()
                if len(parts) >= 2 and parts[1] != "sp.":
                    records = self.optimum_by_species_name.get(f"{parts[0]} {parts[1]}", [])
                    match_type = "name_binomial" if records else "none"

        if not records:
            return {field: "" for field in self.fields}

        temps = [float(record["temp"]) for record in records]
        ogt = statistics.median(temps)
        return {
            "ogt": f"{ogt:.1f}",
            "regime": regime_of(ogt),
            "ogt_match_type": match_type,
            "ogt_taxid": match_taxid,
            "ogt_species_id": ";".join(sorted({str(record["species_id"]) for record in records if record["species_id"]})),
            "ogt_species": ";".join(sorted({str(record["species"]) for record in records if record["species"]})),
            "ogt_raw_temps": ";".join(str(record["raw_temp"]) for record in records),
            "ogt_sources": ";".join(sorted({str(record["source"]) for record in records if record["source"]})),
            "ogt_types": ";".join(sorted({str(record["type"]) for record in records if record["type"]})),
            "ogt_source_ids": ";".join(str(record["source_id"]) for record in records if record["source_id"]),
            "ogt_parse_modes": ";".join(sorted({str(record["parse_mode"]) for record in records})),
            "ogt_has_range": "yes" if any(record["is_range"] for record in records) else "no",
            "ogt_range_count": str(sum(1 for record in records if record["is_range"])),
            "ogt_row_count": str(len(records)),
        }


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
