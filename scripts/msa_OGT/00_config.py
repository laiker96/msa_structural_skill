"""
Shared configuration and metadata helpers for the msa_OGT pipeline.

This pipeline treats the UniRef search result as the sequence universe and
joins OGT metadata by taxid before downstream filtering/alignment. OGTFinder
temperature ranges are preserved in the output and converted to midpoints only
for numeric regime assignment.
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import statistics
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
QUERY_FASTA = PROJECT_ROOT / "test" / "photoHymenobact.fa"

UNIREF_DIR = Path(os.environ.get("UNIREF_DIR", Path.home() / "databases"))
DATABASES = {
    "uniref90": {
        "fasta": UNIREF_DIR / "uniref90" / "uniref90.fasta.gz",
        "mmseqs_db": UNIREF_DIR / "uniref90" / "uniref90_db",
    },
    "uniref50": {
        "fasta": UNIREF_DIR / "uniref50" / "uniref50.fasta.gz",
        "mmseqs_db": UNIREF_DIR / "uniref50" / "uniref50_db",
    },
}
ACTIVE_DB_NAME = os.environ.get("UNIREF_DB", "uniref90")
if ACTIVE_DB_NAME not in DATABASES:
    raise ValueError(f"Unknown UNIREF_DB={ACTIVE_DB_NAME}; choose from {sorted(DATABASES)}")
ACTIVE_DB_FASTA = DATABASES[ACTIVE_DB_NAME]["fasta"]
ACTIVE_DB_MMSEQS = DATABASES[ACTIVE_DB_NAME]["mmseqs_db"]

OUTPUT_DIR = Path(os.environ.get("MSA_OGT_OUT_DIR", PROJECT_ROOT / "results" / "msa_OGT"))
INTER_DIR = OUTPUT_DIR

REF_DIR = PROJECT_ROOT / "test"
UNIPROT_CLASSI_CPD = REF_DIR / "uniprot_classI_cpd.fa"
CLASSII_KNOWN = REF_DIR / "classII_known.fa"
OTHER_KNOWN = REF_DIR / "other_known.fa"

OGTFINDER_TSV = PROJECT_ROOT / "data" / "growth_temp_dataset_OGTFinder.tsv"

# Broad audit search. Tight production thresholds should be applied later from
# the metadata table rather than during search.
MMSEQS_SENSITIVITY = float(os.environ.get("MSA_OGT_MMSEQS_S", "7.5"))
MMSEQS_EVALUE = os.environ.get("MSA_OGT_MMSEQS_E", "1e-5")
MMSEQS_MIN_SEQ_ID = float(os.environ.get("MSA_OGT_MIN_SEQ_ID", "0.05"))
MMSEQS_COVERAGE = float(os.environ.get("MSA_OGT_COVERAGE", "0.50"))
MMSEQS_MAX_SEQS = int(os.environ.get("MSA_OGT_MAX_SEQS", "50000"))

# Post-search defaults for the class-I candidate set. These do not affect the
# raw annotated hit table.
MIN_LENGTH = int(os.environ.get("MSA_OGT_MIN_LENGTH", "350"))
MAX_LENGTH = int(os.environ.get("MSA_OGT_MAX_LENGTH", "600"))
POST_MIN_QCOV = float(os.environ.get("MSA_OGT_POST_MIN_QCOV", "0.50"))
POST_MIN_IDENTITY = float(os.environ.get("MSA_OGT_POST_MIN_IDENTITY", "5.0"))
POST_MAX_IDENTITY = float(os.environ.get("MSA_OGT_POST_MAX_IDENTITY", "98.0"))

N_THREADS = int(os.environ.get("MSA_OGT_THREADS", "8"))
DEFAULT_PRECLUSTER_IDENTITY = "0" if ACTIVE_DB_NAME == "uniref50" else "0.95"
PRECLUSTER_IDENTITY = float(os.environ.get("MSA_OGT_PRECLUSTER_ID", DEFAULT_PRECLUSTER_IDENTITY))
CLUSTER_COV = float(os.environ.get("MSA_OGT_CLUSTER_COV", "0.80"))
TARGET_REPSET = int(os.environ.get("MSA_OGT_TARGET_REPSET", "800"))
REPSET_MESO_QUOTA = int(os.environ.get("MSA_OGT_MESO_QUOTA", "250"))
REPSET_THERMO_QUOTA = int(os.environ.get("MSA_OGT_THERMO_QUOTA", "200"))
REPSET_LOWID_CLUSTER = float(os.environ.get("MSA_OGT_TOPUP_CLUSTER_ID", "0.50"))
REPSET_MESO_LOWID_CLUSTER = float(os.environ.get("MSA_OGT_MESO_CLUSTER_ID", "0.60"))
REPSET_THERMO_LOWID_CLUSTER = float(os.environ.get("MSA_OGT_THERMO_CLUSTER_ID", "0.60"))

OGT_PSYCHRO_MAX = 20.0
OGT_THERMO_MIN = 45.0


def resolve_bin(name: str) -> str:
    candidate = Path(sys.executable).parent / name
    if candidate.exists():
        return str(candidate)
    project_env = PROJECT_ROOT / "envs" / "ankros" / "bin" / name
    if project_env.exists():
        return str(project_env)
    return name


def read_fasta(path: Path):
    entries = []
    with open(path) as f:
        hdr, parts = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if hdr is not None:
                    entries.append((hdr.split()[0], hdr, "".join(parts)))
                hdr, parts = line[1:], []
            elif hdr is not None:
                parts.append(line)
        if hdr is not None:
            entries.append((hdr.split()[0], hdr, "".join(parts)))
    return entries


def write_fasta(entries, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for _, hdr, seq in entries:
            f.write(f">{hdr}\n")
            seq = seq.replace("-", "").replace("*", "")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def open_fasta(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def extract_organism(header: str) -> str:
    if "Tax=" in header:
        tax = header.split("Tax=", 1)[1]
        if " TaxID=" in tax:
            tax = tax.split(" TaxID=", 1)[0]
        return tax.strip()
    if "[" in header:
        return header.split("[")[-1].rstrip("]").strip()
    return "Unknown"


def extract_genus(header: str) -> str:
    org = extract_organism(header)
    return org.split()[0] if org else "Unknown"


def extract_taxid(header: str) -> str:
    m = re.search(r"\bTaxID=(\d+)\b", header)
    return m.group(1) if m else ""


def extract_acc(sid: str) -> str:
    m = re.match(r"UniRef\d+_([A-Z0-9]+)", sid)
    if m:
        return m.group(1)
    return sid.split("|")[-1] if "|" in sid else sid


def parse_temp_value(raw: str):
    """Return (numeric_temp, parse_mode, is_range)."""
    s = (raw or "").strip().replace("−", "-")
    if not s:
        return None, "missing", False
    range_match = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*-\s*(-?\d+(?:\.\d+)?)\s*", s)
    if range_match:
        lo = float(range_match.group(1))
        hi = float(range_match.group(2))
        return (lo + hi) / 2.0, "range_midpoint", True
    point_match = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*", s)
    if point_match:
        return float(point_match.group(1)), "point", False
    nums = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", s)]
    if nums:
        return statistics.median(nums), "numeric_median_fallback", len(nums) > 1
    return None, "unparsed", False


def regime_of(ogt):
    if ogt is None:
        return ""
    if ogt < OGT_PSYCHRO_MAX:
        return "psychro"
    if ogt >= OGT_THERMO_MIN:
        return "thermo"
    return "meso"


class OGTFinderLookup:
    """Taxid-first OGTFinder lookup with raw temperature provenance."""

    fields = [
        "ogt", "regime", "ogt_match_type", "ogt_taxid", "ogt_species_id",
        "ogt_species", "ogt_raw_temps", "ogt_sources", "ogt_types",
        "ogt_source_ids", "ogt_parse_modes", "ogt_has_range", "ogt_range_count",
        "ogt_row_count",
    ]

    def __init__(self, path: Path = OGTFINDER_TSV):
        self.path = path
        self.taxid_to_species_id = {}
        self.species_by_id = {}
        self.optimum_by_taxid = {}
        self.optimum_by_species_id = {}
        self.optimum_by_species_name = {}
        self._load()

    def _load(self):
        if not self.path.exists():
            raise FileNotFoundError(f"OGTFinder table not found: {self.path}")
        with open(self.path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                taxid = (row.get("ncbiTaxID_new") or "").strip()
                species_id = (row.get("species_id") or "").strip()
                species = " ".join((row.get("species") or "").split())
                if taxid and species_id:
                    self.taxid_to_species_id.setdefault(taxid, species_id)
                if species_id and species:
                    self.species_by_id.setdefault(species_id, species)
                if (row.get("Type") or "").strip().lower() != "optimum":
                    continue
                temp, parse_mode, is_range = parse_temp_value(row.get("Temp", ""))
                if temp is None:
                    continue
                rec = {
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
                    self.optimum_by_taxid.setdefault(taxid, []).append(rec)
                if species_id:
                    self.optimum_by_species_id.setdefault(species_id, []).append(rec)
                if species:
                    self.optimum_by_species_name.setdefault(species, []).append(rec)

    def lookup(self, taxid: str, organism: str = ""):
        taxid = (taxid or "").strip()
        organism = " ".join((organism or "").split())
        records = []
        match_type = "none"
        match_taxid = ""
        species_id = ""

        if taxid and taxid in self.optimum_by_taxid:
            records = self.optimum_by_taxid[taxid]
            match_type = "exact_taxid"
            match_taxid = taxid
            species_id = self.taxid_to_species_id.get(taxid, taxid)
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
                    binomial = f"{parts[0]} {parts[1]}"
                    records = self.optimum_by_species_name.get(binomial, [])
                    match_type = "name_binomial" if records else "none"

        if not records:
            return {field: "" for field in self.fields}

        temps = [rec["temp"] for rec in records]
        ogt = statistics.median(temps)
        species_ids = sorted({rec["species_id"] for rec in records if rec["species_id"]})
        species_names = sorted({rec["species"] for rec in records if rec["species"]})
        return {
            "ogt": f"{ogt:.1f}",
            "regime": regime_of(ogt),
            "ogt_match_type": match_type,
            "ogt_taxid": match_taxid,
            "ogt_species_id": ";".join(species_ids),
            "ogt_species": ";".join(species_names),
            "ogt_raw_temps": ";".join(rec["raw_temp"] for rec in records),
            "ogt_sources": ";".join(sorted({rec["source"] for rec in records if rec["source"]})),
            "ogt_types": ";".join(sorted({rec["type"] for rec in records if rec["type"]})),
            "ogt_source_ids": ";".join(rec["source_id"] for rec in records if rec["source_id"]),
            "ogt_parse_modes": ";".join(sorted({rec["parse_mode"] for rec in records})),
            "ogt_has_range": "yes" if any(rec["is_range"] for rec in records) else "no",
            "ogt_range_count": str(sum(1 for rec in records if rec["is_range"])),
            "ogt_row_count": str(len(records)),
        }
