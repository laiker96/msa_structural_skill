#!/usr/bin/env python3
"""Build a compact taxid-keyed OGT summary from the raw OGTFinder table."""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

cfg = import_module("00_config")


FIELDS = [
    "taxid",
    "name",
    "rank",
    "species_taxid",
    "species_name",
    "ogt",
    "regime",
    "ogt_match_type",
    "ogt_row_count",
    "ogt_temp_min",
    "ogt_temp_max",
    "ogt_raw_temps",
    "ogt_sources",
    "ogt_source_ids",
    "ogt_parse_modes",
    "ogt_has_range",
    "ogt_range_count",
    "ogt_source_taxids",
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
]


def clean_id(value: str) -> str:
    text = (value or "").strip()
    return "" if text.lower() == "nan" else text


def unique(values: list[str]) -> list[str]:
    return cfg.unique_preserve_order([value for value in values if value])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Raw OGTFinder TSV.")
    parser.add_argument("--output", type=Path, default=cfg.DEFAULT_OGT_SUMMARY_TSV)
    return parser.parse_args()


def summarize_records(taxid: str, info: dict[str, str], records: list[dict[str, object]], match_type: str) -> dict[str, str]:
    temps = [float(record["temp"]) for record in records]
    ogt = statistics.median(temps)
    species_ids = unique([str(record["species_id"]) for record in records] or [clean_id(info.get("species_id", ""))])
    species_names = unique([str(record["species"]) for record in records] or [" ".join((info.get("species") or "").split())])
    return {
        "taxid": taxid,
        "name": " ".join((info.get("name_new") or "").split()),
        "rank": clean_id(info.get("rank_new", "")),
        "species_taxid": clean_id(info.get("species_id", "")),
        "species_name": " ".join((info.get("species") or "").split()),
        "ogt": f"{ogt:.1f}",
        "regime": cfg.regime_of(ogt),
        "ogt_match_type": match_type,
        "ogt_row_count": str(len(records)),
        "ogt_temp_min": f"{min(temps):.1f}",
        "ogt_temp_max": f"{max(temps):.1f}",
        "ogt_raw_temps": ";".join(str(record["raw_temp"]) for record in records),
        "ogt_sources": ";".join(sorted(set(str(record["source"]) for record in records if record["source"]))),
        "ogt_source_ids": ";".join(str(record["source_id"]) for record in records if record["source_id"]),
        "ogt_parse_modes": ";".join(sorted(set(str(record["parse_mode"]) for record in records))),
        "ogt_has_range": "yes" if any(record["is_range"] for record in records) else "no",
        "ogt_range_count": str(sum(1 for record in records if record["is_range"])),
        "ogt_source_taxids": ";".join(unique([str(record["taxid"]) for record in records])),
        "superkingdom": clean_id(info.get("superkingdom", "")),
        "phylum": clean_id(info.get("phylum", "")),
        "class": clean_id(info.get("class", "")),
        "order": clean_id(info.get("order", "")),
        "family": clean_id(info.get("family", "")),
        "genus": clean_id(info.get("genus", "")),
        # These are not written, but help keep the transformation auditable if
        # this function is reused interactively.
        "_ogt_species_id": ";".join(species_ids),
        "_ogt_species": ";".join(species_names),
    }


def main() -> None:
    args = parse_args()
    taxid_to_info: dict[str, dict[str, str]] = {}
    taxid_to_species_id: dict[str, str] = {}
    optimum_by_taxid: dict[str, list[dict[str, object]]] = defaultdict(list)
    optimum_by_species_id: dict[str, list[dict[str, object]]] = defaultdict(list)

    with args.input.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxid = clean_id(row.get("ncbiTaxID_new", ""))
            species_id = clean_id(row.get("species_id", ""))
            species = " ".join((row.get("species") or "").split())
            if not taxid:
                continue
            taxid_to_info.setdefault(taxid, row)
            if species_id:
                taxid_to_species_id.setdefault(taxid, species_id)
            if (row.get("Type") or "").strip().lower() != "optimum":
                continue
            temp, parse_mode, is_range = cfg.parse_temp_value(row.get("Temp", ""))
            if temp is None:
                continue
            record = {
                "taxid": taxid,
                "species_id": species_id,
                "species": species,
                "temp": temp,
                "raw_temp": row.get("Temp", ""),
                "source": row.get("Source", ""),
                "source_id": row.get("Source_ID", ""),
                "parse_mode": parse_mode,
                "is_range": is_range,
            }
            optimum_by_taxid[taxid].append(record)
            if species_id:
                optimum_by_species_id[species_id].append(record)

    rows = []
    for taxid in sorted(taxid_to_info, key=lambda value: int(value)):
        species_id = taxid_to_species_id.get(taxid, "")
        if taxid in optimum_by_taxid:
            records = optimum_by_taxid[taxid]
            match_type = "exact_taxid"
        elif species_id in optimum_by_species_id:
            records = optimum_by_species_id[species_id]
            match_type = "species_taxid"
        else:
            continue
        rows.append(summarize_records(taxid, taxid_to_info[taxid], records, match_type))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})
    print(f"Saved {len(rows)} taxid OGT summaries: {args.output}")


if __name__ == "__main__":
    main()
