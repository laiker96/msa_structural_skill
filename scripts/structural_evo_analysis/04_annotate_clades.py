#!/usr/bin/env python3
"""Step 04: summarize and annotate tree clades from user-provided metadata."""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=cfg.OUTPUT_DIR / "tree" / "query_msa.treefile")
    parser.add_argument("--metadata", type=Path, required=True, help="TSV/CSV with tip IDs and trait values.")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "metadata_clades")
    parser.add_argument("--id-column", default="id")
    parser.add_argument("--trait-column", required=True)
    parser.add_argument("--delimiter", default="\t", help="Metadata delimiter. Default: tab.")
    parser.add_argument("--trait-type", choices=["continuous", "categorical"], default="continuous")
    parser.add_argument("--low-threshold", type=float, default=None)
    parser.add_argument("--high-threshold", type=float, default=None)
    parser.add_argument("--low-label", default="low")
    parser.add_argument("--mid-label", default="mid")
    parser.add_argument("--high-label", default="high")
    parser.add_argument("--min-labelled", type=int, default=5)
    parser.add_argument("--min-fraction", type=float, default=0.80)
    parser.add_argument(
        "--match-accession",
        action="store_true",
        help="Also match tree tips against normalized accession values from the metadata ID column.",
    )
    return parser.parse_args()


def parse_float(text: str) -> float | None:
    try:
        value = float((text or "").strip())
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def category_for(value: str, args: argparse.Namespace) -> str:
    if args.trait_type == "categorical":
        return value.strip()
    numeric = parse_float(value)
    if numeric is None or args.low_threshold is None or args.high_threshold is None:
        return ""
    if numeric <= args.low_threshold:
        return args.low_label
    if numeric >= args.high_threshold:
        return args.high_label
    return args.mid_label


def load_metadata(args: argparse.Namespace) -> tuple[dict[str, dict[str, str]], list[str]]:
    with args.metadata.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=args.delimiter)
        fields = reader.fieldnames or []
        missing = [name for name in [args.id_column, args.trait_column] if name not in fields]
        if missing:
            raise SystemExit(f"ERROR: metadata missing columns {missing}; available columns: {fields}")
        by_id: dict[str, dict[str, str]] = {}
        for row in reader:
            sid = row.get(args.id_column, "").strip()
            if not sid:
                continue
            row["_trait_value"] = row.get(args.trait_column, "").strip()
            row["_trait_category"] = category_for(row["_trait_value"], args)
            by_id[sid] = row
            if args.match_accession:
                by_id.setdefault(cfg.normalize_accession(sid), row)
        return by_id, fields


def terminal_names(clade) -> list[str]:
    return [tip.name for tip in clade.get_terminals() if tip.name]


def summarize_values(rows: list[dict[str, str]]) -> dict[str, str]:
    values = [parse_float(row.get("_trait_value", "")) for row in rows]
    clean = sorted(value for value in values if value is not None)
    if not clean:
        return {"trait_mean": "", "trait_median": "", "trait_min": "", "trait_max": ""}
    return {
        "trait_mean": f"{sum(clean) / len(clean):.6g}",
        "trait_median": f"{statistics.median(clean):.6g}",
        "trait_min": f"{clean[0]:.6g}",
        "trait_max": f"{clean[-1]:.6g}",
    }


def clade_rows(tree, metadata: dict[str, dict[str, str]], args: argparse.Namespace) -> list[dict[str, object]]:
    rows = []
    clade_index = 0
    for clade in tree.find_clades(order="level"):
        tips = terminal_names(clade)
        if len(tips) <= 1:
            continue
        labelled_rows = [metadata[tip] for tip in tips if tip in metadata and metadata[tip].get("_trait_value")]
        categories = [row.get("_trait_category", "") for row in labelled_rows if row.get("_trait_category", "")]
        counts = Counter(categories)
        top_category, top_count = ("", 0)
        if counts:
            top_category, top_count = counts.most_common(1)[0]
        labelled = len(labelled_rows)
        fraction = top_count / len(categories) if categories else 0.0
        called = (
            bool(top_category)
            and labelled >= args.min_labelled
            and fraction >= args.min_fraction
        )
        clade_index += 1
        clade_id = f"clade_{clade_index:04d}"
        row = {
            "clade_id": clade_id,
            "called_category": top_category if called else "",
            "top_category": top_category,
            "top_category_fraction": f"{fraction:.4f}" if categories else "",
            "n_tips": len(tips),
            "n_labelled": labelled,
            "n_categorized": len(categories),
            "category_counts": ";".join(f"{key}:{counts[key]}" for key in sorted(counts)),
            "tips": ",".join(tips),
        }
        row.update(summarize_values(labelled_rows))
        rows.append(row)
        if called:
            clade.name = f"{clade_id}_{top_category}"
    return rows


def write_tip_metadata(path: Path, tree, metadata: dict[str, dict[str, str]]) -> None:
    rows = []
    for tip in terminal_names(tree.root):
        row = metadata.get(tip, {})
        rows.append({
            "id": tip,
            "trait_value": row.get("_trait_value", ""),
            "trait_category": row.get("_trait_category", ""),
            "metadata_found": "yes" if row else "no",
        })
    cfg.write_tsv(path, rows, ["id", "trait_value", "trait_category", "metadata_found"])


def main() -> None:
    args = parse_args()
    if not args.tree.exists():
        raise SystemExit(f"ERROR: tree not found: {args.tree}")
    if args.trait_type == "continuous" and (args.low_threshold is None) != (args.high_threshold is None):
        raise SystemExit("ERROR: provide both --low-threshold and --high-threshold, or neither.")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    metadata, _fields = load_metadata(args)
    tree = Phylo.read(args.tree, "newick")
    rows = clade_rows(tree, metadata, args)
    fields = [
        "clade_id", "called_category", "top_category", "top_category_fraction",
        "n_tips", "n_labelled", "n_categorized", "trait_mean", "trait_median",
        "trait_min", "trait_max", "category_counts", "tips",
    ]
    cfg.write_tsv(args.out_dir / "clade_annotations.tsv", rows, fields)
    called_rows = [row for row in rows if row.get("called_category")]
    cfg.write_tsv(args.out_dir / "called_clades.tsv", called_rows, fields)
    write_tip_metadata(args.out_dir / "tip_metadata.tsv", tree, metadata)
    Phylo.write(tree, args.out_dir / "annotated_clades.nwk", "newick")
    print(f"Saved: {args.out_dir / 'clade_annotations.tsv'} ({len(rows)} clades)")
    print(f"Saved: {args.out_dir / 'called_clades.tsv'} ({len(called_rows)} called clades)")
    print(f"Saved: {args.out_dir / 'annotated_clades.nwk'}")


if __name__ == "__main__":
    main()
