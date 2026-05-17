#!/usr/bin/env python3
"""Step 09: summarize structure scores by OGT clade or representative MSA group."""
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=cfg.OUTPUT_DIR / "repset_metadata.tsv")
    parser.add_argument("--scores", type=Path, default=cfg.OUTPUT_DIR / "structure_scores" / "global_scores.tsv")
    parser.add_argument("--called-clades", type=Path, default=None)
    parser.add_argument("--query-pdb", type=Path, default=cfg.STRUCTURE_DIR / "query.pdb")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "vulnerability")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: str) -> float | None:
    try:
        return float(value) if value != "" else None
    except ValueError:
        return None


def mean_text(values: list[float]) -> str:
    return f"{statistics.mean(values):.6g}" if values else ""


def group_members(metadata_rows: list[dict[str, str]], called_clades: list[dict[str, str]]) -> dict[str, set[str]]:
    if called_clades:
        groups = {}
        for row in called_clades:
            label = row.get("clade_id", "")
            category = row.get("called_category", "")
            group = f"{label}_{category}" if category else label
            tips = {tip for tip in row.get("tips", "").split(",") if tip}
            if group and tips:
                groups[group] = tips
        if groups:
            return groups
    return {"all_representatives": {row["id"] for row in metadata_rows if row.get("id") and row.get("selected_for_msa", "yes") == "yes"}}


def score_id_for_metadata(row: dict[str, str]) -> str:
    return row.get("accession", "") or cfg.normalize_accession(row.get("id", ""))


def summarize_group(group: str, members: set[str], metadata_by_id: dict[str, dict[str, str]], scores_by_id: dict[str, dict[str, str]]) -> dict[str, str]:
    score_rows = []
    for sid in members:
        meta = metadata_by_id.get(sid, {})
        score_id = score_id_for_metadata(meta)
        if score_id in scores_by_id:
            score_rows.append(scores_by_id[score_id])
    camsol = [value for row in score_rows if (value := safe_float(row.get("global_structural_solubility_score", ""))) is not None]
    camsol_v2 = [value for row in score_rows if (value := safe_float(row.get("v2_structural_solubility_score", ""))) is not None]
    a3d = [value for row in score_rows if (value := safe_float(row.get("positive_aggrescan3d_burden", ""))) is not None]
    return {
        "group": group,
        "n_msa_members": str(len(members)),
        "n_scored_structures": str(len(score_rows)),
        "mean_global_structural_solubility_score": mean_text(camsol),
        "mean_v2_structural_solubility_score": mean_text(camsol_v2),
        "mean_positive_aggrescan3d_burden": mean_text(a3d),
    }


def main() -> None:
    args = parse_args()
    metadata_rows = read_tsv(args.metadata)
    score_rows = read_tsv(args.scores)
    called_clades = read_tsv(args.called_clades) if args.called_clades else []
    metadata_by_id = {row["id"]: row for row in metadata_rows if row.get("id")}
    scores_by_id = {row["structure_id"]: row for row in score_rows if row.get("structure_id")}

    groups = group_members(metadata_rows, called_clades)
    rows = [summarize_group(group, members, metadata_by_id, scores_by_id) for group, members in sorted(groups.items())]

    query_stem = args.query_pdb.stem
    query_rows = [
        row for row in score_rows
        if row.get("structure_id") == query_stem or Path(row.get("source_path", "")).name == args.query_pdb.name
    ]
    if query_rows:
        query = query_rows[0]
        for row in rows:
            for source, target in [
                ("global_structural_solubility_score", "query_global_structural_solubility_score"),
                ("v2_structural_solubility_score", "query_v2_structural_solubility_score"),
                ("positive_aggrescan3d_burden", "query_positive_aggrescan3d_burden"),
            ]:
                row[target] = query.get(source, "")
    fields = [
        "group", "n_msa_members", "n_scored_structures",
        "mean_global_structural_solubility_score", "query_global_structural_solubility_score",
        "mean_v2_structural_solubility_score", "query_v2_structural_solubility_score",
        "mean_positive_aggrescan3d_burden", "query_positive_aggrescan3d_burden",
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cfg.write_tsv(args.out_dir / "group_score_summary.tsv", rows, fields)
    member_rows = [
        {"group": group, "id": sid}
        for group, members in sorted(groups.items())
        for sid in sorted(members)
    ]
    cfg.write_tsv(args.out_dir / "group_members.tsv", member_rows, ["group", "id"])
    print(f"Saved: {args.out_dir / 'group_score_summary.tsv'} ({len(rows)} groups)")


if __name__ == "__main__":
    main()
