#!/usr/bin/env python3
"""Step 02: align the query plus retained homologs with MAFFT."""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=cfg.OUTPUT_DIR / "repset.fa")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR)
    parser.add_argument("--alignment", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=cfg.N_THREADS)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def query_column_map(entries: list[tuple[str, str, str]]) -> list[dict[str, object]]:
    if not entries:
        raise SystemExit("ERROR: alignment is empty")
    query_id, _header, query_aln = entries[0]
    qpos = 0
    rows = []
    for aln_col, aa in enumerate(query_aln):
        if aa.upper() in cfg.AA_SET:
            qpos += 1
            rows.append({
                "alignment_col": aln_col,
                "query_id": query_id,
                "query_pos": qpos,
                "query_aa": aa.upper(),
                "is_query_gap": "no",
            })
        else:
            rows.append({
                "alignment_col": aln_col,
                "query_id": query_id,
                "query_pos": "",
                "query_aa": "",
                "is_query_gap": "yes",
            })
    return rows


def write_summary(path: Path, entries: list[tuple[str, str, str]], alignment: Path) -> None:
    width = len(entries[0][2]) if entries else 0
    nongap = sum(sum(aa.upper() in cfg.AA_SET for aa in seq) for _sid, _hdr, seq in entries)
    occupancy = nongap / (len(entries) * width) if entries and width else 0.0
    rows = [{
        "n_sequences": len(entries),
        "alignment_width": width,
        "mean_occupancy": f"{occupancy:.4f}",
        "alignment": str(alignment),
    }]
    cfg.write_tsv(path, rows, ["n_sequences", "alignment_width", "mean_occupancy", "alignment"])


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"ERROR: input FASTA not found: {args.input}")
    out_alignment = args.alignment or args.out_dir / "repset_aligned.fa"
    if out_alignment.exists() and not args.force:
        print(f"Using existing alignment: {out_alignment}")
    else:
        mafft = cfg.resolve_bin("mafft")
        out_alignment.parent.mkdir(parents=True, exist_ok=True)
        cmd = [mafft, "--auto", "--thread", str(args.threads), str(args.input)]
        print("Running: " + " ".join(cmd))
        with out_alignment.open("w") as out:
            completed = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        if completed.returncode != 0:
            print(completed.stderr[-2000:])
            raise SystemExit(f"MAFFT failed with exit code {completed.returncode}")
    entries = cfg.read_fasta(out_alignment)
    map_rows = query_column_map(entries)
    map_path = args.out_dir / "query_column_map.tsv"
    fields = ["alignment_col", "query_id", "query_pos", "query_aa", "is_query_gap"]
    cfg.write_tsv(map_path, map_rows, fields)
    write_summary(args.out_dir / "alignment_summary.tsv", entries, out_alignment)
    print(f"Saved: {out_alignment}")
    print(f"Saved: {map_path}")


if __name__ == "__main__":
    main()
