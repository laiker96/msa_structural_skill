#!/usr/bin/env python3
"""
Step 09: concatenate HMM-anchored domains with the refined linker window.

Default splice:
  - HMM/QC query columns for ANKros 1-119
  - refined linker-window alignment for ANKros 120-215
  - HMM/QC query columns for ANKros 216-437

Using the wider refined window preserves linker boundary context and avoids
drawing an arbitrary line through insert-rich boundary blocks.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

QUERY_TAG = "photoHymenobact"
HMM_ALN = cfg.INTER_DIR / "repset_hmmalign_qc.fa"
LINKER_DIR = cfg.INTER_DIR / "linker_refined"
LINKER_ALN = LINKER_DIR / "linker_refined.fa"
LINKER_MAP = LINKER_DIR / "linker_refined_column_map.tsv"
OUT_FA = cfg.INTER_DIR / "repset_hmmalign_linker_refined.fa"
OUT_MAP = cfg.INTER_DIR / "repset_hmmalign_linker_refined_column_map.tsv"
OUT_SUMMARY = cfg.INTER_DIR / "repset_hmmalign_linker_refined_summary.tsv"

DOMAINS = {
    "antenna": (1, 130),
    "linker": (131, 205),
    "catalytic": (206, 437),
}


def read_fasta(path: Path):
    entries = []
    hdr, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    entries.append((hdr.split()[0], hdr, "".join(parts)))
                hdr, parts = line[1:].strip(), []
            elif line:
                parts.append(line.strip())
    if hdr is not None:
        entries.append((hdr.split()[0], hdr, "".join(parts)))
    return entries


def write_fasta(entries, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for _sid, hdr, seq in entries:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def find_query(entries):
    for sid, _hdr, seq in entries:
        if QUERY_TAG in sid:
            return sid, seq
    raise SystemExit(f"ERROR: query tag {QUERY_TAG!r} not found")


def query_pos_by_col(query_aln):
    pos_by_col = {}
    pos = 0
    for col, aa in enumerate(query_aln):
        if aa != "-":
            pos += 1
            pos_by_col[col] = pos
    return pos_by_col


def domain_of(qpos):
    if qpos == "":
        return "insert"
    qpos = int(qpos)
    for name, (start, end) in DOMAINS.items():
        if start <= qpos <= end:
            return name
    if qpos < DOMAINS["antenna"][0]:
        return "n_terminal"
    if qpos > DOMAINS["catalytic"][1]:
        return "c_terminal"
    return "interdomain_flank"


def parse_source_range(text):
    if not text:
        return None
    if "-" in text:
        left, right = text.split("-", 1)
        return int(left), int(right)
    value = int(text)
    return value, value


def read_linker_map(path: Path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def select_linker_columns(linker_map, splice_start, splice_end, full_window):
    selected = []
    for idx, row in enumerate(linker_map):
        qpos = row.get("qpos", "")
        if qpos:
            if splice_start <= int(qpos) <= splice_end:
                selected.append(idx)
            continue
        source_range = parse_source_range(row.get("source_range", ""))
        if full_window or (
            source_range
            and splice_start <= source_range[0]
            and source_range[1] <= splice_end
        ):
            selected.append(idx)
    return selected


def load_entries(hmm_path, linker_path):
    hmm_entries = read_fasta(hmm_path)
    linker_entries = read_fasta(linker_path)
    hmm_by_id = {sid: (hdr, seq) for sid, hdr, seq in hmm_entries}
    linker_by_id = {sid: (hdr, seq) for sid, hdr, seq in linker_entries}
    missing = [sid for sid in linker_by_id if sid not in hmm_by_id]
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} linker ids absent from HMM alignment")
    return hmm_entries, linker_entries, hmm_by_id, linker_by_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmm-aln", default=str(HMM_ALN))
    parser.add_argument("--linker-aln", default=str(LINKER_ALN))
    parser.add_argument("--linker-map", default=str(LINKER_MAP))
    parser.add_argument("--out-fa", default=str(OUT_FA))
    parser.add_argument("--out-map", default=str(OUT_MAP))
    parser.add_argument("--out-summary", default=str(OUT_SUMMARY))
    parser.add_argument("--splice-start", type=int, default=120)
    parser.add_argument("--splice-end", type=int, default=215)
    args = parser.parse_args()

    hmm_path = Path(args.hmm_aln)
    linker_path = Path(args.linker_aln)
    linker_map_path = Path(args.linker_map)
    for path in [hmm_path, linker_path, linker_map_path]:
        if not path.exists():
            raise SystemExit(f"ERROR: required input not found: {path}")

    hmm_entries, linker_entries, hmm_by_id, linker_by_id = load_entries(hmm_path, linker_path)
    _query_sid, query = find_query(hmm_entries)
    hmm_qpos_by_col = query_pos_by_col(query)
    hmm_left_cols = [
        col for col, qpos in hmm_qpos_by_col.items()
        if qpos < args.splice_start
    ]
    hmm_right_cols = [
        col for col, qpos in hmm_qpos_by_col.items()
        if qpos > args.splice_end
    ]

    linker_map = read_linker_map(linker_map_path)
    source_ranges = [
        parse_source_range(row.get("source_range", ""))
        for row in linker_map
        if row.get("source_range", "")
    ]
    window_start = min(start for start, _end in source_ranges)
    window_end = max(end for _start, end in source_ranges)
    full_window = args.splice_start <= window_start and args.splice_end >= window_end
    linker_cols = select_linker_columns(
        linker_map, args.splice_start, args.splice_end, full_window)

    out_entries = []
    for sid, linker_hdr, linker_seq in linker_entries:
        hmm_hdr, hmm_seq = hmm_by_id[sid]
        left = "".join(hmm_seq[col] for col in hmm_left_cols)
        middle = "".join(linker_seq[col] for col in linker_cols)
        right = "".join(hmm_seq[col] for col in hmm_right_cols)
        out_entries.append((sid, linker_hdr or hmm_hdr, left + middle + right))

    out_fa = Path(args.out_fa)
    out_map = Path(args.out_map)
    out_summary = Path(args.out_summary)
    write_fasta(out_entries, out_fa)

    map_rows = []
    for col in hmm_left_cols:
        qpos = hmm_qpos_by_col[col]
        map_rows.append({
            "source": "hmm_qc",
            "source_col": col,
            "source_block": "",
            "source_range": "",
            "qpos": qpos,
            "region": domain_of(qpos),
        })
    for col in linker_cols:
        row = linker_map[col]
        qpos = row.get("qpos", "")
        map_rows.append({
            "source": "linker_refined",
            "source_col": col,
            "source_block": row.get("block_label", ""),
            "source_range": row.get("source_range", ""),
            "qpos": qpos,
            "region": domain_of(qpos),
        })
    for col in hmm_right_cols:
        qpos = hmm_qpos_by_col[col]
        map_rows.append({
            "source": "hmm_qc",
            "source_col": col,
            "source_block": "",
            "source_range": "",
            "qpos": qpos,
            "region": domain_of(qpos),
        })

    out_map.parent.mkdir(parents=True, exist_ok=True)
    with open(out_map, "w", newline="") as f:
        fields = ["out_col", "source", "source_col", "source_block", "source_range", "qpos", "region"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for out_col, row in enumerate(map_rows):
            out = dict(row)
            out["out_col"] = out_col
            writer.writerow(out)

    query_out = next(seq for sid, _hdr, seq in out_entries if QUERY_TAG in sid)
    region_counts = {name: 0 for name in ["antenna", "linker", "catalytic", "insert"]}
    for row in map_rows:
        region = row["region"]
        region_counts[region] = region_counts.get(region, 0) + 1
    summary = {
        "n_sequences": len(out_entries),
        "alignment_width": len(out_entries[0][2]) if out_entries else 0,
        "query_residues": sum(aa != "-" for aa in query_out),
        "splice_start": args.splice_start,
        "splice_end": args.splice_end,
        "hmm_left_columns": len(hmm_left_cols),
        "linker_refined_columns": len(linker_cols),
        "hmm_right_columns": len(hmm_right_cols),
        "antenna_columns": region_counts.get("antenna", 0),
        "linker_columns": region_counts.get("linker", 0),
        "catalytic_columns": region_counts.get("catalytic", 0),
        "insert_columns": region_counts.get("insert", 0),
        "hmm_alignment": str(hmm_path),
        "linker_alignment": str(linker_path),
        "linker_column_map": str(linker_map_path),
    }
    with open(out_summary, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary)

    print(f"Saved: {out_fa} ({summary['n_sequences']} seqs x {summary['alignment_width']} cols)")
    print(f"Saved: {out_map}")
    print(f"Saved: {out_summary}")


if __name__ == "__main__":
    main()
