#!/usr/bin/env python3
"""
Step 08: anchor-guided linker realignment.

The full-sequence HMM alignment is kept as the ANKros coordinate anchor. This
step extracts a linker-centered window, keeps conserved query-position anchor
islands fixed, realigns the variable intervals with MAFFT E-INS-i, and writes a
compact linker-focused MSA plus metrics.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
QUERY_TAG = "photoHymenobact"

INPUT_ALN = cfg.INTER_DIR / "repset_hmmalign_qc.fa"
META_TSV = cfg.INTER_DIR / "repset_metadata_qc.tsv"
OUT_DIR = cfg.INTER_DIR / "linker_refined"

DEFAULT_WINDOW = (120, 215)
DEFAULT_LINKER = (131, 205)
DEFAULT_ANCHORS = "auto"

MAFFT = cfg.resolve_bin("mafft")


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


def write_aligned_fasta(entries, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for _sid, hdr, seq in entries:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def normalise_entries(entries):
    return [(sid, hdr, seq.upper().replace(".", "-")) for sid, hdr, seq in entries]


def find_query(entries):
    for i, (sid, hdr, _seq) in enumerate(entries):
        if QUERY_TAG in sid or QUERY_TAG in hdr:
            return i
    sys.exit(f"ERROR: query tag {QUERY_TAG!r} not found in alignment")


def query_maps(query_aln):
    pos_to_col = {}
    col_to_pos = {}
    qpos = 0
    for col, aa in enumerate(query_aln):
        if aa != "-":
            qpos += 1
            pos_to_col[qpos] = col
            col_to_pos[col] = qpos
    return pos_to_col, col_to_pos


def raw_span_for_block(pos_to_col, lo, hi, window_end):
    raw_start = pos_to_col[lo]
    if hi < window_end:
        raw_end = pos_to_col[hi + 1] - 1
    else:
        raw_end = pos_to_col[hi]
    return raw_start, raw_end


def parse_range_list(text):
    ranges = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            left, right = item.split("-", 1)
            lo, hi = int(left), int(right)
        else:
            lo = hi = int(item)
        if lo > hi:
            lo, hi = hi, lo
        ranges.append((lo, hi))
    if not ranges:
        return []
    ranges.sort()
    merged = [ranges[0]]
    for lo, hi in ranges[1:]:
        prev_lo, prev_hi = merged[-1]
        if lo <= prev_hi + 1:
            merged[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged.append((lo, hi))
    return merged


def load_scope_ids(scope):
    if scope == "all":
        return None
    if not META_TSV.exists():
        sys.exit(f"ERROR: {META_TSV} not found; cannot use --scope {scope}")
    keep = set()
    with open(META_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            regime = row.get("regime", "")
            if scope == "annotated" and regime:
                keep.add(row["id"])
            elif scope == regime:
                keep.add(row["id"])
    return keep


def filter_entries(entries, scope):
    keep = load_scope_ids(scope)
    if keep is None:
        return entries
    query = entries[find_query(entries)]
    filtered = [entry for entry in entries if entry[0] in keep]
    if not any(entry[0] == query[0] for entry in filtered):
        filtered.insert(0, query)
    return filtered


def column_metrics(entries, pos_to_col, start, end):
    seqs = [seq for _sid, _hdr, seq in entries]
    rows = []
    for qpos in range(start, end + 1):
        col = pos_to_col.get(qpos)
        if col is None:
            continue
        residues = [
            seq[col]
            for seq in seqs
            if col < len(seq) and seq[col] in AA_ALPHABET
        ]
        occupancy = len(residues) / len(seqs) if seqs else 0.0
        if residues:
            top_aa, top_n = Counter(residues).most_common(1)[0]
            conservation = top_n / len(residues)
        else:
            top_aa, conservation = "-", 0.0
        rows.append({
            "qpos": qpos,
            "raw_col": col,
            "top_aa": top_aa,
            "occupancy": occupancy,
            "conservation": conservation,
            "score": occupancy * conservation,
            "n_residues": len(residues),
        })
    return rows


def discover_anchors(metrics, min_occ, min_cons, min_len, max_gap):
    selected = [
        row for row in metrics
        if row["occupancy"] >= min_occ and row["conservation"] >= min_cons
    ]
    if not selected:
        return []
    groups = []
    current = [selected[0]]
    for row in selected[1:]:
        if row["qpos"] <= current[-1]["qpos"] + max_gap + 1:
            current.append(row)
        else:
            groups.append(current)
            current = [row]
    groups.append(current)
    return [
        (group[0]["qpos"], group[-1]["qpos"])
        for group in groups
        if group[-1]["qpos"] - group[0]["qpos"] + 1 >= min_len
    ]


def clamp_anchors(anchors, start, end, pos_to_col):
    clamped = []
    for lo, hi in anchors:
        lo = max(lo, start)
        hi = min(hi, end)
        if lo > hi:
            continue
        missing = [pos for pos in range(lo, hi + 1) if pos not in pos_to_col]
        if missing:
            sys.exit(f"ERROR: anchor {lo}-{hi} contains absent query positions: {missing}")
        clamped.append((lo, hi))
    return parse_range_list(",".join(
        f"{lo}-{hi}" if lo != hi else str(lo)
        for lo, hi in clamped
    ))


def ungap(seq):
    return "".join(aa for aa in seq if aa in AA_ALPHABET)


def extract_raw_interval(seq, pos_to_col, start, end):
    if start > end:
        return ""
    return seq[pos_to_col[start]:pos_to_col[end] + 1]


def write_temp_fasta(entries, path):
    with open(path, "w") as f:
        for sid, _hdr, seq in entries:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def read_aligned_by_sid(path):
    return {sid: seq for sid, _hdr, seq in normalise_entries(read_fasta(path))}


def run_mafft_segment(segment_entries, workdir, label):
    nonempty = [(sid, hdr, seq) for sid, hdr, seq in segment_entries if seq]
    if not nonempty:
        return {sid: "" for sid, _hdr, _seq in segment_entries}, 0
    if len(nonempty) == 1:
        width = len(nonempty[0][2])
        return {
            sid: seq if seq else "-" * width
            for sid, _hdr, seq in segment_entries
        }, width

    in_fa = workdir / f"{label}.in.fa"
    out_fa = workdir / f"{label}.mafft.fa"
    write_temp_fasta(nonempty, in_fa)
    cmd = [
        MAFFT, "--genafpair", "--maxiterate", "1000",
        "--thread", str(cfg.N_THREADS), "--quiet", str(in_fa),
    ]
    with open(out_fa, "w") as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(res.stderr[-1500:])
        sys.exit(f"ERROR: MAFFT failed for {label} (exit {res.returncode})")
    aligned_nonempty = read_aligned_by_sid(out_fa)
    width = len(next(iter(aligned_nonempty.values()))) if aligned_nonempty else 0
    return {
        sid: aligned_nonempty.get(sid, "-" * width)
        for sid, _hdr, _seq in segment_entries
    }, width


def percentile(sorted_values, frac):
    if not sorted_values:
        return 0
    idx = round((len(sorted_values) - 1) * frac)
    return sorted_values[idx]


def split_core_by_length(segment_entries, query_sid, iqr_factor, min_core):
    lengths = {sid: len(seq) for sid, _hdr, seq in segment_entries}
    nonzero = sorted(length for length in lengths.values() if length > 0)
    if len(nonzero) < max(min_core, 4):
        return list(segment_entries), [], length_summary(nonzero)

    q1 = percentile(nonzero, 0.25)
    q3 = percentile(nonzero, 0.75)
    iqr = max(1, q3 - q1)
    low = q1 - iqr_factor * iqr
    high = q3 + iqr_factor * iqr

    core, outliers = [], []
    for entry in segment_entries:
        sid = entry[0]
        length = lengths[sid]
        if sid == query_sid or (length > 0 and low <= length <= high):
            core.append(entry)
        else:
            outliers.append(entry)

    if len([entry for entry in core if entry[2]]) < min_core:
        return list(segment_entries), [], length_summary(nonzero)
    summary = length_summary(sorted(lengths[entry[0]] for entry in core if lengths[entry[0]] > 0))
    summary.update({"outlier_low": low, "outlier_high": high})
    return core, outliers, summary


def length_summary(lengths):
    return {
        "core_len_min": min(lengths) if lengths else 0,
        "core_len_median": statistics.median(lengths) if lengths else 0,
        "core_len_max": max(lengths) if lengths else 0,
        "outlier_low": "",
        "outlier_high": "",
    }


def add_outliers_keeplength(core_aligned, outlier_entries, workdir, label):
    width = len(next(iter(core_aligned.values()))) if core_aligned else 0
    if not outlier_entries or width == 0:
        return dict(core_aligned), width

    nonempty_outliers = [(sid, hdr, seq) for sid, hdr, seq in outlier_entries if seq]
    out = dict(core_aligned)
    if not nonempty_outliers:
        out.update({sid: "-" * width for sid, _hdr, _seq in outlier_entries})
        return out, width

    core_fa = workdir / f"{label}.core_aln.fa"
    add_fa = workdir / f"{label}.outliers.fa"
    out_fa = workdir / f"{label}.projected.fa"
    write_temp_fasta([(sid, sid, seq) for sid, seq in core_aligned.items()], core_fa)
    write_temp_fasta(nonempty_outliers, add_fa)
    cmd = [
        MAFFT, "--add", str(add_fa), "--keeplength",
        "--thread", str(cfg.N_THREADS), "--quiet", str(core_fa),
    ]
    with open(out_fa, "w") as handle:
        res = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(res.stderr[-1500:])
        sys.exit(f"ERROR: MAFFT --add --keeplength failed for {label} (exit {res.returncode})")
    projected = read_aligned_by_sid(out_fa)
    projected.update({sid: "-" * width for sid, _hdr, seq in outlier_entries if not seq})
    return projected, width


def run_mafft_segment_robust(segment_entries, workdir, label, query_sid, iqr_factor, min_core):
    core_entries, outlier_entries, summary = split_core_by_length(
        segment_entries, query_sid, iqr_factor, min_core)
    core_aligned, width = run_mafft_segment(core_entries, workdir, f"{label}.core")
    aligned, width = add_outliers_keeplength(core_aligned, outlier_entries, workdir, label)
    for sid, _hdr, _seq in segment_entries:
        aligned.setdefault(sid, "-" * width)
    summary.update({
        "n_core": len(core_entries),
        "n_outliers": len(outlier_entries),
        "robust": 1,
    })
    return aligned, width, summary


def fixed_anchor_block(entries, pos_to_col, col_to_pos, lo, hi, window_end):
    raw_start, raw_end = raw_span_for_block(pos_to_col, lo, hi, window_end)
    cols = list(range(raw_start, raw_end + 1))
    block = {
        sid: "".join(seq[col] if col < len(seq) else "-" for col in cols)
        for sid, _hdr, seq in entries
    }
    colmap = []
    for raw_col in cols:
        colmap.append({
            "block_type": "anchor",
            "block_label": f"anchor_{lo}_{hi}",
            "qpos": col_to_pos.get(raw_col, ""),
            "raw_col": raw_col,
            "source_range": f"{lo}-{hi}",
        })
    return block, colmap


def variable_block(entries, pos_to_col, lo, hi, workdir, label, args):
    raw_start, raw_end = raw_span_for_block(pos_to_col, lo, hi, args.window_end)
    raw_segments = []
    lengths = []
    for sid, hdr, seq in entries:
        segment = ungap(seq[raw_start:raw_end + 1])
        raw_segments.append((sid, hdr, segment))
        lengths.append(len(segment))

    query_sid = entries[find_query(entries)][0]
    if args.robust:
        aligned, width, robust_summary = run_mafft_segment_robust(
            raw_segments, workdir, label, query_sid, args.robust_iqr, args.robust_min_core)
    else:
        aligned, width = run_mafft_segment(raw_segments, workdir, label)
        robust_summary = {
            "n_core": len(raw_segments),
            "n_outliers": 0,
            "core_len_min": "",
            "core_len_median": "",
            "core_len_max": "",
            "outlier_low": "",
            "outlier_high": "",
            "robust": 0,
        }

    query_aln = aligned[query_sid]
    qpos = lo - 1
    colmap = []
    for aa in query_aln:
        if aa == "-":
            mapped_qpos = ""
            raw_col = ""
        else:
            qpos += 1
            mapped_qpos = qpos
            raw_col = pos_to_col.get(qpos, "")
        colmap.append({
            "block_type": "variable",
            "block_label": label,
            "qpos": mapped_qpos,
            "raw_col": raw_col,
            "source_range": f"{lo}-{hi}",
        })

    stats = {
        "block_label": label,
        "block_type": "variable",
        "qpos_start": lo,
        "qpos_end": hi,
        "n_sequences": len(entries),
        "n_nonempty": sum(1 for length in lengths if length > 0),
        "raw_len_min": min(lengths) if lengths else 0,
        "raw_len_median": statistics.median(lengths) if lengths else 0,
        "raw_len_max": max(lengths) if lengths else 0,
        "aligned_width": width,
    }
    stats.update(robust_summary)
    return aligned, colmap, stats


def build_blocks(window_start, window_end, anchors):
    blocks = []
    cursor = window_start
    for lo, hi in anchors:
        if cursor < lo:
            blocks.append(("variable", cursor, lo - 1))
        blocks.append(("anchor", lo, hi))
        cursor = hi + 1
    if cursor <= window_end:
        blocks.append(("variable", cursor, window_end))
    return blocks


def alignment_summary(entries, seq_parts, linker_range, window_start):
    start, end = linker_range
    full = [(sid, hdr, "".join(seq_parts[sid])) for sid, hdr, _seq in entries]
    query_sid = entries[find_query(entries)][0]
    query = next(seq for sid, _hdr, seq in full if sid == query_sid)
    qpos_by_col = {}
    qpos = window_start - 1
    for col, aa in enumerate(query):
        if aa != "-":
            qpos += 1
            qpos_by_col[col] = qpos
    linker_cols = [col for col, pos in qpos_by_col.items() if start <= pos <= end]
    if not linker_cols:
        return {}
    occupancies = []
    for col in linker_cols:
        occupancies.append(sum(seq[col] != "-" for _sid, _hdr, seq in full) / len(full))
    return {
        "n_sequences": len(full),
        "alignment_width": len(query),
        "linker_query_columns": len(linker_cols),
        "linker_column_occupancy_min": min(occupancies),
        "linker_column_occupancy_median": statistics.median(occupancies),
        "linker_column_occupancy_mean": statistics.mean(occupancies),
    }


def write_sequence_metrics(entries, seq_parts, pos_to_col, colmap_rows, args, out_path):
    linker_cols = {
        idx for idx, row in enumerate(colmap_rows)
        if str(row.get("qpos", "")).isdigit()
        and DEFAULT_LINKER[0] <= int(row["qpos"]) <= DEFAULT_LINKER[1]
    }
    rows = []
    for sid, hdr, input_seq in entries:
        output_seq = "".join(seq_parts[sid])
        input_window = ungap(extract_raw_interval(
            input_seq, pos_to_col, args.window_start, args.window_end))
        input_linker = ungap(extract_raw_interval(
            input_seq, pos_to_col, DEFAULT_LINKER[0], DEFAULT_LINKER[1]))
        output_residues = sum(aa in AA_ALPHABET for aa in output_seq)
        output_query_linker = sum(
            output_seq[col] in AA_ALPHABET
            for col in linker_cols
            if col < len(output_seq)
        )
        rows.append({
            "id": sid,
            "header": hdr,
            "input_window_residues": len(input_window),
            "output_window_residues": output_residues,
            "window_residue_delta_input_minus_output": len(input_window) - output_residues,
            "window_residue_retention": output_residues / len(input_window) if input_window else 0,
            "input_linker_residues": len(input_linker),
            "output_query_linker_residues": output_query_linker,
        })

    fields = [
        "id", "header", "input_window_residues", "output_window_residues",
        "window_residue_delta_input_minus_output", "window_residue_retention",
        "input_linker_residues", "output_query_linker_residues",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["window_residue_retention"] = f"{out['window_residue_retention']:.4f}"
            writer.writerow(out)


def write_outputs(entries, seq_parts, colmap_rows, block_stats, anchor_metrics, anchors, args, out_dir, pos_to_col):
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = out_dir / "linker_refined.fa"
    out_colmap = out_dir / "linker_refined_column_map.tsv"
    out_blocks = out_dir / "linker_refined_block_stats.tsv"
    out_anchors = out_dir / "linker_refined_anchor_metrics.tsv"
    out_summary = out_dir / "linker_refined_summary.tsv"
    out_sequence_metrics = out_dir / "linker_refined_sequence_metrics.tsv"

    out_entries = [
        (sid, hdr, "".join(seq_parts[sid]))
        for sid, hdr, _seq in entries
    ]
    write_aligned_fasta(out_entries, out_fasta)

    with open(out_colmap, "w", newline="") as f:
        cols = ["out_col", "block_type", "block_label", "qpos", "raw_col", "source_range"]
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, row in enumerate(colmap_rows):
            out = dict(row)
            out["out_col"] = i
            writer.writerow(out)

    block_cols = [
        "block_label", "block_type", "qpos_start", "qpos_end",
        "n_sequences", "n_nonempty", "raw_len_min", "raw_len_median",
        "raw_len_max", "aligned_width", "robust", "n_core", "n_outliers",
        "core_len_min", "core_len_median", "core_len_max",
        "outlier_low", "outlier_high",
    ]
    with open(out_blocks, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=block_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in block_stats:
            writer.writerow({col: row.get(col, "") for col in block_cols})

    anchor_pos = {pos for lo, hi in anchors for pos in range(lo, hi + 1)}
    anchor_cols = [
        "qpos", "raw_col", "top_aa", "occupancy", "conservation",
        "score", "n_residues", "is_anchor",
    ]
    with open(out_anchors, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=anchor_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in anchor_metrics:
            out = dict(row)
            for key in ["occupancy", "conservation", "score"]:
                out[key] = f"{out[key]:.4f}"
            out["is_anchor"] = "1" if row["qpos"] in anchor_pos else "0"
            writer.writerow(out)

    write_sequence_metrics(entries, seq_parts, pos_to_col, colmap_rows, args, out_sequence_metrics)

    summary = alignment_summary(entries, seq_parts, DEFAULT_LINKER, args.window_start)
    summary.update({
        "input": str(args.input),
        "scope": args.scope,
        "window_start": args.window_start,
        "window_end": args.window_end,
        "linker_start": DEFAULT_LINKER[0],
        "linker_end": DEFAULT_LINKER[1],
        "anchors": ",".join(f"{lo}-{hi}" if lo != hi else str(lo) for lo, hi in anchors),
        "robust": int(args.robust),
    })
    with open(out_summary, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            key: f"{value:.4f}" if isinstance(value, float) else value
            for key, value in summary.items()
        })

    manifest = out_dir / "README.md"
    manifest.write_text(
        "# Linker Refined Alignment\n\n"
        f"Input: `{args.input}`\n\n"
        f"Window: ANKros {args.window_start}-{args.window_end}\n\n"
        f"Linker: ANKros {DEFAULT_LINKER[0]}-{DEFAULT_LINKER[1]}\n\n"
        f"Scope: `{args.scope}`\n\n"
        "Anchors: "
        + ", ".join(f"{lo}-{hi}" if lo != hi else str(lo) for lo, hi in anchors)
        + "\n\nVariable intervals were ungapped per sequence and realigned with "
        "`mafft --genafpair --maxiterate 1000`.\n\n"
        f"Robust mode: `{args.robust}`"
        + (f" (IQR factor {args.robust_iqr}, min core {args.robust_min_core}).\n"
           if args.robust else ".\n")
    )
    return out_fasta, out_colmap, out_blocks, out_anchors, out_summary, out_sequence_metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(INPUT_ALN),
                        help="QC-passed full HMM alignment.")
    parser.add_argument("--scope", default="all",
                        choices=["all", "annotated", "psychro", "meso", "thermo"],
                        help="Subset to refine. The ANKros query is always retained.")
    parser.add_argument("--window-start", type=int, default=DEFAULT_WINDOW[0])
    parser.add_argument("--window-end", type=int, default=DEFAULT_WINDOW[1])
    parser.add_argument("--anchors", default=DEFAULT_ANCHORS,
                        help='Comma-separated ANKros anchor ranges, or "auto".')
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--robust", dest="robust", action="store_true", default=False,
                        help="Align central length core, then add outliers with MAFFT --add --keeplength.")
    parser.add_argument("--no-robust", dest="robust", action="store_false",
                        help="Align each full variable block in one MAFFT run. This is the default.")
    parser.add_argument("--robust-iqr", type=float, default=0.5)
    parser.add_argument("--robust-min-core", type=int, default=30)
    parser.add_argument("--auto-min-occ", type=float, default=0.70)
    parser.add_argument("--auto-min-cons", type=float, default=0.35)
    parser.add_argument("--auto-min-len", type=int, default=2)
    parser.add_argument("--auto-max-gap", type=int, default=2)
    parser.add_argument("--keep-temp", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out_dir) if args.out_dir else (
        OUT_DIR if args.scope == "all" else cfg.INTER_DIR / f"linker_refined_{args.scope}"
    )
    if not input_path.exists():
        sys.exit(f"ERROR: input alignment not found: {input_path}")
    if args.window_start > args.window_end:
        sys.exit("ERROR: --window-start must be <= --window-end")

    print(f"Loading {input_path}...")
    entries = normalise_entries(read_fasta(input_path))
    entries = filter_entries(entries, args.scope)
    query_idx = find_query(entries)
    query_aln = entries[query_idx][2]
    pos_to_col, _col_to_pos = query_maps(query_aln)
    missing = [
        pos for pos in range(args.window_start, args.window_end + 1)
        if pos not in pos_to_col
    ]
    if missing:
        sys.exit(f"ERROR: query positions absent from input alignment: {missing[:20]}")
    print(f"  {len(entries)} sequences x {len(query_aln)} input columns")

    metrics = column_metrics(entries, pos_to_col, args.window_start, args.window_end)
    if args.anchors == "auto":
        anchors = discover_anchors(
            metrics, args.auto_min_occ, args.auto_min_cons,
            args.auto_min_len, args.auto_max_gap)
        if not anchors:
            sys.exit("ERROR: --anchors auto found no conserved islands")
    else:
        anchors = parse_range_list(args.anchors)
    anchors = clamp_anchors(anchors, args.window_start, args.window_end, pos_to_col)
    print("  Anchors: " + ", ".join(f"{lo}-{hi}" if lo != hi else str(lo) for lo, hi in anchors))

    blocks = build_blocks(args.window_start, args.window_end, anchors)
    print("  Blocks:")
    for typ, lo, hi in blocks:
        print(f"    {typ:8s} {lo}-{hi}")

    seq_parts = {sid: [] for sid, _hdr, _seq in entries}
    colmap_rows = []
    block_stats = []

    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tmp_linker_", dir=out_dir) as tmp:
        workdir = Path(tmp)
        for idx, (typ, lo, hi) in enumerate(blocks, start=1):
            label = f"{idx:02d}_{typ}_{lo}_{hi}"
            if typ == "anchor":
                block, cmap = fixed_anchor_block(entries, pos_to_col, _col_to_pos, lo, hi, args.window_end)
                for sid, part in block.items():
                    seq_parts[sid].append(part)
                colmap_rows.extend(cmap)
                block_stats.append({
                    "block_label": label,
                    "block_type": "anchor",
                    "qpos_start": lo,
                    "qpos_end": hi,
                    "n_sequences": len(entries),
                    "n_nonempty": len(entries),
                    "raw_len_min": len(cmap),
                    "raw_len_median": len(cmap),
                    "raw_len_max": len(cmap),
                    "aligned_width": len(cmap),
                    "robust": int(args.robust),
                    "n_core": len(entries),
                    "n_outliers": 0,
                })
            else:
                aligned, cmap, stats = variable_block(entries, pos_to_col, lo, hi, workdir, label, args)
                for sid, part in aligned.items():
                    seq_parts[sid].append(part)
                colmap_rows.extend(cmap)
                block_stats.append(stats)
        if args.keep_temp:
            saved = out_dir / "tmp_last"
            if saved.exists():
                shutil.rmtree(saved)
            shutil.copytree(workdir, saved)
            print(f"  Temp files saved to {saved}")

    outputs = write_outputs(
        entries, seq_parts, colmap_rows, block_stats, metrics, anchors, args, out_dir, pos_to_col)
    width = len("".join(next(iter(seq_parts.values()))))
    print(f"\nSaved: {outputs[0]} ({len(entries)} seqs x {width} cols)")
    for path in outputs[1:]:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
