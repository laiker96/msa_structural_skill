#!/usr/bin/env python3
"""
Step 09: domain sequence logos for the refined alignment.

The logos are gap-weighted amino-acid information logos: each column's stack
height is information content multiplied by non-gap occupancy. A small
occupancy track is drawn below each logo row so indel-rich linker columns are
not visually over-interpreted.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)
MAX_BITS = math.log2(len(AA_ORDER))

ALN_FA = cfg.INTER_DIR / "repset_hmmalign_linker_refined.fa"
COL_MAP = cfg.INTER_DIR / "repset_hmmalign_linker_refined_column_map.tsv"
OUT_DIR = cfg.INTER_DIR / "figures" / "domain_logos"

DOMAINS = {
    "antenna": (1, 130),
    "linker": (131, 205),
    "catalytic": (206, 437),
}

AA_COLORS = {}
for aas, color in [
    ("AILMV", "#1976D2"),
    ("KRH", "#D32F2F"),
    ("DE", "#C2185B"),
    ("STNQ", "#388E3C"),
    ("FWY", "#F57C00"),
    ("P", "#FBC02D"),
    ("G", "#6D4C41"),
    ("C", "#7B1FA2"),
]:
    for aa in aas:
        AA_COLORS[aa] = color

FONT = FontProperties(family="DejaVu Sans", weight="bold")
STACK_GAP = 0.035
ANKROS_HIGHLIGHT = "#D81B60"
LETTER_WIDTHS = {
    aa: TextPath((0, 0), aa, size=1, prop=FONT).get_extents().width
    for aa in AA_ORDER
}
MAX_LETTER_WIDTH = max(LETTER_WIDTHS.values())


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


def read_column_map(path: Path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def find_query(entries):
    for sid, hdr, seq in entries:
        if "photoHymenobact" in sid or "photoHymenobact" in hdr:
            return seq
    raise SystemExit("ERROR: ANKros query sequence not found in alignment")


def parse_range(text):
    text = (text or "").strip()
    if not text:
        return None
    if "-" in text:
        left, right = text.split("-", 1)
        return int(left), int(right)
    value = int(text)
    return value, value


def contained_in_domain(source_range, start, end):
    parsed = parse_range(source_range)
    return bool(parsed and start <= parsed[0] and parsed[1] <= end)


def select_domain_columns(colmap, domain, include_contained_inserts=True):
    start, end = DOMAINS[domain]
    cols = []
    for row in colmap:
        out_col = int(row["out_col"])
        if row["region"] == domain:
            cols.append((out_col, row, False))
        elif (
            include_contained_inserts
            and row["region"] == "insert"
            and contained_in_domain(row.get("source_range", ""), start, end)
        ):
            cols.append((out_col, row, True))
    return cols


def column_stats(entries, selected_cols, query_seq):
    rows = []
    nseq = len(entries)
    for domain_col, (aln_col, row, is_insert) in enumerate(selected_cols, start=1):
        residues = [
            seq[aln_col].upper()
            for _sid, _hdr, seq in entries
            if aln_col < len(seq) and seq[aln_col].upper() in AA_SET
        ]
        counts = Counter(residues)
        n_res = sum(counts.values())
        occupancy = n_res / nseq if nseq else 0.0
        probs = {aa: count / n_res for aa, count in counts.items()} if n_res else {}
        entropy = -sum(p * math.log2(p) for p in probs.values()) if probs else 0.0
        information = (MAX_BITS - entropy) * occupancy if probs else 0.0
        heights = {aa: probs.get(aa, 0.0) * information for aa in AA_ORDER}
        top_aa, top_count = ("-", 0)
        if counts:
            top_aa, top_count = counts.most_common(1)[0]
        query_aa = ""
        if aln_col < len(query_seq) and query_seq[aln_col].upper() in AA_SET:
            query_aa = query_seq[aln_col].upper()
        rows.append({
            "domain_col": domain_col,
            "alignment_col": aln_col,
            "qpos": row.get("qpos", ""),
            "region": row.get("region", ""),
            "source": row.get("source", ""),
            "source_block": row.get("source_block", ""),
            "source_range": row.get("source_range", ""),
            "is_insert": "yes" if is_insert or row.get("region") == "insert" else "no",
            "occupancy": occupancy,
            "gap_fraction": 1.0 - occupancy,
            "information_bits_gap_weighted": information,
            "top_aa": top_aa,
            "top_aa_frequency_among_residues": top_count / n_res if n_res else 0.0,
            "n_residues": n_res,
            "query_aa": query_aa,
            "heights": heights,
        })
    return rows


def draw_letter(ax, aa, x, y, height):
    if height <= 0.01:
        return
    text_path = TextPath((0, 0), aa, size=1, prop=FONT)
    bbox = text_path.get_extents()
    if bbox.width == 0 or bbox.height == 0:
        return
    natural_width = LETTER_WIDTHS.get(aa, bbox.width) / MAX_LETTER_WIDTH
    target_width = max(0.18, 0.82 * natural_width)
    sx = target_width / bbox.width
    sy = height / bbox.height
    trans = (
        Affine2D()
        .scale(sx, sy)
        .translate(x + 0.5 - (bbox.x0 + bbox.width / 2) * sx, y - bbox.y0 * sy)
        + ax.transData
    )
    patch = PathPatch(
        text_path,
        transform=trans,
        facecolor=AA_COLORS.get(aa, "#555555"),
        edgecolor="none",
        clip_on=True,
    )
    ax.add_patch(patch)


def label_for_column(row):
    if row["qpos"]:
        return row["qpos"]
    if row["is_insert"] == "yes":
        return "+"
    return ""


def draw_logo(domain, rows, out_png, wrap):
    if not rows:
        return
    chunks = [rows[i:i + wrap] for i in range(0, len(rows), wrap)]
    fig_w = max(12, min(26, wrap * 0.23))
    fig_h = max(3.3, len(chunks) * 2.65)
    fig, axes = plt.subplots(len(chunks), 1, figsize=(fig_w, fig_h), squeeze=False)
    axes = axes.ravel()

    for row_idx, (ax, chunk) in enumerate(zip(axes, chunks), start=1):
        ax.set_xlim(0, len(chunk))
        ax.set_ylim(-0.85, MAX_BITS + 0.15)
        ax.axhline(0, color="#333333", linewidth=0.7)
        ax.set_ylabel("bits")
        ax.set_yticks([0, 1, 2, 3, 4])

        for x, stat in enumerate(chunk):
            occ_h = 0.30 * stat["occupancy"]
            ax.add_patch(Rectangle(
                (x + 0.08, -0.42), 0.84, occ_h,
                facecolor="#9E9E9E", edgecolor="none", alpha=0.75,
            ))
            if stat["is_insert"] == "yes":
                ax.add_patch(Rectangle(
                    (x, -0.82), 1.0, MAX_BITS + 0.97,
                    facecolor="#F5F5F5", edgecolor="none", zorder=-1,
                ))
            if stat["query_aa"]:
                ax.add_patch(Rectangle(
                    (x + 0.14, -0.80), 0.72, 0.26,
                    facecolor="white", edgecolor=ANKROS_HIGHLIGHT,
                    linewidth=0.75, alpha=0.95, zorder=2,
                ))
                ax.text(
                    x + 0.5, -0.67, stat["query_aa"],
                    ha="center", va="center", fontsize=6.5,
                    color=AA_COLORS.get(stat["query_aa"], "#333333"),
                    fontweight="bold", zorder=3,
                )
            elif stat["is_insert"] == "yes":
                ax.text(
                    x + 0.5, -0.67, "+",
                    ha="center", va="center", fontsize=5.5,
                    color="#9E9E9E",
                )

            y = 0.0
            for aa, height in sorted(stat["heights"].items(), key=lambda item: item[1]):
                draw_letter(ax, aa, x, y, height)
                y += height + STACK_GAP

        tick_positions = []
        tick_labels = []
        for x, stat in enumerate(chunk):
            label = label_for_column(stat)
            if not label:
                continue
            if label == "+":
                if x % 10 == 0:
                    tick_positions.append(x + 0.5)
                    tick_labels.append("+")
                continue
            qpos = int(label)
            if qpos % 10 == 0 or x == 0 or x == len(chunk) - 1:
                tick_positions.append(x + 0.5)
                tick_labels.append(str(qpos))
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8, rotation=90)
        ax.text(
            -0.6, -0.67, "ANKros",
            ha="right", va="center", fontsize=8,
            color=ANKROS_HIGHLIGHT, fontweight="bold", clip_on=False,
        )
        ax.set_xlabel("ANKros residue position (+ = contained insert column)")
        ax.set_title(f"{domain.capitalize()} logo, row {row_idx}/{len(chunks)}", loc="left")

    fig.suptitle(
        f"{domain.capitalize()} domain sequence logo "
        "(gap-weighted amino-acid information; grey bars = occupancy)",
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def write_stats(domain, rows, out_tsv):
    fields = [
        "domain", "domain_col", "alignment_col", "qpos", "region", "source",
        "source_block", "source_range", "is_insert", "occupancy",
        "gap_fraction", "information_bits_gap_weighted", "top_aa",
        "top_aa_frequency_among_residues", "n_residues", "query_aa",
    ]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {field: row.get(field, "") for field in fields}
            out["domain"] = domain
            for key in ["occupancy", "gap_fraction", "information_bits_gap_weighted",
                        "top_aa_frequency_among_residues"]:
                out[key] = f"{out[key]:.4f}"
            writer.writerow(out)


def write_summary(summary_rows, out_path):
    fields = [
        "domain", "n_columns", "n_query_columns", "n_insert_columns",
        "median_occupancy", "mean_occupancy", "median_information_bits_gap_weighted",
        "mean_information_bits_gap_weighted", "png", "tsv",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", default=str(ALN_FA))
    parser.add_argument("--column-map", default=str(COL_MAP))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--wrap", type=int, default=80)
    parser.add_argument("--include-contained-inserts", action="store_true",
                        help="Also plot insert columns fully contained inside each domain.")
    args = parser.parse_args()

    alignment = Path(args.alignment)
    colmap_path = Path(args.column_map)
    out_dir = Path(args.out_dir)
    if not alignment.exists():
        raise SystemExit(f"ERROR: alignment not found: {alignment}")
    if not colmap_path.exists():
        raise SystemExit(f"ERROR: column map not found: {colmap_path}")

    entries = read_fasta(alignment)
    colmap = read_column_map(colmap_path)
    if not entries:
        raise SystemExit("ERROR: empty alignment")
    if len(entries[0][2]) != len(colmap):
        raise SystemExit(
            f"ERROR: alignment width {len(entries[0][2])} != column map rows {len(colmap)}")
    query_seq = find_query(entries)

    summary_rows = []
    for domain in ["antenna", "linker", "catalytic"]:
        selected = select_domain_columns(
            colmap, domain, include_contained_inserts=args.include_contained_inserts)
        rows = column_stats(entries, selected, query_seq)
        png = out_dir / f"{domain}_logo.png"
        tsv = out_dir / f"{domain}_logo_columns.tsv"
        draw_logo(domain, rows, png, args.wrap)
        write_stats(domain, rows, tsv)
        occupancies = [row["occupancy"] for row in rows]
        infos = [row["information_bits_gap_weighted"] for row in rows]
        n_insert = sum(1 for row in rows if row["is_insert"] == "yes")
        summary_rows.append({
            "domain": domain,
            "n_columns": len(rows),
            "n_query_columns": len(rows) - n_insert,
            "n_insert_columns": n_insert,
            "median_occupancy": f"{median(occupancies):.4f}",
            "mean_occupancy": f"{mean(occupancies):.4f}",
            "median_information_bits_gap_weighted": f"{median(infos):.4f}",
            "mean_information_bits_gap_weighted": f"{mean(infos):.4f}",
            "png": str(png),
            "tsv": str(tsv),
        })
        print(f"Saved: {png}")
        print(f"Saved: {tsv}")

    summary_path = out_dir / "domain_logo_summary.tsv"
    write_summary(summary_rows, summary_path)
    print(f"Saved: {summary_path}")


def mean(values):
    return sum(values) / len(values) if values else 0.0


def median(values):
    if not values:
        return 0.0
    vals = sorted(values)
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


if __name__ == "__main__":
    main()
