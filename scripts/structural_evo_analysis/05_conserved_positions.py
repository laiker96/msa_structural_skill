#!/usr/bin/env python3
"""Step 05: compute conserved alignment/query positions."""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, default=cfg.OUTPUT_DIR / "repset_aligned.fa")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "conservation")
    parser.add_argument("--min-occupancy", type=float, default=0.70)
    parser.add_argument("--min-top-frequency", type=float, default=0.80)
    parser.add_argument("--include-query-gap-columns", action="store_true")
    parser.add_argument("--no-plot", action="store_true", help="Do not write conservation_plot.png.")
    return parser.parse_args()


def column_stats(entries: list[tuple[str, str, str]], col: int, query_pos: str, query_aa: str) -> dict[str, object]:
    symbols = [seq[col].upper() for _sid, _hdr, seq in entries]
    residues = [aa for aa in symbols if aa in cfg.AA_SET]
    counts = Counter(residues)
    nseq = len(entries)
    nres = len(residues)
    occupancy = nres / nseq if nseq else 0.0
    top_aa, top_count = ("", 0)
    if counts:
        top_aa, top_count = counts.most_common(1)[0]
    top_freq_residues = top_count / nres if nres else 0.0
    top_freq_sequences = top_count / nseq if nseq else 0.0
    probs = [count / nres for count in counts.values()] if nres else []
    entropy = -sum(p * math.log2(p) for p in probs) if probs else 0.0
    norm_entropy = entropy / math.log2(len(cfg.AA_SET)) if probs else 0.0
    conservation = (1.0 - norm_entropy) * occupancy
    return {
        "alignment_col": col,
        "query_pos": query_pos,
        "query_aa": query_aa,
        "top_aa": top_aa,
        "n_sequences": nseq,
        "n_residues": nres,
        "occupancy": f"{occupancy:.4f}",
        "top_aa_count": top_count,
        "top_aa_frequency_among_residues": f"{top_freq_residues:.4f}",
        "top_aa_frequency_among_sequences": f"{top_freq_sequences:.4f}",
        "entropy_bits": f"{entropy:.4f}",
        "conservation": f"{conservation:.4f}",
    }


def write_conservation_plot(rows: list[dict[str, object]], out_path: Path) -> None:
    if not rows:
        return
    os.environ.setdefault("MPLCONFIGDIR", str(out_path.parent / ".matplotlib"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"WARNING: matplotlib not available; skipping conservation plot: {exc}", file=sys.stderr)
        return

    query_rows = [row for row in rows if str(row.get("query_pos", "")).isdigit()]
    plot_rows = query_rows if query_rows else rows
    x_key = "query_pos" if query_rows else "alignment_col"
    x_label = "Query position" if query_rows else "Alignment column"
    x_values = [int(row[x_key]) for row in plot_rows]
    conservation = [float(row["conservation"]) for row in plot_rows]
    occupancy = [float(row["occupancy"]) for row in plot_rows]
    conserved_x = [int(row[x_key]) for row in plot_rows if row.get("is_conserved") == "yes"]
    conserved_y = [float(row["conservation"]) for row in plot_rows if row.get("is_conserved") == "yes"]

    width = max(10.0, min(24.0, len(plot_rows) / 35.0))
    fig, ax = plt.subplots(figsize=(width, 4.8))
    ax.plot(x_values, conservation, color="#2f6f9f", linewidth=1.2, label="Conservation")
    ax.plot(x_values, occupancy, color="#8a8f3a", linewidth=0.9, alpha=0.75, label="Occupancy")
    if conserved_x:
        ax.scatter(conserved_x, conserved_y, s=10, color="#b23a48", label="Conserved calls", zorder=3)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Score")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("MSA conservation")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.alignment.exists():
        raise SystemExit(f"ERROR: alignment not found: {args.alignment}")
    entries = cfg.read_fasta(args.alignment)
    if not entries:
        raise SystemExit("ERROR: empty alignment")
    width = len(entries[0][2])
    if any(len(seq) != width for _sid, _hdr, seq in entries):
        raise SystemExit("ERROR: alignment contains unequal sequence lengths")
    query_seq = entries[0][2]
    query_pos = 0
    rows = []
    conserved = []
    for col, aa in enumerate(query_seq):
        if aa.upper() in cfg.AA_SET:
            query_pos += 1
            qpos_text = str(query_pos)
            query_aa = aa.upper()
        else:
            qpos_text = ""
            query_aa = ""
            if not args.include_query_gap_columns:
                continue
        row = column_stats(entries, col, qpos_text, query_aa)
        is_conserved = (
            float(row["occupancy"]) >= args.min_occupancy
            and float(row["top_aa_frequency_among_residues"]) >= args.min_top_frequency
        )
        row["is_conserved"] = "yes" if is_conserved else "no"
        rows.append(row)
        if is_conserved:
            conserved.append(row)
    fields = [
        "alignment_col", "query_pos", "query_aa", "top_aa", "n_sequences",
        "n_residues", "occupancy", "top_aa_count",
        "top_aa_frequency_among_residues", "top_aa_frequency_among_sequences",
        "entropy_bits", "conservation", "is_conserved",
    ]
    cfg.write_tsv(args.out_dir / "position_conservation.tsv", rows, fields)
    cfg.write_tsv(args.out_dir / "conserved_positions.tsv", conserved, fields)
    if not args.no_plot:
        plot_path = args.out_dir / "conservation_plot.png"
        write_conservation_plot(rows, plot_path)
        if plot_path.exists():
            print(f"Saved: {plot_path}")
    print(f"Saved: {args.out_dir / 'position_conservation.tsv'} ({len(rows)} columns)")
    print(f"Saved: {args.out_dir / 'conserved_positions.tsv'} ({len(conserved)} conserved)")


if __name__ == "__main__":
    main()
