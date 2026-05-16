#!/usr/bin/env python3
"""
Step 10: ANKros-position conservation profile for the refined alignment.

Plots a gap-weighted, entropy-based MSA conservation score on the ANKros query
coordinate frame (1-437). The x-axis is ANKros position; the y-value is derived
from the entire alignment column and stays in [0, 1].
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

ALN_FA = cfg.INTER_DIR / "repset_hmmalign_linker_refined.fa"
COL_MAP = cfg.INTER_DIR / "repset_hmmalign_linker_refined_column_map.tsv"
OUT_DIR = cfg.INTER_DIR / "figures" / "conservation_profile"
AA_COUNT = 20

DOMAINS = {
    "antenna": (1, 130, "#4fc3f7"),
    "linker": (131, 205, "#90a4ae"),
    "catalytic": (206, 437, "#ce93d8"),
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


def read_column_map(path: Path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def find_query(entries):
    for sid, hdr, seq in entries:
        if "photoHymenobact" in sid or "photoHymenobact" in hdr:
            return seq
    raise SystemExit("ERROR: ANKros query sequence not found in alignment")

def compute_profile(entries, colmap):
    nseq = len(entries)
    rows = []
    for row in colmap:
        qpos = row.get("qpos", "").strip()
        if not qpos:
            continue
        aln_col = int(row["out_col"])
        residues = [
            seq[aln_col].upper()
            for _sid, _hdr, seq in entries
            if aln_col < len(seq) and seq[aln_col] not in "-."
        ]
        counts = Counter(residues)
        non_gap = sum(counts.values())
        occupancy = non_gap / nseq if nseq else 0.0
        if counts:
            top_aa, top_n = counts.most_common(1)[0]
            probs = [count / non_gap for count in counts.values()]
            entropy = -sum(p * math.log2(p) for p in probs)
            normalized_entropy = entropy / math.log2(AA_COUNT) if non_gap else 0.0
            conservation = (1.0 - normalized_entropy) * occupancy
        else:
            top_aa, conservation = "-", 0.0
        rows.append({
            "qpos": int(qpos),
            "alignment_col": aln_col,
            "region": row.get("region", ""),
            "source": row.get("source", ""),
            "source_block": row.get("source_block", ""),
            "source_range": row.get("source_range", ""),
            "occupancy": occupancy,
            "conservation": conservation,
            "top_aa": top_aa,
            "n_residues": non_gap,
        })
    return rows


def moving_average(values, window):
    if window <= 1 or not values:
        return list(values)
    half = window // 2
    smoothed = []
    for i in range(len(values)):
        lo = max(0, i - half)
        hi = min(len(values), i + half + 1)
        smoothed.append(sum(values[lo:hi]) / (hi - lo))
    return smoothed


def plot_profile(rows, out_png, out_tsv, smooth_window):
    rows = sorted(rows, key=lambda r: r["qpos"])
    qpos = [row["qpos"] for row in rows]
    conservation = [row["conservation"] for row in rows]
    occupancy = [row["occupancy"] for row in rows]
    smooth = moving_average(conservation, smooth_window)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        fields = ["qpos", "alignment_col", "region", "source", "source_block",
                  "source_range", "occupancy", "conservation", "smoothed_conservation",
                  "top_aa", "n_residues"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row, sm in zip(rows, smooth):
            out = dict(row)
            out["occupancy"] = f"{out['occupancy']:.4f}"
            out["conservation"] = f"{out['conservation']:.4f}"
            out["smoothed_conservation"] = f"{sm:.4f}"
            writer.writerow({field: out.get(field, "") for field in fields})

    fig, ax = plt.subplots(figsize=(16, 4.8))
    ax.plot(qpos, conservation, color="#546E7A", linewidth=1.0, alpha=0.4, label="raw")
    ax.plot(qpos, smooth, color="#1E88E5", linewidth=2.0, label=f"{smooth_window}-res moving avg")
    ax.fill_between(qpos, 0, occupancy, color="#B0BEC5", alpha=0.25, label="occupancy")

    for name, (start, end, color) in DOMAINS.items():
        ax.axvspan(start, end, color=color, alpha=0.12, lw=0)

    ax.set_xlim(1, 437)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("ANKros residue position")
    ax.set_ylabel("MSA conservation")
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6)
    ax.legend(loc="lower left", frameon=False, ncol=3)
    ax.set_title("Refined MSA conservation across ANKros coordinates")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", default=str(ALN_FA))
    parser.add_argument("--column-map", default=str(COL_MAP))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--smooth-window", type=int, default=9)
    args = parser.parse_args()

    aln_path = Path(args.alignment)
    colmap_path = Path(args.column_map)
    out_dir = Path(args.out_dir)
    if not aln_path.exists():
        raise SystemExit(f"ERROR: alignment not found: {aln_path}")
    if not colmap_path.exists():
        raise SystemExit(f"ERROR: column map not found: {colmap_path}")

    entries = read_fasta(aln_path)
    colmap = read_column_map(colmap_path)
    if not entries:
        raise SystemExit("ERROR: empty alignment")
    query_seq = find_query(entries)
    if sum(aa != "-" for aa in query_seq) != 437:
        raise SystemExit("ERROR: query sequence does not span 437 ANKros residues")

    rows = compute_profile(entries, colmap)
    out_png = out_dir / "conservation_profile.png"
    out_tsv = out_dir / "conservation_profile.tsv"
    plot_profile(rows, out_png, out_tsv, args.smooth_window)

    print(f"Saved: {out_png}")
    print(f"Saved: {out_tsv}")


if __name__ == "__main__":
    main()
