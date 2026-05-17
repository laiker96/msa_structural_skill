#!/usr/bin/env python3
"""Step 10: write group consensus tables and sequence-logo PNGs."""
from __future__ import annotations

import argparse
import csv
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
    parser.add_argument("--group-members", type=Path, default=cfg.OUTPUT_DIR / "vulnerability" / "group_members.tsv")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "logos")
    parser.add_argument("--wrap", type=int, default=80)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_groups(path: Path) -> dict[str, set[str]]:
    groups: dict[str, set[str]] = {}
    for row in read_tsv(path):
        groups.setdefault(row["group"], set()).add(row["id"])
    return groups


def column_stats(entries: list[tuple[str, str, str]]) -> list[dict[str, object]]:
    width = len(entries[0][2])
    out = []
    for col in range(width):
        residues = [seq[col].upper() for _sid, _hdr, seq in entries if seq[col].upper() in cfg.AA_SET]
        counts = Counter(residues)
        n = sum(counts.values())
        top_aa, top_count = counts.most_common(1)[0] if counts else ("-", 0)
        probs = {aa: count / n for aa, count in counts.items()} if n else {}
        entropy = -sum(p * math.log2(p) for p in probs.values()) if probs else 0.0
        information = (math.log2(len(cfg.AA_SET)) - entropy) * (n / len(entries)) if entries else 0.0
        out.append({
            "alignment_col": col,
            "top_aa": top_aa,
            "top_frequency": top_count / n if n else 0.0,
            "occupancy": n / len(entries) if entries else 0.0,
            "information_bits": information,
            "counts": counts,
        })
    return out


def draw_logo(group: str, stats: list[dict[str, object]], out_png: Path, wrap: int) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(out_png.parent / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        **{aa: "#2f6fbb" for aa in "AILMVP"},
        **{aa: "#d1495b" for aa in "KRH"},
        **{aa: "#8f2d56" for aa in "DE"},
        **{aa: "#2a9d8f" for aa in "STNQ"},
        **{aa: "#f4a261" for aa in "FWY"},
        "G": "#6d597a",
        "C": "#7b2cbf",
    }
    chunks = [stats[i:i + wrap] for i in range(0, len(stats), wrap)]
    fig, axes = plt.subplots(len(chunks), 1, figsize=(max(12, wrap * 0.22), max(3, len(chunks) * 2.4)), squeeze=False)
    for ax, chunk in zip(axes.ravel(), chunks):
        ax.set_xlim(0, len(chunk))
        ax.set_ylim(0, math.log2(len(cfg.AA_SET)) + 0.2)
        ax.set_ylabel("bits")
        for x, stat in enumerate(chunk):
            counts: Counter = stat["counts"]
            n = sum(counts.values())
            y = 0.0
            for aa, count in sorted(counts.items(), key=lambda item: item[1]):
                height = (count / n) * float(stat["information_bits"]) if n else 0.0
                if height <= 0.03:
                    continue
                ax.text(x + 0.5, y, aa, ha="center", va="bottom", fontsize=8 + height * 6, color=colors.get(aa, "#444444"), fontweight="bold")
                y += height
        ax.set_xticks([])
    fig.suptitle(f"Sequence logo: {group}")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    entries = cfg.read_fasta(args.alignment)
    if not entries:
        raise SystemExit("ERROR: empty alignment")
    groups = load_groups(args.group_members)
    by_id = {sid: (sid, hdr, seq) for sid, hdr, seq in entries}
    rows = []
    for group, ids in sorted(groups.items()):
        selected = [by_id[sid] for sid in ids if sid in by_id]
        if not selected:
            continue
        stats = column_stats(selected)
        consensus = "".join(str(stat["top_aa"]) for stat in stats if stat["top_aa"] != "-")
        rows.append({"group": group, "n_sequences": len(selected), "consensus_sequence": consensus})
        draw_logo(group, stats, args.out_dir / f"{group}.logo.png", args.wrap)
    cfg.write_tsv(args.out_dir / "group_consensus.tsv", rows, ["group", "n_sequences", "consensus_sequence"])
    print(f"Saved logos and consensus for {len(rows)} groups under: {args.out_dir}")


if __name__ == "__main__":
    main()
