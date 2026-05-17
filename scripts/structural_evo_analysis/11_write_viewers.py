#!/usr/bin/env python3
"""Step 11: write self-contained HTML tree and alignment viewers."""
from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

cfg = import_module("00_config")


COLORS = {
    "low": "#2563eb",
    "mid": "#16a34a",
    "high": "#dc2626",
    "psychro": "#2563eb",
    "meso": "#16a34a",
    "thermo": "#dc2626",
    "": "#6b7280",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, default=cfg.OUTPUT_DIR / "repset_aligned.fa")
    parser.add_argument("--tree", type=Path, default=cfg.OUTPUT_DIR / "tree" / "query_msa.treefile")
    parser.add_argument("--metadata", type=Path, default=cfg.OUTPUT_DIR / "repset_metadata.tsv")
    parser.add_argument("--tip-metadata", type=Path, default=None)
    parser.add_argument("--called-clades", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "viewers")
    parser.add_argument("--max-alignment-seqs", type=int, default=300)
    return parser.parse_args()


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metadata_by_id(path: Path) -> dict[str, dict[str, str]]:
    return {row["id"]: row for row in read_tsv(path) if row.get("id")}


def group_by_id(group_members: Path) -> dict[str, str]:
    out = {}
    for row in read_tsv(group_members):
        out[row["id"]] = row["group"]
    return out


def tree_coordinates(tree):
    terminals = tree.get_terminals()
    y = {tip: idx for idx, tip in enumerate(terminals)}
    def assign(clade):
        if clade in y:
            return y[clade]
        vals = [assign(child) for child in clade.clades]
        y[clade] = sum(vals) / len(vals)
        return y[clade]
    assign(tree.root)
    depths = tree.depths()
    max_depth = max(depths.values()) or 1.0
    return depths, y, max_depth, len(terminals)


def write_tree_html(args: argparse.Namespace, meta: dict[str, dict[str, str]], groups: dict[str, str]) -> None:
    if not args.tree.exists():
        return
    tree = Phylo.read(args.tree, "newick")
    depths, yvals, max_depth, n_terms = tree_coordinates(tree)
    width = 1200
    height = max(420, n_terms * 18 + 80)
    left, right, top = 35, 360, 35
    scale_x = (width - left - right) / max_depth
    scale_y = (height - top * 2) / max(n_terms - 1, 1)
    def x(clade): return left + depths[clade] * scale_x
    def y(clade): return top + yvals[clade] * scale_y
    parts = []
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parts.append(f'<path d="M{x(clade):.2f},{y(child):.2f} L{x(child):.2f},{y(child):.2f}" class="branch"/>')
            parts.append(f'<path d="M{x(clade):.2f},{y(clade):.2f} L{x(clade):.2f},{y(child):.2f}" class="branch"/>')
    for tip in tree.get_terminals():
        row = meta.get(tip.name, {})
        category = row.get("trait_category") or row.get("regime", "")
        color = COLORS.get(category, COLORS[""])
        label = row.get("organism") or tip.name
        group = groups.get(tip.name, "")
        title = html.escape(f"{tip.name}\\n{label}\\nOGT={row.get('ogt','')}\\ngroup={group}", quote=True)
        parts.append(f'<circle cx="{x(tip):.2f}" cy="{y(tip):.2f}" r="4" fill="{color}"><title>{title}</title></circle>')
        parts.append(f'<text x="{x(tip)+8:.2f}" y="{y(tip)+4:.2f}" class="tip">{html.escape(label[:90])}</text>')
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>MSA tree</title>
<style>body{{font-family:Arial,sans-serif;margin:0}}.wrap{{padding:16px}}svg{{width:100%;height:auto}}.branch{{stroke:#374151;stroke-width:1;fill:none}}.tip{{font-size:11px;fill:#111827}}</style></head>
<body><div class="wrap"><h1>MSA Tree</h1><svg viewBox="0 0 {width} {height}">{''.join(parts)}</svg></div></body></html>"""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "tree.html").write_text(doc)


def write_alignment_html(args: argparse.Namespace, meta: dict[str, dict[str, str]], groups: dict[str, str]) -> None:
    entries = cfg.read_fasta(args.alignment)[: args.max_alignment_seqs]
    rows = []
    for sid, _hdr, seq in entries:
        row = meta.get(sid, {})
        group = groups.get(sid, "all_representatives")
        label = row.get("organism") or sid
        residues = "".join(f'<span class="aa aa-{html.escape(aa)}">{html.escape(aa)}</span>' for aa in seq)
        rows.append(f'<tr data-group="{html.escape(group)}"><th>{html.escape(label[:60])}<br><small>{html.escape(group)}</small></th><td>{residues}</td></tr>')
    group_names = sorted(set(groups.values()) or {"all_representatives"})
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>MSA alignment</title>
<style>body{{font-family:Arial,sans-serif;margin:0}}.wrap{{padding:16px}}.scroll{{overflow:auto;border:1px solid #d1d5db;max-height:80vh}}table{{border-collapse:collapse}}th{{position:sticky;left:0;background:#fff;text-align:left;font-size:12px;min-width:260px;border-right:1px solid #ddd}}td{{font-family:ui-monospace,monospace;white-space:nowrap}}.aa{{display:inline-block;width:0.72em;text-align:center}}.aa--{{color:#bbb}}select{{margin-bottom:10px}}</style>
<script>function filt(v){{document.querySelectorAll('tbody tr').forEach(r=>r.style.display=(v==='all'||r.dataset.group===v)?'':'none')}}</script></head>
<body><div class="wrap"><h1>MSA Alignment</h1><label>Group <select onchange="filt(this.value)"><option value="all">all</option>{''.join(f'<option>{html.escape(g)}</option>' for g in group_names)}</select></label><div class="scroll"><table><tbody>{''.join(rows)}</tbody></table></div></div></body></html>"""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "alignment.html").write_text(doc)


def main() -> None:
    args = parse_args()
    meta = metadata_by_id(args.metadata)
    if args.tip_metadata and args.tip_metadata.exists():
        for row in read_tsv(args.tip_metadata):
            if row.get("id") in meta:
                meta[row["id"]].update(row)
    groups = group_by_id(args.out_dir.parent / "vulnerability" / "group_members.tsv")
    write_tree_html(args, meta, groups)
    write_alignment_html(args, meta, groups)
    print(f"Saved viewers under: {args.out_dir}")


if __name__ == "__main__":
    main()
