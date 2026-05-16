#!/usr/bin/env python3
"""Step 03: prepare a tree alignment and run IQ-TREE."""
from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, default=cfg.OUTPUT_DIR / "repset_aligned.fa")
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "tree")
    parser.add_argument("--prefix", default="query_msa")
    parser.add_argument("--min-column-occupancy", type=float, default=0.50)
    parser.add_argument("--keep-query-gap-columns", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--fast", action="store_true", help="Use IQ-TREE -fast with LG.")
    parser.add_argument("--model", default="LG+F+R4")
    parser.add_argument("--fast-model", default="LG")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", default=str(cfg.N_THREADS))
    return parser.parse_args()


def selected_columns(entries: list[tuple[str, str, str]], min_occupancy: float, keep_query_gaps: bool) -> list[int]:
    width = len(entries[0][2])
    query_seq = entries[0][2]
    selected = []
    for col in range(width):
        if not keep_query_gaps and query_seq[col].upper() not in cfg.AA_SET:
            continue
        occupancy = sum(seq[col].upper() in cfg.AA_SET for _sid, _hdr, seq in entries) / len(entries)
        if occupancy >= min_occupancy:
            selected.append(col)
    return selected


def subset_alignment(entries: list[tuple[str, str, str]], cols: list[int]) -> list[tuple[str, str, str]]:
    out = []
    for sid, header, seq in entries:
        sub = "".join(seq[col] for col in cols)
        if any(aa.upper() in cfg.AA_SET for aa in sub):
            out.append((sid, header, sub))
    return out


def write_column_map(path: Path, cols: list[int]) -> None:
    rows = [{"tree_col": idx, "source_alignment_col": col} for idx, col in enumerate(cols)]
    cfg.write_tsv(path, rows, ["tree_col", "source_alignment_col"])


def iqtree_version(iqtree: str) -> str:
    try:
        completed = subprocess.run([iqtree, "--version"], capture_output=True, text=True, timeout=10)
        return (completed.stdout or completed.stderr).strip().splitlines()[0]
    except Exception as exc:
        return f"unavailable: {exc}"


def write_manifest(args: argparse.Namespace, cmd: list[str], tree_msa: Path, manifest: Path) -> None:
    body = [
        "# IQ-TREE run manifest",
        f"timestamp_utc: {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"source_alignment: {args.alignment}",
        f"tree_alignment: {tree_msa}",
        f"tree_alignment_sha256: {cfg.sha256_of(tree_msa)}",
        f"min_column_occupancy: {args.min_column_occupancy}",
        f"keep_query_gap_columns: {args.keep_query_gap_columns}",
        f"mode: {'fast' if args.fast else 'full'}",
        f"model: {args.fast_model if args.fast else args.model}",
        f"bootstrap: {0 if args.fast else args.bootstrap}",
        f"seed: {args.seed}",
        f"threads: {args.threads}",
        f"iqtree: {iqtree_version(cfg.resolve_bin('iqtree'))}",
        f"command: {' '.join(cmd)}",
    ]
    manifest.write_text("\n".join(body) + "\n")


def run_iqtree(args: argparse.Namespace, tree_msa: Path) -> Path:
    iqtree = cfg.resolve_bin("iqtree")
    prefix = args.out_dir / args.prefix
    treefile = Path(str(prefix) + ".treefile")
    if treefile.exists() and not args.force:
        print(f"Using existing tree: {treefile}")
        return treefile
    if args.force:
        for path in args.out_dir.glob(args.prefix + ".*"):
            path.unlink()
    if args.fast:
        cmd = [
            iqtree, "-s", str(tree_msa), "-m", args.fast_model, "-fast",
            "-nt", str(args.threads), "-seed", str(args.seed), "-pre", str(prefix),
            "-quiet", "-redo",
        ]
    else:
        cmd = [
            iqtree, "-s", str(tree_msa), "-m", args.model, "-T", str(args.threads),
            "-seed", str(args.seed), "--prefix", str(prefix), "--redo",
        ]
        if args.bootstrap > 0:
            cmd += ["-B", str(args.bootstrap)]
    print("Running: " + " ".join(cmd))
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        print(completed.stdout[-1500:])
        print(completed.stderr[-1500:])
        raise SystemExit(f"IQ-TREE failed with exit code {completed.returncode}")
    if not treefile.exists():
        raise SystemExit(f"ERROR: IQ-TREE finished but did not write {treefile}")
    write_manifest(args, cmd, tree_msa, args.out_dir / "run_manifest.txt")
    shutil.copy2(treefile, args.out_dir.parent / "query_msa.nwk")
    return treefile


def main() -> None:
    args = parse_args()
    if not args.alignment.exists():
        raise SystemExit(f"ERROR: alignment not found: {args.alignment}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    entries = cfg.read_fasta(args.alignment)
    if len(entries) < 4:
        raise SystemExit(f"ERROR: at least four aligned sequences are needed for a tree; got {len(entries)}")
    cols = selected_columns(entries, args.min_column_occupancy, args.keep_query_gap_columns)
    if not cols:
        raise SystemExit("ERROR: no columns passed tree-alignment filters")
    tree_entries = subset_alignment(entries, cols)
    tree_msa = args.out_dir / "tree_alignment.fa"
    cfg.write_fasta(tree_entries, tree_msa)
    write_column_map(args.out_dir / "tree_column_map.tsv", cols)
    summary = [{
        "n_sequences": len(tree_entries),
        "source_width": len(entries[0][2]),
        "tree_width": len(cols),
        "min_column_occupancy": args.min_column_occupancy,
        "tree_alignment": str(tree_msa),
    }]
    cfg.write_tsv(
        args.out_dir / "tree_alignment_summary.tsv",
        summary,
        ["n_sequences", "source_width", "tree_width", "min_column_occupancy", "tree_alignment"],
    )
    print(f"Saved: {tree_msa} ({len(tree_entries)} seqs x {len(cols)} cols)")
    if args.prepare_only:
        return
    treefile = run_iqtree(args, tree_msa)
    print(f"Tree: {treefile}")


if __name__ == "__main__":
    main()
