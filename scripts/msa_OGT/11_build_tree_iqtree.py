#!/usr/bin/env python3
"""
Step 11: IQ-TREE on the antenna + catalytic domains.

The linker is intentionally excluded from the default tree input because it is
indel-rich and alignment-sensitive. This step extracts ANKros query-position
columns for:

  - antenna:   1-130
  - catalytic: 206-437

from the refined full-domain alignment, writes a dedicated tree MSA, and runs
IQ-TREE with a fixed seed/thread count for reproducibility. By default, curated
class II CPD photolyase outgroups from test/outgroup.fa are added with
MAFFT --add --keeplength and passed to IQ-TREE with -o.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

ALN_FA = cfg.INTER_DIR / "repset_hmmalign_linker_refined.fa"
COL_MAP = cfg.INTER_DIR / "repset_hmmalign_linker_refined_column_map.tsv"
TREE_DIR = cfg.INTER_DIR / "tree_antenna_catalytic"
OUTGROUP_FA = cfg.PROJECT_ROOT / "test" / "outgroup.fa"
IQTREE_BIN = cfg.resolve_bin("iqtree")
MAFFT_BIN = cfg.resolve_bin("mafft")

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


def read_column_map(path: Path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iqtree_version() -> str:
    try:
        out = subprocess.run(
            [IQTREE_BIN, "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception as exc:
        return f"(unavailable: {exc})"


def parse_domains(domain_text: str):
    domains = [part.strip() for part in domain_text.split(",") if part.strip()]
    unknown = [name for name in domains if name not in DOMAINS]
    if unknown:
        raise SystemExit(f"ERROR: unknown domain(s): {', '.join(unknown)}")
    if not domains:
        raise SystemExit("ERROR: no domains selected")
    return domains


def select_columns(colmap, domains):
    selected = []
    domain_set = set(domains)
    for row in colmap:
        qpos = (row.get("qpos") or "").strip()
        region = (row.get("region") or "").strip()
        if qpos and region in domain_set:
            selected.append(row)
    return selected


def subset_alignment(entries, selected_rows):
    cols = [int(row["out_col"]) for row in selected_rows]
    out = []
    for sid, hdr, seq in entries:
        missing = [col for col in cols if col >= len(seq)]
        if missing:
            raise SystemExit(f"ERROR: sequence {sid} is shorter than mapped column {missing[0]}")
        sub = "".join(seq[col] for col in cols)
        if set(sub) <= {"-", "."}:
            continue
        out.append((sid, hdr, sub))
    return out


def write_selected_map(selected_rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "tree_col", "source_alignment_col", "source", "source_col",
        "source_block", "source_range", "qpos", "region",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for tree_col, row in enumerate(selected_rows):
            writer.writerow({
                "tree_col": tree_col,
                "source_alignment_col": row.get("out_col", ""),
                "source": row.get("source", ""),
                "source_col": row.get("source_col", ""),
                "source_block": row.get("source_block", ""),
                "source_range": row.get("source_range", ""),
                "qpos": row.get("qpos", ""),
                "region": row.get("region", ""),
            })


def add_outgroup_alignment(ingroup_msa: Path, tree_msa: Path, outgroup_fa: Path,
                           outgroup_ids_path: Path):
    """Add curated outgroups to the prepared tree MSA without changing its columns."""
    if not outgroup_fa.exists():
        shutil.copy2(ingroup_msa, tree_msa)
        outgroup_ids_path.write_text("")
        return []

    og_entries = read_fasta(outgroup_fa)
    if not og_entries:
        shutil.copy2(ingroup_msa, tree_msa)
        outgroup_ids_path.write_text("")
        return []

    cmd = [
        MAFFT_BIN, "--add", str(outgroup_fa), "--keeplength",
        "--thread", str(cfg.N_THREADS), str(ingroup_msa),
    ]
    print("Adding outgroup: " + " ".join(cmd))
    with open(tree_msa, "w") as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(res.stderr[-1500:])
        raise SystemExit(f"MAFFT --add outgroup failed with exit code {res.returncode}")

    outgroup_ids = [sid for sid, _hdr, _seq in og_entries]
    outgroup_ids_path.write_text("\n".join(outgroup_ids) + "\n")
    return outgroup_ids


def write_summary(path: Path, entries, selected_rows, domains, source_alignment, source_map,
                  ingroup_msa, tree_msa, outgroup_fa, outgroup_ids, ingroup_count):
    counts = {name: 0 for name in DOMAINS}
    for row in selected_rows:
        region = row.get("region", "")
        counts[region] = counts.get(region, 0) + 1
    width = len(entries[0][2]) if entries else 0
    nseq = len(entries)
    non_gap = sum(sum(aa not in "-." for aa in seq) for _sid, _hdr, seq in entries)
    occupancy = non_gap / (nseq * width) if nseq and width else 0.0
    summary = {
        "n_sequences": nseq,
        "n_ingroup_sequences": ingroup_count,
        "n_outgroup_sequences": len(outgroup_ids),
        "alignment_width": width,
        "domains": ",".join(domains),
        "antenna_columns": counts.get("antenna", 0),
        "linker_columns": counts.get("linker", 0),
        "catalytic_columns": counts.get("catalytic", 0),
        "mean_occupancy": f"{occupancy:.4f}",
        "source_alignment": str(source_alignment),
        "source_column_map": str(source_map),
        "ingroup_alignment": str(ingroup_msa),
        "outgroup_fasta": str(outgroup_fa) if outgroup_ids else "",
        "outgroup_ids": ",".join(outgroup_ids),
        "tree_alignment": str(tree_msa),
    }
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary)


def prepare_tree_alignment(aln_path: Path, map_path: Path, tree_msa: Path, tree_map: Path,
                           summary_path: Path, domains, outgroup_fa: Path | None,
                           outgroup_ids_path: Path):
    entries = read_fasta(aln_path)
    colmap = read_column_map(map_path)
    if not entries:
        raise SystemExit("ERROR: empty source alignment")
    selected_rows = select_columns(colmap, domains)
    if not selected_rows:
        raise SystemExit("ERROR: selected domains produced zero columns")
    out_entries = subset_alignment(entries, selected_rows)
    if len(out_entries) < 4:
        raise SystemExit(f"ERROR: too few non-empty sequences for tree: {len(out_entries)}")
    ingroup_msa = tree_msa.with_name(f"{tree_msa.stem}_ingroup{tree_msa.suffix}")
    write_fasta(out_entries, ingroup_msa)
    if outgroup_fa is None:
        shutil.copy2(ingroup_msa, tree_msa)
        outgroup_ids_path.write_text("")
        outgroup_ids = []
    else:
        outgroup_ids = add_outgroup_alignment(ingroup_msa, tree_msa, outgroup_fa, outgroup_ids_path)
    final_entries = read_fasta(tree_msa)
    write_selected_map(selected_rows, tree_map)
    write_summary(
        summary_path, final_entries, selected_rows, domains, aln_path, map_path,
        ingroup_msa, tree_msa, outgroup_fa, outgroup_ids, len(out_entries),
    )
    return final_entries, selected_rows, outgroup_ids


def write_manifest(args, cmd, tree_msa, domains, manifest_path: Path, outgroup_ids):
    is_auto = str(args.threads).upper() == "AUTO"
    body = [
        "# IQ-TREE run manifest",
        f"timestamp_utc: {_dt.datetime.now(_dt.UTC).isoformat(timespec='seconds')}",
        f"mode:          {'full' if args.full else 'fast'}",
        f"domains:       {','.join(domains)}",
        f"model:         {args.model if args.full else args.fast_model + ' (-fast)'}",
        f"seed:          {args.seed}",
        f"threads:       {args.threads}"
            + (" [WARNING: AUTO is non-deterministic]" if is_auto else " [deterministic]"),
        f"bootstrap:     {args.bootstrap if args.full else 0} (bcor={args.bcor})",
        f"iqtree:        {iqtree_version()}",
        f"msa:           {tree_msa}",
        f"msa_sha256:    {sha256_of(tree_msa)}",
        f"outgroup:      {','.join(outgroup_ids) if outgroup_ids else '(none; unrooted)'}",
        f"command:       {' '.join(cmd)}",
    ]
    manifest_path.write_text("\n".join(body) + "\n")


def read_manifest(path: Path):
    values = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        if ":" not in line or line.startswith("#"):
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip()
    return values


def run_iqtree(args, tree_msa: Path, prefix: Path, domains, tree_copy: Path,
               manifest_path: Path, outgroup_ids):
    treefile = Path(str(prefix) + ".treefile")
    if treefile.exists() and not args.force:
        current_hash = sha256_of(tree_msa)
        current_outgroup = ",".join(outgroup_ids) if outgroup_ids else "(none; unrooted)"
        manifest = read_manifest(manifest_path)
        if (manifest.get("msa_sha256") == current_hash
                and manifest.get("outgroup") == current_outgroup):
            print(f"[skip] {treefile.name} exists and matches current MSA/outgroup; use --force to rebuild")
        else:
            print(
                f"[skip] {treefile.name} exists; reusing it without rerunning IQ-TREE. "
                "Use --force to rebuild from the current prepared MSA/outgroup."
            )
        if not tree_copy.exists():
            shutil.copy2(treefile, tree_copy)
        return treefile

    if args.force or treefile.exists():
        for path in prefix.parent.glob(prefix.name + ".*"):
            path.unlink()

    outgroup_arg = ["-o", ",".join(outgroup_ids)] if outgroup_ids else []

    if args.full:
        cmd = [
            IQTREE_BIN, "-s", str(tree_msa),
            "-m", args.model,
            "-T", str(args.threads),
            "-seed", str(args.seed),
            "--prefix", str(prefix),
            "--redo",
        ] + outgroup_arg
        if args.bootstrap > 0:
            cmd += ["-B", str(args.bootstrap), "-bcor", str(args.bcor)]
    else:
        cmd = [
            IQTREE_BIN, "-s", str(tree_msa),
            "-m", args.fast_model,
            "-fast",
            "-nt", str(args.threads),
            "-seed", str(args.seed),
            "-pre", str(prefix),
            "-quiet",
            "-redo",
        ] + outgroup_arg

    print("Running: " + " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-1500:])
        print(res.stderr[-1000:])
        raise SystemExit(f"IQ-TREE failed with exit code {res.returncode}")
    if not treefile.exists():
        raise SystemExit(f"ERROR: IQ-TREE finished but did not write {treefile}")

    shutil.copy2(treefile, tree_copy)
    write_manifest(args, cmd, tree_msa, domains, manifest_path, outgroup_ids)
    return treefile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", default=str(ALN_FA))
    parser.add_argument("--column-map", default=str(COL_MAP))
    parser.add_argument("--out-dir", default=str(TREE_DIR))
    parser.add_argument("--domains", default="antenna,catalytic",
                        help="Comma-separated domain names. Default: antenna,catalytic")
    parser.add_argument("--outgroup", default=str(OUTGROUP_FA),
                        help="Outgroup FASTA to add with MAFFT --add. Default: test/outgroup.fa")
    parser.add_argument("--no-outgroup", action="store_true",
                        help="Do not add outgroups; build an unrooted tree.")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Write the tree alignment and metadata but do not run IQ-TREE.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing IQ-TREE outputs for this prefix and rebuild the tree.")
    parser.add_argument("--fast", dest="full", action="store_false",
                        help="Run the quick LG -fast topology instead of the default full ML tree.")
    parser.set_defaults(full=True)
    parser.add_argument("--fast-model", default="LG",
                        help="Model used for fast mode. Default: LG")
    parser.add_argument("--model", default="LG+F+R4",
                        help="Model used for the default full ML mode. Default: LG+F+R4")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--bcor", type=float, default=0.95)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", default=str(cfg.N_THREADS))
    args = parser.parse_args()

    aln_path = Path(args.alignment)
    map_path = Path(args.column_map)
    out_dir = Path(args.out_dir)
    if not aln_path.exists():
        raise SystemExit(f"ERROR: alignment not found: {aln_path}")
    if not map_path.exists():
        raise SystemExit(f"ERROR: column map not found: {map_path}")

    domains = parse_domains(args.domains)
    domain_slug = "_".join(domains)
    tree_msa = out_dir / f"{domain_slug}.fa"
    tree_map = out_dir / f"{domain_slug}_column_map.tsv"
    summary_path = out_dir / f"{domain_slug}_summary.tsv"
    outgroup_ids_path = out_dir / f"{domain_slug}_outgroup_ids.txt"
    prefix = out_dir / f"ankros_{domain_slug}"
    tree_copy = out_dir.parent / f"ankros_{domain_slug}.nw"
    manifest_path = out_dir / "run_manifest.txt"
    treefile = Path(str(prefix) + ".treefile")

    existing_bundle = [tree_msa, tree_map, summary_path, treefile]
    if not args.prepare_only and not args.force and all(path.exists() for path in existing_bundle):
        print(
            f"[skip] Existing tree bundle found in {out_dir}; "
            "use --force to regenerate the tree alignment and rerun IQ-TREE."
        )
        if not tree_copy.exists():
            shutil.copy2(treefile, tree_copy)
            print(f"Copy: {tree_copy}")
        else:
            print(f"Tree: {treefile}")
            print(f"Copy: {tree_copy}")
        print(f"Manifest: {manifest_path}")
        return

    outgroup_fa = None if args.no_outgroup else Path(args.outgroup)
    entries, selected_rows, outgroup_ids = prepare_tree_alignment(
        aln_path, map_path, tree_msa, tree_map, summary_path, domains,
        outgroup_fa, outgroup_ids_path)
    print(f"Saved: {tree_msa} ({len(entries)} seqs x {len(selected_rows)} cols)")
    print(f"Saved: {tree_map}")
    print(f"Saved: {summary_path}")
    if outgroup_ids:
        print(f"Outgroup: {len(outgroup_ids)} sequences ({', '.join(outgroup_ids)})")
        print(f"Saved: {outgroup_ids_path}")
    else:
        print("Outgroup: none (tree will be unrooted)")

    if args.prepare_only:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    treefile = run_iqtree(args, tree_msa, prefix, domains, tree_copy, manifest_path, outgroup_ids)
    print(f"Tree: {treefile}")
    print(f"Copy: {tree_copy}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
