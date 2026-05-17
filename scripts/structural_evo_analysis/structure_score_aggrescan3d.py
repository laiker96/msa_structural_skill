#!/usr/bin/env python3
"""Structural-evolution Aggrescan3D scoring for query and AFDB structures."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

DEFAULT_AFDB_DIR = cfg.STRUCTURE_DIR / "afdb"
DEFAULT_QUERY_PDB = cfg.STRUCTURE_DIR / "query.pdb"
DEFAULT_OUT_DIR = cfg.OUTPUT_DIR / "structure_scores" / "aggrescan3d"
DEFAULT_AGGRESCAN = ROOT / "envs" / "aggrescan3d" / "bin" / "aggrescan"
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
DEFAULT_VALIDATION_STATUSES = {
    "pass_quick_classI_structural_match",
    "review_borderline_not_mismatch",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Aggrescan3D and collect global/per-residue aggregation scores."
    )
    parser.add_argument("--afdb-dir", type=Path, default=DEFAULT_AFDB_DIR)
    parser.add_argument(
        "--query-pdb",
        "--ankros-pdb",
        dest="query_pdb",
        type=Path,
        default=DEFAULT_QUERY_PDB,
        help="Optional local query PDB. --ankros-pdb is a legacy alias.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--aggrescan-bin", type=Path, default=DEFAULT_AGGRESCAN)
    parser.add_argument("--distance", type=int, default=10)
    parser.add_argument(
        "--patch-radius",
        type=float,
        default=6.0,
        help="Side-chain atom contact radius used to cluster positive A3D residues into patches.",
    )
    parser.add_argument("--ph", type=float, default=None)
    parser.add_argument("--chain", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--pattern", default="*.pdb")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run Aggrescan3D even when complete raw outputs already exist.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker processes for per-structure Aggrescan3D runs.",
    )
    parser.add_argument(
        "--structure-validation",
        type=Path,
        help="Optional structural validation TSV; AFDB structures with excluded statuses are skipped.",
    )
    parser.add_argument(
        "--validation-statuses",
        default=",".join(sorted(DEFAULT_VALIDATION_STATUSES)),
        help="Comma-separated validation statuses to include.",
    )
    return parser.parse_args()


def structure_id_from_path(path: Path) -> str:
    if path.name.startswith("AF-"):
        match = re.match(r"AF-(.+)-F\d+-model", path.stem)
        if match:
            return match.group(1)
    return path.stem


def load_validated_accessions(path: Path | None, statuses: str) -> set[str] | None:
    if path is None:
        return None
    allowed_statuses = {status.strip() for status in statuses.split(",") if status.strip()}
    accessions = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession", "status"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: missing structural validation columns in {path}: {sorted(missing)}")
        for row in reader:
            if row.get("status", "") in allowed_statuses and row.get("accession", ""):
                accessions.add(row["accession"])
    return accessions


def input_structures(args: argparse.Namespace) -> list[Path]:
    structures = []
    if args.query_pdb.exists():
        structures.append(args.query_pdb)
    validated = load_validated_accessions(args.structure_validation, args.validation_statuses)
    for pdb in sorted(args.afdb_dir.glob(args.pattern)):
        if validated is not None and structure_id_from_path(pdb) not in validated:
            continue
        structures.append(pdb)
    if args.limit is not None:
        structures = structures[: args.limit]
    return structures


def run_aggrescan(args: argparse.Namespace, pdb: Path, work_dir: Path) -> None:
    if not args.force and (work_dir / "A3D.csv").exists() and (work_dir / "A3D_summary.json").exists():
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    mpl_dir = args.out_dir / ".matplotlib"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(args.aggrescan_bin),
        "-i", str(pdb),
        "-w", str(work_dir),
        "-D", str(args.distance),
        "-v", "1",
    ]
    if args.chain:
        cmd.extend(["-C", args.chain])
    if args.ph is not None:
        cmd.extend(["--ph", str(args.ph)])
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(mpl_dir)
    completed = subprocess.run(
        cmd,
        cwd=str(ROOT),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (work_dir / "aggrescan.log").write_text(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(f"Aggrescan3D failed for {pdb} with exit code {completed.returncode}")


def parse_residue_scores(structure_id: str, source_path: Path, csv_path: Path) -> list[dict]:
    rows = []
    with csv_path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "structure_id": structure_id,
                "source_path": str(source_path),
                "protein_state": row["protein"],
                "chain": row["chain"],
                "resseq": int(row["residue"]),
                "residue_id": f"{row['chain']}:{row['residue']}",
                "one_letter": row["residue_name"],
                "aggrescan3d_score": float(row["score"]),
            })
    return rows


def positive_burden(rows: list[dict]) -> dict:
    scores = [float(row["aggrescan3d_score"]) for row in rows]
    if not scores:
        return {
            "positive_aggrescan3d_burden": float("nan"),
            "positive_aggrescan3d_residues": 0,
            "residue_count": 0,
        }
    positive = [score for score in scores if score > 0.0]
    return {
        "positive_aggrescan3d_burden": sum(positive) / len(scores),
        "positive_aggrescan3d_residues": len(positive),
        "residue_count": len(scores),
    }


def residue_contact_atoms(pdb_path: Path) -> dict[tuple[str, int], tuple]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    atoms_by_residue = {}
    for residue in structure.get_residues():
        hetflag, resseq, _icode = residue.id
        if hetflag.strip():
            continue
        chain = residue.get_parent().id
        sidechain_atoms = [
            atom for atom in residue.get_atoms()
            if atom.get_name().strip() not in BACKBONE_ATOMS
            and (atom.element or "").strip().upper() != "H"
        ]
        if not sidechain_atoms and "CA" in residue:
            sidechain_atoms = [residue["CA"]]
        atoms_by_residue[(chain, resseq)] = tuple(sidechain_atoms)
    return atoms_by_residue


def cluster_positive_patches(
    rows: list[dict],
    scored_pdb: Path,
    patch_radius: float,
) -> dict:
    positive_rows = [row for row in rows if float(row["aggrescan3d_score"]) > 0.0]
    if not positive_rows or not scored_pdb.exists():
        return {
            "positive_patch_count": 0,
            "max_positive_patch_size": 0,
            "max_positive_patch_sum": 0.0,
            "max_positive_patch_mean": 0.0,
            "max_positive_patch_residues": "",
        }

    atoms_by_residue = residue_contact_atoms(scored_pdb)
    indexed_rows = [
        row for row in positive_rows
        if atoms_by_residue.get((row["chain"], int(row["resseq"])))
    ]
    if not indexed_rows:
        return {
            "positive_patch_count": 0,
            "max_positive_patch_size": 0,
            "max_positive_patch_sum": 0.0,
            "max_positive_patch_mean": 0.0,
            "max_positive_patch_residues": "",
        }

    parent = list(range(len(indexed_rows)))

    def find(idx: int) -> int:
        while parent[idx] != idx:
            parent[idx] = parent[parent[idx]]
            idx = parent[idx]
        return idx

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_left] = root_right

    atoms = []
    owner = {}
    for idx, row in enumerate(indexed_rows):
        for atom in atoms_by_residue[(row["chain"], int(row["resseq"]))]:
            owner[id(atom)] = idx
            atoms.append(atom)

    if atoms:
        search = NeighborSearch(atoms)
        for atom_a, atom_b in search.search_all(patch_radius, level="A"):
            left = owner[id(atom_a)]
            right = owner[id(atom_b)]
            if left != right:
                union(left, right)

    groups = {}
    for idx, row in enumerate(indexed_rows):
        groups.setdefault(find(idx), []).append(row)
    patches = list(groups.values())

    def patch_sum(patch: list[dict]) -> float:
        return sum(float(row["aggrescan3d_score"]) for row in patch)

    best = max(patches, key=lambda patch: (patch_sum(patch), len(patch)))
    best_sum = patch_sum(best)
    best_residues = "+".join(
        f"{row['one_letter']}{row['resseq']}"
        for row in sorted(best, key=lambda row: (row["chain"], int(row["resseq"])))
    )
    return {
        "positive_patch_count": len(patches),
        "max_positive_patch_size": len(best),
        "max_positive_patch_sum": best_sum,
        "max_positive_patch_mean": best_sum / len(best),
        "max_positive_patch_residues": best_residues,
    }


def parse_global_scores(
    structure_id: str,
    source_path: Path,
    summary_path: Path,
    residue_rows: list[dict],
    scored_pdb: Path,
    patch_radius: float,
) -> list[dict]:
    summary = json.loads(summary_path.read_text())
    rows = []
    for chain, values in summary.items():
        if chain == "All":
            chain_residues = residue_rows
        else:
            chain_residues = [row for row in residue_rows if row["chain"] == chain]
        patch_metrics = cluster_positive_patches(chain_residues, scored_pdb, patch_radius)
        rows.append({
            "structure_id": structure_id,
            "source_path": str(source_path),
            "chain": chain,
            **positive_burden(chain_residues),
            **patch_metrics,
            "min_aggrescan3d_score": values.get("min_value"),
            "max_aggrescan3d_score": values.get("max_value"),
            "total_aggrescan3d_score": values.get("total_value"),
            "mean_aggrescan3d_score": values.get("avg_value"),
        })
    return rows


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for field in fields:
                value = row.get(field, "")
                if isinstance(value, float):
                    value = f"{value:.6f}"
                out[field] = value
            writer.writerow(out)


def write_readme(out_dir: Path, args: argparse.Namespace) -> None:
    ph_line = f"- pH-specific Aggrescan3D scale: {args.ph}\n" if args.ph is not None else "- pH-specific scale: not used\n"
    text = f"""# Aggrescan3D Aggregability Scores

Generated by `scripts/structural_evo_analysis/structure_score_aggrescan3d.py`.

## Method

The script calls the installed Aggrescan3D CLI at `{args.aggrescan_bin}` and collects:

- `global_scores.tsv`: global/chain-level min, max, total, and mean A3D scores from `A3D_summary.json`.
- `per_residue_scores.tsv`: per-residue A3D scores from `A3D.csv`.
- `scored_pdbs/*.pdb`: Aggrescan3D `output.pdb` files, with A3D residue scores in the B-factor column.
- `raw/<structure_id>/`: raw Aggrescan3D outputs and logs.

Parameters:

- distance cutoff: {args.distance}
{ph_line}- chain filter: {args.chain or "none"}
"""
    (out_dir / "README.md").write_text(text)


def score_one_structure(payload: tuple[argparse.Namespace, Path]) -> tuple[list[dict], list[dict], dict | None]:
    args, pdb = payload
    structure_id = structure_id_from_path(pdb)
    work_dir = args.out_dir / "raw" / structure_id
    scored_pdb_dir = args.out_dir / "scored_pdbs"
    try:
        run_aggrescan(args, pdb, work_dir)
        structure_residue_rows = parse_residue_scores(structure_id, pdb, work_dir / "A3D.csv")
        structure_global_rows = parse_global_scores(
            structure_id,
            pdb,
            work_dir / "A3D_summary.json",
            structure_residue_rows,
            work_dir / "output.pdb",
            args.patch_radius,
        )
        output_pdb = work_dir / "output.pdb"
        if output_pdb.exists():
            shutil.copy2(output_pdb, scored_pdb_dir / f"{structure_id}.aggrescan3d_bfactor.pdb")
        return structure_residue_rows, structure_global_rows, None
    except Exception as exc:
        return [], [], {
            "structure_id": structure_id,
            "source_path": str(pdb),
            "error": str(exc),
        }


def main() -> None:
    args = parse_args()
    if not args.aggrescan_bin.exists():
        raise FileNotFoundError(
            f"Aggrescan3D binary not found at {args.aggrescan_bin}. "
            "Create it with `bash setup_envs.sh --with-aggrescan3d`."
        )

    raw_dir = args.out_dir / "raw"
    scored_pdb_dir = args.out_dir / "scored_pdbs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    scored_pdb_dir.mkdir(parents=True, exist_ok=True)

    structures = input_structures(args)
    residue_rows = []
    global_rows = []
    failures = []
    payloads = [(args, pdb) for pdb in structures]
    if args.workers > 1 and len(payloads) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            results = executor.map(score_one_structure, payloads)
            for structure_residue_rows, structure_global_rows, failure in results:
                residue_rows.extend(structure_residue_rows)
                global_rows.extend(structure_global_rows)
                if failure:
                    failures.append(failure)
    else:
        for payload in payloads:
            structure_residue_rows, structure_global_rows, failure = score_one_structure(payload)
            residue_rows.extend(structure_residue_rows)
            global_rows.extend(structure_global_rows)
            if failure:
                failures.append(failure)

    global_fields = [
        "structure_id", "source_path", "chain", "residue_count",
        "positive_aggrescan3d_burden", "positive_aggrescan3d_residues",
        "positive_patch_count", "max_positive_patch_size", "max_positive_patch_sum",
        "max_positive_patch_mean", "max_positive_patch_residues", "min_aggrescan3d_score",
        "max_aggrescan3d_score", "total_aggrescan3d_score", "mean_aggrescan3d_score",
    ]
    residue_fields = [
        "structure_id", "source_path", "protein_state", "chain", "resseq",
        "residue_id", "one_letter", "aggrescan3d_score",
    ]
    write_tsv(args.out_dir / "global_scores.tsv", global_rows, global_fields)
    write_tsv(args.out_dir / "per_residue_scores.tsv", residue_rows, residue_fields)
    if failures:
        write_tsv(args.out_dir / "failures.tsv", failures, ["structure_id", "source_path", "error"])
    else:
        stale_failures = args.out_dir / "failures.tsv"
        if stale_failures.exists():
            stale_failures.unlink()
    write_readme(args.out_dir, args)
    if failures:
        raise SystemExit(f"Aggrescan3D failed for {len(failures)} structure(s); see {args.out_dir / 'failures.tsv'}")


if __name__ == "__main__":
    main()
