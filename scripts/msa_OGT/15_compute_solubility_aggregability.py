#!/usr/bin/env python3
"""Step 15: compute and merge solubility/Aggrescan3D structure scores."""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AFDB_DIR = ROOT / "structures" / "afdb"
DEFAULT_ANKROS_PDB = (
    ROOT
    / "results"
    / "structural"
    / "docked_holo"
    / "ankros_fad_fmn_donor_holo.pdb"
)
DEFAULT_OUT = ROOT / "results" / "msa_OGT" / "structure_scores"
DEFAULT_CAMSOL_OUT = DEFAULT_OUT / "camsol"
DEFAULT_AGGRESCAN3D_OUT = DEFAULT_OUT / "aggrescan3d"
DEFAULT_AGGRESCAN = ROOT / "envs" / "aggrescan3d" / "bin" / "aggrescan"
DEFAULT_CAMSOL_SCRIPT = ROOT / "scripts" / "msa_OGT" / "structure_score_camsol.py"
DEFAULT_AGGRESCAN_SCRIPT = ROOT / "scripts" / "msa_OGT" / "structure_score_aggrescan3d.py"
DEFAULT_STRUCTURE_VALIDATION = (
    ROOT / "results" / "msa_OGT" / "structure_validation" / "quick_classI_mismatch_review.tsv"
)
DEFAULT_VALIDATION_STATUSES = "pass_quick_classI_structural_match,review_borderline_not_mismatch"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the local CamSol-style and Aggrescan3D scorers, then merge their "
            "per-residue/global outputs for the MSA HTML step."
        )
    )
    parser.add_argument("--afdb-dir", type=Path, default=DEFAULT_AFDB_DIR)
    parser.add_argument("--ankros-pdb", type=Path, default=DEFAULT_ANKROS_PDB)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--camsol-out-dir", type=Path, default=DEFAULT_CAMSOL_OUT)
    parser.add_argument("--aggrescan3d-out-dir", type=Path, default=DEFAULT_AGGRESCAN3D_OUT)
    parser.add_argument("--camsol-script", type=Path, default=DEFAULT_CAMSOL_SCRIPT)
    parser.add_argument("--aggrescan3d-script", type=Path, default=DEFAULT_AGGRESCAN_SCRIPT)
    parser.add_argument("--aggrescan-bin", type=Path, default=DEFAULT_AGGRESCAN)
    parser.add_argument("--pattern", default="*.pdb")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--sasa-points", type=int, default=96)
    parser.add_argument("--aggrescan-distance", type=int, default=10)
    parser.add_argument("--patch-radius", type=float, default=6.0)
    parser.add_argument("--workers", type=int, default=1, help="Parallel worker processes for structure scorers.")
    parser.add_argument("--structure-validation", type=Path, default=DEFAULT_STRUCTURE_VALIDATION)
    parser.add_argument(
        "--validation-statuses",
        default=DEFAULT_VALIDATION_STATUSES,
        help="Comma-separated validation statuses to include from the structural validation TSV.",
    )
    parser.add_argument(
        "--skip-structure-validation-filter",
        action="store_true",
        help="Score all AFDB structures even if the structural validation TSV exists.",
    )
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--skip-compute", action="store_true", help="Only merge existing scorer outputs.")
    parser.add_argument("--skip-camsol", action="store_true")
    parser.add_argument("--skip-aggrescan3d", action="store_true")
    parser.add_argument("--force-aggrescan3d", action="store_true")
    parser.add_argument("--no-pdb", action="store_true", help="Do not write CamSol scored PDBs.")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run_logged(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        cmd,
        cwd=str(ROOT),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(
        "COMMAND: " + " ".join(cmd) + "\n\n" + completed.stdout + f"\nEXIT_STATUS={completed.returncode}\n"
    )
    if completed.returncode != 0:
        raise SystemExit(f"command failed; see {log_path}")


def maybe_run_scorers(args: argparse.Namespace) -> None:
    if args.skip_compute:
        return
    validation_args = []
    if not args.skip_structure_validation_filter:
        if not args.structure_validation.exists():
            raise SystemExit(
                f"missing structural validation table: {args.structure_validation}\n"
                "Run scripts/msa_OGT/14_validate_structural_matches.py first, "
                "or use --skip-structure-validation-filter."
            )
        validation_args = [
            "--structure-validation",
            str(args.structure_validation),
            "--validation-statuses",
            args.validation_statuses,
        ]
    if not args.skip_camsol:
        cmd = [
            str(args.python),
            str(args.camsol_script),
            "--afdb-dir",
            str(args.afdb_dir),
            "--ankros-pdb",
            str(args.ankros_pdb),
            "--out-dir",
            str(args.camsol_out_dir),
            "--pattern",
            args.pattern,
            "--sasa-points",
            str(args.sasa_points),
            "--workers",
            str(args.workers),
            *validation_args,
        ]
        if args.limit is not None:
            cmd.extend(["--limit", str(args.limit)])
        if args.no_pdb:
            cmd.append("--no-pdb")
        run_logged(cmd, args.out_dir / "logs" / "camsol.log")

    if not args.skip_aggrescan3d:
        cmd = [
            str(args.python),
            str(args.aggrescan3d_script),
            "--afdb-dir",
            str(args.afdb_dir),
            "--ankros-pdb",
            str(args.ankros_pdb),
            "--out-dir",
            str(args.aggrescan3d_out_dir),
            "--aggrescan-bin",
            str(args.aggrescan_bin),
            "--pattern",
            args.pattern,
            "--distance",
            str(args.aggrescan_distance),
            "--patch-radius",
            str(args.patch_radius),
            "--workers",
            str(args.workers),
            *validation_args,
        ]
        if args.limit is not None:
            cmd.extend(["--limit", str(args.limit)])
        if args.force_aggrescan3d:
            cmd.append("--force")
        run_logged(cmd, args.out_dir / "logs" / "aggrescan3d.log")


def residue_key(row: dict[str, str]) -> tuple[str, str, int, str]:
    return (
        row.get("structure_id", ""),
        row.get("chain", ""),
        int(row.get("resseq", "")),
        row.get("icode", ""),
    )


def merge_per_residue(camsol_path: Path | None, aggrescan_path: Path | None) -> list[dict[str, object]]:
    rows: dict[tuple[str, str, int, str], dict[str, object]] = {}
    if camsol_path is not None:
        for row in read_tsv(camsol_path):
            key = residue_key(row)
            rows[key] = {
                "structure_id": row.get("structure_id", ""),
                "source_path": row.get("source_path", ""),
                "chain": row.get("chain", ""),
                "resseq": row.get("resseq", ""),
                "icode": row.get("icode", ""),
                "residue_id": row.get("residue_id", ""),
                "one_letter": row.get("one_letter", ""),
                "camsol_intrinsic_solubility_score": row.get("intrinsic_solubility_score", ""),
                "camsol_structural_solubility_score": row.get("structural_solubility_score", ""),
            }

    if aggrescan_path is not None:
        for row in read_tsv(aggrescan_path):
            key = (row.get("structure_id", ""), row.get("chain", ""), int(row.get("resseq", "")), "")
            merged = rows.setdefault(
                key,
                {
                    "structure_id": row.get("structure_id", ""),
                    "source_path": row.get("source_path", ""),
                    "chain": row.get("chain", ""),
                    "resseq": row.get("resseq", ""),
                    "icode": "",
                    "residue_id": row.get("residue_id", ""),
                    "one_letter": row.get("one_letter", ""),
                },
            )
            score = row.get("aggrescan3d_score", "")
            merged["aggrescan3d_score"] = score
            try:
                merged["aggrescan3d_positive"] = f"{max(0.0, float(score)):.6f}"
            except ValueError:
                merged["aggrescan3d_positive"] = ""

    return [rows[key] for key in sorted(rows, key=lambda item: (item[0], item[1], item[2], item[3]))]


def merge_global(camsol_path: Path | None, aggrescan_path: Path | None) -> list[dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    if camsol_path is not None:
        for row in read_tsv(camsol_path):
            sid = row.get("structure_id", "")
            rows[sid] = {
                "structure_id": sid,
                "source_path": row.get("source_path", ""),
                "residue_count": row.get("residue_count", ""),
                "global_intrinsic_solubility_score": row.get("global_intrinsic_solubility_score", ""),
                "global_structural_solubility_score": row.get("global_structural_solubility_score", ""),
                "v2_intrinsic_solubility_score": row.get("v2_intrinsic_solubility_score", ""),
                "v2_structural_solubility_score": row.get("v2_structural_solubility_score", ""),
            }

    if aggrescan_path is not None:
        aggrescan_rows = read_tsv(aggrescan_path)
        preferred = [row for row in aggrescan_rows if row.get("chain") == "All"]
        if not preferred:
            preferred = aggrescan_rows
        for row in preferred:
            sid = row.get("structure_id", "")
            merged = rows.setdefault(sid, {"structure_id": sid, "source_path": row.get("source_path", "")})
            for field in [
                "positive_aggrescan3d_burden",
                "positive_aggrescan3d_residues",
                "positive_patch_count",
                "max_positive_patch_size",
                "max_positive_patch_sum",
                "max_positive_patch_mean",
                "max_positive_patch_residues",
                "min_aggrescan3d_score",
                "max_aggrescan3d_score",
                "mean_aggrescan3d_score",
            ]:
                merged[field] = row.get(field, "")
    return [rows[sid] for sid in sorted(rows)]


def validate_inputs(args: argparse.Namespace) -> None:
    missing = []
    if not args.skip_camsol and not args.camsol_out_dir.joinpath("per_residue_scores.tsv").exists():
        missing.append(args.camsol_out_dir / "per_residue_scores.tsv")
    if not args.skip_camsol and not args.camsol_out_dir.joinpath("global_scores.tsv").exists():
        missing.append(args.camsol_out_dir / "global_scores.tsv")
    if not args.skip_aggrescan3d and not args.aggrescan3d_out_dir.joinpath("per_residue_scores.tsv").exists():
        missing.append(args.aggrescan3d_out_dir / "per_residue_scores.tsv")
    if not args.skip_aggrescan3d and not args.aggrescan3d_out_dir.joinpath("global_scores.tsv").exists():
        missing.append(args.aggrescan3d_out_dir / "global_scores.tsv")
    if missing:
        raise SystemExit("missing scorer output(s):\n" + "\n".join(str(path) for path in missing))


def main() -> None:
    args = parse_args()
    maybe_run_scorers(args)
    validate_inputs(args)

    camsol_per_residue = None if args.skip_camsol else args.camsol_out_dir / "per_residue_scores.tsv"
    aggrescan_per_residue = None if args.skip_aggrescan3d else args.aggrescan3d_out_dir / "per_residue_scores.tsv"
    per_residue_rows = merge_per_residue(camsol_per_residue, aggrescan_per_residue)
    per_residue_fields = [
        "structure_id",
        "source_path",
        "chain",
        "resseq",
        "icode",
        "residue_id",
        "one_letter",
        "camsol_intrinsic_solubility_score",
        "camsol_structural_solubility_score",
        "aggrescan3d_score",
        "aggrescan3d_positive",
    ]
    write_tsv(args.out_dir / "per_residue_scores.tsv", per_residue_rows, per_residue_fields)

    camsol_global = None if args.skip_camsol else args.camsol_out_dir / "global_scores.tsv"
    aggrescan_global = None if args.skip_aggrescan3d else args.aggrescan3d_out_dir / "global_scores.tsv"
    global_rows = merge_global(camsol_global, aggrescan_global)
    global_fields = [
        "structure_id",
        "source_path",
        "residue_count",
        "global_intrinsic_solubility_score",
        "global_structural_solubility_score",
        "v2_intrinsic_solubility_score",
        "v2_structural_solubility_score",
        "positive_aggrescan3d_burden",
        "positive_aggrescan3d_residues",
        "positive_patch_count",
        "max_positive_patch_size",
        "max_positive_patch_sum",
        "max_positive_patch_mean",
        "max_positive_patch_residues",
        "min_aggrescan3d_score",
        "max_aggrescan3d_score",
        "mean_aggrescan3d_score",
    ]
    write_tsv(args.out_dir / "global_scores.tsv", global_rows, global_fields)
    print(f"Saved: {args.out_dir / 'per_residue_scores.tsv'}")
    print(f"Saved: {args.out_dir / 'global_scores.tsv'}")


if __name__ == "__main__":
    main()
