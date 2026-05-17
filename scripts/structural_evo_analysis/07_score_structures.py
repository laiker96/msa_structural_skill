#!/usr/bin/env python3
"""Step 07: run CamSol-style and Aggrescan3D scoring on query/AFDB structures."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--afdb-dir", type=Path, default=cfg.STRUCTURE_DIR / "afdb")
    parser.add_argument(
        "--query-pdb",
        type=Path,
        default=cfg.STRUCTURE_DIR / "query.pdb",
        help="Required local query PDB for vulnerability analysis.",
    )
    parser.add_argument("--out-dir", type=Path, default=cfg.OUTPUT_DIR / "structure_scores")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--pattern", default="*.pdb")
    parser.add_argument("--sasa-points", type=int, default=96)
    parser.add_argument("--aggrescan-distance", type=int, default=10)
    parser.add_argument("--patch-radius", type=float, default=6.0)
    parser.add_argument("--aggrescan-bin", type=Path, default=cfg.PROJECT_ROOT / "envs" / "aggrescan3d" / "bin" / "aggrescan")
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--skip-compute", action="store_true")
    parser.add_argument("--skip-camsol", action="store_true")
    parser.add_argument("--skip-aggrescan3d", action="store_true")
    parser.add_argument("--force-aggrescan3d", action="store_true")
    parser.add_argument("--no-pdb", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.query_pdb.exists():
        raise SystemExit(f"ERROR: query PDB is required for vulnerability analysis: {args.query_pdb}")
    compute_script = Path(__file__).resolve().parent / "structure_score_merge.py"
    camsol_out = args.out_dir / "camsol"
    aggrescan_out = args.out_dir / "aggrescan3d"
    cmd = [
        str(args.python),
        str(compute_script),
        "--afdb-dir", str(args.afdb_dir),
        "--query-pdb", str(args.query_pdb),
        "--out-dir", str(args.out_dir),
        "--camsol-out-dir", str(camsol_out),
        "--aggrescan3d-out-dir", str(aggrescan_out),
        "--aggrescan-bin", str(args.aggrescan_bin),
        "--pattern", args.pattern,
        "--sasa-points", str(args.sasa_points),
        "--aggrescan-distance", str(args.aggrescan_distance),
        "--patch-radius", str(args.patch_radius),
        "--workers", str(args.workers),
        "--skip-structure-validation-filter",
    ]
    if args.limit is not None:
        cmd += ["--limit", str(args.limit)]
    if args.skip_compute:
        cmd.append("--skip-compute")
    if args.skip_camsol:
        cmd.append("--skip-camsol")
    if args.skip_aggrescan3d:
        cmd.append("--skip-aggrescan3d")
    if args.force_aggrescan3d:
        cmd.append("--force-aggrescan3d")
    if args.no_pdb:
        cmd.append("--no-pdb")
    print("Running: " + " ".join(cmd))
    completed = subprocess.run(cmd, cwd=str(cfg.PROJECT_ROOT))
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    print(f"Saved merged scores under: {args.out_dir}")


if __name__ == "__main__":
    main()
