#!/usr/bin/env python3
"""Legacy wrapper for structural-evolution structure score merging."""
from __future__ import annotations

import runpy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "scripts" / "structural_evo_analysis" / "structure_score_merge.py"
DEFAULT_QUERY_PDB = (
    ROOT
    / "results"
    / "structural"
    / "docked_holo"
    / "ankros_fad_fmn_donor_holo.pdb"
)


def has_option(names: tuple[str, ...]) -> bool:
    return any(arg == name or arg.startswith(f"{name}=") for arg in sys.argv[1:] for name in names)


def add_default(names: tuple[str, ...], value: Path) -> None:
    if not has_option(names):
        sys.argv.extend([names[0], str(value)])


def main() -> None:
    default_out = ROOT / "results" / "msa_OGT" / "structure_scores"
    add_default(("--afdb-dir",), ROOT / "structures" / "afdb")
    add_default(("--query-pdb", "--ankros-pdb"), DEFAULT_QUERY_PDB)
    add_default(("--out-dir",), default_out)
    add_default(("--camsol-out-dir",), default_out / "camsol")
    add_default(("--aggrescan3d-out-dir",), default_out / "aggrescan3d")
    add_default(
        ("--camsol-script",),
        ROOT / "scripts" / "structural_evo_analysis" / "structure_score_camsol.py",
    )
    add_default(
        ("--aggrescan3d-script",),
        ROOT / "scripts" / "structural_evo_analysis" / "structure_score_aggrescan3d.py",
    )
    add_default(("--aggrescan-bin",), ROOT / "envs" / "aggrescan3d" / "bin" / "aggrescan")
    add_default(
        ("--structure-validation",),
        ROOT / "results" / "msa_OGT" / "structure_validation" / "quick_classI_mismatch_review.tsv",
    )
    sys.argv[0] = str(TARGET)
    runpy.run_path(str(TARGET), run_name="__main__")


if __name__ == "__main__":
    main()
