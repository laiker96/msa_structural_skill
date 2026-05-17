#!/usr/bin/env python3
"""Legacy wrapper for structural-evolution CamSol-style scoring."""
from __future__ import annotations

import runpy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "scripts" / "structural_evo_analysis" / "structure_score_camsol.py"
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
    add_default(("--afdb-dir",), ROOT / "structures" / "afdb")
    add_default(("--query-pdb", "--ankros-pdb"), DEFAULT_QUERY_PDB)
    add_default(("--out-dir",), ROOT / "results" / "solubility" / "camsol")
    sys.argv[0] = str(TARGET)
    runpy.run_path(str(TARGET), run_name="__main__")


if __name__ == "__main__":
    main()
