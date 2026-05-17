#!/usr/bin/env python3
"""Structural-evolution CamSol-style solubility scoring.

This is a transparent local implementation of the public CamSol workflow
described by Sormanni, Aprile and Vendruscolo (JMB 2015) and the CamSol v2
global-score definition described by Sormanni et al. (Sci Rep 2017).

This file is owned by the structural-evolution pipeline so structure scoring
does not depend on the legacy OGT script directory.

Important limitation: the local papers refer to fitted coefficient tables in
supplementary material and to a web-server implementation. Those proprietary
parameters are not present in this repository, so this module follows the
published equations and uses documented physicochemical scales/weights rather
than claiming identity with the closed CamSol server.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBIO, PDBParser
from Bio.PDB.SASA import ShrakeRupley


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

DEFAULT_AFDB_DIR = cfg.STRUCTURE_DIR / "afdb"
DEFAULT_QUERY_PDB = cfg.STRUCTURE_DIR / "query.pdb"
DEFAULT_OUT_DIR = cfg.OUTPUT_DIR / "structure_scores" / "camsol"
DEFAULT_VALIDATION_STATUSES = {
    "pass_quick_classI_structural_match",
    "review_borderline_not_mismatch",
}

WINDOW = 7
PATCH_RADIUS_A = 8.0
SEQUENCE_EXCLUSION = 3
BACKGROUND_SEED = 1
BACKGROUND_SEQUENCES = 100
BACKGROUND_LENGTH = 500

ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
THREE_TO_ONE = {three: one for one, three in ONE_TO_THREE.items()}

# Tien et al./Miller-style maximum residue ASA values, also used elsewhere in
# this repository for relative SASA estimates.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 223.0, "GLU": 225.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# Scales are normalized in code before entering the linear score. Kyte-Doolittle
# is used with inverted sign because CamSol reports larger values as more
# soluble; Chou-Fasman secondary-structure propensities are used as public
# substitutes for the PDB-derived propensities referenced by the paper.
KYTE_DOOLITTLE = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
HELIX_PROPENSITY = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13,
    "G": 0.57, "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21,
    "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11, "R": 0.98,
    "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69,
}
BETA_PROPENSITY = {
    "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38,
    "G": 0.75, "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30,
    "M": 1.05, "N": 0.89, "P": 0.55, "Q": 1.10, "R": 0.93,
    "S": 0.75, "T": 1.19, "V": 1.70, "W": 1.37, "Y": 1.47,
}

# Local weights chosen to match the CamSol qualitative terms: hydrophobicity and
# beta propensity lower solubility, charge/gatekeeping and polar residues raise
# it, and helix propensity has a weak soluble contribution.
WEIGHTS = {
    "hydrophobicity": 0.95,
    "charge": 0.75,
    "helix": 0.12,
    "beta": 0.25,
    "pattern": 0.18,
    "gatekeeping": 0.18,
}


@dataclass
class ResidueRow:
    structure_id: str
    source_path: Path
    chain: str
    resseq: int
    icode: str
    resname: str
    one: str
    residue: object
    atoms: tuple
    asa: float
    rsa: float
    intrinsic_raw: float = 0.0
    intrinsic_score: float = 0.0
    exposed_intrinsic_score: float = 0.0
    exposure_weight: float = 0.0
    patch_score: float = 0.0
    structural_score: float = 0.0

    @property
    def residue_id(self) -> str:
        suffix = self.icode.strip()
        return f"{self.chain}:{self.resseq}{suffix}" if suffix else f"{self.chain}:{self.resseq}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute CamSol-style global and per-residue solubility scores."
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
    parser.add_argument("--sasa-points", type=int, default=96)
    parser.add_argument("--global-threshold", type=float, default=0.7)
    parser.add_argument("--v2-threshold", type=float, default=1.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--pattern", default="*.pdb")
    parser.add_argument("--no-pdb", action="store_true", help="Do not write B-factor scored PDBs.")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker processes for per-structure scoring.",
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


def charge(one: str) -> int:
    if one in {"D", "E"}:
        return -1
    if one in {"K", "R"}:
        return 1
    return 0


def normalized_hydrophobicity(one: str) -> float:
    return KYTE_DOOLITTLE[one] / 4.5


def normalized_propensity(value: float) -> float:
    return value - 1.0


def residue_raw_score(one: str) -> float:
    return (
        -WEIGHTS["hydrophobicity"] * normalized_hydrophobicity(one)
        + WEIGHTS["charge"] * abs(charge(one))
        + WEIGHTS["helix"] * normalized_propensity(HELIX_PROPENSITY[one])
        - WEIGHTS["beta"] * normalized_propensity(BETA_PROPENSITY[one])
    )


def alternating_pattern_score(seq: str, center: int) -> float:
    """Public approximation of CamSol/Zyggregator pattern term.

    Alternating hydrophobic/hydrophilic residues are made less damaging, while
    contiguous hydrophobic patterning is penalized.
    """
    start = max(0, center - 3)
    end = min(len(seq), center + 4)
    frag = seq[start:end]
    if len(frag) < 3:
        return 0.0
    hydrophobic = [KYTE_DOOLITTLE[one] > 1.0 for one in frag]
    transitions = sum(a != b for a, b in zip(hydrophobic, hydrophobic[1:]))
    hydrophobic_pairs = sum(a and b for a, b in zip(hydrophobic, hydrophobic[1:]))
    return (transitions / (len(frag) - 1)) - (hydrophobic_pairs / (len(frag) - 1))


def gatekeeping_score(seq: str, center: int) -> float:
    """Sequence gatekeeping term from same-sign nearby charges."""
    score = 0.0
    for offset in range(-5, 6):
        idx = center + offset
        if idx < 0 or idx >= len(seq):
            continue
        c = charge(seq[idx])
        if c == 0:
            continue
        distance = abs(offset)
        score += math.exp(-0.45 * distance) * abs(c)
    return score


def intrinsic_profile(seq: str) -> list[float]:
    seq = "".join(one for one in seq.upper() if one in KYTE_DOOLITTLE)
    return zscore_profile(raw_intrinsic_profile(seq))


def raw_intrinsic_profile(seq: str) -> list[float]:
    raw = [residue_raw_score(one) for one in seq]
    profile = []
    for i in range(len(seq)):
        start = max(0, i - 3)
        end = min(len(seq), i + 4)
        window_mean = sum(raw[start:end]) / (end - start)
        profile.append(
            window_mean
            + WEIGHTS["pattern"] * alternating_pattern_score(seq, i)
            + WEIGHTS["gatekeeping"] * gatekeeping_score(seq, i)
        )
    return profile


def zscore_profile(values: list[float]) -> list[float]:
    if not values:
        return []
    mean, sd = background_normalization()
    if sd == 0.0:
        return [0.0 for _ in values]
    return [(value - mean) / sd for value in values]


_BACKGROUND_NORMALIZATION: tuple[float, float] | None = None


def background_normalization() -> tuple[float, float]:
    """Fixed random-polypeptide normalization, analogous to CamSol scaling."""
    global _BACKGROUND_NORMALIZATION
    if _BACKGROUND_NORMALIZATION is not None:
        return _BACKGROUND_NORMALIZATION

    rng = random.Random(BACKGROUND_SEED)
    alphabet = tuple(ONE_TO_THREE)
    values = []
    for _ in range(BACKGROUND_SEQUENCES):
        seq = "".join(rng.choice(alphabet) for _ in range(BACKGROUND_LENGTH))
        values.extend(raw_intrinsic_profile(seq))
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    _BACKGROUND_NORMALIZATION = (mean, math.sqrt(variance))
    return _BACKGROUND_NORMALIZATION


def exposure_weight(rsa: float) -> float:
    x = max(0.0, rsa / 100.0)
    if x < 0.05:
        return 0.0
    return 1.0 / (1.0 + math.exp(-10.0 * (x - 0.3)))


def exposure_rescaled(rsa: float) -> float:
    return 0.25 + 0.75 * max(0.0, min(rsa / 100.0, 1.0))


def ca_atom(row: ResidueRow):
    if "CA" in row.residue:
        return row.residue["CA"]
    for atom in row.atoms:
        return atom
    return None


def load_residues(pdb: Path, sasa_points: int) -> list[ResidueRow]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb.stem, str(pdb))
    ShrakeRupley(n_points=sasa_points).compute(structure, level="R")
    structure_id = structure_id_from_path(pdb)

    rows = []
    for residue in structure.get_residues():
        hetflag, resseq, icode = residue.id
        if hetflag.strip():
            continue
        resname = residue.get_resname().upper()
        one = THREE_TO_ONE.get(resname)
        max_asa = MAX_ASA.get(resname)
        if one is None or max_asa is None:
            continue
        asa = float(getattr(residue, "sasa", 0.0))
        rows.append(
            ResidueRow(
                structure_id=structure_id,
                source_path=pdb,
                chain=residue.get_parent().id,
                resseq=resseq,
                icode=icode,
                resname=resname,
                one=one,
                residue=residue,
                atoms=tuple(residue.get_atoms()),
                asa=asa,
                rsa=100.0 * asa / max_asa,
            )
        )
    return rows


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


def add_intrinsic_scores(rows: list[ResidueRow]) -> None:
    by_chain = defaultdict(list)
    for row in rows:
        by_chain[row.chain].append(row)

    for chain_rows in by_chain.values():
        seq = "".join(row.one for row in chain_rows)
        profile = intrinsic_profile(seq)
        for row, score in zip(chain_rows, profile):
            row.intrinsic_score = score
        for i, row in enumerate(chain_rows):
            start = max(0, i - 3)
            end = min(len(chain_rows), i + 4)
            weights = [exposure_rescaled(r.rsa) for r in chain_rows[start:end]]
            denom = sum(weights)
            row.exposed_intrinsic_score = sum(
                weight * r.intrinsic_score for weight, r in zip(weights, chain_rows[start:end])
            ) / denom


def add_structural_scores(rows: list[ResidueRow]) -> None:
    for row in rows:
        row.exposure_weight = exposure_weight(row.rsa)

    ca_atoms = []
    owner = {}
    for idx, row in enumerate(rows):
        atom = ca_atom(row)
        if atom is None:
            continue
        ca_atoms.append(atom)
        owner[id(atom)] = idx

    neighbor_terms = [[] for _ in rows]
    if ca_atoms:
        search = NeighborSearch(ca_atoms)
        for atom_a, atom_b in search.search_all(PATCH_RADIUS_A, level="A"):
            i = owner[id(atom_a)]
            j = owner[id(atom_b)]
            row_i = rows[i]
            row_j = rows[j]
            if row_i.chain != row_j.chain:
                continue
            if abs(row_i.resseq - row_j.resseq) <= SEQUENCE_EXCLUSION:
                continue
            distance = float(atom_a - atom_b)
            distance_weight = max(1.0 - distance / PATCH_RADIUS_A, 0.0)
            if row_j.exposure_weight > 0.0:
                neighbor_terms[i].append(
                    (distance_weight * row_j.exposure_weight, row_j.exposed_intrinsic_score)
                )
            if row_i.exposure_weight > 0.0:
                neighbor_terms[j].append(
                    (distance_weight * row_i.exposure_weight, row_i.exposed_intrinsic_score)
                )

    for i, row in enumerate(rows):
        weighted = neighbor_terms[i]
        if weighted:
            patch = sum(weight * score for weight, score in weighted) / sum(
                weight for weight, _ in weighted
            )
        else:
            patch = 0.0
        row.patch_score = patch
        row.structural_score = row.exposure_weight * (row.exposed_intrinsic_score + patch)


def global_camsol_v1_score(scores: list[float], threshold: float = 0.7) -> float:
    if not scores:
        return float("nan")
    kept = [score for score in scores if score < -threshold or score > threshold]
    return sum(kept) / len(scores)


def global_camsol_v2_tail_score(scores: list[float], threshold: float = 1.0) -> dict[str, float | int]:
    if not scores:
        return {
            "overall_score": float("nan"),
            "positive_tail": float("nan"),
            "negative_tail": float("nan"),
            "positive_count": 0,
            "negative_count": 0,
        }
    positive = [score - threshold for score in scores if score > threshold]
    negative = [score + threshold for score in scores if score < -threshold]
    return {
        "overall_score": (sum(positive) + sum(negative)) / len(scores),
        "positive_tail": sum(positive) / len(scores),
        "negative_tail": sum(negative) / len(scores),
        "positive_count": len(positive),
        "negative_count": len(negative),
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def score_structure(pdb: Path, sasa_points: int, global_threshold: float, v2_threshold: float) -> tuple[dict, list[ResidueRow]]:
    rows = load_residues(pdb, sasa_points)
    add_intrinsic_scores(rows)
    add_structural_scores(rows)
    intrinsic = [row.intrinsic_score for row in rows]
    structural = [row.structural_score for row in rows]
    intrinsic_v2 = global_camsol_v2_tail_score(intrinsic, v2_threshold)
    structural_v2 = global_camsol_v2_tail_score(structural, v2_threshold)
    global_row = {
        "structure_id": structure_id_from_path(pdb),
        "source_path": str(pdb),
        "residue_count": len(rows),
        "chains": ",".join(sorted({row.chain for row in rows})),
        "mean_rsa": mean([row.rsa for row in rows]),
        "mean_intrinsic_score": mean(intrinsic),
        "mean_structural_score": mean(structural),
        "min_intrinsic_score": min(intrinsic) if intrinsic else float("nan"),
        "min_structural_score": min(structural) if structural else float("nan"),
        "max_intrinsic_score": max(intrinsic) if intrinsic else float("nan"),
        "max_structural_score": max(structural) if structural else float("nan"),
        "global_intrinsic_solubility_score": global_camsol_v1_score(intrinsic, global_threshold),
        "global_structural_solubility_score": global_camsol_v1_score(structural, global_threshold),
        "v2_intrinsic_solubility_score": intrinsic_v2["overall_score"],
        "v2_structural_solubility_score": structural_v2["overall_score"],
        "v2_intrinsic_positive_tail": intrinsic_v2["positive_tail"],
        "v2_intrinsic_negative_tail": intrinsic_v2["negative_tail"],
        "v2_structural_positive_tail": structural_v2["positive_tail"],
        "v2_structural_negative_tail": structural_v2["negative_tail"],
        "v2_intrinsic_positive_residues": intrinsic_v2["positive_count"],
        "v2_intrinsic_negative_residues": intrinsic_v2["negative_count"],
        "v2_structural_positive_residues": structural_v2["positive_count"],
        "v2_structural_negative_residues": structural_v2["negative_count"],
    }
    return global_row, rows


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


def residue_dict(row: ResidueRow) -> dict:
    return {
        "structure_id": row.structure_id,
        "source_path": str(row.source_path),
        "chain": row.chain,
        "resseq": row.resseq,
        "icode": row.icode.strip(),
        "residue_id": row.residue_id,
        "resname": row.resname,
        "one_letter": row.one,
        "asa": row.asa,
        "rsa": row.rsa,
        "exposure_weight": row.exposure_weight,
        "intrinsic_solubility_score": row.intrinsic_score,
        "exposure_weighted_intrinsic_score": row.exposed_intrinsic_score,
        "surface_patch_score": row.patch_score,
        "structural_solubility_score": row.structural_score,
    }


def write_bfactor_pdb(path: Path, rows: list[ResidueRow]) -> None:
    if not rows:
        return
    for row in rows:
        for atom in row.atoms:
            atom.set_bfactor(row.structural_score)
    structure = rows[0].residue.get_parent().get_parent().get_parent()
    tmp_path = path.with_suffix(".tmp.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(tmp_path))
    remark = (
        "REMARK   1 B-FACTORS STORE CAMSOL-STYLE STRUCTURAL SOLUBILITY SCORE\n"
        "REMARK   1 POSITIVE VALUES ARE MORE SOLUBLE; NEGATIVE VALUES ARE LESS SOLUBLE\n"
        "REMARK   1 LOCAL PUBLIC-EQUATION APPROXIMATION; NOT CLOSED CAMSOL SERVER OUTPUT\n"
    )
    path.write_text(remark + tmp_path.read_text())
    tmp_path.unlink()


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


def score_structure_worker(payload: tuple[Path, int, float, float, bool, Path]) -> tuple[dict, list[dict]]:
    pdb, sasa_points, global_threshold, v2_threshold, no_pdb, scored_pdb_dir = payload
    global_row, rows = score_structure(pdb, sasa_points, global_threshold, v2_threshold)
    residue_rows = [residue_dict(row) for row in rows]
    if not no_pdb:
        write_bfactor_pdb(scored_pdb_dir / f"{global_row['structure_id']}.camsol_structural_bfactor.pdb", rows)
    return global_row, residue_rows


def write_readme(out_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# CamSol-Style Solubility Scores

Generated by `scripts/structural_evo_analysis/structure_score_camsol.py`.

## Method

- Intrinsic profile follows the public CamSol workflow in Sormanni et al. 2015: amino-acid physicochemical score, seven-residue smoothing, local hydrophobic/hydrophilic pattern correction, and charge gatekeeping.
- Structural correction follows the same paper: residues are weighted by relative SASA, projected onto surface-exposed positions, and smoothed over an 8 A C-alpha patch while excluding sequence neighbors within +/-3 residues.
- Intrinsic scores are normalized against deterministic random polypeptides generated with seed {BACKGROUND_SEED}, {BACKGROUND_SEQUENCES} sequences, and length {BACKGROUND_LENGTH}.
- `global_intrinsic_solubility_score` and `global_structural_solubility_score` use the 2015 thresholded score with threshold +/-{args.global_threshold}.
- `v2_*_solubility_score` uses the CamSol v2 tail-score structure described in Sormanni et al. 2017 with threshold +/-{args.v2_threshold}.

## Limitation

The fitted CamSol coefficient table referenced by the papers is not present in `info/`, and the official CamSol server implementation is not included. These scores are reproducible local approximations using public equations and public amino-acid scales. Treat them as rank/order and hotspot scores, not absolute experimental solubility.
"""
    (out_dir / "README.md").write_text(text)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scored_pdb_dir = args.out_dir / "scored_pdbs"
    if not args.no_pdb:
        scored_pdb_dir.mkdir(parents=True, exist_ok=True)

    structures = input_structures(args)
    global_rows = []
    residue_rows = []
    payloads = [
        (pdb, args.sasa_points, args.global_threshold, args.v2_threshold, args.no_pdb, scored_pdb_dir)
        for pdb in structures
    ]
    if args.workers > 1 and len(payloads) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for global_row, rows in executor.map(score_structure_worker, payloads):
                global_rows.append(global_row)
                residue_rows.extend(rows)
    else:
        for payload in payloads:
            global_row, rows = score_structure_worker(payload)
            global_rows.append(global_row)
            residue_rows.extend(rows)

    global_fields = [
        "structure_id", "source_path", "residue_count", "chains", "mean_rsa",
        "mean_intrinsic_score", "mean_structural_score", "min_intrinsic_score",
        "min_structural_score", "max_intrinsic_score", "max_structural_score",
        "global_intrinsic_solubility_score", "global_structural_solubility_score",
        "v2_intrinsic_solubility_score", "v2_structural_solubility_score",
        "v2_intrinsic_positive_tail", "v2_intrinsic_negative_tail",
        "v2_structural_positive_tail", "v2_structural_negative_tail",
        "v2_intrinsic_positive_residues", "v2_intrinsic_negative_residues",
        "v2_structural_positive_residues", "v2_structural_negative_residues",
    ]
    residue_fields = [
        "structure_id", "source_path", "chain", "resseq", "icode", "residue_id",
        "resname", "one_letter", "asa", "rsa", "exposure_weight",
        "intrinsic_solubility_score", "exposure_weighted_intrinsic_score",
        "surface_patch_score", "structural_solubility_score",
    ]
    write_tsv(args.out_dir / "global_scores.tsv", global_rows, global_fields)
    write_tsv(args.out_dir / "per_residue_scores.tsv", residue_rows, residue_fields)
    write_readme(args.out_dir, args)


if __name__ == "__main__":
    main()
