#!/usr/bin/env python3
"""Step 16: ANKros-positioned consensus logos for thermal regime clades."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import tempfile
import warnings
from collections import Counter
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "results" / "msa_OGT" / "figures" / "thermal_clade_consensus_logos_representatives"
os.environ.setdefault("MPLCONFIGDIR", str(DEFAULT_OUT / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner
from Bio.PDB import NeighborSearch, PDBParser, ShrakeRupley
from Bio.PDB.DSSP import DSSP
import biotite.structure as bstruc
import biotite.structure.io as bstrucio
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D


DEFAULT_ALIGNMENT = ROOT / "results" / "msa_OGT" / "repset_hmmalign_linker_refined.fa"
DEFAULT_COLMAP = ROOT / "results" / "msa_OGT" / "repset_hmmalign_linker_refined_column_map.tsv"
DEFAULT_CLADES = ROOT / "results" / "msa_OGT" / "regime_clades" / "regime_clades.tsv"
DEFAULT_CONTACTS = ROOT / "results" / "structural" / "cofactor_contacts.tsv"
DEFAULT_RESIDUE_ANNOTATIONS = (
    ROOT
    / "results"
    / "structural"
    / "docked_holo"
    / "ankros_residue_annotations.tsv"
)
DEFAULT_POCKET = ROOT / "results" / "structural" / "conservation_per_pocket.tsv"
DEFAULT_ET = ROOT / "results" / "structural" / "et_chain.tsv"
DEFAULT_AFDB_DIR = ROOT / "structures" / "afdb"
DEFAULT_ANKROS_STRUCTURE = (
    ROOT
    / "results"
    / "structural"
    / "docked_holo"
    / "ankros_fad_fmn_donor_holo.pdb"
)
DEFAULT_MERGED_RESIDUE_SCORES = ROOT / "results" / "msa_OGT" / "structure_scores" / "per_residue_scores.tsv"
DEFAULT_CAMSOL_RESIDUES = ROOT / "results" / "solubility" / "camsol" / "per_residue_scores.tsv"
DEFAULT_AGGRESCAN3D_RESIDUES = ROOT / "results" / "aggregability" / "aggrescan3d" / "per_residue_scores.tsv"
ANKROS_ID = "photoHymenobact"
ANKROS_SCORE_IDS = ("photoHymenobact", "ankros_fad_fmn_donor_holo", "ankros_cpd")

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)
MAX_BITS = math.log2(len(AA_ORDER))

DOMAINS = {
    "antenna": (1, 130),
    "linker": (131, 205),
    "catalytic": (206, 437),
}

CATEGORY_COLORS = {
    "surface": "#2F80ED",
    "core": "#555555",
    "functional": "#D1495B",
    "insert": "#9CA3AF",
    "unknown": "#C7C7C7",
}

SECONDARY_STRUCTURE_LABELS = {
    "H": "helix",
    "E": "strand",
    "C": "coil",
}

SECONDARY_STRUCTURE_COLORS = {
    "H": "#D1495B",
    "E": "#2F80ED",
    "C": "#D1D5DB",
}

REGIME_COLORS = {
    "psychro": "#2F80ED",
    "meso": "#2E7D32",
    "thermo": "#D1495B",
}

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

MAX_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}

LIGAND_RESN = {"FAD", "FMN", "MHF", "HDF", "8HDF", "TCP", "CPD", "DNA", "DA", "DT", "DG", "DC"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

# Kyte-Doolittle hydrophobicity. Higher values are more hydrophobic.
HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

# Side-chain/full residue volume scale used only for relative bar lengths.
AA_VOLUME = {
    "G": 60.1,
    "A": 88.6,
    "S": 89.0,
    "C": 108.5,
    "D": 111.1,
    "P": 112.7,
    "N": 114.1,
    "T": 116.1,
    "E": 138.4,
    "V": 140.0,
    "Q": 143.8,
    "H": 153.2,
    "M": 162.9,
    "I": 166.7,
    "L": 166.7,
    "K": 168.6,
    "R": 173.4,
    "F": 189.9,
    "Y": 193.6,
    "W": 227.8,
}

HYDRO_CMAP = LinearSegmentedColormap.from_list(
    "hydropathy",
    ["#2166AC", "#F7F7F7", "#B2182B"],
)
HYDRO_NORM = Normalize(vmin=min(HYDROPATHY.values()), vmax=max(HYDROPATHY.values()))

FONT = FontProperties(family="DejaVu Sans", weight="bold")
LETTER_WIDTHS = {
    aa: TextPath((0, 0), aa, size=1, prop=FONT).get_extents().width
    for aa in AA_ORDER
}
MAX_LETTER_WIDTH = max(LETTER_WIDTHS.values())
MIN_VOLUME = min(AA_VOLUME.values())
MAX_VOLUME = max(AA_VOLUME.values())


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    entries = []
    header = None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    entries.append((header.split()[0], header, "".join(parts)))
                header = line[1:].strip()
                parts = []
            elif line:
                parts.append(line.strip())
    if header is not None:
        entries.append((header.split()[0], header, "".join(parts)))
    return entries


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


def parse_float(text: object) -> float | None:
    try:
        if text == "":
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def format_optional_float(value: float | None, ndigits: int = 4) -> str:
    return f"{value:.{ndigits}f}" if value is not None else ""


def blend_hex(left: str, right: str, fraction: float) -> str:
    fraction = max(0.0, min(1.0, fraction))
    left_rgb = tuple(int(left[idx : idx + 2], 16) for idx in (1, 3, 5))
    right_rgb = tuple(int(right[idx : idx + 2], 16) for idx in (1, 3, 5))
    rgb = tuple(round(left_rgb[idx] + (right_rgb[idx] - left_rgb[idx]) * fraction) for idx in range(3))
    return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"


def structure_id_for_alignment_id(sid: str) -> str:
    if sid == ANKROS_ID or sid.startswith(ANKROS_ID):
        return ANKROS_SCORE_IDS[0]
    if sid.startswith("UniRef90_"):
        return sid.removeprefix("UniRef90_")
    return sid


def load_residue_score_cache(
    camsol_path: Path,
    aggrescan3d_path: Path,
    merged_path: Path | None = None,
) -> tuple[dict[str, dict[int, dict[str, float]]], list[str]]:
    scores: dict[str, dict[int, dict[str, float]]] = {}
    warnings: list[str] = []

    def alias_ankros_scores() -> None:
        for score_id in ANKROS_SCORE_IDS:
            if score_id in scores:
                ankros_scores = scores[score_id]
                break
        else:
            return
        for score_id in ANKROS_SCORE_IDS:
            scores.setdefault(score_id, ankros_scores)

    if merged_path is not None and merged_path.exists():
        for row in read_tsv(merged_path):
            structure_id = row.get("structure_id", "")
            if not structure_id or row.get("chain", "") not in {"", "A"}:
                continue
            try:
                resseq = int(row.get("resseq", ""))
            except ValueError:
                continue
            residue_scores = scores.setdefault(structure_id, {}).setdefault(resseq, {})
            intrinsic = parse_float(row.get("camsol_intrinsic_solubility_score", ""))
            structural = parse_float(row.get("camsol_structural_solubility_score", ""))
            aggrescan = parse_float(row.get("aggrescan3d_score", ""))
            if intrinsic is not None:
                residue_scores["camsol_intrinsic"] = intrinsic
            if structural is not None:
                residue_scores["camsol_structural"] = structural
            if aggrescan is not None:
                residue_scores["aggrescan3d"] = aggrescan
                residue_scores["aggrescan3d_positive"] = max(0.0, aggrescan)
        alias_ankros_scores()
        return scores, warnings

    if camsol_path.exists():
        for row in read_tsv(camsol_path):
            structure_id = row.get("structure_id", "")
            if not structure_id or row.get("chain", "") not in {"", "A"}:
                continue
            try:
                resseq = int(row.get("resseq", ""))
            except ValueError:
                continue
            residue_scores = scores.setdefault(structure_id, {}).setdefault(resseq, {})
            intrinsic = parse_float(row.get("intrinsic_solubility_score", ""))
            structural = parse_float(row.get("structural_solubility_score", ""))
            if intrinsic is not None:
                residue_scores["camsol_intrinsic"] = intrinsic
            if structural is not None:
                residue_scores["camsol_structural"] = structural
    else:
        warnings.append(f"missing CamSol residue scores: {camsol_path}")

    if aggrescan3d_path.exists():
        for row in read_tsv(aggrescan3d_path):
            structure_id = row.get("structure_id", "")
            if not structure_id or row.get("protein_state", "") != "folded":
                continue
            if row.get("chain", "") not in {"", "A"}:
                continue
            try:
                resseq = int(row.get("resseq", ""))
            except ValueError:
                continue
            score = parse_float(row.get("aggrescan3d_score", ""))
            if score is None:
                continue
            residue_scores = scores.setdefault(structure_id, {}).setdefault(resseq, {})
            residue_scores["aggrescan3d"] = score
            residue_scores["aggrescan3d_positive"] = max(0.0, score)
    else:
        warnings.append(f"missing Aggrescan3D residue scores: {aggrescan3d_path}")

    alias_ankros_scores()
    return scores, warnings


def split_tips(text: str) -> list[str]:
    return [part.strip() for part in (text or "").split(",") if part.strip()]


def parse_range(text: str) -> tuple[int, int] | None:
    text = (text or "").strip()
    if not text:
        return None
    if "-" in text:
        left, right = text.split("-", 1)
        return int(left), int(right)
    value = int(text)
    return value, value


def overlaps_domain(source_range: str, start: int, end: int) -> bool:
    parsed = parse_range(source_range)
    return bool(parsed and parsed[0] <= end and start <= parsed[1])


def find_query(entries: list[tuple[str, str, str]]) -> tuple[str, str, str]:
    for sid, header, seq in entries:
        if ANKROS_ID in sid or ANKROS_ID in header:
            return sid, header, seq
    raise SystemExit(f"ERROR: ANKros query '{ANKROS_ID}' not found in alignment")


def alignment_sequence_positions(seq: str) -> dict[int, int]:
    positions = {}
    seq_pos = 0
    for aln_col, aa in enumerate(seq):
        if aa.upper() in AA_SET:
            seq_pos += 1
            positions[aln_col] = seq_pos
    return positions


def afdb_path_for_sid(sid: str, afdb_dir: Path) -> Path | None:
    if not sid.startswith("UniRef90_"):
        return None
    accession = sid.removeprefix("UniRef90_")
    path = afdb_dir / f"AF-{accession}-F1-model_v6.pdb"
    return path if path.exists() else None


def insert_run_lengths(colmap: list[dict[str, str]]) -> dict[int, int]:
    lengths: dict[int, int] = {}
    run: list[int] = []
    for row in colmap:
        out_col = int(row["out_col"])
        if row.get("region") == "insert":
            run.append(out_col)
            continue
        if run:
            for col in run:
                lengths[col] = len(run)
            run = []
    if run:
        for col in run:
            lengths[col] = len(run)
    return lengths


def select_domain_columns(
    colmap: list[dict[str, str]],
    domain: str,
    max_insert_run: int,
) -> list[tuple[int, dict[str, str], bool, int]]:
    start, end = DOMAINS[domain]
    run_lengths = insert_run_lengths(colmap)
    selected = []
    for row in colmap:
        out_col = int(row["out_col"])
        qpos = row.get("qpos", "")
        if qpos and start <= int(qpos) <= end:
            selected.append((out_col, row, False, 0))
            continue
        if row.get("region") != "insert":
            continue
        run_len = run_lengths.get(out_col, 0)
        if run_len <= max_insert_run and overlaps_domain(row.get("source_range", ""), start, end):
            selected.append((out_col, row, True, run_len))
    return selected


def select_all_columns(
    colmap: list[dict[str, str]],
    max_insert_run: int,
) -> list[tuple[int, dict[str, str], bool, int]]:
    start = min(start for start, _end in DOMAINS.values())
    end = max(end for _start, end in DOMAINS.values())
    run_lengths = insert_run_lengths(colmap)
    selected = []
    for row in colmap:
        out_col = int(row["out_col"])
        qpos = row.get("qpos", "")
        if qpos and start <= int(qpos) <= end:
            selected.append((out_col, row, False, 0))
            continue
        if row.get("region") != "insert":
            continue
        run_len = run_lengths.get(out_col, 0)
        if run_len <= max_insert_run and overlaps_domain(row.get("source_range", ""), start, end):
            selected.append((out_col, row, True, run_len))
    return selected


def choose_representative_entry(
    selected_entries: list[tuple[str, str, str]],
    selected_cols: list[tuple[int, dict[str, str], bool, int]],
    scored_structure_ids: set[str] | None = None,
) -> tuple[str, str, str] | None:
    if not selected_entries:
        return None
    scored_structure_ids = scored_structure_ids or set()
    aln_cols = [aln_col for aln_col, _row, _is_insert, _insert_run_len in selected_cols]
    ranked = []
    for sid, header, seq in selected_entries:
        n_residues = sum(1 for aln_col in aln_cols if aln_col < len(seq) and seq[aln_col].upper() in AA_SET)
        score_priority = 0 if structure_id_for_alignment_id(sid) in scored_structure_ids else 1
        ranked.append((score_priority, -n_residues, sid, header, seq))
    _score_priority, _neg_residue_count, sid, header, seq = sorted(ranked)[0]
    return sid, header, seq


def aa_counts(entries: list[tuple[str, str, str]], aln_col: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for _sid, _header, seq in entries:
        if aln_col < len(seq):
            aa = seq[aln_col].upper()
            if aa in AA_SET:
                counts[aa] += 1
    return counts


def gap_count(entries: list[tuple[str, str, str]], aln_col: int) -> int:
    gaps = 0
    for _sid, _header, seq in entries:
        if aln_col >= len(seq) or seq[aln_col].upper() not in AA_SET:
            gaps += 1
    return gaps


def top_residue(counts: Counter[str], gaps: int) -> tuple[str, int]:
    if not counts and not gaps:
        return "-", 0
    ranked = [("-", gaps), *counts.items()]
    return sorted(
        ranked,
        key=lambda item: (-item[1], -1 if item[0] == "-" else AA_ORDER.index(item[0])),
    )[0]


def information_bits(counts: Counter[str], nseq: int) -> tuple[float, float, float]:
    n_res = sum(counts.values())
    if not nseq or not n_res:
        return 0.0, 0.0, 1.0
    occupancy = n_res / nseq
    probs = [count / n_res for count in counts.values()]
    entropy = -sum(p * math.log2(p) for p in probs)
    return (MAX_BITS - entropy) * occupancy, occupancy, 1.0 - occupancy


def parse_resnum(text: str) -> int | None:
    match = re.search(r"(\d+)", text or "")
    return int(match.group(1)) if match else None


def protein_residues_from_structure(path: Path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, path)
    ShrakeRupley().compute(structure, level="R")
    residues = []
    for model in structure:
        for chain in model:
            chain_residues = [
                residue
                for residue in chain
                if residue.id[0] == " " and THREE_TO_ONE.get(residue.get_resname().strip().upper()) in AA_SET
            ]
            if len(chain_residues) > len(residues):
                residues = chain_residues
        break
    return structure, residues


def structure_ligand_atoms(structure):
    atoms = []
    for residue in structure.get_residues():
        resname = residue.get_resname().strip().upper()
        if residue.id[0] != " " and resname in LIGAND_RESN:
            atoms.extend(atom for atom in residue.get_atoms() if atom.element != "H")
    return atoms


def residue_sidechain_contact_counts(
    residues: list,
    cutoff_a: float,
) -> tuple[dict[int, int], dict[int, str]]:
    """Count unique residue-residue side-chain heavy-atom contacts."""
    atoms = []
    owner_by_atom_id: dict[int, int] = {}
    residue_label_by_pos: dict[int, str] = {}
    contacts: dict[int, set[int]] = {idx: set() for idx in range(1, len(residues) + 1)}

    for seq_pos, residue in enumerate(residues, start=1):
        chain_id = residue.get_parent().id
        resseq = residue.id[1]
        aa = THREE_TO_ONE.get(residue.get_resname().strip().upper(), "")
        residue_label_by_pos[seq_pos] = f"{aa}{resseq}{chain_id}"
        for atom in residue.get_atoms():
            atom_name = atom.get_name().strip()
            element = (atom.element or "").strip().upper()
            if atom_name in BACKBONE_ATOMS or element == "H":
                continue
            owner_by_atom_id[id(atom)] = seq_pos
            atoms.append(atom)

    if atoms:
        search = NeighborSearch(atoms)
        for atom_a, atom_b in search.search_all(cutoff_a, level="A"):
            left = owner_by_atom_id[id(atom_a)]
            right = owner_by_atom_id[id(atom_b)]
            if left == right:
                continue
            contacts[left].add(right)
            contacts[right].add(left)

    counts = {seq_pos: len(partners) for seq_pos, partners in contacts.items()}
    partner_labels = {
        seq_pos: ";".join(residue_label_by_pos[pos] for pos in sorted(partners))
        for seq_pos, partners in contacts.items()
    }
    return counts, partner_labels


def psea_secondary_structure_from_pdb(path: Path) -> dict[tuple[str, int, str], str]:
    """Assign H/E/C with Biotite's P-SEA implementation, keyed by chain/residue."""
    atom_array = bstrucio.load_structure(str(path))
    if type(atom_array).__name__ == "AtomArrayStack":
        atom_array = atom_array[0]
    protein = atom_array[bstruc.filter_amino_acids(atom_array)]
    by_residue: dict[tuple[str, int, str], str] = {}
    state_map = {"a": "H", "b": "E", "c": "C"}
    for chain_id in sorted(set(str(chain) for chain in protein.chain_id.tolist())):
        chain_atoms = protein[protein.chain_id == chain_id]
        ca_atoms = chain_atoms[chain_atoms.atom_name == "CA"]
        if len(ca_atoms) == 0:
            continue
        try:
            states = bstruc.annotate_sse(ca_atoms)
        except Exception:
            continue
        for atom, state in zip(ca_atoms, states):
            key = (str(atom.chain_id), int(atom.res_id), str(atom.ins_code).strip())
            by_residue[key] = state_map.get(str(state), "C")
    return by_residue


def find_mkdssp_executable() -> str | None:
    env_bin = Path(sys.executable).resolve().parent / "mkdssp"
    if env_bin.exists():
        return str(env_bin)
    return None


def dssp_secondary_structure_from_pdb(path: Path) -> dict[tuple[str, int, str], str]:
    """Assign collapsed H/E/C states with mkdssp, keyed by chain/residue."""
    mkdssp = find_mkdssp_executable()
    if mkdssp is None:
        return {}
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(path.stem, path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            dssp = DSSP(structure[0], str(path), dssp=mkdssp)
    except Exception:
        try:
            with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=True) as handle:
                with path.open() as pdb_handle:
                    for line in pdb_handle:
                        if line.startswith(("ATOM  ", "HETATM", "TER", "END")):
                            handle.write(line)
                handle.flush()
                structure = parser.get_structure(path.stem, handle.name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    dssp = DSSP(structure[0], handle.name, dssp=mkdssp)
        except Exception:
            return {}
    state_map = {
        "H": "H",
        "G": "H",
        "I": "H",
        "E": "E",
        "B": "E",
        "T": "C",
        "S": "C",
        "-": "C",
    }
    by_residue: dict[tuple[str, int, str], str] = {}
    for (chain_id, residue_id), values in dssp.property_dict.items():
        hetflag, resseq, icode = residue_id
        if hetflag != " ":
            continue
        by_residue[(str(chain_id), int(resseq), str(icode).strip())] = state_map.get(str(values[2]), "C")
    return by_residue


def assign_secondary_structure(path: Path, residues: list) -> dict[int, str]:
    by_residue = dssp_secondary_structure_from_pdb(path) or psea_secondary_structure_from_pdb(path)
    secondary_structure: dict[int, str] = {}
    for seq_pos, residue in enumerate(residues, start=1):
        chain_id = residue.get_parent().id
        resseq = residue.id[1]
        icode = str(residue.id[2]).strip()
        secondary_structure[seq_pos] = by_residue.get((chain_id, resseq, icode), "C")
    return secondary_structure


def structure_features_from_pdb(
    path: Path,
    exposed_rsa_pct: float,
    ligand_cutoff_a: float,
    sidechain_contact_cutoff_a: float,
    include_ligands: bool,
) -> tuple[dict[int, dict[str, object]], dict[int, set[str]]]:
    structure, residues = protein_residues_from_structure(path)
    ligand_atoms = structure_ligand_atoms(structure) if include_ligands else []
    sidechain_counts, sidechain_partners = residue_sidechain_contact_counts(residues, sidechain_contact_cutoff_a)
    secondary_structure = assign_secondary_structure(path, residues)
    features: dict[int, dict[str, object]] = {}
    evidence: dict[int, set[str]] = {}

    for seq_pos, residue in enumerate(residues, start=1):
        aa = THREE_TO_ONE.get(residue.get_resname().strip().upper(), "")
        if aa not in AA_SET:
            continue
        sasa = float(getattr(residue, "sasa", 0.0) or 0.0)
        rsa_pct = 100.0 * sasa / MAX_ASA.get(aa, 200.0)
        category = "surface" if rsa_pct >= exposed_rsa_pct else "core"
        labels: set[str] = set()
        if ligand_atoms:
            for atom in residue.get_atoms():
                if atom.element == "H":
                    continue
                if any(atom - ligand_atom <= ligand_cutoff_a for ligand_atom in ligand_atoms):
                    category = "functional"
                    labels.add("structure_ligand_contact")
                    break
        features[seq_pos] = {
            "aa": aa,
            "category": category,
            "sasa_a2": sasa,
            "rsa_pct": rsa_pct,
            "sidechain_contact_count": sidechain_counts.get(seq_pos, 0),
            "sidechain_contact_residues": sidechain_partners.get(seq_pos, ""),
            "secondary_structure": secondary_structure.get(seq_pos, "C"),
        }
        if labels:
            evidence[seq_pos] = labels
    return features, evidence


def map_entry_positions_to_structure(
    alignment_seq: str,
    structure_features: dict[int, dict[str, object]],
) -> dict[int, int]:
    """Map ungapped alignment sequence positions to structure sequence positions."""
    entry_seq = "".join(aa.upper() for aa in alignment_seq if aa.upper() in AA_SET)
    if not entry_seq or not structure_features:
        return {}

    structure_items = sorted(structure_features.items())
    structure_seq = "".join(str(row["aa"]) for _seq_pos, row in structure_items)
    if not structure_seq:
        return {}

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(entry_seq, structure_seq)[0]

    mapping: dict[int, int] = {}
    for entry_span, structure_span in zip(alignment.aligned[0], alignment.aligned[1]):
        entry_start, entry_end = int(entry_span[0]), int(entry_span[1])
        structure_start, structure_end = int(structure_span[0]), int(structure_span[1])
        block_len = min(entry_end - entry_start, structure_end - structure_start)
        for offset in range(block_len):
            entry_pos = entry_start + offset + 1
            structure_index = structure_start + offset
            if entry_seq[entry_pos - 1] != structure_seq[structure_index]:
                continue
            structure_pos = structure_items[structure_index][0]
            mapping[entry_pos] = structure_pos
    return mapping


def load_ankros_structure_features(
    structure_path: Path,
    exposed_rsa_pct: float,
    ligand_cutoff_a: float,
    sidechain_contact_cutoff_a: float,
    extra_catalytic_sites: set[int],
    extra_evidence: dict[int, set[str]],
) -> tuple[dict[int, dict[str, object]], dict[int, set[str]]]:
    features, evidence = structure_features_from_pdb(
        structure_path,
        exposed_rsa_pct,
        ligand_cutoff_a,
        sidechain_contact_cutoff_a,
        include_ligands=True,
    )
    for qpos in extra_catalytic_sites:
        row = features.setdefault(qpos, {"aa": "", "category": "unknown", "sasa_a2": "", "rsa_pct": ""})
        row["category"] = "functional"
        evidence.setdefault(qpos, set()).update(extra_evidence.get(qpos, {"mapped_catalytic_site"}))
    for qpos, labels in evidence.items():
        if qpos in features:
            features[qpos]["functional_evidence"] = "; ".join(sorted(labels))
    return features, evidence


def load_one_afdb_feature(
    job: tuple[str, str, Path, float, float, float],
) -> tuple[str, str, dict[int, dict[str, object]], dict[int, int], str] | None:
    sid, seq, path, exposed_rsa_pct, ligand_cutoff_a, sidechain_contact_cutoff_a = job
    try:
        features, _evidence = structure_features_from_pdb(
            path,
            exposed_rsa_pct,
            ligand_cutoff_a,
            sidechain_contact_cutoff_a,
            include_ligands=False,
        )
    except Exception as exc:
        return sid, str(path), {}, {}, str(exc)
    return sid, str(path), features, map_entry_positions_to_structure(seq, features), ""


def load_afdb_feature_cache(
    entries: list[tuple[str, str, str]],
    afdb_dir: Path,
    exposed_rsa_pct: float,
    ligand_cutoff_a: float,
    sidechain_contact_cutoff_a: float,
    workers: int = 1,
) -> tuple[dict[str, dict[int, dict[str, object]]], dict[str, str], dict[str, dict[int, int]]]:
    jobs: list[tuple[str, str, Path, float, float, float]] = []
    for sid, _header, seq in entries:
        path = afdb_path_for_sid(sid, afdb_dir)
        if path is not None:
            jobs.append((sid, seq, path, exposed_rsa_pct, ligand_cutoff_a, sidechain_contact_cutoff_a))

    feature_by_sid: dict[str, dict[int, dict[str, object]]] = {}
    structure_path_by_sid: dict[str, str] = {}
    structure_pos_by_entry_pos_by_sid: dict[str, dict[int, int]] = {}
    failures = 0

    def store_result(
        result: tuple[str, str, dict[int, dict[str, object]], dict[int, int], str] | None,
    ) -> None:
        nonlocal failures
        if result is None:
            failures += 1
            return
        sid, path, features, position_map, error = result
        if error:
            failures += 1
            return
        feature_by_sid[sid] = features
        structure_path_by_sid[sid] = path
        structure_pos_by_entry_pos_by_sid[sid] = position_map

    if jobs:
        worker_count = max(1, int(workers))
        print(f"Loading AFDB structure features for {len(jobs)} entries with {worker_count} worker(s)...", flush=True)
        if worker_count == 1:
            for index, job in enumerate(jobs, start=1):
                store_result(load_one_afdb_feature(job))
                if index % 25 == 0 or index == len(jobs):
                    print(f"  processed {index}/{len(jobs)} AFDB structures", flush=True)
        else:
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                future_to_sid = {executor.submit(load_one_afdb_feature, job): job[0] for job in jobs}
                for index, future in enumerate(as_completed(future_to_sid), start=1):
                    try:
                        store_result(future.result())
                    except Exception:
                        failures += 1
                    if index % 25 == 0 or index == len(jobs):
                        print(f"  processed {index}/{len(jobs)} AFDB structures", flush=True)
        if failures:
            print(f"WARNING: skipped {failures} AFDB structures that could not be parsed or mapped", flush=True)
    return feature_by_sid, structure_path_by_sid, structure_pos_by_entry_pos_by_sid


def load_catalytic_sites(
    contacts_path: Path,
    residue_annotations_path: Path,
    pocket_path: Path,
    et_path: Path,
    contact_cutoff: float,
    et_fad_cutoff: float,
) -> tuple[set[int], dict[int, set[str]]]:
    sites: set[int] = set()
    evidence: dict[int, set[str]] = {}

    def add(resnum: int | None, label: str) -> None:
        if resnum is None:
            return
        sites.add(resnum)
        evidence.setdefault(resnum, set()).add(label)

    if residue_annotations_path.exists():
        for row in read_tsv(residue_annotations_path):
            try:
                resnum = int(row.get("resi", ""))
                dist = float(row.get("min_distance_A", ""))
            except (TypeError, ValueError):
                continue
            annotation = row.get("annotation", "").strip()
            if not annotation:
                continue
            if annotation == "fad_contacts_4A":
                label = "FAD"
            elif annotation == "fmn_antenna_pocket_contacts_4A":
                label = "FMN"
            elif annotation == "cpd_contacts_4A":
                label = "CPD"
            elif annotation == "dna_damaged_strand_D_contacts_4A":
                label = "DNA damaged strand D"
            elif annotation == "dna_complementary_strand_E_contacts_4A":
                label = "DNA complementary strand E"
            else:
                label = annotation
            add(resnum, f"{label} {dist:.3f} A")

    if contacts_path.exists():
        for row in read_tsv(contacts_path):
            try:
                dist = float(row["distance"])
            except (KeyError, TypeError, ValueError):
                continue
            if dist <= contact_cutoff:
                add(parse_resnum(row.get("resnum", "")), row.get("cofactor", "ligand").lower())

    if pocket_path.exists():
        for row in read_tsv(pocket_path):
            add(parse_resnum(row.get("ank_resi", "")), f"antenna_{row.get('candidate', '').lower()}")

    if et_path.exists():
        for row in read_tsv(et_path):
            try:
                dist = float(row["edge_to_edge_dist"])
            except (KeyError, TypeError, ValueError):
                continue
            left = row.get("residue1", "")
            right = row.get("residue2", "")
            if right == "FAD" and dist <= et_fad_cutoff:
                add(parse_resnum(left), "near_fad_et")
            if row.get("et_candidate", "").lower() == "yes":
                add(parse_resnum(left), "et_candidate")
                add(parse_resnum(right), "et_candidate")

    return sites, evidence


def majority_category(counter: Counter[str], default: str = "unknown") -> str:
    if not counter:
        return default
    priority = {"functional": 0, "surface": 1, "core": 2, "unknown": 3, "insert": 4}
    return sorted(counter.items(), key=lambda item: (-item[1], priority.get(item[0], 99)))[0][0]


def clade_structural_category(
    top_aa: str,
    qpos: str,
    selected_entries: list[tuple[str, str, str]],
    aln_col: int,
    seq_pos_by_sid: dict[str, dict[int, int]],
    afdb_features: dict[str, dict[int, dict[str, object]]],
    structure_pos_by_entry_pos_by_sid: dict[str, dict[int, int]],
    ankros_features: dict[int, dict[str, object]],
) -> tuple[str, int, int, float | None, float | None]:
    if top_aa not in AA_SET:
        return "unknown", 0, 0, None, None
    forced_functional = False
    if qpos:
        ankros_row = ankros_features.get(int(qpos), {})
        if ankros_row.get("category") == "functional":
            forced_functional = True

    votes: Counter[str] = Counter()
    n_structures = 0
    n_matching_structures = 0
    sasa_values: list[float] = []
    rsa_values: list[float] = []
    for sid, _header, seq in selected_entries:
        if aln_col >= len(seq) or seq[aln_col].upper() != top_aa:
            continue
        seq_pos = seq_pos_by_sid.get(sid, {}).get(aln_col)
        if seq_pos is None:
            continue
        structure_pos = structure_pos_by_entry_pos_by_sid.get(sid, {}).get(seq_pos)
        if structure_pos is None:
            continue
        feature = afdb_features.get(sid, {}).get(structure_pos)
        if not feature:
            continue
        n_structures += 1
        if feature.get("aa") != top_aa:
            continue
        n_matching_structures += 1
        try:
            sasa_values.append(float(feature.get("sasa_a2") or 0.0))
            rsa_values.append(float(feature.get("rsa_pct") or 0.0))
        except (TypeError, ValueError):
            pass
        category = str(feature.get("category") or "unknown")
        if category == "functional":
            category = "surface" if float(feature.get("rsa_pct") or 0.0) >= 20.0 else "core"
        votes[category] += 1
    category = "functional" if forced_functional else majority_category(votes)
    mean_sasa = sum(sasa_values) / len(sasa_values) if sasa_values else None
    mean_rsa = sum(rsa_values) / len(rsa_values) if rsa_values else None
    return category, n_structures, n_matching_structures, mean_sasa, mean_rsa


def hydropathy_color(aa: str) -> tuple[float, float, float, float]:
    value = HYDROPATHY.get(aa, 0.0)
    return HYDRO_CMAP(HYDRO_NORM(value))


def volume_fraction(aa: str) -> float:
    value = AA_VOLUME.get(aa, MIN_VOLUME)
    return (value - MIN_VOLUME) / (MAX_VOLUME - MIN_VOLUME)


def draw_letter(ax, aa: str, x: float, y: float, height: float, color: str, width_scale: float = 0.76) -> None:
    if aa not in AA_SET or height <= 0.015:
        return
    text_path = TextPath((0, 0), aa, size=1, prop=FONT)
    bbox = text_path.get_extents()
    if bbox.width == 0 or bbox.height == 0:
        return
    sx = width_scale / bbox.width
    sy = height / bbox.height
    trans = (
        Affine2D()
        .scale(sx, sy)
        .translate(x + 0.5 - (bbox.x0 + bbox.width / 2) * sx, y - bbox.y0 * sy)
        + ax.transData
    )
    patch = matplotlib.patches.PathPatch(
        MplPath(text_path.vertices, text_path.codes),
        transform=trans,
        facecolor=color,
        edgecolor="none",
        clip_on=True,
    )
    ax.add_patch(patch)


def add_size_bar(
    ax,
    aa: str,
    x: float,
    y: float,
    max_height: float,
    width: float = 0.52,
    x_center: float | None = None,
) -> None:
    center = x + 0.5 if x_center is None else x + x_center
    left = center - width / 2
    ax.add_patch(
        Rectangle(
            (left, y),
            width,
            max_height,
            facecolor="#F2F2F2",
            edgecolor="#8A8A8A",
            linewidth=0.25,
        )
    )
    if aa not in AA_SET:
        return
    height = 0.10 * max_height + 0.90 * max_height * volume_fraction(aa)
    ax.add_patch(
        Rectangle(
            (left, y),
            width,
            height,
            facecolor=hydropathy_color(aa),
            edgecolor="#333333",
            linewidth=0.35,
        )
    )


def add_ankros_track(
    ax,
    query_aa: str,
    category: str,
    x: float,
    y: float,
    height: float,
    width_scale: float = 0.62,
) -> None:
    if query_aa not in AA_SET:
        return
    color = CATEGORY_COLORS.get(category, CATEGORY_COLORS["unknown"])
    draw_letter(ax, query_aa, x, y - height / 2.0, height, color, width_scale=width_scale)


def x_ticks_for_chunk(chunk: list[dict[str, object]]) -> tuple[list[float], list[str]]:
    positions: list[float] = []
    labels: list[str] = []
    for x, stat in enumerate(chunk):
        qpos = str(stat.get("qpos", ""))
        if qpos:
            pos_int = int(qpos)
            if pos_int % 10 == 0 or x == 0 or x == len(chunk) - 1:
                positions.append(x + 0.5)
                labels.append(qpos)
        elif x == 0 or x == len(chunk) - 1 or x % 10 == 0:
            positions.append(x + 0.5)
            labels.append("+")
    return positions, labels


def add_legends(fig, axes) -> None:
    patches = [
        Patch(facecolor=CATEGORY_COLORS["surface"], label="SAS exposed"),
        Patch(facecolor=CATEGORY_COLORS["core"], label="SAS non-exposed/core"),
        Patch(facecolor=CATEGORY_COLORS["functional"], label="Functional/ligand/antenna"),
        Patch(facecolor=CATEGORY_COLORS["insert"], label="Insert/no ANKros position"),
        Patch(facecolor="#F2F2F2", edgecolor="#8A8A8A", label="AA volume bar slot"),
    ]
    fig.legend(handles=patches, loc="upper center", bbox_to_anchor=(0.56, 0.948), ncol=5, fontsize=8, frameon=False)
    cax = fig.add_axes([0.76, 0.026, 0.20, 0.016])
    sm = ScalarMappable(norm=HYDRO_NORM, cmap=HYDRO_CMAP)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("")
    cb.ax.tick_params(labelsize=6)
    fig.text(
        0.86,
        0.047,
        "vertical bar height = AA volume; bar color = hydrophobicity",
        ha="center",
        va="bottom",
        fontsize=7,
    )


def draw_clade_logo(
    domain: str,
    clade_id: str,
    regime: str,
    nseq: int,
    stats: list[dict[str, object]],
    out_png: Path,
    wrap: int,
) -> None:
    if not stats:
        return
    chunks = [stats[i : i + wrap] for i in range(0, len(stats), wrap)]
    fig_w = max(16.0, min(38.0, wrap * 0.34))
    fig_h = max(5.2, len(chunks) * 3.8)
    fig, axes = plt.subplots(len(chunks), 1, figsize=(fig_w, fig_h), squeeze=False)
    axes = axes.ravel()

    for idx, (ax, chunk) in enumerate(zip(axes, chunks), start=1):
        ax.set_xlim(0, wrap)
        ax.set_ylim(-1.55, MAX_BITS + 2.45)
        ax.axhline(0, color="#333333", linewidth=0.7)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_ylabel("bits")

        for x, stat in enumerate(chunk):
            aa = str(stat["top_aa"])
            category = str(stat["category"])
            height = float(stat["information_bits"])
            if stat["is_insert"] == "yes":
                ax.add_patch(
                    Rectangle((x, -1.55), 1.0, MAX_BITS + 3.70, facecolor="#F4F4F5", edgecolor="none", zorder=-10)
                )
            draw_letter(ax, aa, x, 0.0, height, CATEGORY_COLORS.get(category, CATEGORY_COLORS["unknown"]))
            add_size_bar(ax, aa, x, -1.34, max_height=1.08, width=0.30)
            add_size_bar(
                ax,
                str(stat.get("query_aa", "")),
                x,
                MAX_BITS + 0.64,
                max_height=1.10,
                width=0.16,
                x_center=0.84,
            )
            add_ankros_track(
                ax,
                str(stat.get("query_aa", "")),
                str(stat.get("query_category", "unknown")),
                x,
                MAX_BITS + 1.30,
                height=1.32,
                width_scale=0.36,
            )

        tick_pos, tick_lab = x_ticks_for_chunk(chunk)
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lab, fontsize=7, rotation=90)
        ax.text(-0.35, MAX_BITS + 1.30, "ANKros", ha="right", va="center", fontsize=8, clip_on=False)
        ax.text(-0.35, MAX_BITS + 0.94, "ANKros\nvolume", ha="right", va="center", fontsize=7, clip_on=False)
        ax.text(-0.35, -0.80, "Consensus\nvolume", ha="right", va="center", fontsize=7, clip_on=False)
        ax.set_title(f"{domain} {clade_id} ({regime}, n={nseq}) row {idx}/{len(chunks)}", loc="left", fontsize=10)

    fig.suptitle(f"ANKros vs {clade_id} consensus: {domain} (n={nseq})", fontsize=12, y=0.985)
    add_legends(fig, axes)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.77, bottom=0.14, hspace=0.92)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def draw_domain_overview(
    domain: str,
    clade_stats: list[tuple[str, str, int, list[dict[str, object]]]],
    out_png: Path,
    wrap: int,
) -> None:
    nonempty = [(cid, regime, nseq, stats) for cid, regime, nseq, stats in clade_stats if stats]
    if not nonempty:
        return
    n_chunks = max(math.ceil(len(stats) / wrap) for _cid, _regime, _nseq, stats in nonempty)
    fig_w = max(16.0, min(38.0, wrap * 0.34))
    fig_h = max(6.0, len(nonempty) * n_chunks * 2.25 + 2.4)
    fig, axes = plt.subplots(len(nonempty) * n_chunks, 1, figsize=(fig_w, fig_h), squeeze=False)
    axes = axes.ravel()

    axis_idx = 0
    for cid, regime, nseq, stats in nonempty:
        chunks = [stats[i : i + wrap] for i in range(0, len(stats), wrap)]
        for chunk_idx in range(n_chunks):
            ax = axes[axis_idx]
            axis_idx += 1
            if chunk_idx >= len(chunks):
                ax.axis("off")
                continue
            chunk = chunks[chunk_idx]
            ax.set_xlim(0, wrap)
            ax.set_ylim(-1.45, MAX_BITS + 2.35)
            ax.axhline(0, color="#333333", linewidth=0.45)
            ax.set_yticks([])
            label = f"{cid}\n{regime}, n={nseq}"
            if chunk_idx:
                label = ""
            ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9)

            for x, stat in enumerate(chunk):
                aa = str(stat["top_aa"])
                category = str(stat["category"])
                if stat["is_insert"] == "yes":
                    ax.add_patch(
                        Rectangle((x, -1.45), 1.0, MAX_BITS + 3.45, facecolor="#F4F4F5", edgecolor="none", zorder=-10)
                    )
                draw_letter(
                    ax,
                    aa,
                    x,
                    0.0,
                    float(stat["information_bits"]),
                    CATEGORY_COLORS.get(category, CATEGORY_COLORS["unknown"]),
                    width_scale=0.66,
                )
                add_size_bar(ax, aa, x, -1.24, max_height=0.98, width=0.28)
                add_size_bar(
                    ax,
                    str(stat.get("query_aa", "")),
                    x,
                    MAX_BITS + 0.62,
                    max_height=1.00,
                    width=0.14,
                    x_center=0.84,
                )
                add_ankros_track(
                    ax,
                    str(stat.get("query_aa", "")),
                    str(stat.get("query_category", "unknown")),
                    x,
                    MAX_BITS + 1.26,
                    height=1.18,
                    width_scale=0.34,
                )

            tick_pos, tick_lab = x_ticks_for_chunk(chunk)
            if axis_idx == len(axes) or chunk_idx == n_chunks - 1:
                ax.set_xticks(tick_pos)
                ax.set_xticklabels(tick_lab, fontsize=6, rotation=90)
            else:
                ax.set_xticks([])

    for ax in axes[axis_idx:]:
        ax.axis("off")

    fig.suptitle(f"ANKros vs clade consensus logos: {domain}", fontsize=12, y=0.985)
    add_legends(fig, axes)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.84, bottom=0.10, hspace=0.82)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def metric_value(aa: str, metric: str) -> float | None:
    if aa not in AA_SET:
        return None
    if metric == "hydropathy":
        return HYDROPATHY[aa]
    if metric == "volume":
        return AA_VOLUME[aa]
    raise ValueError(f"unknown metric: {metric}")


def unweighted_metric_value(aa: str, metric: str, _rsa_pct: object = "") -> float | None:
    return metric_value(aa, metric)


def smooth_series(points: list[tuple[int, float]], window: int) -> tuple[list[int], list[float]]:
    if window <= 1 or len(points) <= 2:
        return [x for x, _y in points], [y for _x, y in points]
    half = window // 2
    xs: list[int] = []
    ys: list[float] = []
    for idx, (x, _y) in enumerate(points):
        left = max(0, idx - half)
        right = min(len(points), idx + half + 1)
        vals = [value for _pos, value in points[left:right]]
        xs.append(x)
        ys.append(sum(vals) / len(vals))
    return xs, ys


def draw_metric_tracks(
    domain: str,
    clade_stats: list[tuple[str, str, int, list[dict[str, object]]]],
    out_png: Path,
    smooth_window: int,
) -> None:
    nonempty = [(cid, regime, nseq, stats) for cid, regime, nseq, stats in clade_stats if stats]
    if not nonempty:
        return
    start, end = DOMAINS[domain]
    fig, axes = plt.subplots(2, 1, figsize=(18.0, 7.2), sharex=True)
    metric_defs = [
        ("hydropathy", "Kyte-Doolittle hydrophobicity", (-4.8, 4.8)),
        ("volume", "AA volume", (50.0, 235.0)),
    ]

    query_by_qpos: dict[int, tuple[str, object]] = {}
    for _cid, _regime, _nseq, stats in nonempty:
        for stat in stats:
            qpos = str(stat.get("qpos", ""))
            query_aa = str(stat.get("query_aa", ""))
            if qpos and query_aa in AA_SET:
                query_by_qpos.setdefault(int(qpos), (query_aa, stat.get("query_rsa_pct", "")))

    handles = []
    labels = []
    for ax, (metric, ylabel, ylim) in zip(axes, metric_defs):
        ax.set_xlim(start, end)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylabel)
        ax.axhline(0, color="#8A8A8A", linewidth=0.7, linestyle="--" if metric == "hydropathy" else "-")
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.4)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.4)

        query_points = sorted(
            (qpos, value)
            for qpos, (aa, rsa_pct) in query_by_qpos.items()
            if (value := unweighted_metric_value(aa, metric, rsa_pct)) is not None
        )
        qx, qy = smooth_series(query_points, smooth_window)
        query_line = ax.plot(qx, qy, color="#111111", linewidth=1.9, label="ANKros")[0]
        if metric == "hydropathy":
            handles.append(query_line)
            labels.append("ANKros")

        by_regime_qpos: dict[str, dict[int, list[float]]] = {}
        clade_counts: Counter[str] = Counter()
        for cid, regime, _nseq, stats in nonempty:
            color = REGIME_COLORS.get(regime, "#6B7280")
            points = []
            for stat in stats:
                qpos = str(stat.get("qpos", ""))
                aa = str(stat.get("top_aa", ""))
                if not qpos:
                    continue
                value = unweighted_metric_value(aa, metric, stat.get("consensus_rsa_pct", ""))
                if value is None:
                    continue
                pos = int(qpos)
                points.append((pos, value))
                by_regime_qpos.setdefault(regime, {}).setdefault(pos, []).append(value)
            if not points:
                continue
            clade_counts[regime] += 1
            xs, ys = smooth_series(sorted(points), smooth_window)
            ax.plot(xs, ys, color=color, linewidth=0.8, alpha=0.28)

        for regime in sorted(by_regime_qpos, key=lambda item: ["psychro", "meso", "thermo"].index(item) if item in {"psychro", "meso", "thermo"} else 99):
            qpos_values = by_regime_qpos[regime]
            mean_points = sorted(
                (qpos, sum(values) / len(values))
                for qpos, values in qpos_values.items()
            )
            xs, ys = smooth_series(mean_points, smooth_window)
            label = f"{regime} mean ({clade_counts[regime]} clades)"
            line = ax.plot(
                xs,
                ys,
                color=REGIME_COLORS.get(regime, "#6B7280"),
                linewidth=2.2,
                label=label,
            )[0]
            if metric == "hydropathy":
                handles.append(line)
                labels.append(label)

    axes[-1].set_xlabel("ANKros position")
    axes[-1].set_xticks(list(range(start if start % 10 == 0 else start + (10 - start % 10), end + 1, 10)))
    fig.suptitle(f"ANKros-frame residue hydrophobicity and volume: {domain}", fontsize=13, y=0.975)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.925), ncol=min(4, len(handles)), frameon=False)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.84, bottom=0.10, hspace=0.22)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def html_float(value: object) -> float | None:
    try:
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def domain_for_qpos(qpos: object) -> str:
    qpos_text = str(qpos or "")
    if not qpos_text:
        return ""
    pos = int(qpos_text)
    for domain, (start, end) in DOMAINS.items():
        if start <= pos <= end:
            return domain
    return ""


def domain_start_class(qpos: object) -> str:
    qpos_text = str(qpos or "")
    if not qpos_text:
        return ""
    pos = int(qpos_text)
    starts = {start for start, _end in DOMAINS.values()}
    return " domain-start" if pos in starts else ""


def html_metric_value(aa: str, metric: str, _rsa_pct: object = "") -> float | None:
    return metric_value(aa, metric)


def html_heatmap_score_cell(
    value: float | None,
    metric: str,
    title: str,
    qpos: object = "",
    data_attrs: str = "",
) -> str:
    start_class = domain_start_class(qpos)
    if value is None:
        return f'<div class="cell score-cell heatmap-cell empty{start_class}" {data_attrs} title="{escape(title, quote=True)}"></div>'
    if metric == "camsol":
        fraction = min(1.0, abs(value) / 3.5)
        color = blend_hex("#F7F7F7", "#2E7D32" if value >= 0 else "#B2182B", fraction)
    elif metric == "aggrescan3d":
        fraction = min(1.0, abs(value) / 4.0)
        color = blend_hex("#F7F7F7", "#B2182B" if value >= 0 else "#2166AC", fraction)
    else:
        raise ValueError(f"unknown heatmap metric: {metric}")
    return (
        f'<div class="cell score-cell heatmap-cell {metric}-heatmap{start_class}" '
        f'{data_attrs} style="background:{color}" title="{escape(title, quote=True)}"></div>'
    )


def html_contact_cell(
    count: int | None,
    max_count: int,
    title: str,
    qpos: object = "",
    data_attrs: str = "",
) -> str:
    start_class = domain_start_class(qpos)
    if count is None:
        return f'<div class="cell score-cell contact-cell empty{start_class}" {data_attrs} title="{escape(title, quote=True)}"></div>'
    fraction = min(1.0, max(0.0, float(count)) / max(1.0, float(max_count)))
    color = blend_hex("#F7F7F7", "#5B5F97", fraction)
    return (
        f'<div class="cell score-cell contact-cell{start_class}" '
        f'{data_attrs} style="background:{color}" title="{escape(title, quote=True)}">'
        f'<span>{escape(str(count))}</span>'
        "</div>"
    )


def html_secondary_structure_cell(ss: object, title: str, qpos: object = "", data_attrs: str = "") -> str:
    state = str(ss or "")
    if state not in SECONDARY_STRUCTURE_LABELS:
        state = ""
    start_class = domain_start_class(qpos)
    if not state:
        return f'<div class="cell ss-cell empty{start_class}" {data_attrs} title="{escape(title, quote=True)}"></div>'
    color = SECONDARY_STRUCTURE_COLORS[state]
    label = SECONDARY_STRUCTURE_LABELS[state]
    return (
        f'<div class="cell ss-cell ss-{state.lower()}{start_class}" {data_attrs} '
        f'style="background:{color}" title="{escape(title, quote=True)}">'
        f'<span>{escape(state)}</span>'
        f'<b>{escape(label)}</b>'
        "</div>"
    )


def html_category_class(category: object) -> str:
    category_text = str(category or "unknown")
    if category_text not in CATEGORY_COLORS:
        category_text = "unknown"
    return f"cat-{category_text}"


def html_letter_cell(
    aa: str,
    letter_category: object,
    column_category: object,
    title: str,
    font_size_px: float,
    qpos: object = "",
    data_attrs: str = "",
) -> str:
    display = aa if aa in AA_SET else "-"
    letter_color = CATEGORY_COLORS.get(str(letter_category or "unknown"), CATEGORY_COLORS["unknown"])
    classes = f"cell aa-cell {html_category_class(column_category)}{domain_start_class(qpos)}"
    extra = " gap" if aa not in AA_SET else ""
    return (
        f'<div class="{classes}{extra}" {data_attrs} title="{escape(title, quote=True)}">'
        f'<span style="color:{letter_color};font-size:{font_size_px:.1f}px">{escape(display)}</span>'
        "</div>"
    )


def html_score_cell(value: float | None, metric: str, title: str, qpos: object = "", data_attrs: str = "") -> str:
    start_class = domain_start_class(qpos)
    if value is None:
        return f'<div class="cell score-cell empty{start_class}" {data_attrs} title="{escape(title, quote=True)}"></div>'
    if metric == "hydropathy":
        fraction = min(1.0, abs(value) / max(abs(min(HYDROPATHY.values())), abs(max(HYDROPATHY.values()))))
        height_pct = 48.0 * fraction
        if value >= 0:
            style = f"height:{height_pct:.2f}%;bottom:50%;"
            cls = "hydro-positive"
        else:
            style = f"height:{height_pct:.2f}%;top:50%;"
            cls = "hydro-negative"
        return (
            f'<div class="cell score-cell hydro-score{start_class}" {data_attrs} title="{escape(title, quote=True)}">'
            f'<span class="{cls}" style="{style}"></span>'
            "</div>"
        )
    if metric == "volume":
        height_pct = 100.0 * min(1.0, max(0.0, value) / MAX_VOLUME)
        return (
            f'<div class="cell score-cell volume-score{start_class}" {data_attrs} title="{escape(title, quote=True)}">'
            f'<span style="height:{height_pct:.2f}%"></span>'
            "</div>"
        )
    if metric == "camsol":
        return html_heatmap_score_cell(value, metric, title, qpos=qpos, data_attrs=data_attrs)
    if metric == "aggrescan3d":
        return html_heatmap_score_cell(value, metric, title, qpos=qpos, data_attrs=data_attrs)
    raise ValueError(f"unknown metric: {metric}")


def html_conservation_cell(stat: dict[str, object], data_attrs: str = "") -> str:
    qpos = str(stat.get("qpos", ""))
    top_aa = str(stat.get("top_aa", ""))
    conservation = html_float(stat.get("conservation", "")) or 0.0
    occupancy = html_float(stat.get("occupancy", "")) or 0.0
    info = html_float(stat.get("information_bits", "")) or 0.0
    raw_height_scale = min(1.0, max(0.0, conservation))
    title = (
        f"Total MSA conservation at {qpos or 'insert'}: {conservation:.3f}; "
        f"top={top_aa or '-'}; occupancy={occupancy:.3f}; information_bits={info:.3f}"
    )
    label = top_aa if top_aa in AA_SET else "-"
    height_scale = max(0.28, raw_height_scale) if label in AA_SET else raw_height_scale
    extra = " gap" if top_aa not in AA_SET else ""
    return (
        f'<div class="cell conservation-cell{domain_start_class(qpos)}{extra}" {data_attrs} title="{escape(title, quote=True)}">'
        f'<span class="conservation-logo" style="transform:scaleX(0.82) scaleY({height_scale:.4f})">{escape(label)}</span>'
        "</div>"
    )


def html_position_label(stat: dict[str, object]) -> str:
    qpos = str(stat.get("qpos", ""))
    if qpos:
        return qpos if int(qpos) % 10 == 0 else ""
    return "+"


def html_row(label: str, cells_html: str, classes: str = "", data_attrs: str = "") -> str:
    return (
        f'<div class="row-label {classes}" {data_attrs}>{label}</div>'
        f'<div class="row-grid {classes}" {data_attrs}>{cells_html}</div>'
    )


def filter_empty_html_columns(
    clade_stats: list[tuple[str, str, int, list[dict[str, object]]]],
) -> list[tuple[str, str, int, list[dict[str, object]]]]:
    keep = html_nonempty_column_indices(clade_stats)
    return [
        (cid, regime, nseq, [stats[idx] for idx in keep])
        for cid, regime, nseq, stats in clade_stats
    ]


def html_nonempty_column_indices(
    clade_stats: list[tuple[str, str, int, list[dict[str, object]]]],
) -> list[int]:
    if not clade_stats:
        return []
    n_cols = min(len(stats) for _cid, _regime, _nseq, stats in clade_stats)
    keep: list[int] = []
    for idx in range(n_cols):
        if any(
            str(stats[idx].get("query_aa", "")) in AA_SET or str(stats[idx].get("top_aa", "")) in AA_SET
            for _cid, _regime, _nseq, stats in clade_stats
        ):
            keep.append(idx)
    return keep


def filter_html_columns_by_indices(
    rows: list[dict[str, object]],
    keep: list[int],
) -> list[dict[str, object]]:
    return [rows[idx] for idx in keep if idx < len(rows)]


def total_msa_conservation_stats(
    entries: list[tuple[str, str, str]],
    selected_cols: list[tuple[int, dict[str, str], bool, int]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    nseq = len(entries)
    for domain_col, (aln_col, colrow, is_insert, insert_run_len) in enumerate(selected_cols, start=1):
        counts = aa_counts(entries, aln_col)
        gaps = gap_count(entries, aln_col)
        n_residues = sum(counts.values())
        if counts:
            top_aa, top_count = sorted(counts.items(), key=lambda item: (-item[1], AA_ORDER.index(item[0])))[0]
        else:
            top_aa, top_count = "-", 0
        info, occupancy, gap_fraction = information_bits(counts, nseq)
        conservation = (top_count / n_residues) * occupancy if top_aa in AA_SET and n_residues else 0.0
        rows.append(
            {
                "domain_col": domain_col,
                "alignment_col": aln_col,
                "qpos": colrow.get("qpos", ""),
                "is_insert": "yes" if is_insert else "no",
                "insert_run_len": insert_run_len if is_insert else "",
                "source_block": colrow.get("source_block", ""),
                "source_range": colrow.get("source_range", ""),
                "n_sequences": nseq,
                "n_residues": n_residues,
                "n_gaps": gaps,
                "occupancy": f"{occupancy:.4f}",
                "gap_fraction": f"{gap_fraction:.4f}",
                "information_bits": f"{info:.4f}" if top_aa in AA_SET else "0.0000",
                "top_aa": top_aa,
                "top_aa_count": top_count,
                "top_aa_frequency_among_residues": (
                    f"{(top_count / n_residues):.4f}" if top_aa in AA_SET and n_residues else "0.0000"
                ),
                "top_symbol_frequency_among_sequences": f"{(top_count / nseq):.4f}" if nseq else "0.0000",
                "conservation": f"{conservation:.4f}",
            }
        )
    return rows


def write_interactive_html(
    out_html: Path,
    clade_stats: list[tuple[str, str, int, list[dict[str, object]]]],
    total_conservation_stats: list[dict[str, object]],
    default_visible_per_regime: int,
) -> None:
    nonempty = [(cid, regime, nseq, stats) for cid, regime, nseq, stats in clade_stats if stats]
    if not nonempty:
        return
    regime_order = {"psychro": 0, "meso": 1, "thermo": 2}
    nonempty = sorted(nonempty, key=lambda item: (regime_order.get(item[1], 99), item[0]))
    reference_stats = nonempty[0][3]
    n_cols = len(reference_stats)
    total_conservation_stats = total_conservation_stats[:n_cols]
    max_per_regime = Counter(regime for _cid, regime, _nseq, _stats in nonempty)
    max_visible = max(max_per_regime.values()) if max_per_regime else 1
    default_visible = max(1, min(default_visible_per_regime, max_visible))

    domain_cells = []
    current_domain = ""
    for idx, stat in enumerate(reference_stats):
        qpos = str(stat.get("qpos", ""))
        domain = domain_for_qpos(qpos)
        label = ""
        if domain and domain != current_domain:
            label = domain
            current_domain = domain
        domain_cells.append(
            f'<div class="cell domain-cell{domain_start_class(qpos)}" data-col="{idx}">'
            f'<span>{escape(label)}</span>'
            "</div>"
        )
    position_cells = "".join(
        f'<div class="cell pos-cell {html_category_class(stat.get("query_category", "unknown"))}{domain_start_class(stat.get("qpos", ""))}" data-col="{idx}">{escape(html_position_label(stat))}</div>'
        for idx, stat in enumerate(reference_stats)
    )
    rows_html = [
        html_row(
            "Total MSA<br><span>conservation</span>",
            "".join(
                html_conservation_cell(stat, data_attrs=f'data-col="{idx}"')
                for idx, stat in enumerate(total_conservation_stats)
            ),
            classes="conservation-row",
        ),
        html_row("Domain", "".join(domain_cells), classes="domain-row"),
        html_row("Position", position_cells, classes="position-row"),
    ]

    ankros_letters = []
    ankros_secondary_structure = []
    ankros_hydro = []
    ankros_volume = []
    ankros_contacts = []
    ankros_annotations = []
    ankros_camsol = []
    ankros_a3d = []
    max_ankros_contacts = max(
        [
            int(html_float(stat.get("query_sidechain_contact_count", "")) or 0)
            for stat in reference_stats
            if str(stat.get("qpos", ""))
        ]
        or [1]
    )
    for idx, stat in enumerate(reference_stats):
        aa = str(stat.get("query_aa", ""))
        category = stat.get("query_category", "unknown")
        qpos = str(stat.get("qpos", ""))
        title = f"ANKros {aa or '-'} at {qpos or 'insert'}; category={category}; RSA={stat.get('query_rsa_pct', '')}"
        query_evidence = str(stat.get("query_functional_evidence", ""))
        if query_evidence:
            title += f"; holo contacts={query_evidence}"
        col_attrs = f'data-col="{idx}"'
        ankros_letters.append(html_letter_cell(aa, category, category, title, 18.0, qpos=qpos, data_attrs=col_attrs))
        query_ss = str(stat.get("query_secondary_structure", ""))
        ankros_secondary_structure.append(
            html_secondary_structure_cell(
                query_ss,
                (
                    "ANKros secondary structure "
                    f"at {qpos or 'insert'}: {SECONDARY_STRUCTURE_LABELS.get(query_ss, 'unmapped')}"
                ),
                qpos=qpos,
                data_attrs=col_attrs,
            )
        )
        hyd_value = html_metric_value(aa, "hydropathy", stat.get("query_rsa_pct", ""))
        vol_value = html_metric_value(aa, "volume", stat.get("query_rsa_pct", ""))
        ankros_hydro.append(html_score_cell(hyd_value, "hydropathy", f"ANKros hydropathy: {hyd_value}", qpos=qpos, data_attrs=col_attrs))
        ankros_volume.append(html_score_cell(vol_value, "volume", f"ANKros volume: {vol_value}", qpos=qpos, data_attrs=col_attrs))
        contact_value = html_float(stat.get("query_sidechain_contact_count", ""))
        contact_count = int(contact_value) if contact_value is not None else None
        contact_partners = str(stat.get("query_sidechain_contact_residues", ""))
        ankros_contacts.append(
            html_contact_cell(
                contact_count,
                max_ankros_contacts,
                (
                    "ANKros side-chain contacts within cutoff: "
                    f"{contact_count if contact_count is not None else 'unmapped'}; "
                    f"partners: {contact_partners or 'none'}"
                ),
                qpos=qpos,
                data_attrs=col_attrs,
            )
        )
        annotation_labels = [part.strip() for part in query_evidence.split(";") if part.strip()]
        ankros_annotations.append(
            html_contact_cell(
                len(annotation_labels) if annotation_labels else None,
                5,
                (
                    "Measured holo/receptor contacts to FAD, FMN, CPD, or DNA: "
                    f"{'; '.join(annotation_labels) if annotation_labels else 'none'}"
                ),
                qpos=qpos,
                data_attrs=col_attrs,
            )
        )
        camsol_value = html_float(stat.get("ankros_camsol_structural", ""))
        a3d_value = html_float(stat.get("ankros_aggrescan3d", ""))
        ankros_camsol.append(
            html_score_cell(
                camsol_value,
                "camsol",
                f"ANKros CamSol structural solubility: {camsol_value}",
                qpos=qpos,
                data_attrs=col_attrs,
            )
        )
        ankros_a3d.append(
            html_score_cell(
                a3d_value,
                "aggrescan3d",
                f"ANKros Aggrescan3D score: {a3d_value}",
                qpos=qpos,
                data_attrs=col_attrs,
            )
        )
    rows_html.append(html_row("ANKros", "".join(ankros_letters), classes="ankros-row"))
    rows_html.append(html_row("SS", "".join(ankros_secondary_structure), classes="score-row ankros-row secondary-structure-row"))
    rows_html.append(html_row("hydro", "".join(ankros_hydro), classes="score-row ankros-row"))
    rows_html.append(html_row("volume", "".join(ankros_volume), classes="score-row ankros-row"))
    rows_html.append(html_row("SC contacts", "".join(ankros_contacts), classes="score-row ankros-row"))
    rows_html.append(html_row("holo contacts", "".join(ankros_annotations), classes="score-row ankros-row"))
    rows_html.append(html_row("CamSol", "".join(ankros_camsol), classes="score-row ankros-row"))
    rows_html.append(html_row("A3D", "".join(ankros_a3d), classes="score-row ankros-row"))

    rank_by_regime: Counter[str] = Counter()
    for clade_id, regime, nseq, stats in nonempty:
        rank_by_regime[regime] += 1
        rank = rank_by_regime[regime]
        data_attrs = f'data-regime="{escape(regime, quote=True)}" data-rank="{rank}"'
        sequence_cells = []
        representative_secondary_structure = []
        representative_camsol = []
        representative_a3d = []
        hydro_cells = []
        volume_cells = []
        for idx, stat in enumerate(stats):
            aa = str(stat.get("top_aa", ""))
            representative_id = str(stat.get("representative_id", ""))
            qcat = stat.get("query_category", "unknown")
            category = stat.get("category", "unknown")
            info = html_float(stat.get("information_bits", "")) or 0.0
            font_size = 10.0 + 15.0 * min(1.0, info / MAX_BITS)
            qpos = str(stat.get("qpos", ""))
            title = (
                f"{clade_id} {aa or '-'} at {qpos or 'insert'}; "
                f"n={nseq}; category={category}; freq={stat.get('top_symbol_frequency_among_sequences', '')}; "
                f"info={stat.get('information_bits', '')}; RSA={stat.get('consensus_rsa_pct', '')}"
            )
            col_attrs = f'data-col="{idx}" data-has-residue="{1 if aa in AA_SET else 0}"'
            sequence_cells.append(html_letter_cell(aa, category, qcat, title, font_size, qpos=qpos, data_attrs=col_attrs))
            rep_ss = str(stat.get("representative_secondary_structure", ""))
            rep_pos = stat.get("representative_structure_pos", "")
            representative_secondary_structure.append(
                html_secondary_structure_cell(
                    rep_ss,
                    (
                        f"{representative_id or clade_id} secondary structure "
                        f"at structure residue {rep_pos or 'unmapped'}: "
                        f"{SECONDARY_STRUCTURE_LABELS.get(rep_ss, 'unmapped')}"
                    ),
                    qpos=qpos,
                    data_attrs=f'data-col="{idx}"',
                )
            )
            rep_camsol_value = html_float(stat.get("representative_camsol_structural", ""))
            rep_a3d_value = html_float(stat.get("representative_aggrescan3d", ""))
            rep_pos = stat.get("representative_structure_pos", "")
            representative_camsol.append(
                html_score_cell(
                    rep_camsol_value,
                    "camsol",
                    f"{representative_id or clade_id} CamSol structural solubility at structure residue {rep_pos or 'unmapped'}: {rep_camsol_value}",
                    qpos=qpos,
                    data_attrs=f'data-col="{idx}"',
                )
            )
            representative_a3d.append(
                html_score_cell(
                    rep_a3d_value,
                    "aggrescan3d",
                    f"{representative_id or clade_id} Aggrescan3D at structure residue {rep_pos or 'unmapped'}: {rep_a3d_value}",
                    qpos=qpos,
                    data_attrs=f'data-col="{idx}"',
                )
            )
            hyd_value = html_metric_value(aa, "hydropathy", stat.get("consensus_rsa_pct", ""))
            vol_value = html_metric_value(aa, "volume", stat.get("consensus_rsa_pct", ""))
            hydro_cells.append(html_score_cell(hyd_value, "hydropathy", f"{clade_id} hydropathy: {hyd_value}", qpos=qpos, data_attrs=f'data-col="{idx}"'))
            volume_cells.append(html_score_cell(vol_value, "volume", f"{clade_id} volume: {vol_value}", qpos=qpos, data_attrs=f'data-col="{idx}"'))
        label = f"{escape(clade_id)}<br><span>{escape(regime)}, n={nseq}</span>"
        representative_id_label = escape(str(stats[0].get("representative_id", ""))) if stats else ""
        row_classes = f"clade-row regime-{escape(regime)}"
        rows_html.append(html_row(label, "".join(sequence_cells), classes=row_classes, data_attrs=data_attrs))
        rows_html.append(html_row(f"rep SS<br><span>{representative_id_label}</span>", "".join(representative_secondary_structure), classes=f"{row_classes} representative-row score-row secondary-structure-row", data_attrs=data_attrs))
        rows_html.append(html_row(f"rep CamSol<br><span>{representative_id_label}</span>", "".join(representative_camsol), classes=f"{row_classes} representative-row score-row", data_attrs=data_attrs))
        rows_html.append(html_row(f"rep A3D<br><span>{representative_id_label}</span>", "".join(representative_a3d), classes=f"{row_classes} representative-row score-row", data_attrs=data_attrs))
        rows_html.append(html_row("hydro", "".join(hydro_cells), classes=f"{row_classes} score-row", data_attrs=data_attrs))
        rows_html.append(html_row("volume", "".join(volume_cells), classes=f"{row_classes} score-row", data_attrs=data_attrs))

    regime_controls = "".join(
        f'<label><input type="checkbox" class="regime-toggle" value="{escape(regime)}" checked> {escape(regime)} ({count})</label>'
        for regime, count in sorted(max_per_regime.items(), key=lambda item: regime_order.get(item[0], 99))
    )
    colors_json = json.dumps({key: value for key, value in CATEGORY_COLORS.items()})
    column_meta_json = json.dumps(
        [
            {
                "qpos": str(stat.get("qpos", "")),
                "isInsert": str(stat.get("is_insert", "")) == "yes",
            }
            for stat in reference_stats
        ]
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ANKros thermal clade logos</title>
<style>
:root {{
  --cell-w: 22px;
  --visible-cols: {n_cols};
  --cell-h: 34px;
  --score-h: 34px;
  --label-w: 190px;
  --blue: #2F80ED;
  --green: #2E7D32;
  --red: #D1495B;
}}
body {{ margin: 0; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #1f2937; background: #ffffff; }}
header {{ position: sticky; top: 0; z-index: 20; background: #ffffff; border-bottom: 1px solid #d1d5db; padding: 12px 18px; }}
h1 {{ margin: 0 0 10px; font-size: 20px; font-weight: 650; }}
.controls {{ display: flex; flex-wrap: wrap; gap: 14px; align-items: center; font-size: 13px; }}
.controls input[type="number"] {{ width: 64px; }}
.legend {{ display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin-top: 8px; font-size: 12px; }}
.swatch {{ display: inline-block; width: 18px; height: 12px; border: 1px solid #9ca3af; vertical-align: -2px; }}
.viewer {{ padding: 16px 18px 28px; }}
.sequence-scroll {{ overflow-x: auto; border: 1px solid #d1d5db; }}
.track-table {{ display: grid; grid-template-columns: var(--label-w) max-content; width: max-content; min-width: 100%; }}
.row-label {{ position: sticky; left: 0; z-index: 5; width: var(--label-w); box-sizing: border-box; padding: 4px 10px; background: #ffffff; border-right: 1px solid #d1d5db; border-bottom: 1px solid #e5e7eb; font-size: 13px; text-align: right; line-height: 1.1; }}
.row-label span {{ color: #6b7280; font-size: 12px; }}
.row-grid {{ display: grid; grid-template-columns: repeat(var(--visible-cols), var(--cell-w)); border-bottom: 1px solid #e5e7eb; }}
.cell {{ width: var(--cell-w); box-sizing: border-box; border-right: 1px solid #f3f4f6; text-align: center; position: relative; overflow: hidden; }}
.aa-cell {{ height: var(--cell-h); line-height: var(--cell-h); font-family: "DejaVu Sans Mono", ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 800; }}
.aa-cell span {{ display: inline-block; width: 100%; transform: scaleX(0.82); transform-origin: center; }}
.gap span {{ color: #9ca3af !important; font-weight: 600; }}
.domain-cell {{ height: 24px; line-height: 24px; overflow: visible; background: #ffffff; }}
.domain-cell span {{ position: absolute; left: 5px; top: 1px; z-index: 2; font-size: 12px; font-weight: 700; color: #111827; text-transform: uppercase; letter-spacing: 0.04em; }}
.pos-cell {{ height: 18px; line-height: 18px; font-size: 10px; color: #374151; }}
.conservation-cell {{ height: 34px; background: #ffffff; }}
.conservation-logo {{ position: absolute; left: 0; right: 0; bottom: -2px; display: block; height: 32px; line-height: 32px; font-family: "DejaVu Sans Mono", ui-monospace, SFMono-Regular, Menlo, monospace; color: #111827; font-weight: 900; font-size: 32px; transform-origin: bottom center; }}
.score-cell {{ height: var(--score-h); background: #fafafa; }}
.hydro-score::after {{ content: ""; position: absolute; left: 0; right: 0; top: 50%; border-top: 1px solid #9ca3af; }}
.hydro-score span, .volume-score span {{ position: absolute; left: 4px; right: 4px; display: block; border: 1px solid rgba(31,41,55,0.45); }}
.hydro-positive {{ background: #B2182B; }}
.hydro-negative {{ background: #2166AC; }}
.volume-score span {{ bottom: 1px; background: #64748B; }}
.heatmap-cell {{ border-right-color: rgba(31,41,55,0.16); }}
.contact-cell {{ border-right-color: rgba(31,41,55,0.16); }}
.contact-cell span {{ display: block; color: #111827; font-size: 9px; line-height: var(--score-h); opacity: 0.82; }}
.ss-cell {{ height: var(--score-h); border-right-color: rgba(31,41,55,0.16); color: #ffffff; }}
.ss-cell span {{ display: block; font-size: 11px; line-height: 15px; font-weight: 800; margin-top: 2px; }}
.ss-cell b {{ display: block; font-size: 7px; line-height: 8px; font-weight: 600; opacity: 0.92; }}
.ss-cell.empty {{ background: #fafafa; }}
.cat-surface {{ background: rgba(47, 128, 237, 0.14); }}
.cat-core {{ background: rgba(85, 85, 85, 0.13); }}
.cat-functional {{ background: rgba(209, 73, 91, 0.20); }}
.cat-insert {{ background: rgba(156, 163, 175, 0.18); }}
.cat-unknown {{ background: rgba(199, 199, 199, 0.10); }}
.domain-start {{ border-left: 9px solid #ffffff; box-shadow: inset 1px 0 0 #9ca3af; }}
.score-row .row-label, .score-row.row-label {{ color: #6b7280; font-size: 11px; }}
.representative-row .row-label, .representative-row.row-label {{ color: #4b5563; font-size: 11px; }}
.ankros-row .row-label {{ font-weight: 700; }}
.regime-psychro .row-label {{ border-left: 4px solid var(--blue); }}
.regime-meso .row-label {{ border-left: 4px solid var(--green); }}
.regime-thermo .row-label {{ border-left: 4px solid var(--red); }}
.hidden {{ display: none; }}
.col-hidden {{ display: none; }}
</style>
</head>
<body>
<header>
  <h1>ANKros thermal clade consensus logos</h1>
  <div class="controls">
    <label>Visible clades per regime
      <input id="visiblePerRegime" type="number" min="1" max="{max_visible}" value="{default_visible}">
    </label>
    {regime_controls}
    <span>Columns: <span id="visibleColumnCount">{n_cols}</span> / {n_cols}</span>
  </div>
  <div class="legend">
    <span><span class="swatch cat-surface"></span> surface/SAS exposed</span>
    <span><span class="swatch cat-core"></span> buried/core</span>
    <span><span class="swatch cat-functional"></span> functional/ligand/antenna/ET</span>
    <span><span class="swatch cat-insert"></span> insert/no ANKros position</span>
    <span>Letter color = residue/consensus structural class; column shading = ANKros structural class.</span>
    <span>SS uses mkdssp/DSSP when available, with Biotite P-SEA fallback: H = helix, E = strand, C = coil.</span>
    <span>Hydropathy uses Kyte-Doolittle values; volume uses the AA volume scale.</span>
    <span>SC contacts counts unique ANKros residue side-chain heavy-atom contacts within the configured cutoff.</span>
    <span>Holo annotations use measured FAD/FMN/CPD/DNA distances from the FAD/FMN holo-builder QC table.</span>
    <span>CamSol heatmap uses ANKros structural solubility: green = more soluble, red = less soluble.</span>
    <span>A3D heatmap: red = aggregation-prone, blue = low aggregation score.</span>
    <span>Rep tracks show secondary structure, CamSol, and A3D scores for the best-covered representative structure in each displayed clade.</span>
  </div>
</header>
<main class="viewer">
  <div class="sequence-scroll">
    <div class="track-table">
      {''.join(rows_html)}
    </div>
  </div>
</main>
<script>
const categoryColors = {colors_json};
const columnMeta = {column_meta_json};
function updateRows() {{
  const visible = Number(document.getElementById('visiblePerRegime').value || 1);
  const enabled = new Set(Array.from(document.querySelectorAll('.regime-toggle:checked')).map(el => el.value));
  document.querySelectorAll('.clade-row').forEach(el => {{
    const regime = el.dataset.regime;
    const rank = Number(el.dataset.rank || 0);
    el.classList.toggle('hidden', !enabled.has(regime) || rank > visible);
  }});
  updateColumns();
}}
function updateColumns() {{
  const visibleRows = Array.from(document.querySelectorAll('.row-grid.clade-row:not(.score-row):not(.hidden)'));
  const keep = columnMeta.map((meta, idx) => {{
    if (!meta.isInsert || meta.qpos) {{
      return true;
    }}
    return visibleRows.some(row => {{
      const cell = row.querySelector(`[data-col="${{idx}}"]`);
      return cell && cell.dataset.hasResidue === '1';
    }});
  }});
  const visibleCount = keep.filter(Boolean).length;
  document.documentElement.style.setProperty('--visible-cols', String(Math.max(1, visibleCount)));
  document.getElementById('visibleColumnCount').textContent = String(visibleCount);
  document.querySelectorAll('.row-grid').forEach(row => {{
    Array.from(row.children).forEach((cell, idx) => {{
      cell.classList.toggle('col-hidden', !keep[idx]);
    }});
  }});
}}
document.getElementById('visiblePerRegime').addEventListener('input', updateRows);
document.querySelectorAll('.regime-toggle').forEach(el => el.addEventListener('change', updateRows));
updateRows();
</script>
</body>
</html>
"""
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html)


def build_column_stats(
    domain: str,
    clade_id: str,
    regime: str,
    selected_entries: list[tuple[str, str, str]],
    selected_cols: list[tuple[int, dict[str, str], bool, int]],
    query_seq: str,
    seq_pos_by_sid: dict[str, dict[int, int]],
    afdb_features: dict[str, dict[int, dict[str, object]]],
    structure_pos_by_entry_pos_by_sid: dict[str, dict[int, int]],
    ankros_features: dict[int, dict[str, object]],
    residue_score_cache: dict[str, dict[int, dict[str, float]]],
    representative_entry: tuple[str, str, str] | None,
    include_empty_inserts: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    nseq = len(selected_entries)
    representative_id = representative_entry[0] if representative_entry else ""
    representative_seq = representative_entry[2] if representative_entry else ""
    representative_structure_id = structure_id_for_alignment_id(representative_id) if representative_id else ""
    for domain_col, (aln_col, colrow, is_insert, insert_run_len) in enumerate(selected_cols, start=1):
        counts = aa_counts(selected_entries, aln_col)
        gaps = gap_count(selected_entries, aln_col)
        top_aa, top_count = top_residue(counts, gaps)
        if is_insert and not include_empty_inserts and top_aa == "-":
            continue
        info, occupancy, gap_fraction = information_bits(counts, nseq)
        query_aa = ""
        if aln_col < len(query_seq) and query_seq[aln_col].upper() in AA_SET:
            query_aa = query_seq[aln_col].upper()
        qpos = colrow.get("qpos", "")
        (
            category,
            n_structures_for_top_aa,
            n_matching_structure_residues,
            consensus_sasa,
            consensus_rsa,
        ) = clade_structural_category(
            top_aa,
            qpos,
            selected_entries,
            aln_col,
            seq_pos_by_sid,
            afdb_features,
            structure_pos_by_entry_pos_by_sid,
            ankros_features,
        )
        if is_insert and category == "unknown":
            category = "insert"
        query_category = "insert" if is_insert and not qpos else "unknown"
        query_sasa = ""
        query_rsa = ""
        query_sidechain_contact_count = ""
        query_sidechain_contact_residues = ""
        query_secondary_structure = ""
        query_functional_evidence = ""
        if qpos:
            query_feature = ankros_features.get(int(qpos), {})
            query_category = str(query_feature.get("category") or "unknown")
            query_sasa = query_feature.get("sasa_a2", "")
            query_rsa = query_feature.get("rsa_pct", "")
            query_sidechain_contact_count = query_feature.get("sidechain_contact_count", "")
            query_sidechain_contact_residues = query_feature.get("sidechain_contact_residues", "")
            query_secondary_structure = query_feature.get("secondary_structure", "")
            query_functional_evidence = query_feature.get("functional_evidence", "")
        query_scores = residue_score_cache.get(ANKROS_SCORE_IDS[0], {}).get(int(qpos), {}) if qpos else {}
        representative_aa = ""
        representative_structure_pos = ""
        representative_scores: dict[str, float] = {}
        representative_secondary_structure = ""
        if representative_seq and aln_col < len(representative_seq):
            representative_symbol = representative_seq[aln_col].upper()
            representative_aa = representative_symbol if representative_symbol in AA_SET else "-"
            if representative_aa in AA_SET:
                representative_entry_pos = seq_pos_by_sid.get(representative_id, {}).get(aln_col)
                if representative_entry_pos is not None:
                    representative_structure_pos = structure_pos_by_entry_pos_by_sid.get(representative_id, {}).get(
                        representative_entry_pos,
                        "",
                    )
                if isinstance(representative_structure_pos, int):
                    representative_scores = residue_score_cache.get(representative_structure_id, {}).get(
                        representative_structure_pos,
                        {},
                    )
                    representative_secondary_structure = str(
                        afdb_features.get(representative_id, {})
                        .get(representative_structure_pos, {})
                        .get("secondary_structure", "")
                    )
        rows.append(
            {
                "clade_id": clade_id,
                "regime": regime,
                "domain": domain,
                "domain_col": domain_col,
                "alignment_col": aln_col,
                "qpos": qpos,
                "is_insert": "yes" if is_insert else "no",
                "insert_run_len": insert_run_len if is_insert else "",
                "source_block": colrow.get("source_block", ""),
                "source_range": colrow.get("source_range", ""),
                "n_sequences": nseq,
                "n_residues": sum(counts.values()),
                "n_gaps": gaps,
                "n_structures_for_top_aa": n_structures_for_top_aa,
                "n_matching_structure_residues": n_matching_structure_residues,
                "consensus_sasa_a2": f"{consensus_sasa:.2f}" if consensus_sasa is not None else "",
                "consensus_rsa_pct": f"{consensus_rsa:.2f}" if consensus_rsa is not None else "",
                "occupancy": f"{occupancy:.4f}",
                "gap_fraction": f"{gap_fraction:.4f}",
                "information_bits": f"{info:.4f}" if top_aa in AA_SET else "0.0000",
                "top_aa": top_aa,
                "top_aa_count": top_count,
                "top_aa_frequency_among_residues": (
                    f"{(top_count / sum(counts.values())):.4f}" if top_aa in AA_SET and counts else "0.0000"
                ),
                "top_symbol_frequency_among_sequences": f"{(top_count / nseq):.4f}" if nseq else "0.0000",
                "representative_id": representative_id,
                "representative_structure_id": representative_structure_id,
                "representative_structure_pos": representative_structure_pos,
                "representative_aa": representative_aa,
                "representative_secondary_structure": representative_secondary_structure,
                "representative_matches_consensus": (
                    "yes" if representative_aa in AA_SET and representative_aa == top_aa else "no"
                ),
                "query_aa": query_aa,
                "query_secondary_structure": query_secondary_structure,
                "query_functional_evidence": query_functional_evidence,
                "query_matches_consensus": "yes" if query_aa and query_aa == top_aa else "no",
                "query_category": query_category,
                "query_sasa_a2": f"{float(query_sasa):.2f}" if isinstance(query_sasa, float) else query_sasa,
                "query_rsa_pct": f"{float(query_rsa):.2f}" if isinstance(query_rsa, float) else query_rsa,
                "query_aa_volume": f"{AA_VOLUME.get(query_aa, 0.0):.1f}" if query_aa in AA_SET else "",
                "query_aa_hydropathy": f"{HYDROPATHY.get(query_aa, 0.0):.1f}" if query_aa in AA_SET else "",
                "query_sidechain_contact_count": query_sidechain_contact_count,
                "query_sidechain_contact_residues": query_sidechain_contact_residues,
                "ankros_camsol_intrinsic": format_optional_float(query_scores.get("camsol_intrinsic")),
                "ankros_camsol_structural": format_optional_float(query_scores.get("camsol_structural")),
                "ankros_aggrescan3d": format_optional_float(query_scores.get("aggrescan3d")),
                "ankros_aggrescan3d_positive": format_optional_float(query_scores.get("aggrescan3d_positive")),
                "representative_camsol_intrinsic": format_optional_float(representative_scores.get("camsol_intrinsic")),
                "representative_camsol_structural": format_optional_float(representative_scores.get("camsol_structural")),
                "representative_aggrescan3d": format_optional_float(representative_scores.get("aggrescan3d")),
                "representative_aggrescan3d_positive": format_optional_float(representative_scores.get("aggrescan3d_positive")),
                "category": category,
                "aa_volume": f"{AA_VOLUME.get(top_aa, 0.0):.1f}" if top_aa in AA_SET else "",
                "aa_hydropathy": f"{HYDROPATHY.get(top_aa, 0.0):.1f}" if top_aa in AA_SET else "",
            }
        )
    return rows


def validate_inputs(paths: list[tuple[Path, str]]) -> None:
    missing = [f"{label}: {path}" for path, label in paths if not path.exists()]
    if missing:
        raise SystemExit("ERROR: required input(s) missing:\n" + "\n".join(missing))


def select_representative_clades(
    clade_rows: list[dict[str, str]],
    regimes: set[str],
    tip_set: str,
    entry_by_id: dict[str, tuple[str, str, str]],
    max_per_regime: int,
    min_sequences: int,
) -> list[dict[str, str]]:
    filtered = [row for row in clade_rows if not regimes or row.get("regime", "") in regimes]
    if max_per_regime == 0:
        return filtered

    by_regime: dict[str, list[tuple[int, int, str, dict[str, str]]]] = {}
    for row in filtered:
        tips = split_tips(row.get(f"{tip_set}_tips", ""))
        n_available = sum(1 for tip in tips if tip in entry_by_id)
        if n_available < min_sequences:
            continue
        try:
            n_total = int(row.get("n_total_tips", "0") or 0)
        except ValueError:
            n_total = 0
        by_regime.setdefault(row.get("regime", ""), []).append(
            (-n_available, -n_total, row.get("clade_id", ""), row)
        )

    selected: list[dict[str, str]] = []
    for regime in sorted(by_regime):
        ranked = sorted(by_regime[regime])
        selected.extend(item[-1] for item in ranked[:max_per_regime])
    return selected


def merge_clades_by_regime(
    clade_rows: list[dict[str, str]],
    regimes: set[str],
    tip_set: str,
    entry_by_id: dict[str, tuple[str, str, str]],
    min_sequences: int,
) -> list[dict[str, str]]:
    tips_by_regime: dict[str, set[str]] = defaultdict(set)
    source_count_by_regime: Counter[str] = Counter()
    missing_by_regime: Counter[str] = Counter()
    for row in clade_rows:
        regime = row.get("regime", "")
        if regimes and regime not in regimes:
            continue
        source_count_by_regime[regime] += 1
        for tip in split_tips(row.get(f"{tip_set}_tips", "")):
            if tip in entry_by_id:
                tips_by_regime[regime].add(tip)
            else:
                missing_by_regime[regime] += 1

    merged_rows: list[dict[str, str]] = []
    regime_order = {"psychro": 0, "meso": 1, "thermo": 2}
    for regime in sorted(tips_by_regime, key=lambda item: (regime_order.get(item, 99), item)):
        tips = sorted(tips_by_regime[regime])
        if len(tips) < min_sequences:
            continue
        merged_rows.append(
            {
                "clade_id": f"{regime}_merged",
                "regime": regime,
                f"{tip_set}_tips": ",".join(tips),
                "n_total_tips": str(len(tips)),
                "n_labelled": str(len(tips)) if tip_set == "labelled" else "",
                "n_unlabelled": str(len(tips)) if tip_set == "unlabelled" else "",
                "source_clades": str(source_count_by_regime[regime]),
                "missing_tips": str(missing_by_regime[regime]),
            }
        )
    return merged_rows


def numeric_value(row: dict[str, object], field: str) -> float | None:
    try:
        value = row.get(field, "")
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def summarize_core_metrics(all_stats: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    query_by_qpos: dict[str, dict[str, object]] = {}
    for stat in all_stats:
        qpos = str(stat.get("qpos", ""))
        if qpos and qpos not in query_by_qpos:
            query_by_qpos[qpos] = stat

    def summarize(label: str, regime: str, stats: list[dict[str, object]], aa_field: str, category_field: str) -> None:
        values = []
        for stat in stats:
            aa = str(stat.get(aa_field, ""))
            if aa not in AA_SET or str(stat.get(category_field, "")) != "core":
                continue
            hydro = HYDROPATHY.get(aa)
            volume = AA_VOLUME.get(aa)
            if hydro is None or volume is None:
                continue
            values.append((stat, hydro, volume))
        if not values:
            rows.append(
                {
                    "group_id": label,
                    "regime": regime,
                    "basis": category_field,
                    "n_core_observations": 0,
                    "n_unique_qpos": 0,
                    "mean_hydropathy": "",
                    "mean_volume": "",
                    "mean_sasa_a2": "",
                    "mean_rsa_pct": "",
                }
            )
            return
        sasa_field = "query_sasa_a2" if aa_field == "query_aa" else "consensus_sasa_a2"
        rsa_field = "query_rsa_pct" if aa_field == "query_aa" else "consensus_rsa_pct"
        sasa_values = [value for stat, _hydro, _volume in values if (value := numeric_value(stat, sasa_field)) is not None]
        rsa_values = [value for stat, _hydro, _volume in values if (value := numeric_value(stat, rsa_field)) is not None]
        rows.append(
            {
                "group_id": label,
                "regime": regime,
                "basis": category_field,
                "n_core_observations": len(values),
                "n_unique_qpos": len({str(stat.get("qpos", "")) for stat, _hydro, _volume in values if stat.get("qpos", "")}),
                "mean_hydropathy": f"{sum(hydro for _stat, hydro, _volume in values) / len(values):.4f}",
                "mean_volume": f"{sum(volume for _stat, _hydro, volume in values) / len(values):.4f}",
                "mean_sasa_a2": f"{sum(sasa_values) / len(sasa_values):.4f}" if sasa_values else "",
                "mean_rsa_pct": f"{sum(rsa_values) / len(rsa_values):.4f}" if rsa_values else "",
            }
        )

    summarize("ANKros", "ankros", list(query_by_qpos.values()), "query_aa", "query_category")

    by_regime: dict[str, list[dict[str, object]]] = defaultdict(list)
    by_clade: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for stat in all_stats:
        regime = str(stat.get("regime", ""))
        clade_id = str(stat.get("clade_id", ""))
        by_regime[regime].append(stat)
        by_clade[(regime, clade_id)].append(stat)

    regime_order = {"psychro": 0, "meso": 1, "thermo": 2}
    for regime in sorted(by_regime, key=lambda item: (regime_order.get(item, 99), item)):
        summarize(f"{regime}_consensus", regime, by_regime[regime], "top_aa", "category")
    for regime, clade_id in sorted(by_clade, key=lambda item: (regime_order.get(item[0], 99), item[1])):
        summarize(clade_id, regime, by_clade[(regime, clade_id)], "top_aa", "category")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, default=DEFAULT_ALIGNMENT)
    parser.add_argument("--column-map", type=Path, default=DEFAULT_COLMAP)
    parser.add_argument("--clades", type=Path, default=DEFAULT_CLADES)
    parser.add_argument("--cofactor-contacts", type=Path, default=DEFAULT_CONTACTS)
    parser.add_argument("--residue-annotations", type=Path, default=DEFAULT_RESIDUE_ANNOTATIONS)
    parser.add_argument("--pocket-conservation", type=Path, default=DEFAULT_POCKET)
    parser.add_argument("--et-chain", type=Path, default=DEFAULT_ET)
    parser.add_argument("--afdb-dir", type=Path, default=DEFAULT_AFDB_DIR)
    parser.add_argument("--ankros-structure", type=Path, default=DEFAULT_ANKROS_STRUCTURE)
    parser.add_argument(
        "--merged-residue-scores",
        type=Path,
        default=DEFAULT_MERGED_RESIDUE_SCORES,
        help=(
            "Merged per-residue score table from 15_compute_solubility_aggregability.py. "
            "If present, this is used before the separate CamSol/Aggrescan3D tables."
        ),
    )
    parser.add_argument("--camsol-residue-scores", type=Path, default=DEFAULT_CAMSOL_RESIDUES)
    parser.add_argument("--aggrescan3d-residue-scores", type=Path, default=DEFAULT_AGGRESCAN3D_RESIDUES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--html", action="store_true", help="Write an interactive HTML view instead of PNG figures.")
    parser.add_argument("--html-output", type=Path, default=None, help="Output path for --html mode.")
    parser.add_argument(
        "--html-default-visible",
        type=int,
        default=2,
        help="Initial number of clades shown per regime in the HTML control.",
    )
    parser.add_argument("--tip-set", choices=["labelled", "all", "unlabelled"], default="labelled")
    parser.add_argument("--domains", default="antenna,linker,catalytic")
    parser.add_argument("--regimes", default="psychro,meso,thermo")
    parser.add_argument(
        "--max-clades-per-regime",
        type=int,
        default=2,
        help="Select this many largest clades per regime; use 0 to include all clades.",
    )
    parser.add_argument(
        "--merge-regime-clades",
        action="store_true",
        help=(
            "Build one consensus/logo per regime from the union of all MSA sequences in that regime. "
            "When set, --max-clades-per-regime is ignored."
        ),
    )
    parser.add_argument("--min-sequences", type=int, default=2)
    parser.add_argument("--wrap", type=int, default=70)
    parser.add_argument("--max-insert-run", type=int, default=20)
    parser.add_argument("--exposed-rsa-pct", type=float, default=20.0)
    parser.add_argument("--contact-cutoff-a", type=float, default=4.0)
    parser.add_argument(
        "--sidechain-contact-cutoff-a",
        type=float,
        default=4.5,
        help="Heavy-atom side-chain distance cutoff for the ANKros SC contacts HTML row.",
    )
    parser.add_argument("--et-fad-cutoff-a", type=float, default=6.0)
    parser.add_argument("--include-empty-inserts", action="store_true")
    parser.add_argument("--skip-per-clade", action="store_true")
    parser.add_argument("--skip-domain-overview", action="store_true")
    parser.add_argument("--skip-metric-tracks", action="store_true")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for per-AFDB-structure feature extraction.",
    )
    parser.add_argument(
        "--metric-smooth-window",
        type=int,
        default=5,
        help="Centered rolling window, in ANKros positions, for hydrophobicity/volume line tracks; use 1 for raw values.",
    )
    args = parser.parse_args()

    os.environ.setdefault("MPLCONFIGDIR", str(args.out_dir / ".matplotlib"))
    validate_inputs(
        [
            (args.alignment, "alignment"),
            (args.column_map, "column map"),
            (args.clades, "clade table"),
            (args.ankros_structure, "ANKros structure"),
            (args.afdb_dir, "AFDB structure directory"),
            (args.residue_annotations, "holo residue annotations"),
        ]
    )
    if args.min_sequences < 1:
        raise SystemExit("ERROR: --min-sequences must be >= 1")
    if args.wrap < 10:
        raise SystemExit("ERROR: --wrap must be >= 10")
    if args.max_insert_run < 0:
        raise SystemExit("ERROR: --max-insert-run must be >= 0")
    if args.max_clades_per_regime < 0:
        raise SystemExit("ERROR: --max-clades-per-regime must be >= 0")
    if args.metric_smooth_window < 1:
        raise SystemExit("ERROR: --metric-smooth-window must be >= 1")
    if args.html_default_visible < 1:
        raise SystemExit("ERROR: --html-default-visible must be >= 1")
    if args.sidechain_contact_cutoff_a <= 0:
        raise SystemExit("ERROR: --sidechain-contact-cutoff-a must be > 0")

    domains = [part.strip() for part in args.domains.split(",") if part.strip()]
    unknown_domains = [domain for domain in domains if domain not in DOMAINS]
    if unknown_domains:
        raise SystemExit(f"ERROR: unknown domain(s): {', '.join(unknown_domains)}")
    regimes = {part.strip() for part in args.regimes.split(",") if part.strip()}

    entries = read_fasta(args.alignment)
    if not entries:
        raise SystemExit("ERROR: empty alignment")
    _query_id, _query_header, query_seq = find_query(entries)
    colmap = read_tsv(args.column_map)
    if len(query_seq) != len(colmap):
        raise SystemExit(f"ERROR: query width {len(query_seq)} != column map rows {len(colmap)}")

    entry_by_id = {sid: (sid, header, seq) for sid, header, seq in entries}
    seq_pos_by_sid = {sid: alignment_sequence_positions(seq) for sid, _header, seq in entries}
    extra_catalytic_sites, catalytic_evidence = load_catalytic_sites(
        args.cofactor_contacts,
        args.residue_annotations,
        args.pocket_conservation,
        args.et_chain,
        args.contact_cutoff_a,
        args.et_fad_cutoff_a,
    )
    ankros_features, ankros_structure_evidence = load_ankros_structure_features(
        args.ankros_structure,
        args.exposed_rsa_pct,
        args.contact_cutoff_a,
        args.sidechain_contact_cutoff_a,
        extra_catalytic_sites,
        catalytic_evidence,
    )
    residue_score_cache, score_warnings = load_residue_score_cache(
        args.camsol_residue_scores,
        args.aggrescan3d_residue_scores,
        args.merged_residue_scores,
    )
    scored_structure_ids = set(residue_score_cache)
    for warning in score_warnings:
        print(f"WARNING: {warning}")
    catalytic_evidence = ankros_structure_evidence
    source_clade_rows = read_tsv(args.clades)
    if args.merge_regime_clades:
        clade_rows = merge_clades_by_regime(
            source_clade_rows,
            regimes,
            args.tip_set,
            entry_by_id,
            args.min_sequences,
        )
    else:
        clade_rows = select_representative_clades(
            source_clade_rows,
            regimes,
            args.tip_set,
            entry_by_id,
            args.max_clades_per_regime,
            args.min_sequences,
        )
    selected_sids = {
        tip
        for clade in clade_rows
        for tip in split_tips(clade.get(f"{args.tip_set}_tips", ""))
        if tip in entry_by_id
    }
    selected_entries_for_structures = [entry_by_id[sid] for sid in sorted(selected_sids)]
    afdb_features, afdb_structure_paths, structure_pos_by_entry_pos_by_sid = load_afdb_feature_cache(
        selected_entries_for_structures,
        args.afdb_dir,
        args.exposed_rsa_pct,
        args.contact_cutoff_a,
        args.sidechain_contact_cutoff_a,
        args.workers,
    )

    all_stats: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    domain_grouped: dict[str, list[tuple[str, str, int, list[dict[str, object]]]]] = {domain: [] for domain in domains}
    html_grouped: list[tuple[str, str, int, list[dict[str, object]]]] = []
    selected_cols_by_domain = {
        domain: select_domain_columns(colmap, domain, args.max_insert_run)
        for domain in domains
    }
    selected_cols_for_html = select_all_columns(colmap, args.max_insert_run) if args.html else []
    total_conservation_for_html = (
        total_msa_conservation_stats(entries, selected_cols_for_html)
        if args.html
        else []
    )
    html_out = args.html_output or (args.out_dir / "thermal_clade_consensus_logos.html")

    for clade in clade_rows:
        regime = clade.get("regime", "")
        if regimes and regime not in regimes:
            continue
        clade_id = clade.get("clade_id", "")
        tip_text = clade.get(f"{args.tip_set}_tips", "")
        tips = split_tips(tip_text)
        selected_entries = [entry_by_id[tip] for tip in tips if tip in entry_by_id]
        missing = len([tip for tip in tips if tip not in entry_by_id])
        if len(selected_entries) < args.min_sequences:
            run_rows.append(
                {
                    "clade_id": clade_id,
                    "regime": regime,
                    "status": "skipped_too_few_sequences",
                    "tip_set": args.tip_set,
                    "n_sequences": len(selected_entries),
                    "n_missing_tips": missing,
                    "source_clades": clade.get("source_clades", ""),
                    "domain": "",
                    "n_columns": "",
                    "logo_png": "",
                }
            )
            continue

        if args.html:
            representative_entry = choose_representative_entry(
                selected_entries,
                selected_cols_for_html,
                scored_structure_ids,
            )
            stats = build_column_stats(
                "all",
                clade_id,
                regime,
                selected_entries,
                selected_cols_for_html,
                query_seq,
                seq_pos_by_sid,
                afdb_features,
                structure_pos_by_entry_pos_by_sid,
                ankros_features,
                residue_score_cache,
                representative_entry,
                True,
            )
            all_stats.extend(stats)
            html_grouped.append((clade_id, regime, len(selected_entries), stats))
            run_rows.append(
                {
                    "clade_id": clade_id,
                    "regime": regime,
                    "status": "ok",
                    "tip_set": args.tip_set,
                    "n_sequences": len(selected_entries),
                    "n_missing_tips": missing,
                    "source_clades": clade.get("source_clades", ""),
                    "domain": "all",
                    "n_columns": len(stats),
                    "logo_png": str(html_out),
                }
            )
            continue

        for domain in domains:
            representative_entry = choose_representative_entry(
                selected_entries,
                selected_cols_by_domain[domain],
                scored_structure_ids,
            )
            stats = build_column_stats(
                domain,
                clade_id,
                regime,
                selected_entries,
                selected_cols_by_domain[domain],
                query_seq,
                seq_pos_by_sid,
                afdb_features,
                structure_pos_by_entry_pos_by_sid,
                ankros_features,
                residue_score_cache,
                representative_entry,
                args.include_empty_inserts,
            )
            all_stats.extend(stats)
            domain_grouped[domain].append((clade_id, regime, len(selected_entries), stats))
            png = args.out_dir / "by_clade" / f"{clade_id}_{regime}" / f"{domain}_logo.png"
            if not args.skip_per_clade:
                draw_clade_logo(domain, clade_id, regime, len(selected_entries), stats, png, args.wrap)
            run_rows.append(
                {
                    "clade_id": clade_id,
                    "regime": regime,
                    "status": "ok",
                    "tip_set": args.tip_set,
                    "n_sequences": len(selected_entries),
                    "n_missing_tips": missing,
                    "source_clades": clade.get("source_clades", ""),
                    "domain": domain,
                    "n_columns": len(stats),
                    "logo_png": str(png) if not args.skip_per_clade else "",
                }
            )

    if args.html:
        keep_html_columns = html_nonempty_column_indices(html_grouped)
        html_grouped = [
            (cid, regime, nseq, [stats[idx] for idx in keep_html_columns])
            for cid, regime, nseq, stats in html_grouped
        ]
        total_conservation_for_html = filter_html_columns_by_indices(total_conservation_for_html, keep_html_columns)
        all_stats = [stat for _cid, _regime, _nseq, stats in html_grouped for stat in stats]
        n_columns_by_clade = {cid: len(stats) for cid, _regime, _nseq, stats in html_grouped}
        for row in run_rows:
            clade_id = str(row.get("clade_id", ""))
            if clade_id in n_columns_by_clade:
                row["n_columns"] = n_columns_by_clade[clade_id]
        write_interactive_html(html_out, html_grouped, total_conservation_for_html, args.html_default_visible)
    elif not args.skip_domain_overview:
        for domain, grouped in domain_grouped.items():
            draw_domain_overview(domain, grouped, args.out_dir / "by_domain" / f"{domain}_clades.png", args.wrap)
    if not args.html and not args.skip_metric_tracks:
        for domain, grouped in domain_grouped.items():
            draw_metric_tracks(
                domain,
                grouped,
                args.out_dir / "by_domain" / f"{domain}_hydropathy_volume_tracks.png",
                args.metric_smooth_window,
            )

    stats_fields = [
        "clade_id",
        "regime",
        "domain",
        "domain_col",
        "alignment_col",
        "qpos",
        "is_insert",
        "insert_run_len",
        "source_block",
        "source_range",
        "n_sequences",
        "n_residues",
        "n_gaps",
        "n_structures_for_top_aa",
        "n_matching_structure_residues",
        "consensus_sasa_a2",
        "consensus_rsa_pct",
        "occupancy",
        "gap_fraction",
        "information_bits",
        "top_aa",
        "top_aa_count",
        "top_aa_frequency_among_residues",
        "top_symbol_frequency_among_sequences",
        "representative_id",
        "representative_structure_id",
        "representative_structure_pos",
        "representative_aa",
        "representative_secondary_structure",
        "representative_matches_consensus",
        "query_aa",
        "query_secondary_structure",
        "query_functional_evidence",
        "query_matches_consensus",
        "query_category",
        "query_sasa_a2",
        "query_rsa_pct",
        "query_aa_volume",
        "query_aa_hydropathy",
        "query_sidechain_contact_count",
        "query_sidechain_contact_residues",
        "ankros_camsol_intrinsic",
        "ankros_camsol_structural",
        "ankros_aggrescan3d",
        "ankros_aggrescan3d_positive",
        "representative_camsol_intrinsic",
        "representative_camsol_structural",
        "representative_aggrescan3d",
        "representative_aggrescan3d_positive",
        "category",
        "aa_volume",
        "aa_hydropathy",
    ]
    run_fields = [
        "clade_id",
        "regime",
        "status",
        "tip_set",
        "n_sequences",
        "n_missing_tips",
        "source_clades",
        "domain",
        "n_columns",
        "logo_png",
    ]
    write_tsv(args.out_dir / "column_stats.tsv", all_stats, stats_fields)
    if args.html:
        write_tsv(
            args.out_dir / "total_msa_conservation.tsv",
            total_conservation_for_html,
            [
                "domain_col",
                "alignment_col",
                "qpos",
                "is_insert",
                "insert_run_len",
                "source_block",
                "source_range",
                "n_sequences",
                "n_residues",
                "n_gaps",
                "occupancy",
                "gap_fraction",
                "information_bits",
                "top_aa",
                "top_aa_count",
                "top_aa_frequency_among_residues",
                "top_symbol_frequency_among_sequences",
                "conservation",
            ],
        )
    write_tsv(args.out_dir / "run_summary.tsv", run_rows, run_fields)
    write_tsv(
        args.out_dir / "core_metric_summary.tsv",
        summarize_core_metrics(all_stats),
        [
            "group_id",
            "regime",
            "basis",
            "n_core_observations",
            "n_unique_qpos",
            "mean_hydropathy",
            "mean_volume",
            "mean_sasa_a2",
            "mean_rsa_pct",
        ],
    )

    evidence_rows = [
        {"qpos": resnum, "evidence": ";".join(sorted(labels))}
        for resnum, labels in sorted(catalytic_evidence.items())
    ]
    write_tsv(args.out_dir / "catalytic_site_evidence.tsv", evidence_rows, ["qpos", "evidence"])
    structure_rows = [
        {
            "id": sid,
            "structure_pdb": path,
            "n_featured_residues": len(afdb_features.get(sid, {})),
            "n_mapped_sequence_positions": len(structure_pos_by_entry_pos_by_sid.get(sid, {})),
        }
        for sid, path in sorted(afdb_structure_paths.items())
    ]
    write_tsv(
        args.out_dir / "afdb_structure_support.tsv",
        structure_rows,
        ["id", "structure_pdb", "n_featured_residues", "n_mapped_sequence_positions"],
    )

    n_ok = sum(1 for row in run_rows if row.get("status") == "ok")
    n_skipped = sum(1 for row in run_rows if row.get("status") != "ok")
    print(f"Saved domain/clade logo rows: {n_ok}")
    if n_skipped:
        print(f"Skipped clades/domains: {n_skipped}")
    print(f"Saved: {args.out_dir / 'column_stats.tsv'}")
    print(f"Saved: {args.out_dir / 'run_summary.tsv'}")
    print(f"Saved: {args.out_dir / 'core_metric_summary.tsv'}")
    if args.html:
        print(f"Saved: {html_out}")


if __name__ == "__main__":
    main()
