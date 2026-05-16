#!/usr/bin/env python3
"""Step 17: classify homologues by extension of the ANKros 211-220 helix.

The classifier is intentionally local and auditable:

* ANKros qpos 211-220 define the core helix.
* ANKros qpos upstream of that core, plus intervening insertion columns, define
  the candidate extension region.
* Secondary structure is assigned from the available structure model, using
  mkdssp/DSSP when available and Biotite P-SEA as a fallback.

This does not claim evolutionary causality. It classifies AFDB-backed
representatives that can be mapped back to the refined MSA.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "results" / "msa_OGT" / "helix_extension_classifier"
os.environ.setdefault("MPLCONFIGDIR", str(DEFAULT_OUT / ".matplotlib"))

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import biotite.structure as bstruc
import biotite.structure.io as bstrucio


DEFAULT_ALIGNMENT = ROOT / "results" / "msa_OGT" / "repset_hmmalign_linker_refined.fa"
DEFAULT_COLMAP = ROOT / "results" / "msa_OGT" / "repset_hmmalign_linker_refined_column_map.tsv"
DEFAULT_METADATA = ROOT / "results" / "msa_OGT" / "regime_clades" / "regime_clades_tip_metadata.tsv"
DEFAULT_AFDB_DIR = ROOT / "structures" / "afdb"
DEFAULT_ANKROS_STRUCTURE = (
    ROOT
    / "results"
    / "structural"
    / "docked_holo"
    / "ankros_fad_fmn_donor_holo.pdb"
)
ANKROS_ID = "photoHymenobact"
EXTENSION_PARTIAL_CALLS = {"partial_extension"}

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)
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
    "SEC": "U",
    "PYL": "O",
}


@dataclass(frozen=True)
class ResidueFeature:
    seq_pos: int
    aa: str
    chain: str
    resseq: int
    icode: str
    ss: str
    plddt: float | None


@dataclass(frozen=True)
class LocalResidue:
    entry_pos: int
    structure_pos: int
    aln_col: int
    qpos: int | None
    region: str
    aa: str
    ss: str
    plddt: float | None


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    entries: list[tuple[str, str, str]] = []
    header: str | None = None
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


def parse_range(text: str) -> tuple[int, int]:
    if "-" in text:
        left, right = text.split("-", 1)
        return int(left), int(right)
    value = int(text)
    return value, value


def parse_float(text: str) -> float | None:
    try:
        if text == "":
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def fmt_float(value: float | None, ndigits: int = 3) -> str:
    return "" if value is None else f"{value:.{ndigits}f}"


def extension_signal(row: dict[str, object], min_full_run: int, min_partial_run: int) -> str:
    call_text = str(row.get("call", "") or "")
    try:
        extension_run_len = int(row.get("extension_run_len", 0) or 0)
    except (TypeError, ValueError):
        extension_run_len = 0
    if extension_run_len >= min_full_run:
        return "positive"
    if extension_run_len >= min_partial_run or call_text in EXTENSION_PARTIAL_CALLS:
        return "partial"
    if call_text in {"missing_structure", "structure_error", "missing_core_region"}:
        return "missing"
    return "negative"


def is_extension_positive(row: dict[str, object]) -> bool:
    return str(row.get("extension_signal", "")) == "positive"


def is_extension_partial(row: dict[str, object]) -> bool:
    return str(row.get("extension_signal", "")) == "partial"


def legacy_extension_signal(call: object) -> str:
    call_text = str(call or "")
    if call_text == "extended_helix":
        return "positive"
    if call_text in EXTENSION_PARTIAL_CALLS:
        return "partial"
    if call_text in {"missing_structure", "structure_error", "missing_core_region"}:
        return "missing"
    return "negative"


def alignment_sequence_positions(seq: str) -> dict[int, int]:
    positions: dict[int, int] = {}
    seq_pos = 0
    for aln_col, aa in enumerate(seq):
        if aa.upper() in AA_SET:
            seq_pos += 1
            positions[aln_col] = seq_pos
    return positions


def afdb_path_for_sid(sid: str, afdb_dir: Path) -> Path | None:
    if sid == ANKROS_ID or sid.startswith(ANKROS_ID):
        return None
    if sid.startswith("UniRef90_"):
        accession = sid.removeprefix("UniRef90_")
        path = afdb_dir / f"AF-{accession}-F1-model_v6.pdb"
        return path if path.exists() else None
    return None


def structure_path_for_sid(sid: str, afdb_dir: Path, ankros_structure: Path) -> Path | None:
    if sid == ANKROS_ID or sid.startswith(ANKROS_ID):
        return ankros_structure if ankros_structure.exists() else None
    return afdb_path_for_sid(sid, afdb_dir)


def protein_residues_from_structure(path: Path) -> list:
    parser = PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        structure = parser.get_structure(path.stem, path)
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                aa = THREE_TO_ONE.get(residue.get_resname().strip().upper())
                if aa in AA_SET and "CA" in residue:
                    residues.append(residue)
        break
    return residues


def find_mkdssp_executable() -> str | None:
    env_bin = Path(sys.executable).resolve().parent / "mkdssp"
    if env_bin.exists():
        return str(env_bin)
    return None


def dssp_secondary_structure_from_pdb(path: Path) -> dict[tuple[str, int, str], str]:
    mkdssp = find_mkdssp_executable()
    if mkdssp is None:
        return {}
    parser = PDBParser(QUIET=True)
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

    by_residue: dict[tuple[str, int, str], str] = {}
    for (chain_id, residue_id), values in dssp.property_dict.items():
        hetflag, resseq, icode = residue_id
        if hetflag != " ":
            continue
        by_residue[(str(chain_id), int(resseq), str(icode).strip())] = state_map.get(str(values[2]), "C")
    return by_residue


def psea_secondary_structure_from_pdb(path: Path) -> dict[tuple[str, int, str], str]:
    state_map = {"a": "H", "b": "E", "c": "C"}
    by_residue: dict[tuple[str, int, str], str] = {}
    atom_array = bstrucio.load_structure(str(path))
    if type(atom_array).__name__ == "AtomArrayStack":
        atom_array = atom_array[0]
    protein = atom_array[bstruc.filter_amino_acids(atom_array)]
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


def load_structure_features(path: Path) -> dict[int, ResidueFeature]:
    residues = protein_residues_from_structure(path)
    ss_by_residue = dssp_secondary_structure_from_pdb(path) or psea_secondary_structure_from_pdb(path)
    raw_plddt = [float(residue["CA"].get_bfactor()) for residue in residues if "CA" in residue]
    bfactor_is_plddt = bool(raw_plddt and max(raw_plddt) > 1.0)
    features: dict[int, ResidueFeature] = {}
    for seq_pos, residue in enumerate(residues, start=1):
        chain = str(residue.get_parent().id)
        resseq = int(residue.id[1])
        icode = str(residue.id[2]).strip()
        aa = THREE_TO_ONE.get(residue.get_resname().strip().upper(), "")
        plddt = float(residue["CA"].get_bfactor()) if bfactor_is_plddt and "CA" in residue else None
        features[seq_pos] = ResidueFeature(
            seq_pos=seq_pos,
            aa=aa,
            chain=chain,
            resseq=resseq,
            icode=icode,
            ss=ss_by_residue.get((chain, resseq, icode), "C"),
            plddt=plddt,
        )
    return features


def map_entry_positions_to_structure(
    alignment_seq: str,
    structure_features: dict[int, ResidueFeature],
) -> tuple[dict[int, int], float]:
    entry_seq = "".join(aa.upper() for aa in alignment_seq if aa.upper() in AA_SET)
    if not entry_seq or not structure_features:
        return {}, 0.0

    structure_items = sorted(structure_features.items())
    structure_seq = "".join(feature.aa for _seq_pos, feature in structure_items)
    if not structure_seq:
        return {}, 0.0

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(entry_seq, structure_seq)[0]

    mapping: dict[int, int] = {}
    matches = 0
    aligned_pairs = 0
    for entry_span, structure_span in zip(alignment.aligned[0], alignment.aligned[1]):
        entry_start, entry_end = int(entry_span[0]), int(entry_span[1])
        structure_start, structure_end = int(structure_span[0]), int(structure_span[1])
        block_len = min(entry_end - entry_start, structure_end - structure_start)
        for offset in range(block_len):
            entry_pos = entry_start + offset + 1
            structure_index = structure_start + offset
            aligned_pairs += 1
            if entry_seq[entry_pos - 1] != structure_seq[structure_index]:
                continue
            matches += 1
            mapping[entry_pos] = structure_items[structure_index][0]
    denominator = max(1, min(len(entry_seq), len(structure_seq)))
    identity = matches / denominator
    if aligned_pairs == 0:
        identity = 0.0
    return mapping, identity


def columns_for_qpos_range(colmap: list[dict[str, str]], qstart: int, qend: int) -> set[int]:
    qpos_cols = [int(row["out_col"]) for row in colmap if row.get("qpos") and qstart <= int(row["qpos"]) <= qend]
    if not qpos_cols:
        return set()
    first_col = min(qpos_cols)
    last_col = max(qpos_cols)
    return {int(row["out_col"]) for row in colmap if first_col <= int(row["out_col"]) <= last_col}


def local_residues_for_region(
    seq: str,
    seqpos_by_col: dict[int, int],
    entry_to_structure: dict[int, int],
    structure_features: dict[int, ResidueFeature],
    colmap_by_col: dict[int, dict[str, str]],
    cols: set[int],
    region: str,
) -> list[LocalResidue]:
    residues: dict[int, LocalResidue] = {}
    for aln_col in sorted(cols):
        if aln_col >= len(seq):
            continue
        entry_pos = seqpos_by_col.get(aln_col)
        if entry_pos is None:
            continue
        structure_pos = entry_to_structure.get(entry_pos)
        if structure_pos is None:
            continue
        feature = structure_features.get(structure_pos)
        if feature is None:
            continue
        row = colmap_by_col.get(aln_col, {})
        qpos = int(row["qpos"]) if row.get("qpos") else None
        residues.setdefault(
            entry_pos,
            LocalResidue(
                entry_pos=entry_pos,
                structure_pos=structure_pos,
                aln_col=aln_col,
                qpos=qpos,
                region=region,
                aa=feature.aa,
                ss=feature.ss,
                plddt=feature.plddt,
            ),
        )
    return sorted(residues.values(), key=lambda residue: residue.entry_pos)


def helix_runs(residues: list[LocalResidue]) -> list[list[LocalResidue]]:
    runs: list[list[LocalResidue]] = []
    current: list[LocalResidue] = []
    previous_pos: int | None = None
    for residue in residues:
        if residue.ss == "H" and (previous_pos is None or residue.entry_pos == previous_pos + 1):
            current.append(residue)
        else:
            if current:
                runs.append(current)
            current = [residue] if residue.ss == "H" else []
        previous_pos = residue.entry_pos
    if current:
        runs.append(current)
    return runs


def nearest_upstream_run(
    upstream_runs: list[list[LocalResidue]],
    core_residues: list[LocalResidue],
    max_gap: int,
) -> list[LocalResidue]:
    if not upstream_runs:
        return []
    core_anchor_positions = [residue.entry_pos for residue in core_residues if residue.ss == "H"] or [
        residue.entry_pos for residue in core_residues
    ]
    if not core_anchor_positions:
        return max(upstream_runs, key=len)
    core_start = min(core_anchor_positions)
    candidates = [run for run in upstream_runs if core_start - run[-1].entry_pos - 1 <= max_gap]
    if not candidates:
        return []
    return max(candidates, key=lambda run: (len(run), run[-1].entry_pos))


def mean(values: list[float | None]) -> float | None:
    usable = [value for value in values if value is not None]
    return statistics.mean(usable) if usable else None


def summarize_region(residues: list[LocalResidue]) -> dict[str, object]:
    runs = helix_runs(residues)
    helix_residues = [residue for residue in residues if residue.ss == "H"]
    qpos_values = [residue.qpos for residue in residues if residue.qpos is not None]
    return {
        "mapped_residues": len(residues),
        "helix_residues": len(helix_residues),
        "helix_fraction": len(helix_residues) / len(residues) if residues else 0.0,
        "max_helix_run": max((len(run) for run in runs), default=0),
        "proline_count": sum(1 for residue in residues if residue.aa == "P"),
        "proline_fraction": sum(1 for residue in residues if residue.aa == "P") / len(residues) if residues else 0.0,
        "mean_plddt": mean([residue.plddt for residue in residues]),
        "qpos_span": f"{min(qpos_values)}-{max(qpos_values)}" if qpos_values else "",
        "entry_pos_span": f"{residues[0].entry_pos}-{residues[-1].entry_pos}" if residues else "",
        "structure_pos_span": f"{residues[0].structure_pos}-{residues[-1].structure_pos}" if residues else "",
    }


def classify_entry(
    seq: str,
    structure_features: dict[int, ResidueFeature],
    entry_to_structure: dict[int, int],
    mapping_identity: float,
    colmap_by_col: dict[int, dict[str, str]],
    upstream_cols: set[int],
    core_cols: set[int],
    args: argparse.Namespace,
) -> dict[str, object]:
    seqpos_by_col = alignment_sequence_positions(seq)
    upstream = local_residues_for_region(
        seq, seqpos_by_col, entry_to_structure, structure_features, colmap_by_col, upstream_cols, "upstream"
    )
    core = local_residues_for_region(
        seq, seqpos_by_col, entry_to_structure, structure_features, colmap_by_col, core_cols, "core"
    )
    upstream_summary = summarize_region(upstream)
    core_summary = summarize_region(core)
    upstream_runs = helix_runs(upstream)
    core_present = (
        int(core_summary["mapped_residues"]) >= args.min_core_mapped
        and int(core_summary["helix_residues"]) >= args.min_core_helix_residues
        and float(core_summary["helix_fraction"]) >= args.min_core_fraction
    )
    extension_run = nearest_upstream_run(upstream_runs, core, args.max_extension_gap)
    extension_len = len(extension_run)
    if not core:
        call = "missing_core_region"
    elif mapping_identity < args.min_mapping_identity:
        call = "low_mapping_identity"
    elif core_present and extension_len >= args.min_extension_run:
        call = "extended_helix"
    elif core_present and extension_len >= args.min_partial_extension_run:
        call = "partial_extension"
    elif core_present:
        call = "core_only"
    else:
        call = "no_core_helix"

    extension_qpos = [residue.qpos for residue in extension_run if residue.qpos is not None]
    extension_plddt = mean([residue.plddt for residue in extension_run])
    local_plddt = mean([residue.plddt for residue in upstream + core])
    quality_flags: list[str] = []
    if local_plddt is not None and local_plddt < args.min_local_plddt:
        quality_flags.append("low_local_plddt")
    if int(upstream_summary["mapped_residues"]) < args.min_upstream_mapped:
        quality_flags.append("low_upstream_mapping")
    if mapping_identity < args.min_mapping_identity:
        quality_flags.append("low_mapping_identity")

    return {
        "call": call,
        "quality_flags": ";".join(quality_flags),
        "mapping_identity": fmt_float(mapping_identity, 4),
        "core_mapped_residues": core_summary["mapped_residues"],
        "core_helix_residues": core_summary["helix_residues"],
        "core_helix_fraction": fmt_float(float(core_summary["helix_fraction"]), 3),
        "core_max_helix_run": core_summary["max_helix_run"],
        "core_qpos_span": core_summary["qpos_span"],
        "core_entry_pos_span": core_summary["entry_pos_span"],
        "core_structure_pos_span": core_summary["structure_pos_span"],
        "core_mean_plddt": fmt_float(core_summary["mean_plddt"], 2),
        "upstream_mapped_residues": upstream_summary["mapped_residues"],
        "upstream_helix_residues": upstream_summary["helix_residues"],
        "upstream_helix_fraction": fmt_float(float(upstream_summary["helix_fraction"]), 3),
        "upstream_max_helix_run": upstream_summary["max_helix_run"],
        "upstream_qpos_span": upstream_summary["qpos_span"],
        "upstream_entry_pos_span": upstream_summary["entry_pos_span"],
        "upstream_structure_pos_span": upstream_summary["structure_pos_span"],
        "upstream_mean_plddt": fmt_float(upstream_summary["mean_plddt"], 2),
        "upstream_proline_count": upstream_summary["proline_count"],
        "upstream_proline_fraction": fmt_float(float(upstream_summary["proline_fraction"]), 3),
        "extension_run_len": extension_len,
        "extension_qpos_span": f"{min(extension_qpos)}-{max(extension_qpos)}" if extension_qpos else "",
        "extension_entry_pos_span": f"{extension_run[0].entry_pos}-{extension_run[-1].entry_pos}" if extension_run else "",
        "extension_structure_pos_span": f"{extension_run[0].structure_pos}-{extension_run[-1].structure_pos}" if extension_run else "",
        "extension_sequence": "".join(residue.aa for residue in extension_run),
        "extension_mean_plddt": fmt_float(extension_plddt, 2),
        "local_mean_plddt": fmt_float(local_plddt, 2),
    }


def summarize_by_group(rows: list[dict[str, object]], group_field: str) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = str(row.get(group_field, "") or "unknown")
        grouped[key].append(row)
    summary_rows: list[dict[str, object]] = []
    for key, group in sorted(grouped.items()):
        calls = Counter(str(row.get("call", "")) for row in group)
        callable_rows = [row for row in group if row.get("call") not in {"missing_structure", "structure_error"}]
        callable_calls = Counter(str(row.get("call", "")) for row in callable_rows)
        n_callable = len(callable_rows)
        extended = callable_calls.get("extended_helix", 0)
        no_core = callable_calls.get("no_core_helix", 0)
        partial = callable_calls.get("partial_extension", 0)
        core_only = callable_calls.get("core_only", 0)
        positive = sum(1 for row in callable_rows if is_extension_positive(row))
        signal_partial = sum(1 for row in callable_rows if is_extension_partial(row))
        summary_rows.append(
            {
                group_field: key,
                "n_total": len(group),
                "n_callable": n_callable,
                "n_extension_positive": positive,
                "n_extension_partial_signal": signal_partial,
                "n_extended_helix": extended,
                "n_partial_extension": partial,
                "n_core_only": core_only,
                "n_no_core_helix": no_core,
                "n_missing_core_region": callable_calls.get("missing_core_region", 0),
                "extension_positive_fraction_callable": fmt_float(positive / n_callable if n_callable else None, 3),
                "extension_positive_or_partial_fraction_callable": fmt_float(
                    (positive + signal_partial) / n_callable if n_callable else None, 3
                ),
                "extended_fraction_callable": fmt_float(extended / n_callable if n_callable else None, 3),
                "extended_or_partial_fraction_callable": fmt_float(
                    (extended + partial) / n_callable if n_callable else None, 3
                ),
                "majority_call": calls.most_common(1)[0][0] if calls else "",
                "call_counts": ";".join(f"{call}:{count}" for call, count in sorted(calls.items())),
            }
        )
    return summary_rows


def classify_one_entry(
    item: tuple[str, str, str],
    meta: dict[str, str],
    afdb_dir: Path,
    ankros_structure: Path,
    colmap_by_col: dict[int, dict[str, str]],
    upstream_cols: set[int],
    core_cols: set[int],
    args: argparse.Namespace,
    structure_cache: dict[Path, dict[int, ResidueFeature]] | None = None,
) -> dict[str, object]:
    sid, _header, seq = item
    row: dict[str, object] = {
        "id": sid,
        "accession": meta.get("accession", ""),
        "organism": meta.get("organism", ""),
        "regime": meta.get("threshold_regime") or meta.get("regime", ""),
        "ogt": meta.get("ogt", ""),
        "regime_clade_id": meta.get("regime_clade_id") or meta.get("regime_clade_context_id", ""),
    }
    structure_path = structure_path_for_sid(sid, afdb_dir, ankros_structure)
    row["structure_path"] = (
        str(structure_path.relative_to(ROOT))
        if structure_path and structure_path.is_relative_to(ROOT)
        else str(structure_path or "")
    )
    if structure_path is None:
        row["call"] = "missing_structure"
        row["extension_signal"] = legacy_extension_signal(row["call"])
        return row
    try:
        if structure_cache is not None:
            if structure_path not in structure_cache:
                structure_cache[structure_path] = load_structure_features(structure_path)
            features = structure_cache[structure_path]
        else:
            features = load_structure_features(structure_path)
        entry_to_structure, mapping_identity = map_entry_positions_to_structure(seq, features)
        row.update(
            classify_entry(
                seq,
                features,
                entry_to_structure,
                mapping_identity,
                colmap_by_col,
                upstream_cols,
                core_cols,
                args,
            )
        )
    except Exception as exc:
        row["call"] = "structure_error"
        row["quality_flags"] = type(exc).__name__
    row["extension_signal"] = extension_signal(row, args.min_extension_run, args.min_partial_extension_run)
    return row


def classify_one_entry_worker(payload: tuple) -> dict[str, object]:
    return classify_one_entry(*payload, structure_cache=None)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, default=DEFAULT_ALIGNMENT)
    parser.add_argument("--column-map", type=Path, default=DEFAULT_COLMAP)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--afdb-dir", type=Path, default=DEFAULT_AFDB_DIR)
    parser.add_argument("--ankros-structure", type=Path, default=DEFAULT_ANKROS_STRUCTURE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--upstream-qpos", default="192-210")
    parser.add_argument("--core-qpos", default="211-220")
    parser.add_argument("--min-core-mapped", type=int, default=6)
    parser.add_argument("--min-core-helix-residues", type=int, default=5)
    parser.add_argument("--min-core-fraction", type=float, default=0.50)
    parser.add_argument("--min-upstream-mapped", type=int, default=8)
    parser.add_argument("--min-extension-run", type=int, default=6)
    parser.add_argument("--min-partial-extension-run", type=int, default=3)
    parser.add_argument("--max-extension-gap", type=int, default=5)
    parser.add_argument("--min-local-plddt", type=float, default=70.0)
    parser.add_argument("--min-mapping-identity", type=float, default=0.95)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker processes for per-entry structure classification. Default preserves serial behavior.",
    )
    parser.add_argument("--limit", type=int, default=0, help="Debug limit on alignment entries")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    entries = read_fasta(args.alignment)
    if args.limit:
        entries = entries[: args.limit]
    metadata = {row.get("id", ""): row for row in read_tsv(args.metadata)}
    colmap = read_tsv(args.column_map)
    colmap_by_col = {int(row["out_col"]): row for row in colmap}
    upstream_start, upstream_end = parse_range(args.upstream_qpos)
    core_start, core_end = parse_range(args.core_qpos)
    upstream_cols = columns_for_qpos_range(colmap, upstream_start, upstream_end)
    core_cols = columns_for_qpos_range(colmap, core_start, core_end)
    if not upstream_cols:
        raise SystemExit(f"ERROR: no alignment columns found for upstream qpos {args.upstream_qpos}")
    if not core_cols:
        raise SystemExit(f"ERROR: no alignment columns found for core qpos {args.core_qpos}")

    structure_cache: dict[Path, dict[int, ResidueFeature]] = {}
    if args.workers > 1:
        payloads = [
            (
                item,
                metadata.get(item[0], {}),
                args.afdb_dir,
                args.ankros_structure,
                colmap_by_col,
                upstream_cols,
                core_cols,
                args,
            )
            for item in entries
        ]
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            rows = list(executor.map(classify_one_entry_worker, payloads))
    else:
        rows = [
            classify_one_entry(
                item,
                metadata.get(item[0], {}),
                args.afdb_dir,
                args.ankros_structure,
                colmap_by_col,
                upstream_cols,
                core_cols,
                args,
                structure_cache,
            )
            for item in entries
        ]

    fields = [
        "id",
        "accession",
        "organism",
        "regime",
        "ogt",
        "regime_clade_id",
        "call",
        "extension_signal",
        "quality_flags",
        "mapping_identity",
        "structure_path",
        "core_mapped_residues",
        "core_helix_residues",
        "core_helix_fraction",
        "core_max_helix_run",
        "core_qpos_span",
        "core_entry_pos_span",
        "core_structure_pos_span",
        "core_mean_plddt",
        "upstream_mapped_residues",
        "upstream_helix_residues",
        "upstream_helix_fraction",
        "upstream_max_helix_run",
        "upstream_qpos_span",
        "upstream_entry_pos_span",
        "upstream_structure_pos_span",
        "upstream_mean_plddt",
        "upstream_proline_count",
        "upstream_proline_fraction",
        "extension_run_len",
        "extension_qpos_span",
        "extension_entry_pos_span",
        "extension_structure_pos_span",
        "extension_sequence",
        "extension_mean_plddt",
        "local_mean_plddt",
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / "helix_extension_calls.tsv", rows, fields)

    regime_rows = summarize_by_group(rows, "regime")
    write_tsv(
        args.out_dir / "regime_summary.tsv",
        regime_rows,
        [
            "regime",
            "n_total",
            "n_callable",
            "n_extension_positive",
            "n_extension_partial_signal",
            "n_extended_helix",
            "n_partial_extension",
            "n_core_only",
            "n_no_core_helix",
            "n_missing_core_region",
            "extension_positive_fraction_callable",
            "extension_positive_or_partial_fraction_callable",
            "extended_fraction_callable",
            "extended_or_partial_fraction_callable",
            "majority_call",
            "call_counts",
        ],
    )

    clade_rows = summarize_by_group(rows, "regime_clade_id")
    write_tsv(
        args.out_dir / "clade_summary.tsv",
        clade_rows,
        [
            "regime_clade_id",
            "n_total",
            "n_callable",
            "n_extension_positive",
            "n_extension_partial_signal",
            "n_extended_helix",
            "n_partial_extension",
            "n_core_only",
            "n_no_core_helix",
            "n_missing_core_region",
            "extension_positive_fraction_callable",
            "extension_positive_or_partial_fraction_callable",
            "extended_fraction_callable",
            "extended_or_partial_fraction_callable",
            "majority_call",
            "call_counts",
        ],
    )

    manifest = [
        "Step 17 helix extension classifier",
        f"alignment={args.alignment}",
        f"column_map={args.column_map}",
        f"metadata={args.metadata}",
        f"afdb_dir={args.afdb_dir}",
        f"ankros_structure={args.ankros_structure}",
        f"upstream_qpos={args.upstream_qpos}",
        f"core_qpos={args.core_qpos}",
        f"upstream_alignment_columns={min(upstream_cols)}-{max(upstream_cols)} ({len(upstream_cols)} columns)",
        f"core_alignment_columns={min(core_cols)}-{max(core_cols)} ({len(core_cols)} columns)",
        f"min_core_mapped={args.min_core_mapped}",
        f"min_core_helix_residues={args.min_core_helix_residues}",
        f"min_core_fraction={args.min_core_fraction}",
        f"min_extension_run={args.min_extension_run}",
        f"min_partial_extension_run={args.min_partial_extension_run}",
        f"max_extension_gap={args.max_extension_gap}",
        f"min_local_plddt={args.min_local_plddt}",
        f"min_mapping_identity={args.min_mapping_identity}",
        f"entries={len(entries)}",
        f"structures_loaded={len(structure_cache) if args.workers <= 1 else len({row.get('structure_path', '') for row in rows if row.get('structure_path')})}",
        f"workers={args.workers}",
        f"calls={dict(Counter(str(row.get('call', '')) for row in rows))}",
        "extension_signal_positive=extension_run_len >= min_extension_run",
    ]
    (args.out_dir / "run_manifest.txt").write_text("\n".join(manifest) + "\n")
    print(f"Wrote {args.out_dir / 'helix_extension_calls.tsv'}")
    print(f"Wrote {args.out_dir / 'regime_summary.tsv'}")
    print(f"Wrote {args.out_dir / 'clade_summary.tsv'}")


if __name__ == "__main__":
    main()
