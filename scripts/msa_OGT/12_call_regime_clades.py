#!/usr/bin/env python3
"""
Step 12: call thermal-regime clades on the final OGT-aware MSA tree.

The clade call uses only labelled descendants for the regime decision. Tips
without an OGT/regime are retained as phylogenetic context, but they do not
dilute the target-regime fraction and are not assigned to regime clades unless
--assign-unlabelled is used.
"""
from __future__ import annotations

import argparse
import copy
import csv
import math
import sys
from collections import Counter
from pathlib import Path

from Bio import Phylo
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
import numpy as np

from importlib import import_module

MSA_OGT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MSA_OGT_DIR))
cfg = import_module("00_config")


TREE = cfg.INTER_DIR / "tree_antenna_catalytic" / "ankros_antenna_catalytic.treefile"
METADATA = cfg.INTER_DIR / "repset_metadata_qc.tsv"
OUTGROUP_IDS = cfg.INTER_DIR / "tree_antenna_catalytic" / "antenna_catalytic_outgroup_ids.txt"
OUT_DIR = cfg.INTER_DIR / "regime_clades"
ALN_FA = cfg.INTER_DIR / "repset_hmmalign_linker_refined.fa"
COL_MAP = cfg.INTER_DIR / "repset_hmmalign_linker_refined_column_map.tsv"
AFDB_DIR = cfg.PROJECT_ROOT / "structures" / "afdb"

REGIME_ORDER = {"psychro": 0, "meso": 1, "thermo": 2}
AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
DOMAINS = {
    "antenna": tuple(range(1, 131)),
    "linker": tuple(range(131, 206)),
    "catalytic": tuple(range(206, 438)),
}
GEOMETRY_MIN_COMMON = {
    "antenna": 70,
    "linker": 25,
    "catalytic": 140,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_float(text: str | None):
    text = (text or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def regime_from_ogt(ogt: float | None, psychro_max: float, thermo_min: float) -> str:
    if ogt is None:
        return ""
    if ogt < psychro_max:
        return "psychro"
    if ogt >= thermo_min:
        return "thermo"
    return "meso"


def load_metadata(path: Path, psychro_max: float, thermo_min: float):
    rows = []
    by_id = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        for row in reader:
            sid = row.get("id", "")
            ogt = parse_float(row.get("ogt"))
            threshold_regime = regime_from_ogt(ogt, psychro_max, thermo_min)
            row["_threshold_regime"] = threshold_regime
            rows.append(row)
            if sid:
                by_id[sid] = row
    return rows, by_id, fields


def load_outgroup_ids(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    sid = ""
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if sid:
                    seqs[sid] = "".join(parts)
                sid = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if sid:
        seqs[sid] = "".join(parts)
    return seqs


def read_colmap(path: Path) -> dict[int, int]:
    qcol = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("qpos"):
                qcol[int(row["qpos"])] = int(row["out_col"])
    return qcol


def sequence_positions(alignment_seq: str) -> dict[int, int]:
    out = {}
    seq_pos = 0
    for aln_col, aa in enumerate(alignment_seq):
        if aa.upper() in AA_SET:
            seq_pos += 1
            out[aln_col] = seq_pos
    return out


def afdb_path_for_sid(sid: str, afdb_dir: Path) -> Path | None:
    if not sid.startswith("UniRef90_"):
        return None
    accession = sid.removeprefix("UniRef90_")
    path = afdb_dir / f"AF-{accession}-F1-model_v6.pdb"
    return path if path.exists() else None


def structure_items(path: Path) -> list[dict[str, object]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, path)
    residues = []
    for model in structure:
        for chain in model:
            chain_residues = [
                residue for residue in chain
                if residue.id[0] == " "
                and THREE_TO_ONE.get(residue.get_resname().strip().upper()) in AA_SET
                and "CA" in residue
            ]
            if len(chain_residues) > len(residues):
                residues = chain_residues
        break
    items = []
    for seq_pos, residue in enumerate(residues, start=1):
        ca = residue["CA"]
        items.append({
            "seq_pos": seq_pos,
            "aa": THREE_TO_ONE[residue.get_resname().strip().upper()],
            "plddt": float(ca.get_bfactor()),
            "xyz": np.array(ca.get_coord(), dtype=float),
        })
    return items


def map_entry_positions_to_structure(alignment_seq: str, items: list[dict[str, object]]) -> dict[int, int]:
    entry_seq = "".join(aa.upper() for aa in alignment_seq if aa.upper() in AA_SET)
    structure_seq = "".join(str(item["aa"]) for item in items)
    if not entry_seq or not structure_seq:
        return {}
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(entry_seq, structure_seq)[0]
    out = {}
    for entry_span, structure_span in zip(alignment.aligned[0], alignment.aligned[1]):
        entry_start, entry_end = int(entry_span[0]), int(entry_span[1])
        structure_start, structure_end = int(structure_span[0]), int(structure_span[1])
        for offset in range(min(entry_end - entry_start, structure_end - structure_start)):
            entry_pos = entry_start + offset + 1
            structure_idx = structure_start + offset
            if entry_seq[entry_pos - 1] == structure_seq[structure_idx]:
                out[entry_pos] = int(items[structure_idx]["seq_pos"])
    return out


def kabsch_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    left_centroid = left.mean(axis=0)
    right_centroid = right.mean(axis=0)
    left_c = left - left_centroid
    right_c = right - right_centroid
    u, _s, vt = np.linalg.svd(left_c.T @ right_c)
    rot = u @ vt
    if np.linalg.det(rot) < 0.0:
        vt[-1] *= -1.0
        rot = u @ vt
    aligned = left_c @ rot + right_centroid
    return float(np.sqrt(np.mean(np.sum((aligned - right) ** 2, axis=1))))


def summarize_float(values: list[float]) -> tuple[str, str]:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return "", ""
    mean = sum(clean) / len(clean)
    mid = len(clean) // 2
    median = clean[mid] if len(clean) % 2 else (clean[mid - 1] + clean[mid]) / 2.0
    return f"{mean:.4f}", f"{median:.4f}"


def build_domain_points_loader(
    alignment_path: Path,
    colmap_path: Path,
    afdb_dir: Path,
    min_plddt: float,
):
    seq_by_id = read_fasta(alignment_path)
    qcol = read_colmap(colmap_path)
    seq_pos_by_id = {sid: sequence_positions(seq) for sid, seq in seq_by_id.items()}
    cache: dict[str, dict[str, dict[int, np.ndarray]]] = {}

    def load(sid: str) -> dict[str, dict[int, np.ndarray]]:
        if sid in cache:
            return cache[sid]
        empty = {domain: {} for domain in DOMAINS}
        seq = seq_by_id.get(sid)
        path = afdb_path_for_sid(sid, afdb_dir)
        if not seq or path is None:
            cache[sid] = empty
            return empty
        try:
            items = structure_items(path)
            entry_to_structure = map_entry_positions_to_structure(seq, items)
        except Exception:
            cache[sid] = empty
            return empty
        item_by_pos = {int(item["seq_pos"]): item for item in items}
        points = {domain: {} for domain in DOMAINS}
        for domain, qpos_list in DOMAINS.items():
            for qpos in qpos_list:
                aln_col = qcol.get(qpos)
                if aln_col is None or aln_col >= len(seq):
                    continue
                aa = seq[aln_col].upper()
                if aa not in AA_SET:
                    continue
                entry_pos = seq_pos_by_id.get(sid, {}).get(aln_col)
                structure_pos = entry_to_structure.get(entry_pos) if entry_pos is not None else None
                item = item_by_pos.get(structure_pos) if structure_pos is not None else None
                if not item or item["aa"] != aa or float(item["plddt"]) < min_plddt:
                    continue
                points[domain][qpos] = item["xyz"]
        cache[sid] = points
        return points

    return load


def compute_clade_domain_geometry(clades, alignment_path: Path, colmap_path: Path, afdb_dir: Path, min_plddt: float):
    if not alignment_path.exists() or not colmap_path.exists() or not afdb_dir.exists():
        return {}
    load_points = build_domain_points_loader(alignment_path, colmap_path, afdb_dir, min_plddt)
    geometry_by_clade = {}
    for clade in clades:
        rows = {}
        for domain, min_common in GEOMETRY_MIN_COMMON.items():
            tip_points = []
            for sid in clade["tips"]:
                points = load_points(sid).get(domain, {})
                if len(points) >= min_common:
                    tip_points.append((sid, points))
            rmsds = []
            common_counts = []
            for i, (_left_sid, left_points) in enumerate(tip_points):
                for _right_sid, right_points in tip_points[:i]:
                    common = sorted(set(left_points) & set(right_points))
                    if len(common) < min_common:
                        continue
                    left_arr = np.vstack([left_points[qpos] for qpos in common])
                    right_arr = np.vstack([right_points[qpos] for qpos in common])
                    rmsds.append(kabsch_rmsd(left_arr, right_arr))
                    common_counts.append(float(len(common)))
            mean_rmsd, median_rmsd = summarize_float(rmsds)
            mean_common, _median_common = summarize_float(common_counts)
            rows[f"{domain}_geometry_n_structures"] = str(len(tip_points))
            rows[f"{domain}_geometry_n_pairs"] = str(len(rmsds))
            rows[f"{domain}_mean_ca_rmsd_a"] = mean_rmsd
            rows[f"{domain}_median_ca_rmsd_a"] = median_rmsd
            rows[f"{domain}_mean_common_qpos"] = mean_common
        geometry_by_clade[clade["clade_id"]] = rows
    return geometry_by_clade


def terminal_names(tree) -> set[str]:
    return {leaf.name for leaf in tree.get_terminals() if leaf.name}


def root_on_outgroups(tree, outgroup_ids: list[str]) -> list[str]:
    present = sorted(set(outgroup_ids) & terminal_names(tree))
    if not present:
        return []
    tree.root_with_outgroup(*present)
    return present


def remove_outgroups(tree, outgroup_ids: list[str]) -> None:
    for oid in outgroup_ids:
        if oid not in terminal_names(tree):
            continue
        try:
            tree.prune(oid)
        except ValueError:
            pass


def root_and_remove_outgroups(tree, outgroup_ids: list[str]) -> list[str]:
    present = root_on_outgroups(tree, outgroup_ids)
    for oid in present:
        try:
            tree.prune(oid)
        except ValueError:
            pass
    return present


def support_value(node):
    if node.confidence is None:
        return None
    try:
        return float(node.confidence)
    except (TypeError, ValueError):
        return None


def annotate_counts(tree, regime_by_id: dict[str, str]):
    anno = {}

    def walk(node):
        if node.is_terminal():
            sid = node.name or ""
            regime = regime_by_id.get(sid, "")
            counts = Counter({regime: 1}) if regime else Counter()
            unlabelled = 0 if regime else 1
            total = 1
        else:
            counts = Counter()
            unlabelled = 0
            total = 0
            for child in node.clades:
                child_counts, child_unlabelled, child_total = walk(child)
                counts.update(child_counts)
                unlabelled += child_unlabelled
                total += child_total
        anno[id(node)] = {
            "counts": counts,
            "n_labelled": sum(counts.values()),
            "n_unlabelled": unlabelled,
            "n_total": total,
        }
        return counts, unlabelled, total

    walk(tree.root)
    return anno


def node_call(node, info, args):
    if node.is_terminal():
        return None
    support = support_value(node)
    if support is None:
        if not args.allow_missing_support:
            return None
    elif support < args.min_support:
        return None

    n_labelled = info["n_labelled"]
    n_total = info["n_total"]
    if n_labelled < args.min_labelled or n_total < args.min_total_tips:
        return None
    if n_total and info["n_unlabelled"] / n_total > args.max_unlabelled_fraction:
        return None

    counts = info["counts"]
    if not counts:
        return None
    regime, target_n = counts.most_common(1)[0]
    if regime not in args.regimes:
        return None
    target_fraction = target_n / n_labelled
    if target_fraction + 1e-12 < args.min_target_fraction:
        return None
    return regime, target_n, target_fraction


def find_called_clades(tree, anno, args):
    found = []

    def walk(node):
        info = anno[id(node)]
        call = node_call(node, info, args)
        if call:
            regime, target_n, target_fraction = call
            tips = sorted(leaf.name for leaf in node.get_terminals() if leaf.name)
            found.append({
                "node": node,
                "regime": regime,
                "target_n": target_n,
                "target_fraction": target_fraction,
                "tips": tips,
                **info,
            })
            return
        for child in node.clades:
            if not child.is_terminal():
                walk(child)

    walk(tree.root)
    return found


def assign_clade_ids(clades):
    clades.sort(key=lambda c: (
        REGIME_ORDER.get(c["regime"], 99),
        -c["n_labelled"],
        -c["n_total"],
        -(support_value(c["node"]) or -1),
        c["tips"][0] if c["tips"] else "",
    ))
    seen = Counter()
    for clade in clades:
        regime = clade["regime"]
        seen[regime] += 1
        clade["clade_id"] = f"{regime}_{seen[regime]:03d}"


def geometry_fields() -> list[str]:
    fields = []
    for domain in DOMAINS:
        fields.extend([
            f"{domain}_geometry_n_structures",
            f"{domain}_geometry_n_pairs",
            f"{domain}_mean_ca_rmsd_a",
            f"{domain}_median_ca_rmsd_a",
            f"{domain}_mean_common_qpos",
        ])
    return fields


def write_clades(path: Path, clades, regime_by_id, geometry_by_clade=None):
    geometry_by_clade = geometry_by_clade or {}
    fields = [
        "clade_id", "regime", "support", "branch_length",
        "n_total_tips", "n_labelled", "n_unlabelled",
        "n_target", "target_fraction",
        "n_psychro", "n_meso", "n_thermo",
        "labelled_tips", "unlabelled_tips", "all_tips",
    ] + geometry_fields()
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for clade in clades:
            counts = clade["counts"]
            labelled_tips = [sid for sid in clade["tips"] if regime_by_id.get(sid)]
            unlabelled_tips = [sid for sid in clade["tips"] if not regime_by_id.get(sid)]
            out = {
                "clade_id": clade["clade_id"],
                "regime": clade["regime"],
                "support": "" if support_value(clade["node"]) is None else f"{support_value(clade['node']):.6g}",
                "branch_length": "" if clade["node"].branch_length is None else f"{clade['node'].branch_length:.10g}",
                "n_total_tips": clade["n_total"],
                "n_labelled": clade["n_labelled"],
                "n_unlabelled": clade["n_unlabelled"],
                "n_target": clade["target_n"],
                "target_fraction": f"{clade['target_fraction']:.6f}",
                "n_psychro": counts.get("psychro", 0),
                "n_meso": counts.get("meso", 0),
                "n_thermo": counts.get("thermo", 0),
                "labelled_tips": ",".join(labelled_tips),
                "unlabelled_tips": ",".join(unlabelled_tips),
                "all_tips": ",".join(clade["tips"]),
            }
            out.update(geometry_by_clade.get(clade["clade_id"], {}))
            writer.writerow(out)


def write_tip_metadata(path: Path, rows, original_fields, clades, regime_by_id, args):
    context_by_tip = {}
    assigned_by_tip = {}
    clade_by_id = {clade["clade_id"]: clade for clade in clades}
    for clade in clades:
        cid = clade["clade_id"]
        for sid in clade["tips"]:
            context_by_tip[sid] = cid
            if args.assign_unlabelled or regime_by_id.get(sid):
                assigned_by_tip[sid] = cid

    extra_fields = [
        "threshold_regime", "regime_clade_id", "regime_clade",
        "regime_clade_support", "regime_clade_n_labelled",
        "regime_clade_n_total_tips", "regime_clade_target_fraction",
        "regime_clade_context_id", "eligible_for_clade_call",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(original_fields) + extra_fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            sid = row.get("id", "")
            cid = assigned_by_tip.get(sid, "")
            context_id = context_by_tip.get(sid, "")
            clade = clade_by_id.get(cid) if cid else None
            out = {field: row.get(field, "") for field in original_fields}
            out.update({
                "threshold_regime": row.get("_threshold_regime", ""),
                "regime_clade_id": cid,
                "regime_clade": clade["regime"] if clade else "",
                "regime_clade_support": (
                    "" if not clade or support_value(clade["node"]) is None
                    else f"{support_value(clade['node']):.6g}"
                ),
                "regime_clade_n_labelled": clade["n_labelled"] if clade else "",
                "regime_clade_n_total_tips": clade["n_total"] if clade else "",
                "regime_clade_target_fraction": (
                    f"{clade['target_fraction']:.6f}" if clade else ""
                ),
                "regime_clade_context_id": context_id,
                "eligible_for_clade_call": "yes" if regime_by_id.get(sid) else "no",
            })
            writer.writerow(out)


def write_summary(path: Path, args, tree_path: Path, metadata_path: Path, outgroups, clades, rows, regime_by_id):
    regime_totals = Counter(regime_by_id.values())
    clade_totals = Counter(clade["regime"] for clade in clades)
    assigned_labelled = set()
    assigned_total = set()
    for clade in clades:
        assigned_total.update(clade["tips"])
        assigned_labelled.update(sid for sid in clade["tips"] if regime_by_id.get(sid))
    body = [
        "# Regime clade calling summary",
        f"tree: {tree_path}",
        f"metadata: {metadata_path}",
        f"rooted_outgroups_removed: {','.join(outgroups) if outgroups else '(none)'}",
        f"psychro_max: {args.psychro_max}",
        f"thermo_min: {args.thermo_min}",
        f"min_support: {args.min_support}",
        f"allow_missing_support: {args.allow_missing_support}",
        f"min_labelled: {args.min_labelled}",
        f"min_total_tips: {args.min_total_tips}",
        f"min_target_fraction: {args.min_target_fraction}",
        f"max_unlabelled_fraction: {args.max_unlabelled_fraction}",
        f"assign_unlabelled: {args.assign_unlabelled}",
        f"domain_geometry: {'disabled' if args.skip_domain_geometry else 'enabled'}",
        f"geometry_min_plddt: {args.geometry_min_plddt}",
        f"metadata_rows: {len(rows)}",
        f"labelled_total: {sum(regime_totals.values())}",
        f"unlabelled_total: {len(rows) - sum(regime_totals.values())}",
        f"labelled_psychro: {regime_totals.get('psychro', 0)}",
        f"labelled_meso: {regime_totals.get('meso', 0)}",
        f"labelled_thermo: {regime_totals.get('thermo', 0)}",
        f"clades_total: {len(clades)}",
        f"clades_psychro: {clade_totals.get('psychro', 0)}",
        f"clades_meso: {clade_totals.get('meso', 0)}",
        f"clades_thermo: {clade_totals.get('thermo', 0)}",
        f"labelled_tips_in_called_clades: {len(assigned_labelled)}",
        f"all_tips_in_called_clade_contexts: {len(assigned_total)}",
    ]
    path.write_text("\n".join(body) + "\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", default=str(TREE), help="Newick tree; default is final IQ-TREE .treefile")
    parser.add_argument("--metadata", default=str(METADATA), help="Final repset metadata TSV")
    parser.add_argument("--outgroup-ids", default=str(OUTGROUP_IDS), help="One outgroup tip id per line")
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--prefix", default="regime_clades")
    parser.add_argument("--alignment", default=str(ALN_FA), help="Refined full-domain MSA for AFDB domain geometry")
    parser.add_argument("--column-map", default=str(COL_MAP), help="Refined MSA ANKros column map for AFDB domain geometry")
    parser.add_argument("--afdb-dir", default=str(AFDB_DIR), help="Local AFDB PDB directory for domain geometry summaries")
    parser.add_argument("--psychro-max", type=float, default=cfg.OGT_PSYCHRO_MAX)
    parser.add_argument("--thermo-min", type=float, default=cfg.OGT_THERMO_MIN)
    parser.add_argument("--min-support", type=float, default=90.0)
    parser.add_argument("--allow-missing-support", dest="allow_missing_support", action="store_true",
                        help="Keep clade candidates when branch support is absent; this is the default.")
    parser.add_argument("--require-support", dest="allow_missing_support", action="store_false",
                        help="Require branch support and apply --min-support.")
    parser.set_defaults(allow_missing_support=True)
    parser.add_argument("--min-labelled", type=int, default=3)
    parser.add_argument("--min-total-tips", type=int, default=3)
    parser.add_argument("--min-target-fraction", type=float, default=1.0,
                        help="1.0 requires no labelled regime mixing; lower values enable majority calls")
    parser.add_argument("--max-unlabelled-fraction", type=float, default=1.0,
                        help="Maximum unlabelled fraction inside a called clade; default disables this filter")
    parser.add_argument("--regimes", default="psychro,meso,thermo",
                        help="Comma-separated regimes eligible for clade calls")
    parser.add_argument("--assign-unlabelled", action="store_true",
                        help="Assign unlabelled tips inside called clades in the tip metadata")
    parser.add_argument("--skip-domain-geometry", action="store_true",
                        help="Do not append within-clade AFDB domain RMSD summaries to the clade table")
    parser.add_argument("--geometry-min-plddt", type=float, default=70.0,
                        help="Minimum AFDB CA pLDDT used for within-clade domain RMSD summaries")
    args = parser.parse_args()
    args.regimes = {r.strip() for r in args.regimes.split(",") if r.strip()}
    unknown = args.regimes - set(REGIME_ORDER)
    if unknown:
        raise SystemExit(f"ERROR: unknown regime(s): {', '.join(sorted(unknown))}")
    if args.psychro_max >= args.thermo_min:
        raise SystemExit("ERROR: --psychro-max must be smaller than --thermo-min")
    if not 0 < args.min_target_fraction <= 1:
        raise SystemExit("ERROR: --min-target-fraction must be in (0, 1]")
    if not 0 <= args.max_unlabelled_fraction <= 1:
        raise SystemExit("ERROR: --max-unlabelled-fraction must be in [0, 1]")
    if args.geometry_min_plddt < 0:
        raise SystemExit("ERROR: --geometry-min-plddt must be >= 0")
    return args


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    metadata_path = Path(args.metadata)
    outgroup_path = Path(args.outgroup_ids)
    out_dir = Path(args.out_dir)

    for path, label in [(tree_path, "tree"), (metadata_path, "metadata")]:
        if not path.exists():
            raise SystemExit(f"ERROR: {label} not found: {path}")

    rows, meta_by_id, original_fields = load_metadata(metadata_path, args.psychro_max, args.thermo_min)
    regime_by_id = {
        sid: row["_threshold_regime"]
        for sid, row in meta_by_id.items()
        if row.get("_threshold_regime")
    }

    tree = Phylo.read(str(tree_path), "newick")
    outgroup_ids = load_outgroup_ids(outgroup_path)
    removed_outgroups = root_on_outgroups(tree, outgroup_ids)
    rooted_with_outgroups = copy.deepcopy(tree)
    remove_outgroups(tree, removed_outgroups)

    tree_tips = terminal_names(tree)
    metadata_tips = set(meta_by_id)
    missing_from_tree = sorted(metadata_tips - tree_tips)
    if missing_from_tree:
        print(f"[warn] {len(missing_from_tree)} metadata ids are absent from the tree", file=sys.stderr)

    anno = annotate_counts(tree, regime_by_id)
    clades = find_called_clades(tree, anno, args)
    assign_clade_ids(clades)
    geometry_by_clade = {}
    if not args.skip_domain_geometry:
        geometry_by_clade = compute_clade_domain_geometry(
            clades,
            Path(args.alignment),
            Path(args.column_map),
            Path(args.afdb_dir),
            args.geometry_min_plddt,
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    clades_tsv = out_dir / f"{args.prefix}.tsv"
    tip_tsv = out_dir / f"{args.prefix}_tip_metadata.tsv"
    summary_txt = out_dir / f"{args.prefix}_summary.txt"
    rooted_with_outgroups_tree = out_dir / f"{args.prefix}_rooted_with_outgroups.nwk"
    rooted_tree = out_dir / f"{args.prefix}_rooted_ingroup.nwk"

    Phylo.write(rooted_with_outgroups, str(rooted_with_outgroups_tree), "newick")
    write_clades(clades_tsv, clades, regime_by_id, geometry_by_clade)
    write_tip_metadata(tip_tsv, rows, original_fields, clades, regime_by_id, args)
    write_summary(summary_txt, args, tree_path, metadata_path, removed_outgroups, clades, rows, regime_by_id)
    Phylo.write(tree, str(rooted_tree), "newick")

    clade_totals = Counter(clade["regime"] for clade in clades)
    print(f"Called {len(clades)} clades "
          f"(psychro={clade_totals.get('psychro', 0)}, "
          f"meso={clade_totals.get('meso', 0)}, "
          f"thermo={clade_totals.get('thermo', 0)})")
    print(f"Saved: {clades_tsv}")
    print(f"Saved: {tip_tsv}")
    print(f"Saved: {summary_txt}")
    print(f"Saved: {rooted_with_outgroups_tree}")
    print(f"Saved: {rooted_tree}")


if __name__ == "__main__":
    main()
