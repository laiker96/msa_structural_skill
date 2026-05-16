#!/usr/bin/env python3
"""
Step 06: iterative alignment QC and domain-level metrics.

This step does not realign. It scores the HMM-aligned representative set in the
ANKros coordinate frame, removes clear sequence outliers iteratively, and writes
global plus per-domain coverage/identity metrics.
"""
from __future__ import annotations

import csv
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


QUERY_ID = "photoHymenobact"
DOMAINS = {
    "antenna": (1, 130),
    "linker": (131, 205),
    "catalytic": (206, 437),
}

MATCH_ALN = cfg.INTER_DIR / "repset_hmmalign_matchcols.fa"
FULL_ALN = cfg.INTER_DIR / "repset_hmmalign.fa"
RAW_HMM_ALN = cfg.INTER_DIR / "hmm" / "repset_hmmalign.afa"
REPSET_TSV = cfg.INTER_DIR / "repset_metadata.tsv"

QC_METRICS = cfg.INTER_DIR / "alignment_qc_metrics.tsv"
QC_DOMAIN_METRICS = cfg.INTER_DIR / "alignment_qc_domain_metrics.tsv"
QC_ITERATIONS = cfg.INTER_DIR / "alignment_qc_iterations.tsv"
QC_REJECTED = cfg.INTER_DIR / "alignment_qc_rejected.tsv"
QC_CLEAN_MATCH = cfg.INTER_DIR / "repset_hmmalign_matchcols_qc.fa"
QC_CLEAN_FULL = cfg.INTER_DIR / "repset_hmmalign_qc.fa"
QC_CLEAN_META = cfg.INTER_DIR / "repset_metadata_qc.tsv"

MIN_GLOBAL_COV = float(os.environ.get("MSA_OGT_QC_MIN_GLOBAL_COV", "0.75"))
MIN_ANTENNA_COV = float(os.environ.get("MSA_OGT_QC_MIN_ANTENNA_COV", "0.70"))
MIN_LINKER_COV = float(os.environ.get("MSA_OGT_QC_MIN_LINKER_COV", "0.35"))
MIN_CATALYTIC_COV = float(os.environ.get("MSA_OGT_QC_MIN_CATALYTIC_COV", "0.80"))
MAX_LOW_OCC_RES_FRAC = float(os.environ.get("MSA_OGT_QC_MAX_LOW_OCC_RES_FRAC", "0.35"))
LOW_OCC_COL_THRESHOLD = float(os.environ.get("MSA_OGT_QC_LOW_OCC_COL_THRESHOLD", "0.10"))
MAX_ITER = int(os.environ.get("MSA_OGT_QC_MAX_ITER", "5"))


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
        for sid, hdr, seq in entries:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def find_query(entries):
    for sid, _hdr, seq in entries:
        if sid == QUERY_ID or QUERY_ID in sid:
            return sid, seq
    raise SystemExit("ERROR: query sequence not found in alignment")


def query_position_map(query_aln):
    pos_by_col = {}
    pos = 0
    for col, aa in enumerate(query_aln):
        if aa != "-":
            pos += 1
            pos_by_col[col] = pos
    return pos_by_col


def domain_columns(pos_by_col):
    out = {}
    for name, (start, end) in DOMAINS.items():
        out[name] = [col for col, pos in pos_by_col.items() if start <= pos <= end]
    out["global"] = sorted(pos_by_col)
    return out


def safe_div(num, den):
    return num / den if den else 0.0


def load_metadata():
    if not REPSET_TSV.exists():
        return {}, []
    with open(REPSET_TSV, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        return {row["id"]: row for row in rows}, reader.fieldnames or []


def load_insert_counts():
    if not RAW_HMM_ALN.exists():
        return defaultdict(int), defaultdict(lambda: defaultdict(int))
    raw = read_fasta(RAW_HMM_ALN)
    query_sid, query = find_query([(sid, hdr, seq.upper().replace(".", "-")) for sid, hdr, seq in raw])
    # Insert columns in AFA output are lowercase in at least one sequence.
    insert_cols = {
        i
        for i, col in enumerate(zip(*(seq for _sid, seq in [(sid, seq) for sid, _hdr, seq in raw])))
        if any(char.islower() for char in col)
    }
    insert_total = defaultdict(int)
    insert_domain = defaultdict(lambda: defaultdict(int))
    prev_query_pos = None
    query_pos_by_raw_col = {}
    pos = 0
    for col, aa in enumerate(query):
        if aa != "-":
            pos += 1
            prev_query_pos = pos
        query_pos_by_raw_col[col] = prev_query_pos

    for sid, _hdr, seq in raw:
        for col in insert_cols:
            aa = seq[col]
            if not aa.isalpha():
                continue
            insert_total[sid] += 1
            qpos = query_pos_by_raw_col.get(col)
            if qpos is None:
                domain = "n_terminal_insert"
            else:
                domain = next(
                    (name for name, (start, end) in DOMAINS.items() if start <= qpos <= end),
                    "c_terminal_insert",
                )
            insert_domain[sid][domain] += 1
    return insert_total, insert_domain


def compute_metrics(entries, active_ids, meta_by_id, insert_total, insert_domain):
    query_sid, query = find_query(entries)
    pos_by_col = query_position_map(query)
    cols_by_domain = domain_columns(pos_by_col)
    active_entries = [(sid, hdr, seq) for sid, hdr, seq in entries if sid in active_ids]
    active_cols = list(zip(*(seq for _sid, _hdr, seq in active_entries)))
    occupancy = [
        safe_div(sum(aa != "-" for aa in col), len(active_entries))
        for col in active_cols
    ]
    low_occ_cols = {i for i, occ in enumerate(occupancy) if occ < LOW_OCC_COL_THRESHOLD}

    metrics = {}
    domain_rows = {}
    for sid, hdr, seq in entries:
        row = meta_by_id.get(sid, {})
        clean_len = sum(aa != "-" for aa in seq)
        low_occ_res = sum(1 for i, aa in enumerate(seq) if aa != "-" and i in low_occ_cols)
        record = {
            "id": sid,
            "header": hdr,
            "active": "yes" if sid in active_ids else "no",
            "source": row.get("source", ""),
            "regime": row.get("regime", ""),
            "ogt": row.get("ogt", ""),
            "pident": row.get("pident", ""),
            "class_pident": row.get("class_pident", ""),
            "class_bits": row.get("class_bits", ""),
            "aligned_residues": str(clean_len),
            "insert_residues": str(insert_total.get(sid, 0)),
            "low_occupancy_residues": str(low_occ_res),
            "low_occupancy_residue_fraction": f"{safe_div(low_occ_res, clean_len):.4f}",
        }
        per_domain = {}
        for domain, cols in cols_by_domain.items():
            q_res_cols = [col for col in cols if query[col] != "-"]
            present = [col for col in q_res_cols if seq[col] != "-"]
            identical = [
                col for col in present
                if seq[col].upper() == query[col].upper()
            ]
            coverage = safe_div(len(present), len(q_res_cols))
            identity = safe_div(len(identical), len(present))
            gap_fraction = 1.0 - coverage
            insert_res = insert_total.get(sid, 0) if domain == "global" else insert_domain[sid].get(domain, 0)
            per_domain[domain] = {
                "coverage": coverage,
                "identity_to_query": identity,
                "gap_fraction": gap_fraction,
                "non_gap_query_columns": len(present),
                "query_columns": len(q_res_cols),
                "insert_residues": insert_res,
            }
            prefix = f"{domain}_"
            record[prefix + "coverage"] = f"{coverage:.4f}"
            record[prefix + "identity_to_query"] = f"{identity:.4f}"
            record[prefix + "gap_fraction"] = f"{gap_fraction:.4f}"
            record[prefix + "insert_residues"] = str(insert_res)
        metrics[sid] = record
        domain_rows[sid] = per_domain
    return metrics, domain_rows, occupancy


def rejection_reasons(record):
    if record["id"] == QUERY_ID:
        return []
    checks = [
        ("global_coverage", float(record["global_coverage"]), MIN_GLOBAL_COV),
        ("antenna_coverage", float(record["antenna_coverage"]), MIN_ANTENNA_COV),
        ("linker_coverage", float(record["linker_coverage"]), MIN_LINKER_COV),
        ("catalytic_coverage", float(record["catalytic_coverage"]), MIN_CATALYTIC_COV),
    ]
    reasons = [
        f"{name}={value:.3f}<{threshold:.3f}"
        for name, value, threshold in checks
        if value < threshold
    ]
    low_occ_frac = float(record["low_occupancy_residue_fraction"])
    if low_occ_frac > MAX_LOW_OCC_RES_FRAC:
        reasons.append(
            f"low_occupancy_residue_fraction={low_occ_frac:.3f}>{MAX_LOW_OCC_RES_FRAC:.3f}"
        )
    return reasons


def main():
    if not MATCH_ALN.exists() or not FULL_ALN.exists():
        sys.exit("ERROR: run 05_align_hmm.py first")

    match_entries = read_fasta(MATCH_ALN)
    full_entries = read_fasta(FULL_ALN)
    entries = full_entries
    meta_by_id, meta_fields = load_metadata()
    insert_total, insert_domain = load_insert_counts()

    active_ids = {sid for sid, _hdr, _seq in entries}
    query_sid, _query = find_query(entries)
    active_ids.add(query_sid)
    rejected = []
    iteration_rows = []
    final_metrics = {}
    final_domain_rows = {}

    for iteration in range(1, MAX_ITER + 1):
        metrics, domain_rows, occupancy = compute_metrics(
            entries, active_ids, meta_by_id, insert_total, insert_domain
        )
        to_remove = []
        for sid in sorted(active_ids):
            reasons = rejection_reasons(metrics[sid])
            if reasons:
                to_remove.append((sid, reasons))
        query_cols = sorted(query_position_map(find_query(entries)[1]))
        query_occupancy = [occupancy[col] for col in query_cols]
        iteration_rows.append({
            "iteration": str(iteration),
            "active_before": str(len(active_ids)),
            "removed": str(len(to_remove)),
            "mean_full_column_occupancy": f"{statistics.mean(occupancy):.4f}",
            "median_full_column_occupancy": f"{statistics.median(occupancy):.4f}",
            "mean_query_column_occupancy": f"{statistics.mean(query_occupancy):.4f}",
            "median_query_column_occupancy": f"{statistics.median(query_occupancy):.4f}",
        })
        if not to_remove:
            final_metrics, final_domain_rows = metrics, domain_rows
            break
        for sid, reasons in to_remove:
            active_ids.remove(sid)
            rejected.append({
                "id": sid,
                "iteration": str(iteration),
                "reasons": ";".join(reasons),
            })
    else:
        final_metrics, final_domain_rows, _occupancy = compute_metrics(
            entries, active_ids, meta_by_id, insert_total, insert_domain
        )

    metric_fields = [
        "id", "header", "active", "source", "regime", "ogt",
        "pident", "class_pident", "class_bits",
        "aligned_residues", "insert_residues",
        "low_occupancy_residues", "low_occupancy_residue_fraction",
    ]
    for domain in ["global", "antenna", "linker", "catalytic"]:
        metric_fields.extend([
            f"{domain}_coverage",
            f"{domain}_identity_to_query",
            f"{domain}_gap_fraction",
            f"{domain}_insert_residues",
        ])
    with open(QC_METRICS, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=metric_fields, delimiter="\t")
        writer.writeheader()
        for sid, _hdr, _seq in entries:
            writer.writerow({field: final_metrics[sid].get(field, "") for field in metric_fields})

    domain_fields = [
        "id", "active", "domain", "query_start", "query_end",
        "coverage", "identity_to_query", "gap_fraction",
        "non_gap_query_columns", "query_columns", "insert_residues",
    ]
    with open(QC_DOMAIN_METRICS, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=domain_fields, delimiter="\t")
        writer.writeheader()
        for sid, _hdr, _seq in entries:
            for domain in ["antenna", "linker", "catalytic"]:
                start, end = DOMAINS[domain]
                vals = final_domain_rows[sid][domain]
                writer.writerow({
                    "id": sid,
                    "active": "yes" if sid in active_ids else "no",
                    "domain": domain,
                    "query_start": str(start),
                    "query_end": str(end),
                    "coverage": f"{vals['coverage']:.4f}",
                    "identity_to_query": f"{vals['identity_to_query']:.4f}",
                    "gap_fraction": f"{vals['gap_fraction']:.4f}",
                    "non_gap_query_columns": str(vals["non_gap_query_columns"]),
                    "query_columns": str(vals["query_columns"]),
                    "insert_residues": str(vals["insert_residues"]),
                })

    with open(QC_ITERATIONS, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "iteration", "active_before", "removed",
                "mean_full_column_occupancy", "median_full_column_occupancy",
                "mean_query_column_occupancy", "median_query_column_occupancy",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(iteration_rows)

    with open(QC_REJECTED, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "iteration", "reasons"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rejected)

    active_match = [(sid, hdr, seq) for sid, hdr, seq in match_entries if sid in active_ids]
    active_full = [(sid, hdr, seq) for sid, hdr, seq in full_entries if sid in active_ids]
    write_fasta(active_match, QC_CLEAN_MATCH)
    write_fasta(active_full, QC_CLEAN_FULL)

    if meta_fields:
        with open(QC_CLEAN_META, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=meta_fields, delimiter="\t")
            writer.writeheader()
            for sid, _hdr, _seq in active_match:
                if sid in meta_by_id:
                    writer.writerow(meta_by_id[sid])

    print(f"Input alignment: {len(entries)} sequences")
    print(f"QC kept: {len(active_ids)}; rejected: {len(rejected)}")
    print(
        "Thresholds: "
        f"global>={MIN_GLOBAL_COV}, antenna>={MIN_ANTENNA_COV}, "
        f"linker>={MIN_LINKER_COV}, catalytic>={MIN_CATALYTIC_COV}, "
        f"low-occ-res-frac<={MAX_LOW_OCC_RES_FRAC}"
    )
    print(f"Saved: {QC_METRICS}")
    print(f"Saved: {QC_DOMAIN_METRICS}")
    print(f"Saved: {QC_ITERATIONS}")
    print(f"Saved: {QC_REJECTED}")
    print(f"Saved: {QC_CLEAN_MATCH}")
    print(f"Saved: {QC_CLEAN_FULL}")
    if meta_fields:
        print(f"Saved: {QC_CLEAN_META}")


if __name__ == "__main__":
    main()
