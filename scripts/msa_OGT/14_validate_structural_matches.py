#!/usr/bin/env python3
"""Step 14: validate AFDB structures against the ANKros class-I fold."""
from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_QUERY = ROOT / "results" / "structural" / "docked_holo" / "ankros_fad_fmn_donor_holo.pdb"
DEFAULT_AFDB_DIR = ROOT / "structures" / "afdb"
DEFAULT_OUT_DIR = ROOT / "results" / "msa_OGT" / "structure_validation"
DEFAULT_RAW = DEFAULT_OUT_DIR / "ankros_vs_afdb_foldseek.tsv"
DEFAULT_SUMMARY = DEFAULT_OUT_DIR / "ankros_vs_afdb_summary.tsv"
DEFAULT_REVIEW = DEFAULT_OUT_DIR / "quick_classI_mismatch_review.tsv"
DEFAULT_MANIFEST = DEFAULT_OUT_DIR / "run_manifest.txt"
FORMAT_FIELDS = [
    "query", "target", "qlen", "tlen", "alnlen", "qcov", "tcov", "pident",
    "fident", "evalue", "bits", "prob", "lddt", "qtmscore", "ttmscore",
    "alntmscore", "rmsd",
]
FLOAT_FIELDS = {
    "qcov", "tcov", "pident", "fident", "evalue", "bits", "prob", "lddt",
    "qtmscore", "ttmscore", "alntmscore", "rmsd",
}
INT_FIELDS = {"qlen", "tlen", "alnlen"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", type=Path, default=DEFAULT_QUERY)
    parser.add_argument("--afdb-dir", type=Path, default=DEFAULT_AFDB_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--foldseek", type=Path, default=ROOT / "envs" / "ankros" / "bin" / "foldseek")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Re-run Foldseek even if the raw TSV exists.")
    parser.add_argument("--hard-min-qcov", type=float, default=0.90)
    parser.add_argument("--hard-min-tcov", type=float, default=0.80)
    parser.add_argument("--hard-min-qtm", type=float, default=0.80)
    parser.add_argument("--hard-min-alntm", type=float, default=0.80)
    parser.add_argument("--hard-min-prob", type=float, default=0.95)
    parser.add_argument("--review-min-qcov", type=float, default=0.95)
    parser.add_argument("--review-min-tcov", type=float, default=0.90)
    parser.add_argument("--review-min-qtm", type=float, default=0.88)
    parser.add_argument("--review-min-lddt", type=float, default=0.72)
    parser.add_argument("--review-min-pident", type=float, default=20.0)
    return parser.parse_args()


def accession_from_target(target: str) -> str:
    if target.startswith("AF-") and "-F1-model" in target:
        return target.removeprefix("AF-").split("-F1-model", 1)[0]
    return target


def foldseek_version(foldseek: Path) -> str:
    try:
        return subprocess.check_output([str(foldseek), "version"], text=True).strip()
    except Exception:
        return "unknown"


def run_foldseek(args: argparse.Namespace, raw_path: Path) -> None:
    if not args.query.exists():
        raise SystemExit(f"ERROR: query structure not found: {args.query}")
    if not args.afdb_dir.exists():
        raise SystemExit(f"ERROR: AFDB directory not found: {args.afdb_dir}")
    if not args.foldseek.exists():
        raise SystemExit(f"ERROR: Foldseek binary not found: {args.foldseek}")
    if raw_path.exists() and not args.force:
        print(f"Using existing Foldseek TSV: {raw_path}; pass --force to rebuild")
        return

    tmp_dir = args.out_dir / "tmp"
    if tmp_dir.exists():
        import shutil
        shutil.rmtree(tmp_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(args.foldseek), "easy-search",
        str(args.query),
        str(args.afdb_dir),
        str(raw_path),
        str(tmp_dir),
        "--alignment-type", "1",
        "--exhaustive-search", "1",
        "--format-mode", "4",
        "--format-output", ",".join(FORMAT_FIELDS),
        "--threads", str(args.threads),
    ]
    print("Running: " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(ROOT), check=True)


def read_foldseek(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in FORMAT_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"ERROR: missing Foldseek columns in {path}: {missing}")
        for row in reader:
            parsed: dict[str, object] = dict(row)
            for field in INT_FIELDS:
                parsed[field] = int(float(str(parsed[field])))
            for field in FLOAT_FIELDS:
                parsed[field] = float(str(parsed[field]))
            parsed["accession"] = accession_from_target(str(parsed["target"]))
            rows.append(parsed)
    return rows


def quantile(vals: list[float], pct: float) -> float:
    vals = sorted(vals)
    return vals[round((len(vals) - 1) * pct / 100.0)]


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    metrics = ["qcov", "tcov", "qtmscore", "ttmscore", "alntmscore", "lddt", "rmsd", "pident", "prob"]
    fields = ["metric", "n", "min", "p01", "p05", "median", "p95", "max"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for metric in metrics:
            vals = [float(row[metric]) for row in rows]
            writer.writerow({
                "metric": metric,
                "n": len(vals),
                "min": f"{min(vals):.6g}",
                "p01": f"{quantile(vals, 1):.6g}",
                "p05": f"{quantile(vals, 5):.6g}",
                "median": f"{statistics.median(vals):.6g}",
                "p95": f"{quantile(vals, 95):.6g}",
                "max": f"{max(vals):.6g}",
            })


def classify(row: dict[str, object], args: argparse.Namespace) -> tuple[str, list[str]]:
    hard = []
    if float(row["qcov"]) < args.hard_min_qcov:
        hard.append(f"qcov_lt_{args.hard_min_qcov:g}")
    if float(row["tcov"]) < args.hard_min_tcov:
        hard.append(f"tcov_lt_{args.hard_min_tcov:g}")
    if float(row["qtmscore"]) < args.hard_min_qtm:
        hard.append(f"qtmscore_lt_{args.hard_min_qtm:g}")
    if float(row["alntmscore"]) < args.hard_min_alntm:
        hard.append(f"alntmscore_lt_{args.hard_min_alntm:g}")
    if float(row["prob"]) < args.hard_min_prob:
        hard.append(f"prob_lt_{args.hard_min_prob:g}")
    if hard:
        return "fail_possible_mismatch", hard

    review = []
    if float(row["qcov"]) < args.review_min_qcov:
        review.append(f"qcov_lt_{args.review_min_qcov:g}")
    if float(row["tcov"]) < args.review_min_tcov:
        review.append(f"tcov_lt_{args.review_min_tcov:g}")
    if float(row["qtmscore"]) < args.review_min_qtm:
        review.append(f"qtmscore_lt_{args.review_min_qtm:g}")
    if float(row["lddt"]) < args.review_min_lddt:
        review.append(f"lddt_lt_{args.review_min_lddt:g}")
    if float(row["pident"]) < args.review_min_pident:
        review.append(f"pident_lt_{args.review_min_pident:g}")
    if review:
        return "review_borderline_not_mismatch", review
    return "pass_quick_classI_structural_match", []


def write_review(path: Path, rows: list[dict[str, object]], args: argparse.Namespace) -> Counter:
    fields = [
        "status", "reasons", "target", "accession", "qlen", "tlen", "alnlen",
        "qcov", "tcov", "qtmscore", "ttmscore", "alntmscore", "lddt", "rmsd",
        "pident", "prob",
    ]
    out_rows = []
    for row in rows:
        status, reasons = classify(row, args)
        out_rows.append({
            "status": status,
            "reasons": ";".join(reasons),
            **{field: row[field] for field in fields if field not in {"status", "reasons"}},
        })
    out_rows.sort(key=lambda row: (
        row["status"] != "fail_possible_mismatch",
        row["status"] != "review_borderline_not_mismatch",
        float(row["qtmscore"]),
        float(row["qcov"]),
        float(row["tcov"]),
    ))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in out_rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return Counter(str(row["status"]) for row in out_rows)


def write_manifest(path: Path, args: argparse.Namespace, counts: Counter, target_count: int) -> None:
    body = [
        "Step 14 structural validation: ANKros-vs-AFDB quick class-I structural mismatch scan",
        f"foldseek={foldseek_version(args.foldseek)}",
        f"query={args.query}",
        f"afdb_dir={args.afdb_dir}",
        f"target_pdbs={target_count}",
        "alignment_type=1 (TM-align mode)",
        "exhaustive_search=1",
        f"threads={args.threads}",
        (
            "hard_fail_thresholds="
            f"qcov<{args.hard_min_qcov};tcov<{args.hard_min_tcov};"
            f"qtmscore<{args.hard_min_qtm};alntmscore<{args.hard_min_alntm};prob<{args.hard_min_prob}"
        ),
        (
            "soft_review_thresholds="
            f"qcov<{args.review_min_qcov};tcov<{args.review_min_tcov};"
            f"qtmscore<{args.review_min_qtm};lddt<{args.review_min_lddt};pident<{args.review_min_pident}"
        ),
        f"hard_fail_possible_mismatch={counts.get('fail_possible_mismatch', 0)}",
        f"review_borderline_not_mismatch={counts.get('review_borderline_not_mismatch', 0)}",
        f"pass_quick_classI_structural_match={counts.get('pass_quick_classI_structural_match', 0)}",
        f"raw_output={DEFAULT_RAW.name}",
        f"summary={DEFAULT_SUMMARY.name}",
        f"review_table={DEFAULT_REVIEW.name}",
    ]
    path.write_text("\n".join(body) + "\n")


def main() -> None:
    args = parse_args()
    raw_path = args.out_dir / DEFAULT_RAW.name
    summary_path = args.out_dir / DEFAULT_SUMMARY.name
    review_path = args.out_dir / DEFAULT_REVIEW.name
    manifest_path = args.out_dir / DEFAULT_MANIFEST.name
    run_foldseek(args, raw_path)
    rows = read_foldseek(raw_path)
    target_count = len(list(args.afdb_dir.glob("*.pdb")))
    if len(rows) != target_count:
        print(f"WARNING: Foldseek returned {len(rows)} rows for {target_count} AFDB PDBs", file=sys.stderr)
    write_summary(summary_path, rows)
    counts = write_review(review_path, rows, args)
    write_manifest(manifest_path, args, counts, target_count)
    print(f"Validated {len(rows)} Foldseek hits against {target_count} AFDB PDBs")
    print(f"Status counts: {dict(counts)}")
    print(f"Saved: {raw_path}")
    print(f"Saved: {summary_path}")
    print(f"Saved: {review_path}")
    print(f"Saved: {manifest_path}")


if __name__ == "__main__":
    main()
