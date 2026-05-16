#!/usr/bin/env python3
"""
Step 01: broad MMseqs2 search against UniRef.

The search is intentionally permissive. All returned hits keep pident, qcov,
tcov, evalue, bits, and the full UniRef header so later steps can audit the
effect of stricter thresholds after taxid/OGTFinder annotation.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")


MMSEQS = cfg.resolve_bin("mmseqs")


def create_mmseqs_db():
    if Path(f"{cfg.ACTIVE_DB_MMSEQS}.dbtype").exists():
        print(f"MMseqs2 DB exists: {cfg.ACTIVE_DB_MMSEQS}")
        return
    if not cfg.ACTIVE_DB_FASTA.exists():
        sys.exit(f"ERROR: UniRef FASTA not found: {cfg.ACTIVE_DB_FASTA}")
    print(f"Creating MMseqs2 DB from {cfg.ACTIVE_DB_FASTA}")
    subprocess.run([MMSEQS, "createdb", str(cfg.ACTIVE_DB_FASTA), str(cfg.ACTIVE_DB_MMSEQS)],
                   check=True)


def main():
    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = cfg.INTER_DIR / "mmseqs_search_results.tsv"
    if out_tsv.exists():
        n = sum(1 for _ in open(out_tsv))
        print(f"Search results already exist: {out_tsv} ({n} rows)")
        print("Delete that file to re-run the search.")
        return

    if shutil.which(MMSEQS) is None and not Path(MMSEQS).exists():
        sys.exit("ERROR: mmseqs not found on PATH")

    create_mmseqs_db()

    work = cfg.INTER_DIR / "mmseqs_work"
    tmp = work / "tmp"
    query_db = work / "query_db"
    result_db = work / "search_result"
    work.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)

    print("Creating query DB")
    subprocess.run([MMSEQS, "createdb", str(cfg.QUERY_FASTA), str(query_db)], check=True)

    print(
        "Searching "
        f"{cfg.ACTIVE_DB_NAME}: -s {cfg.MMSEQS_SENSITIVITY}, "
        f"-e {cfg.MMSEQS_EVALUE}, --min-seq-id {cfg.MMSEQS_MIN_SEQ_ID}, "
        f"-c {cfg.MMSEQS_COVERAGE}, --max-seqs {cfg.MMSEQS_MAX_SEQS}"
    )
    subprocess.run([
        MMSEQS, "search",
        str(query_db), str(cfg.ACTIVE_DB_MMSEQS), str(result_db), str(tmp),
        "-s", str(cfg.MMSEQS_SENSITIVITY),
        "-e", str(cfg.MMSEQS_EVALUE),
        "--min-seq-id", str(cfg.MMSEQS_MIN_SEQ_ID),
        "-c", str(cfg.MMSEQS_COVERAGE),
        "--cov-mode", "1",
        "--max-seqs", str(cfg.MMSEQS_MAX_SEQS),
        "--threads", str(cfg.N_THREADS),
    ], check=True)

    print("Converting results to TSV")
    subprocess.run([
        MMSEQS, "convertalis",
        str(query_db), str(cfg.ACTIVE_DB_MMSEQS), str(result_db), str(out_tsv),
        "--format-output",
        "target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,"
        "evalue,bits,qcov,tcov,theader",
    ], check=True)

    n = sum(1 for _ in open(out_tsv))
    print(f"Saved: {out_tsv} ({n} rows)")
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
