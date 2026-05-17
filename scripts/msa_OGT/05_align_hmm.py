#!/usr/bin/env python3
"""
Step 05: profile-HMM alignment of the OGT-annotated master set.

The sequence universe comes from step 04. Sampling uses the OGT regime labels
and preserves the selected sequence metadata, including OGTFinder range flags.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
cfg = import_module("00_config")

MMSEQS = cfg.resolve_bin("mmseqs")
MAFFT = cfg.resolve_bin("mafft")
HMMBUILD = cfg.resolve_bin("hmmbuild")
HMMALIGN = cfg.resolve_bin("hmmalign")

MASTER_FA = cfg.INTER_DIR / "master_homologs.fa"
MASTER_TSV = cfg.INTER_DIR / "master_homologs.tsv"
PRECLUSTER_FA = cfg.INTER_DIR / "master_precluster.fa"
PRECLUSTER_MEMBERS = cfg.INTER_DIR / "master_precluster_members.tsv"
REPSET_FA = cfg.INTER_DIR / "repset.fa"
REPSET_TSV = cfg.INTER_DIR / "repset_metadata.tsv"
ALIGNED_FA = cfg.INTER_DIR / "repset_hmmalign.fa"
ALIGNED_MATCH_FA = cfg.INTER_DIR / "repset_hmmalign_matchcols.fa"

HMM_DIR = cfg.INTER_DIR / "hmm"
SEED_FA = cfg.PROJECT_ROOT / "test" / "uniprot_classI_cpd.fa"
SEED_ALN = HMM_DIR / "classI_seed_mafft.fa"
PROFILE_HMM = HMM_DIR / "classI_profile.hmm"
HMM_ALIGN_RAW = HMM_DIR / "repset_hmmalign.afa"


def write_aligned_fasta(entries, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for _sid, hdr, seq in entries:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def require_bins(names):
    missing = []
    for name in names:
        resolved = cfg.resolve_bin(name)
        if shutil.which(resolved) is None and not Path(resolved).exists():
            missing.append(name)
    if missing:
        sys.exit(f"ERROR: missing required alignment binaries: {', '.join(missing)}")


def load_metadata():
    with open(MASTER_TSV, newline="") as f:
        return {row["id"]: row for row in csv.DictReader(f, delimiter="\t")}


def as_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def representative_score(entry, metadata):
    sid = entry[0]
    row = metadata.get(sid, {})
    match_rank = {
        "exact_taxid": 5,
        "species_taxid": 4,
        "name_exact": 3,
        "name_binomial": 2,
    }.get(row.get("ogt_match_type", ""), 0)
    has_ogt = 1 if row.get("ogt") else 0
    point_ogt = 1 if row.get("ogt") and row.get("ogt_has_range") == "no" else 0
    return (
        has_ogt,
        match_rank,
        point_ogt,
        as_float(row.get("class_bits")),
        as_float(row.get("qcov")),
        as_float(row.get("bits")),
        -abs(as_float(row.get("length"), len(entry[2])) - len(entry[2])),
        sid,
    )


def parse_cluster_members(cluster_tsv, entries):
    entry_by_id = {sid: (sid, hdr, seq) for sid, hdr, seq in entries}
    clusters = {}
    assigned = set()
    if cluster_tsv.exists():
        with open(cluster_tsv) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                rep, member = parts[0], parts[1]
                if member in entry_by_id:
                    clusters.setdefault(rep, []).append(entry_by_id[member])
                    assigned.add(member)
    for sid, entry in entry_by_id.items():
        if sid not in assigned:
            clusters.setdefault(sid, []).append(entry)
    return clusters


def choose_cluster_representatives(
    entries, identity, metadata, work_name, membership_out=None, sensitivity=None,
):
    if not entries:
        return []

    work = cfg.INTER_DIR / work_name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    in_fa = work / "input.fa"
    cfg.write_fasta(entries, in_fa)
    prefix = work / "clust"
    tmp = work / "tmp"
    tmp.mkdir()
    cmd = [
        MMSEQS, "easy-cluster", str(in_fa), str(prefix), str(tmp),
        "--min-seq-id", str(identity),
        "-c", str(cfg.CLUSTER_COV),
        "--cluster-reassign", "1",
        "--threads", str(cfg.N_THREADS),
        "-v", "1",
    ]
    if sensitivity is not None:
        cmd.extend(["-s", str(sensitivity)])
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    clusters = parse_cluster_members(Path(str(prefix) + "_cluster.tsv"), entries)
    selected = []
    rows = []
    for cluster_id, members in sorted(clusters.items()):
        best = max(members, key=lambda entry: representative_score(entry, metadata))
        selected.append(best)
        for sid, _hdr, _seq in sorted(members, key=lambda x: x[0]):
            rows.append({
                "cluster_id": cluster_id,
                "selected_id": best[0],
                "member_id": sid,
                "cluster_size": str(len(members)),
            })
    if membership_out:
        with open(membership_out, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["cluster_id", "selected_id", "member_id", "cluster_size"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
    shutil.rmtree(work, ignore_errors=True)
    return selected


def precluster_master(metadata):
    entries = cfg.read_fasta(MASTER_FA)
    if cfg.PRECLUSTER_IDENTITY <= 0:
        print(f"Precluster: skipped for {cfg.ACTIVE_DB_NAME}")
        return entries
    if PRECLUSTER_FA.exists() and PRECLUSTER_FA.stat().st_mtime >= MASTER_FA.stat().st_mtime:
        print(f"Precluster: using existing {PRECLUSTER_FA.name}")
        return cfg.read_fasta(PRECLUSTER_FA)
    selected = choose_cluster_representatives(
        entries, cfg.PRECLUSTER_IDENTITY, metadata,
        "mmseqs_precluster", PRECLUSTER_MEMBERS,
    )
    cfg.write_fasta(selected, PRECLUSTER_FA)
    print(
        f"Precluster: {len(entries)} -> {len(selected)} "
        f"via MMseqs @{cfg.PRECLUSTER_IDENTITY}, metadata-selected reps"
    )
    return selected


def lowid_cluster(entries, identity, target_n, work_name, metadata):
    if len(entries) <= target_n:
        return list(entries)
    reps = choose_cluster_representatives(
        entries, identity, metadata, work_name,
        cfg.INTER_DIR / f"{work_name}_members.tsv",
        sensitivity=6.0,
    )
    reps = sorted(reps, key=lambda entry: representative_score(entry, metadata), reverse=True)
    if len(reps) > target_n:
        reps = reps[:target_n]
    return reps


def build_repset():
    metadata = load_metadata()
    entries = precluster_master(metadata)
    by_regime = {"psychro": [], "meso": [], "thermo": []}
    unlabelled = []
    for entry in entries:
        regime = metadata.get(entry[0], {}).get("regime", "")
        if regime in by_regime:
            by_regime[regime].append(entry)
        else:
            unlabelled.append(entry)

    keep_psychro = by_regime["psychro"]
    keep_meso = lowid_cluster(
        by_regime["meso"], cfg.REPSET_MESO_LOWID_CLUSTER,
        cfg.REPSET_MESO_QUOTA, "mmseqs_meso", metadata,
    )
    keep_thermo = lowid_cluster(
        by_regime["thermo"], cfg.REPSET_THERMO_LOWID_CLUSTER,
        cfg.REPSET_THERMO_QUOTA, "mmseqs_thermo", metadata,
    )
    reps = keep_psychro + keep_meso + keep_thermo

    query = cfg.read_fasta(cfg.QUERY_FASTA)[0]
    if not any(sid == query[0] for sid, _, _ in reps):
        reps.insert(0, query)

    selected_ids = {sid for sid, _, _ in reps}
    unlabelled = [entry for entry in unlabelled if entry[0] not in selected_ids]

    topup_n = max(0, cfg.TARGET_REPSET - len(reps))
    if topup_n:
        reps.extend(lowid_cluster(unlabelled, cfg.REPSET_LOWID_CLUSTER, topup_n, "mmseqs_topup", metadata))

    cfg.write_fasta(reps, REPSET_FA)
    with open(MASTER_TSV, newline="") as f:
        fieldnames = csv.DictReader(f, delimiter="\t").fieldnames or []
    with open(REPSET_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sid, _hdr, _seq in reps:
            writer.writerow(metadata.get(sid, {field: "" for field in fieldnames}) | {"id": sid})

    print(
        f"Repset: {len(reps)} sequences "
        f"(psychro={len(keep_psychro)}, meso={len(keep_meso)}, "
        f"thermo={len(keep_thermo)}, unlabelled_topup={max(0, len(reps) - 1 - len(keep_psychro) - len(keep_meso) - len(keep_thermo))})"
    )
    print(f"Saved: {REPSET_FA}")
    print(f"Saved: {REPSET_TSV}")


def build_seed_alignment():
    HMM_DIR.mkdir(parents=True, exist_ok=True)
    if SEED_ALN.exists() and SEED_ALN.stat().st_mtime >= SEED_FA.stat().st_mtime:
        return
    with open(SEED_ALN, "w") as out:
        subprocess.run([
            MAFFT, "--localpair", "--maxiterate", "1000",
            "--thread", str(cfg.N_THREADS), str(SEED_FA),
        ], stdout=out, stderr=subprocess.PIPE, text=True, check=True)


def build_profile():
    if PROFILE_HMM.exists() and PROFILE_HMM.stat().st_mtime >= SEED_ALN.stat().st_mtime:
        return
    subprocess.run([
        HMMBUILD, "--cpu", str(cfg.N_THREADS),
        "-n", "classI_CPD", str(PROFILE_HMM), str(SEED_ALN),
    ], check=True)


def align_repset():
    with open(HMM_ALIGN_RAW, "w") as out:
        subprocess.run([
            HMMALIGN, "--outformat", "afa", str(PROFILE_HMM), str(REPSET_FA),
        ], stdout=out, stderr=subprocess.PIPE, text=True, check=True)
    entries = cfg.read_fasta(HMM_ALIGN_RAW)
    cleaned = [(sid, hdr, seq.upper().replace(".", "-")) for sid, hdr, seq in entries]
    write_aligned_fasta(cleaned, ALIGNED_FA)
    if cleaned:
        insert_cols = {
            i
            for i, col in enumerate(zip(*(seq for _sid, _hdr, seq in entries)))
            if any(char.islower() for char in col)
        }
        match_cols = [
            i for i in range(len(cleaned[0][2]))
            if i not in insert_cols
        ]
        match_only = [
            (sid, hdr, "".join(seq[i] for i in match_cols))
            for sid, hdr, seq in cleaned
        ]
        write_aligned_fasta(match_only, ALIGNED_MATCH_FA)
    print(f"Saved: {ALIGNED_FA}")
    print(f"Saved: {ALIGNED_MATCH_FA}")


def main():
    if not MASTER_FA.exists() or not MASTER_TSV.exists():
        sys.exit("ERROR: run 04_build_master_set.py first")
    require_bins(["mmseqs", "mafft", "hmmbuild", "hmmalign"])
    build_repset()
    build_seed_alignment()
    build_profile()
    align_repset()


if __name__ == "__main__":
    main()
