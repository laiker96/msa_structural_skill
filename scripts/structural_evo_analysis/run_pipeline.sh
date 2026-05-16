#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <query.fasta> [metadata.tsv] [trait_column]" >&2
  echo "Set SEA_OUT_DIR to change the output directory." >&2
  exit 2
fi

QUERY_FASTA="$1"
METADATA_TSV="${2:-}"
TRAIT_COLUMN="${3:-}"
PYTHON="${PYTHON:-./envs/structural_evo/bin/python}"
OUT_DIR="${SEA_OUT_DIR:-results/structural_evo_analysis}"
TRAIT_TYPE="${SEA_TRAIT_TYPE:-continuous}"
LOW_THRESHOLD="${SEA_LOW_THRESHOLD:-}"
HIGH_THRESHOLD="${SEA_HIGH_THRESHOLD:-}"

"${PYTHON}" scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query "${QUERY_FASTA}" \
  --out-dir "${OUT_DIR}"

"${PYTHON}" scripts/structural_evo_analysis/02_align_mafft.py \
  --out-dir "${OUT_DIR}"

"${PYTHON}" scripts/structural_evo_analysis/03_build_tree_iqtree.py \
  --out-dir "${OUT_DIR}/tree"

if [[ -n "${METADATA_TSV}" && -n "${TRAIT_COLUMN}" ]]; then
  annotate_args=(
    --tree "${OUT_DIR}/tree/query_msa.treefile"
    --metadata "${METADATA_TSV}"
    --trait-column "${TRAIT_COLUMN}"
    --trait-type "${TRAIT_TYPE}"
    --out-dir "${OUT_DIR}/metadata_clades"
  )
  if [[ -n "${LOW_THRESHOLD}" && -n "${HIGH_THRESHOLD}" ]]; then
    annotate_args+=(--low-threshold "${LOW_THRESHOLD}" --high-threshold "${HIGH_THRESHOLD}")
  fi
  "${PYTHON}" scripts/structural_evo_analysis/04_annotate_clades.py \
    "${annotate_args[@]}"
fi

"${PYTHON}" scripts/structural_evo_analysis/05_conserved_positions.py \
  --alignment "${OUT_DIR}/repset_aligned.fa" \
  --out-dir "${OUT_DIR}/conservation"

"${PYTHON}" scripts/structural_evo_analysis/06_download_afdb.py \
  --metadata "${OUT_DIR}/repset_metadata.tsv" \
  --dest "structures/structural_evo_analysis/afdb" \
  --manifest "${OUT_DIR}/afdb_downloads/download_manifest.tsv"

"${PYTHON}" scripts/structural_evo_analysis/07_score_structures.py \
  --afdb-dir "structures/structural_evo_analysis/afdb" \
  --out-dir "${OUT_DIR}/structure_scores"
