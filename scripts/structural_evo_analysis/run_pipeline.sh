#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
ORIGINAL_CWD="$(pwd -P)"

resolve_input_path() {
  local path="$1"
  if [[ "${path}" = /* ]]; then
    printf '%s\n' "${path}"
  elif [[ -e "${ORIGINAL_CWD}/${path}" ]]; then
    printf '%s\n' "${ORIGINAL_CWD}/${path}"
  else
    printf '%s\n' "${ROOT}/${path}"
  fi
}

resolve_output_path() {
  local path="$1"
  if [[ "${path}" = /* ]]; then
    printf '%s\n' "${path}"
  else
    printf '%s\n' "${ROOT}/${path}"
  fi
}

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <query.fasta> [metadata.tsv] [trait_column]" >&2
  echo "Modes: SEA_PIPELINE_MODE=full|msa|conservation (default: full)" >&2
  echo "Defaults: SEA_WORK_DIR=${SEA_WORK_DIR:-${HOME}/structural_evo_analysis}" >&2
  echo "Set SEA_WORK_DIR, SEA_OUT_DIR, or SEA_STRUCTURE_DIR to change run storage." >&2
  echo "Full mode requires query PDB with SEA_QUERY_PDB or ${SEA_STRUCTURE_DIR:-${SEA_WORK_DIR:-${HOME}/structural_evo_analysis}/structures}/query.pdb." >&2
  exit 2
fi

QUERY_FASTA="$(resolve_input_path "$1")"
METADATA_TSV="${2:-}"
if [[ -n "${METADATA_TSV}" ]]; then
  METADATA_TSV="$(resolve_input_path "${METADATA_TSV}")"
fi
TRAIT_COLUMN="${3:-}"
PYTHON="${PYTHON:-${ROOT}/envs/structural_evo/bin/python}"
WORK_DIR="${SEA_WORK_DIR:-${HOME}/structural_evo_analysis}"
OUT_DIR="$(resolve_output_path "${SEA_OUT_DIR:-${WORK_DIR}/results}")"
STRUCTURE_DIR="$(resolve_output_path "${SEA_STRUCTURE_DIR:-${WORK_DIR}/structures}")"
PIPELINE_MODE="${SEA_PIPELINE_MODE:-full}"
case "${PIPELINE_MODE}" in
  full|msa|conservation) ;;
  *)
    echo "ERROR: unsupported SEA_PIPELINE_MODE=${PIPELINE_MODE}; choose full, msa, or conservation." >&2
    exit 2
    ;;
esac
QUERY_PDB=""
if [[ "${PIPELINE_MODE}" == "full" ]]; then
  QUERY_PDB="$(resolve_input_path "${SEA_QUERY_PDB:-${STRUCTURE_DIR}/query.pdb}")"
  if [[ ! -s "${QUERY_PDB}" ]]; then
    echo "ERROR: full vulnerability mode requires query PDB: ${QUERY_PDB}" >&2
    echo "Set SEA_QUERY_PDB, place query.pdb under ${STRUCTURE_DIR}, or use SEA_PIPELINE_MODE=msa/conservation." >&2
    exit 2
  fi
fi
TRAIT_TYPE="${SEA_TRAIT_TYPE:-continuous}"
OGT_AWARE="${SEA_OGT_AWARE:-0}"
LOW_THRESHOLD="${SEA_LOW_THRESHOLD:-20}"
HIGH_THRESHOLD="${SEA_HIGH_THRESHOLD:-45}"
JOIN_OGT=()
if [[ "${OGT_AWARE}" == "1" || "${TRAIT_COLUMN}" == "ogt" ]]; then
  JOIN_OGT=(--join-ogt)
fi

"${PYTHON}" "${SCRIPT_DIR}/01_mmseqs_search.py" \
  --query "${QUERY_FASTA}" \
  --out-dir "${OUT_DIR}" \
  "${JOIN_OGT[@]}"

"${PYTHON}" "${SCRIPT_DIR}/02_align_mafft.py" \
  --out-dir "${OUT_DIR}"

if [[ "${PIPELINE_MODE}" == "msa" ]]; then
  echo "SEA_PIPELINE_MODE=msa complete. Outputs: ${OUT_DIR}/repset.fa and ${OUT_DIR}/repset_aligned.fa"
  exit 0
fi

"${PYTHON}" "${SCRIPT_DIR}/03_build_tree_iqtree.py" \
  --out-dir "${OUT_DIR}/tree"

if [[ "${OGT_AWARE}" == "1" && -z "${METADATA_TSV}" ]]; then
  METADATA_TSV="${OUT_DIR}/repset_metadata.tsv"
  TRAIT_COLUMN="ogt"
fi

CALLED_CLADES=()
TIP_METADATA=()
if [[ -n "${METADATA_TSV}" && -n "${TRAIT_COLUMN}" ]]; then
  annotate_args=(
    --tree "${OUT_DIR}/tree/query_msa.treefile"
    --metadata "${METADATA_TSV}"
    --trait-column "${TRAIT_COLUMN}"
    --trait-type "${TRAIT_TYPE}"
    --out-dir "${OUT_DIR}/metadata_clades"
  )
  annotate_args+=(--low-threshold "${LOW_THRESHOLD}" --high-threshold "${HIGH_THRESHOLD}")
  "${PYTHON}" "${SCRIPT_DIR}/04_annotate_clades.py" \
    "${annotate_args[@]}"
  CALLED_CLADES=(--called-clades "${OUT_DIR}/metadata_clades/called_clades.tsv")
  TIP_METADATA=(--tip-metadata "${OUT_DIR}/metadata_clades/tip_metadata.tsv")
fi

"${PYTHON}" "${SCRIPT_DIR}/05_conserved_positions.py" \
  --alignment "${OUT_DIR}/repset_aligned.fa" \
  --out-dir "${OUT_DIR}/conservation"

if [[ "${PIPELINE_MODE}" == "conservation" ]]; then
  echo "SEA_PIPELINE_MODE=conservation complete. Outputs: ${OUT_DIR}/conservation"
  exit 0
fi

"${PYTHON}" "${SCRIPT_DIR}/06_download_afdb.py" \
  --metadata "${OUT_DIR}/repset_metadata.tsv" \
  --dest "${STRUCTURE_DIR}/afdb" \
  --manifest "${OUT_DIR}/afdb_downloads/download_manifest.tsv"

"${PYTHON}" "${SCRIPT_DIR}/07_score_structures.py" \
  --afdb-dir "${STRUCTURE_DIR}/afdb" \
  --query-pdb "${QUERY_PDB}" \
  --out-dir "${OUT_DIR}/structure_scores"

"${PYTHON}" "${SCRIPT_DIR}/08_vulnerability_analysis.py" \
  --conservation "${OUT_DIR}/conservation/position_conservation.tsv" \
  --scores "${OUT_DIR}/structure_scores/per_residue_scores.tsv" \
  --query-pdb "${QUERY_PDB}" \
  --out-dir "${OUT_DIR}/vulnerability"

"${PYTHON}" "${SCRIPT_DIR}/09_group_score_summary.py" \
  --metadata "${OUT_DIR}/repset_metadata.tsv" \
  --scores "${OUT_DIR}/structure_scores/global_scores.tsv" \
  --query-pdb "${QUERY_PDB}" \
  --out-dir "${OUT_DIR}/vulnerability" \
  "${CALLED_CLADES[@]}"

"${PYTHON}" "${SCRIPT_DIR}/10_sequence_logos.py" \
  --alignment "${OUT_DIR}/repset_aligned.fa" \
  --group-members "${OUT_DIR}/vulnerability/group_members.tsv" \
  --out-dir "${OUT_DIR}/logos"

"${PYTHON}" "${SCRIPT_DIR}/11_write_viewers.py" \
  --alignment "${OUT_DIR}/repset_aligned.fa" \
  --tree "${OUT_DIR}/tree/query_msa.treefile" \
  --metadata "${OUT_DIR}/repset_metadata.tsv" \
  --out-dir "${OUT_DIR}/viewers" \
  "${TIP_METADATA[@]}"
