#!/usr/bin/env bash
# Reproducible setup for the structural-evolution skill pipeline.
#
# Creates a minimal local conda environment at envs/structural_evo from
# config/environment.yml. Optional Aggrescan3D support is installed in a
# separate env because its conda package comes from a dedicated channel.
#
# Usage:
#   bash setup_envs.sh
#   bash setup_envs.sh --with-aggrescan3d
#   bash setup_envs.sh --skip-main --with-aggrescan3d
#
# Logs go to logs/setup_envs.log. Run in tmux if conda solves are slow:
#   tmux new -s structural-evo-setup 'bash setup_envs.sh 2>&1 | tee logs/setup_envs.tmux.log'

set -euo pipefail

THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${THIS}"
CONFIG_DIR="${ROOT}/config"
MAIN_YML="${CONFIG_DIR}/environment.yml"
AGGRESCAN3D_YML="${CONFIG_DIR}/aggrescan3d_environment.yml"
MAIN_ENV="${ROOT}/envs/structural_evo"
AGGRESCAN3D_ENV="${ROOT}/envs/aggrescan3d"
LOG_DIR="${ROOT}/logs"
LOG="${LOG_DIR}/setup_envs.log"

DO_MAIN=1
DO_AGGRESCAN3D=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --skip-main|--skip-structural-evo)
            DO_MAIN=0
            ;;
        --with-aggrescan3d)
            DO_AGGRESCAN3D=1
            ;;
        -h|--help)
            awk '/^set -euo pipefail/{exit} /^# /{sub(/^# /, ""); print} /^#$/{print ""}' "$0"
            exit 0
            ;;
        *)
            echo "unknown flag: $1" >&2
            exit 2
            ;;
    esac
    shift
done

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG}") 2>&1

echo "==[$(date -Is)] setup_envs.sh starting in ${ROOT}"

[[ -f "${MAIN_YML}" ]] || { echo "ERROR: missing environment file: ${MAIN_YML}"; exit 1; }
[[ -f "${AGGRESCAN3D_YML}" ]] || { echo "ERROR: missing environment file: ${AGGRESCAN3D_YML}"; exit 1; }

SOLVER=""
if command -v mamba >/dev/null; then
    SOLVER="mamba"
elif command -v conda >/dev/null; then
    SOLVER="conda"
elif command -v micromamba >/dev/null; then
    SOLVER="micromamba"
else
    echo "ERROR: conda, mamba, or micromamba not on PATH"
    exit 1
fi
echo "    solver: ${SOLVER}"

if [[ "${DO_MAIN}" -eq 1 ]]; then
    echo
    echo "== main env: ${MAIN_ENV} =="
    if [[ -x "${MAIN_ENV}/bin/python" ]]; then
        echo "    updating existing env from ${MAIN_YML}"
        "${SOLVER}" env update -f "${MAIN_YML}" -p "${MAIN_ENV}" --prune -y
    else
        echo "    creating env from ${MAIN_YML}"
        "${SOLVER}" env create -f "${MAIN_YML}" -p "${MAIN_ENV}" -y
    fi
else
    echo "    (--skip-main)"
fi

if [[ "${DO_AGGRESCAN3D}" -eq 1 ]]; then
    echo
    echo "== optional Aggrescan3D env: ${AGGRESCAN3D_ENV} =="
    if [[ -x "${AGGRESCAN3D_ENV}/bin/aggrescan" ]]; then
        echo "    updating existing env from ${AGGRESCAN3D_YML}"
        "${SOLVER}" env update -f "${AGGRESCAN3D_YML}" -p "${AGGRESCAN3D_ENV}" --prune -y
    else
        echo "    creating env from ${AGGRESCAN3D_YML}"
        "${SOLVER}" env create -f "${AGGRESCAN3D_YML}" -p "${AGGRESCAN3D_ENV}" -y
    fi
else
    echo "    (--with-aggrescan3d not requested)"
fi

echo
echo "== smoke tests =="
if [[ "${DO_MAIN}" -eq 1 && -x "${MAIN_ENV}/bin/python" ]]; then
    "${MAIN_ENV}/bin/python" -c "
import sys
import matplotlib
from Bio import Phylo
from Bio.PDB import PDBParser
print('    python:', sys.version.split()[0])
print('    biopython imports: ok')
print('    matplotlib:', matplotlib.__version__)
"
    "${MAIN_ENV}/bin/mmseqs" version
    "${MAIN_ENV}/bin/mafft" --version >/dev/null
    "${MAIN_ENV}/bin/iqtree" --version | head -1
fi

if [[ "${DO_AGGRESCAN3D}" -eq 1 && -x "${AGGRESCAN3D_ENV}/bin/aggrescan" ]]; then
    mkdir -p "${LOG_DIR}/.matplotlib"
    MPLCONFIGDIR="${LOG_DIR}/.matplotlib" "${AGGRESCAN3D_ENV}/bin/aggrescan" --help >/dev/null
    echo "    aggrescan CLI: ok"
fi

echo
echo "==[$(date -Is)] setup_envs.sh done. Log: ${LOG}"
