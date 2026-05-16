#!/usr/bin/env bash
# Reproducible setup for the ANKros structural/MSA/MD environment.
#
# Creates exactly the same software stack on any machine with conda + an
# NVIDIA GPU that has driver >= 525 (covers CUDA 12.x runtime). Tested on
# RTX 3060 (CC 8.6) and GTX 1060 (CC 6.1). Idempotent — re-running is safe.
#
# What it does:
#   1. Extract every AF3 fold from structures/af3/folds_*.zip
#      (photohymenobacter_* drives the unified pipeline; the homolog/oligomer
#      folders feed scripts/experiments/structural_validation/).
#   2. Create envs/ankros from config/environment.yml
#   3. Optionally create envs/amber_cuda and build local Amber/pmemd CUDA
#   4. Optionally create envs/aggrescan3d from config/aggrescan3d_environment.yml
#   5. Optionally clone ThermoMPNN and create envs/thermompnn
#   6. Smoke-test envs (python imports + CLI checks)
#
# Usage:
#   bash setup_envs.sh
#   bash setup_envs.sh --skip-ankros    # optional envs only
#   AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash setup_envs.sh --with-amber-cuda
#   bash setup_envs.sh --with-amber-cuda --amber-cuda-jobs 2
#   bash setup_envs.sh --with-aggrescan3d
#   bash setup_envs.sh --with-thermompnn
#
# Logs go to logs/setup_envs.log. Heavy step is the conda solve for
# envs/ankros (30-90 min). The optional Amber CUDA build is also long-running
# and should be launched inside tmux.

set -euo pipefail

# Pinned versions ─────────────────────────────────────────────────────────────
THERMOMPNN_REPO="https://github.com/Kuhlman-Lab/ThermoMPNN"
THERMOMPNN_COMMIT="2b04fd370e399911b1fa5848112cc9013f084110"  # main, pinned 2026-05-06

# Paths ───────────────────────────────────────────────────────────────────────
THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${THIS}"
CONFIG_DIR="${ROOT}/config"
ANKROS_YML="${CONFIG_DIR}/environment.yml"
AMBER_CUDA_YML="${CONFIG_DIR}/amber_cuda_environment.yml"
AGGRESCAN3D_YML="${CONFIG_DIR}/aggrescan3d_environment.yml"
THERMOMPNN_YML="${CONFIG_DIR}/thermompnn_environment.yml"
ANKROS_ENV="${ROOT}/envs/ankros"
AMBER_CUDA_ENV="${ROOT}/envs/amber_cuda"
AGGRESCAN3D_ENV="${ROOT}/envs/aggrescan3d"
THERMOMPNN_ENV="${ROOT}/envs/thermompnn"
THERMOMPNN_DIR="${ROOT}/external/ThermoMPNN"
AMBER_CUDA_INSTALL="${CONFIG_DIR}/install_amber_cuda.sh"
AF3_DIR="${ROOT}/structures/af3"
LOG_DIR="${ROOT}/logs"
LOG="${LOG_DIR}/setup_envs.log"

# CLI flags ───────────────────────────────────────────────────────────────────
DO_ANKROS=1
DO_AMBER_CUDA=0
DO_AGGRESCAN3D=0
DO_THERMOMPNN=0
AMBER_CUDA_JOBS=""
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --skip-ankros)    DO_ANKROS=0 ;;
        --with-amber-cuda) DO_AMBER_CUDA=1 ;;
        --with-aggrescan3d) DO_AGGRESCAN3D=1 ;;
        --with-thermompnn) DO_THERMOMPNN=1 ;;
        --amber-cuda-jobs)
            shift
            if [[ "$#" -eq 0 ]]; then
                echo "ERROR: --amber-cuda-jobs requires a positive integer value" >&2
                exit 2
            fi
            AMBER_CUDA_JOBS="$1"
            ;;
        --amber-cuda-jobs=*)
            AMBER_CUDA_JOBS="${1#*=}"
            ;;
        -h|--help)
            awk '/^set -euo pipefail/{exit} /^# /{sub(/^# /, ""); print} /^#$/{print ""}' "$0"
            exit 0
            ;;
        *) echo "unknown flag: $1" >&2; exit 2 ;;
    esac
    shift
done

if [[ -n "${AMBER_CUDA_JOBS}" && ! "${AMBER_CUDA_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --amber-cuda-jobs must be a positive integer, got: ${AMBER_CUDA_JOBS}" >&2
    exit 2
fi

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG}") 2>&1

echo "==[$(date -Is)] setup_envs.sh starting in ${ROOT}"

# 0. Sanity ───────────────────────────────────────────────────────────────────
command -v conda >/dev/null || { echo "ERROR: conda not on PATH"; exit 1; }
command -v git   >/dev/null || { echo "ERROR: git not on PATH"; exit 1; }

# Prefer mamba if available — solves are 5-10× faster, identical result given
# fully-pinned envs.
SOLVER="conda"
if command -v mamba >/dev/null; then
    SOLVER="mamba"
fi
echo "    solver: ${SOLVER}"

# Show GPU+driver so reproducibility is self-documenting in the log.
if command -v nvidia-smi >/dev/null; then
    nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total \
               --format=csv,noheader || true
fi

# 1. Extract AF3 CIFs ─────────────────────────────────────────────────────────
# Extract the full archive: photohymenobacter_* drives the unified pipeline,
# and the homolog/oligomer folders (1dnp_*, 1tez_1owl_*, 6kii_*, gst_dimer,
# ldha_*, mbp_*, phr_*, phrb_*, etc.) are inputs for
# scripts/experiments/structural_validation/.
echo
echo "== 1. AF3 zip extraction =="
af3_zip="$(ls "${AF3_DIR}"/folds_*.zip 2>/dev/null | head -1 || true)"
if [[ -z "${af3_zip}" ]]; then
    echo "    no folds_*.zip found in ${AF3_DIR} (skipping; may already be extracted)"
else
    echo "    extracting all entries from $(basename "${af3_zip}") (-n: skip existing)"
    ( cd "${AF3_DIR}" && unzip -nq "$(basename "${af3_zip}")" )
fi

# 2. Create envs/ankros ───────────────────────────────────────────────────────
if [[ "${DO_ANKROS}" -eq 1 ]]; then
    echo
    echo "== 2. envs/ankros (heavy: 30-90 min) =="
    if [[ -x "${ANKROS_ENV}/bin/python" ]]; then
        echo "    already exists at ${ANKROS_ENV} — skipping"
    else
        echo "    creating from ${ANKROS_YML}"
        "${SOLVER}" env create -f "${ANKROS_YML}" -p "${ANKROS_ENV}" -y
    fi
else
    echo "    (--skip-ankros)"
fi

# 3. Optional Amber/pmemd CUDA build ──────────────────────────────────────────
if [[ "${DO_AMBER_CUDA}" -eq 1 ]]; then
    echo
    echo "== 3. optional Amber/pmemd CUDA build =="
    if [[ -d "${AMBER_CUDA_ENV}/conda-meta" ]]; then
        echo "    updating existing Amber CUDA build env at ${AMBER_CUDA_ENV}"
        "${SOLVER}" env update -f "${AMBER_CUDA_YML}" -p "${AMBER_CUDA_ENV}" --prune
    else
        echo "    creating Amber CUDA build env from ${AMBER_CUDA_YML}"
        "${SOLVER}" env create -f "${AMBER_CUDA_YML}" -p "${AMBER_CUDA_ENV}" -y
    fi
    amber_cuda_env=(ANKROS_AMBER_CUDA_ENV_PREFIX="${AMBER_CUDA_ENV}")
    if [[ -n "${AMBER_CUDA_JOBS}" ]]; then
        amber_cuda_env+=(ANKROS_AMBER_CUDA_JOBS="${AMBER_CUDA_JOBS}")
    fi
    env "${amber_cuda_env[@]}" "${AMBER_CUDA_INSTALL}" "${AMBER_CUDA_PACKAGE:-${HOME}/Downloads/pmemd26.tar.bz2}"
else
    echo "    (--with-amber-cuda not requested)"
fi

# 4. Optional Aggrescan3D environment ─────────────────────────────────────────
if [[ "${DO_AGGRESCAN3D}" -eq 1 ]]; then
    echo
    echo "== 4. optional envs/aggrescan3d =="
    if [[ -x "${AGGRESCAN3D_ENV}/bin/aggrescan" ]]; then
        echo "    already exists at ${AGGRESCAN3D_ENV} — skipping"
    else
        echo "    creating from ${AGGRESCAN3D_YML}"
        "${SOLVER}" env create -f "${AGGRESCAN3D_YML}" -p "${AGGRESCAN3D_ENV}" -y
    fi
else
    echo "    (--with-aggrescan3d not requested)"
fi

# 5. Optional ThermoMPNN environment ──────────────────────────────────────────
if [[ "${DO_THERMOMPNN}" -eq 1 ]]; then
    echo
    echo "== 5. optional ThermoMPNN @ ${THERMOMPNN_COMMIT:0:8} =="
    if [[ -d "${THERMOMPNN_DIR}/.git" ]]; then
        have="$(git -C "${THERMOMPNN_DIR}" rev-parse HEAD)"
        if [[ "${have}" == "${THERMOMPNN_COMMIT}" ]]; then
            echo "    already at pinned commit"
        else
            echo "    advancing to pinned commit (current: ${have:0:8})"
            git -C "${THERMOMPNN_DIR}" fetch --depth 1 origin "${THERMOMPNN_COMMIT}"
            git -C "${THERMOMPNN_DIR}" checkout --detach "${THERMOMPNN_COMMIT}"
        fi
    else
        mkdir -p "$(dirname "${THERMOMPNN_DIR}")"
        git clone "${THERMOMPNN_REPO}" "${THERMOMPNN_DIR}"
        git -C "${THERMOMPNN_DIR}" checkout --detach "${THERMOMPNN_COMMIT}"
    fi
    git -C "${THERMOMPNN_DIR}" log -1 --pretty='    HEAD: %H  %s'
    mkdir -p "${THERMOMPNN_DIR}/cache"
    cat > "${THERMOMPNN_DIR}/local.yaml" <<EOF
platform:
  accel: "gpu"
  cache_dir: "${THERMOMPNN_DIR}/cache"
  thermompnn_dir: "${THERMOMPNN_DIR}"

data_loc:
  megascale_csv: ""
  megascale_splits: ""
  megascale_pdbs: ""
  fireprot_csv: ""
  fireprot_splits: ""
  fireprot_pdbs: ""
EOF

    if [[ -x "${THERMOMPNN_ENV}/bin/python" ]]; then
        echo "    already exists at ${THERMOMPNN_ENV} — skipping"
    else
        echo "    creating from ${THERMOMPNN_YML}"
        "${SOLVER}" env create -f "${THERMOMPNN_YML}" -p "${THERMOMPNN_ENV}" -y
    fi
else
    echo "    (--with-thermompnn not requested)"
fi

# 6. Smoke tests ─────────────────────────────────────────────────────────────
echo
echo "== 6. smoke tests =="

if [[ "${DO_ANKROS}" -eq 1 && -x "${ANKROS_ENV}/bin/python" ]]; then
    echo "    -- ankros --"
    "${ANKROS_ENV}/bin/python" -c "
import sys, importlib
mods = ['Bio', 'pymol', 'openmm', 'rdkit', 'numpy', 'scipy', 'matplotlib']
for m in mods:
    importlib.import_module(m)
print('    python:', sys.version.split()[0])
import torch
print('    torch :', torch.__version__, 'cuda:', torch.cuda.is_available())
"
    "${ANKROS_ENV}/bin/mkdssp" --version
    "${ANKROS_ENV}/bin/foldseek" version
fi

if [[ "${DO_AGGRESCAN3D}" -eq 1 && -x "${AGGRESCAN3D_ENV}/bin/aggrescan" ]]; then
    echo "    -- aggrescan3d --"
    mkdir -p "${LOG_DIR}/.matplotlib"
    MPLCONFIGDIR="${LOG_DIR}/.matplotlib" "${AGGRESCAN3D_ENV}/bin/python" -c "
import sys, importlib.metadata
print('    python:', sys.version.split()[0])
print('    Aggrescan3D:', importlib.metadata.version('Aggrescan3D'))
"
    MPLCONFIGDIR="${LOG_DIR}/.matplotlib" "${AGGRESCAN3D_ENV}/bin/aggrescan" --help >/dev/null
    echo "    aggrescan CLI: ok"
fi

if [[ "${DO_THERMOMPNN}" -eq 1 && -x "${THERMOMPNN_ENV}/bin/python" ]]; then
    echo "    -- thermompnn --"
    "${THERMOMPNN_ENV}/bin/python" -c "
import sys
import torch, pandas, omegaconf, Bio
print('    python:', sys.version.split()[0])
print('    torch :', torch.__version__, 'cuda:', torch.cuda.is_available())
"
    test -f "${THERMOMPNN_DIR}/analysis/custom_inference.py"
    test -f "${THERMOMPNN_DIR}/models/thermoMPNN_default.pt" -o -f "${THERMOMPNN_DIR}/models/thermoMPNN_default.ckpt"
    echo "    ThermoMPNN checkout: ok"
fi

echo
echo "==[$(date -Is)] setup_envs.sh done. Log: ${LOG}"
