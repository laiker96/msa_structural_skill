#!/usr/bin/env bash
# Build a local, untracked Amber/pmemd CUDA source archive into envs/amber_cuda.
#
# The Amber source archive is licensed external input and must not be tracked.
# Default input path is ~/Downloads/pmemd26.tar.bz2; override with either:
#
#   AMBER_CUDA_PACKAGE=/path/to/pmemd26.tar.bz2 bash config/install_amber_cuda.sh
#   bash config/install_amber_cuda.sh /path/to/pmemd26.tar.bz2
#
# This is a long build. Run from tmux if executing interactively:
#
#   tmux new -s amber-cuda 'AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash config/install_amber_cuda.sh'

set -euo pipefail

THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${THIS}/.." && pwd)"

PACKAGE="${1:-${AMBER_CUDA_PACKAGE:-${HOME}/Downloads/pmemd26.tar.bz2}}"
package_file="$(basename "${PACKAGE}")"
package_stem="${package_file}"
for suffix in .tar.bz2 .tbz2 .tar.gz .tgz .tar.xz .txz .tar; do
    if [[ "${package_stem}" == *"${suffix}" ]]; then
        package_stem="${package_stem%"${suffix}"}"
        break
    fi
done

AMBER_PACKAGE_NAME="${ANKROS_AMBER_CUDA_PACKAGE_NAME:-${package_stem}}"
ENV_PREFIX="${ANKROS_AMBER_CUDA_ENV_PREFIX:-${ANKROS_ENV_PREFIX:-${ROOT}/envs/amber_cuda}}"
WORK_DIR="${ANKROS_AMBER_CUDA_WORK_DIR:-${ROOT}/external/amber_cuda}"
SRC_DIR="${ANKROS_AMBER_CUDA_SRC_DIR:-${WORK_DIR}/${AMBER_PACKAGE_NAME}_src}"

detected_gpu_targets=""
if command -v nvidia-smi >/dev/null; then
    while IFS= read -r compute_cap; do
        compute_cap="${compute_cap//[[:space:]]/}"
        if [[ "${compute_cap}" =~ ^[0-9]+[.][0-9]+$ ]]; then
            sm="sm_${compute_cap/.}"
            case " ${detected_gpu_targets} " in
                *" ${sm} "*) ;;
                *) detected_gpu_targets="${detected_gpu_targets:+${detected_gpu_targets} }${sm}" ;;
            esac
        fi
    done < <(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null || true)
fi
GPU_TARGETS="${ANKROS_AMBER_CUDA_GPU_TARGETS:-${detected_gpu_targets}}"
GPU_TARGET_LABEL="${GPU_TARGETS// /-}"
ONLY_GPU_TARGETS="${ANKROS_AMBER_CUDA_ONLY_GPU_TARGETS:-1}"
GPU_TARGET_MODE_LABEL=""
GPU_TARGET_MODE="amber-default-plus-overrides"
if [[ -n "${GPU_TARGETS}" && "${ONLY_GPU_TARGETS}" != "0" ]]; then
    GPU_TARGET_MODE_LABEL="-only"
    GPU_TARGET_MODE="only"
fi

BUILD_DIR="${ANKROS_AMBER_CUDA_BUILD_DIR:-${WORK_DIR}/build-${AMBER_PACKAGE_NAME}${GPU_TARGET_LABEL:+-${GPU_TARGET_LABEL}}${GPU_TARGET_MODE_LABEL}}"
INSTALL_PREFIX="${ANKROS_AMBER_CUDA_PREFIX:-${ENV_PREFIX}/opt/${AMBER_PACKAGE_NAME}}"
LOG_DIR="${ROOT}/logs"
LOG="${LOG_DIR}/install_amber_cuda.log"
JOBS="${ANKROS_AMBER_CUDA_JOBS:-2}"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG}") 2>&1

echo "==[$(date -Is)] install_amber_cuda.sh starting"
echo "    root           : ${ROOT}"
echo "    package        : ${PACKAGE}"
echo "    package name   : ${AMBER_PACKAGE_NAME}"
echo "    env prefix     : ${ENV_PREFIX}"
echo "    work dir       : ${WORK_DIR}"
echo "    source dir     : ${SRC_DIR}"
echo "    build dir      : ${BUILD_DIR}"
echo "    install prefix : ${INSTALL_PREFIX}"
echo "    gpu targets    : ${GPU_TARGETS:-Amber default}"
echo "    target mode    : ${GPU_TARGET_MODE}"
echo "    jobs           : ${JOBS}"
echo "    log            : ${LOG}"

if [[ -z "${TMUX:-}" ]]; then
    echo "WARNING: this build may run for a long time; repo policy recommends running it inside tmux."
fi

if [[ ! -f "${PACKAGE}" ]]; then
    echo "ERROR: Amber CUDA package not found: ${PACKAGE}" >&2
    echo "       Set AMBER_CUDA_PACKAGE or pass the path as the first argument." >&2
    exit 2
fi

if [[ ! -d "${ENV_PREFIX}/conda-meta" ]]; then
    echo "ERROR: env prefix does not look usable: ${ENV_PREFIX}" >&2
    echo "       Create it first with: mamba env create -f config/amber_cuda_environment.yml -p ./envs/amber_cuda" >&2
    exit 2
fi

# Prefer compilers, CUDA tools, CMake, and libraries from the target conda env.
export CONDA_PREFIX="${ENV_PREFIX}"
if [[ -d "${ENV_PREFIX}/targets/x86_64-linux" ]]; then
    export CUDA_HOME="${CUDA_HOME:-${ENV_PREFIX}/targets/x86_64-linux}"
else
    export CUDA_HOME="${CUDA_HOME:-${ENV_PREFIX}}"
fi
export CUDA_ROOT="${CUDA_ROOT:-${CUDA_HOME}}"
export PATH="${ENV_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${CUDA_HOME}/lib:${CUDA_HOME}/lib64:${ENV_PREFIX}/lib:${ENV_PREFIX}/lib64:${LD_LIBRARY_PATH:-}"
export CUDA_NVCC_EXECUTABLE="${CUDA_NVCC_EXECUTABLE:-${ENV_PREFIX}/bin/nvcc}"

HOST_CC="${ANKROS_AMBER_CC:-/usr/bin/gcc}"
HOST_CXX="${ANKROS_AMBER_CXX:-/usr/bin/g++}"
HOST_FC="${ANKROS_AMBER_FC:-/usr/bin/gfortran}"
export CC="${HOST_CC}"
export CXX="${HOST_CXX}"
export FC="${HOST_FC}"

for exe in cmake make nvcc; do
    if ! command -v "${exe}" >/dev/null; then
        echo "ERROR: required build tool not found on PATH: ${exe}" >&2
        echo "       Install build tools in envs/amber_cuda or on the host, then re-run." >&2
        exit 2
    fi
done

for exe in "${HOST_CC}" "${HOST_CXX}" "${HOST_FC}"; do
    if [[ ! -x "${exe}" ]]; then
        echo "ERROR: required host compiler not executable: ${exe}" >&2
        echo "       Override with ANKROS_AMBER_CC, ANKROS_AMBER_CXX, or ANKROS_AMBER_FC." >&2
        exit 2
    fi
done

if [[ ! -f "${CUDA_HOME}/include/cuda.h" ]]; then
    echo "ERROR: CUDA headers not found under ${CUDA_HOME}/include" >&2
    echo "       Install cuda-driver-dev and cuda-cudart-dev into envs/amber_cuda, then re-run." >&2
    exit 2
fi

if [[ ! -e "${CUDA_HOME}/lib/libcudart.so" && ! -e "${CUDA_HOME}/lib64/libcudart.so" ]]; then
    echo "ERROR: unversioned libcudart.so not found under ${CUDA_HOME}/lib or ${CUDA_HOME}/lib64" >&2
    echo "       Install cuda-cudart-dev into envs/amber_cuda, then re-run." >&2
    exit 2
fi

if [[ ! -f "${CUDA_HOME}/include/fatbinary_section.h" ]]; then
    echo "ERROR: fatbinary_section.h not found under ${CUDA_HOME}/include" >&2
    echo "       Install cuda-nvcc-dev_linux-64 into envs/amber_cuda, then re-run." >&2
    exit 2
fi

# Conda-forge CUDA places nvvm under the env prefix, while nvcc.profile in
# cuda-nvcc-tools looks for it under targets/x86_64-linux.
if [[ -x "${ENV_PREFIX}/nvvm/bin/cicc" && ! -e "${ENV_PREFIX}/targets/x86_64-linux/nvvm" ]]; then
    ln -s "../../nvvm" "${ENV_PREFIX}/targets/x86_64-linux/nvvm"
fi

echo
echo "== tool versions =="
cmake --version | head -1
make --version | head -1
"${HOST_CC}" --version | head -1
"${HOST_CXX}" --version | head -1
"${HOST_FC}" --version | head -1
echo "host CC : ${HOST_CC}"
echo "host CXX: ${HOST_CXX}"
echo "host FC : ${HOST_FC}"
nvcc --version | tail -1
if command -v nvidia-smi >/dev/null; then
    nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total \
               --format=csv,noheader || true
fi

cuda_version="$(nvcc --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' | head -1)"
cuda_major="${cuda_version%%.*}"
cuda_minor="${cuda_version#*.}"
cuda_minor="${cuda_minor%%.*}"
if [[ -z "${cuda_version}" ]]; then
    echo "ERROR: could not parse CUDA version from nvcc --version" >&2
    exit 2
fi

if (( cuda_major > 12 || (cuda_major == 12 && cuda_minor >= 9) )); then
    echo "ERROR: Amber/pmemd refuses CUDA ${cuda_version}; it requires CUDA >= 7.5 and < 12.9." >&2
    echo "       Use a CUDA 12.8 or older toolkit for this build, or set CUDA_HOME to one before running." >&2
    exit 2
fi

SYSROOT_LIB_DIR="${ENV_PREFIX}/x86_64-conda-linux-gnu/sysroot/lib64"
CMATH_LIB=""
if [[ -e "${SYSROOT_LIB_DIR}/libm.so" ]]; then
    CMATH_LIB="${SYSROOT_LIB_DIR}/libm.so"
elif [[ -e "${SYSROOT_LIB_DIR}/libm.a" ]]; then
    CMATH_LIB="${SYSROOT_LIB_DIR}/libm.a"
fi
CMATH_CMAKE_ARG=()
if [[ -n "${CMATH_LIB}" ]]; then
    CMATH_CMAKE_ARG=(-DCMath_LIBRARIES="${CMATH_LIB}")
fi

GPU_CMAKE_ARG=()
GPU_NVCC_FLAGS=""
if [[ -n "${GPU_TARGETS}" ]]; then
    for sm in ${GPU_TARGETS}; do
        sm="${sm#sm_}"
        if [[ ! "${sm}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: invalid GPU target '${sm}'. Use values like 'sm_61' or '61'." >&2
            exit 2
        fi
        GPU_NVCC_FLAGS="${GPU_NVCC_FLAGS:+${GPU_NVCC_FLAGS};}-gencode;arch=compute_${sm},code=sm_${sm}"
    done
    GPU_NVCC_FLAGS="${GPU_NVCC_FLAGS};-Wno-deprecated-gpu-targets"
    if [[ "${ONLY_GPU_TARGETS}" == "0" ]]; then
        GPU_CMAKE_ARG=(-DCUDA_NVCC_FLAGS="${GPU_NVCC_FLAGS}")
    else
        GPU_CMAKE_ARG=(-DANKROS_CUDA_NVCC_FLAGS:STRING="${GPU_NVCC_FLAGS}")
    fi
fi

echo
echo "== source extraction =="
mkdir -p "${WORK_DIR}"
if [[ -d "${SRC_DIR}" ]]; then
    echo "    reusing existing source tree: ${SRC_DIR}"
else
    echo "    extracting $(basename "${PACKAGE}") into ${WORK_DIR}"
    tar -xjf "${PACKAGE}" -C "${WORK_DIR}"
fi

if [[ ! -f "${SRC_DIR}/CMakeLists.txt" ]]; then
    echo "ERROR: expected ${SRC_DIR}/CMakeLists.txt after extraction" >&2
    exit 2
fi

if [[ -n "${GPU_NVCC_FLAGS}" && "${ONLY_GPU_TARGETS}" != "0" ]]; then
    CUDA_CONFIG="${SRC_DIR}/cmake/CudaConfig.cmake"
    if [[ ! -f "${CUDA_CONFIG}" ]]; then
        echo "ERROR: expected ${CUDA_CONFIG} for CUDA GPU target patching" >&2
        exit 2
    fi
    if grep -q "ANKROS_CUDA_NVCC_FLAGS" "${CUDA_CONFIG}"; then
        echo "    CUDA target override patch already present"
    else
        echo "    patching Amber CUDA target selection to use only: ${GPU_TARGETS}"
        tmp_cuda_config="$(mktemp)"
        awk '
            !patched && $0 ~ /^[[:space:]]*if\(\$\{CUDA_VERSION\} VERSION_EQUAL 7\.5\)/ {
                print "\t\tif(DEFINED ANKROS_CUDA_NVCC_FLAGS AND NOT \"${ANKROS_CUDA_NVCC_FLAGS}\" STREQUAL \"\")"
                print "\t\t\tmessage(STATUS \"Configuring CUDA for ANKros GPU target override: ${ANKROS_CUDA_NVCC_FLAGS}\")"
                print "\t\t\tlist(APPEND CUDA_NVCC_FLAGS ${ANKROS_CUDA_NVCC_FLAGS})"
                sub(/^[[:space:]]*if/, "\t\telseif")
                print
                patched=1
                next
            }
            { print }
            END { if (!patched) exit 42 }
        ' "${CUDA_CONFIG}" > "${tmp_cuda_config}" || {
            status="$?"
            rm -f "${tmp_cuda_config}"
            echo "ERROR: failed to patch ${CUDA_CONFIG} (awk exit ${status})" >&2
            exit 2
        }
        mv "${tmp_cuda_config}" "${CUDA_CONFIG}"
    fi
fi

echo
echo "== configure =="
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"
rm -rf CMakeCache.txt CMakeFiles

cmake "${SRC_DIR}" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_HOME}" \
    -DCUDA_NVCC_EXECUTABLE="${ENV_PREFIX}/bin/nvcc" \
    -DCMAKE_LIBRARY_PATH="${CUDA_HOME}/lib;${CUDA_HOME}/lib64;${ENV_PREFIX}/lib;${SYSROOT_LIB_DIR}" \
    "${CMATH_CMAKE_ARG[@]}" \
    "${GPU_CMAKE_ARG[@]}" \
    -DCOMPILER=MANUAL \
    -DCMAKE_C_COMPILER="${HOST_CC}" \
    -DCMAKE_CXX_COMPILER="${HOST_CXX}" \
    -DCMAKE_Fortran_COMPILER="${HOST_FC}" \
    -DCUDA_HOST_COMPILER="${HOST_CXX}" \
    -DMPI=FALSE \
    -DCUDA=TRUE \
    -DINSTALL_TESTS=FALSE \
    -DDOWNLOAD_MINICONDA=FALSE \
    -DBUILD_PYTHON=FALSE \
    -DBUILD_PERL=FALSE \
    -DBUILD_GUI=FALSE \
    -DPMEMD_ONLY=TRUE \
    -DCHECK_UPDATES=FALSE \
    ${ANKROS_AMBER_CMAKE_ARGS:-}

echo
echo "== build + install =="
make -j "${JOBS}" install

echo
echo "== env wrappers =="
mkdir -p "${ENV_PREFIX}/bin"
wrapper_count=0

write_wrapper() {
    local exe="$1"
    local target="${ENV_PREFIX}/bin/${exe}"

    if [[ -e "${target}" ]] && ! grep -q 'Generated by ANKros Amber CUDA installer' "${target}" 2>/dev/null; then
        echo "ERROR: refusing to overwrite existing non-generated file: ${target}" >&2
        echo "       Move it aside or set ANKROS_AMBER_OVERWRITE_WRAPPERS=1 if this is intentional." >&2
        if [[ "${ANKROS_AMBER_OVERWRITE_WRAPPERS:-0}" != "1" ]]; then
            exit 1
        fi
    fi

    cat > "${target}" <<EOF
#!/usr/bin/env bash
# Generated by ANKros Amber CUDA installer.
export AMBERHOME="${INSTALL_PREFIX}"
export LD_LIBRARY_PATH="\${AMBERHOME}/lib:\${AMBERHOME}/lib64:${ENV_PREFIX}/lib:${ENV_PREFIX}/lib64:\${LD_LIBRARY_PATH:-}"
exec "\${AMBERHOME}/bin/${exe}" "\$@"
EOF
    chmod +x "${target}"
    echo "    ${target} -> ${INSTALL_PREFIX}/bin/${exe}"
    wrapper_count=$((wrapper_count + 1))
}

for exe in pmemd.cuda pmemd.cuda_SPFP pmemd.cuda_DPFP pmemd.cuda.MPI; do
    if [[ -x "${INSTALL_PREFIX}/bin/${exe}" ]]; then
        write_wrapper "${exe}"
    fi
done

if [[ "${wrapper_count}" -eq 0 ]]; then
    echo "ERROR: no pmemd executables found under ${INSTALL_PREFIX}/bin" >&2
    exit 1
fi

echo
echo "== smoke check =="
if [[ -x "${ENV_PREFIX}/bin/pmemd.cuda" ]]; then
    "${ENV_PREFIX}/bin/pmemd.cuda" -h >/dev/null 2>&1 || true
    echo "    pmemd.cuda wrapper exists"
else
    echo "    pmemd.cuda not installed; available wrappers:"
    find "${ENV_PREFIX}/bin" -maxdepth 1 -type f -name 'pmemd*' -printf '    %f\n' | sort
fi

echo
echo "==[$(date -Is)] install_amber_cuda.sh done"
