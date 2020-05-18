#!/usr/bin/env bash
# Copyright (C) 2016-2019 The ESPResSo project
# Copyright (C) 2014 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

abort() {
    echo "An error occurred. Exiting..." >&2
    echo "Command that failed: ${BASH_COMMAND}" >&2
    exit 1
}

trap 'abort' 0
set -e

# HELPER FUNCTIONS

# output value of env variables
outp() {
    for p in ${*}; do
        echo "  ${p}=${!p}"
    done
}

# start a block
start() {
    echo "=================================================="
    echo "START ${1}"
    echo "=================================================="
}

# end a block
end() {
    echo "=================================================="
    echo "END ${1}"
    echo "=================================================="
}

# set a default value to empty environment variables
# cast boolean values to true/false
set_default_value() {
    if [ "${#}" != 2 ]; then
        echo "set_default_value() takes 2 arguments (varname, default), got ${#}"
        exit 1
    fi
    local -r varname="${1}"
    local -r default="${2}"
    local -r varname_alphabet=$(echo "${varname}" | tr -d '[:alnum:]_')
    if [ ! -z "${varname_alphabet}" ]; then
        echo "variable name '${varname}' contains unauthorized symbols"
        exit 1
    fi
    local -r value="${!varname}"
    if [ "${default}" = true ] || [ "${default}" = false ]; then
        # cast boolean values to true/false
        local -r val=$(echo "${value}" | tr '[:upper:]' '[:lower:]')
        if [ "${val}" = false ] || [ "${val}" = "off" ] || [ "${val}" = 0 ] || [ "${val}" = "no" ]; then
            eval "${varname}=false"
        elif [ "${val}" = true ] || [ "${val}" = "on" ] || [ "${val}" = 1 ] || [ "${val}" = "yes" ]; then
            eval "${varname}=true"
        elif [ -z "${val}" ]; then
            eval "${varname}='${default}'"
        else
            echo "Cannot interpret '${value}' as a true/false value in variable '${varname}'"
            exit 1
        fi
    elif [ -z "${value}" ]; then
        eval "${varname}='${default}'"
    fi
}


# handle environment variables
set_default_value srcdir "$(pwd)"
set_default_value cmake_params ""
set_default_value with_fftw true
set_default_value with_coverage false
set_default_value with_ubsan false
set_default_value with_asan false
set_default_value with_static_analysis false
set_default_value myconfig "default"
set_default_value build_procs 2
set_default_value check_procs ${build_procs}
set_default_value check_odd_only false
set_default_value check_gpu_only false
set_default_value check_skip_long false
set_default_value make_check_unit_tests true
set_default_value make_check_python true
set_default_value make_check_tutorials false
set_default_value make_check_samples false
set_default_value make_check_benchmarks false
set_default_value with_cuda false
set_default_value with_cuda_compiler "nvcc"
set_default_value build_type "RelWithAssert"
set_default_value with_ccache false
set_default_value with_scafacos true
set_default_value test_timeout 300
set_default_value hide_gpu false

if [ "${make_check_unit_tests}" = true ] || [ "${make_check_python}" = true ] || [ "${make_check_tutorials}" = true ] || [ "${make_check_samples}" = true ] || [ "${make_check_benchmarks}" = true ]; then
    run_checks=true
else
    run_checks=false
fi

if [ "${with_coverage}" = true ]; then
    build_type="Coverage"
    bash <(curl -s https://codecov.io/env) 1>/dev/null 2>&1
fi

cmake_params="-DCMAKE_BUILD_TYPE=${build_type} -DWARNINGS_ARE_ERRORS=ON -DTEST_NP:INT=${check_procs} ${cmake_params}"
cmake_params="${cmake_params} -DCMAKE_INSTALL_PREFIX=/tmp/espresso-unit-tests"
cmake_params="${cmake_params} -DTEST_TIMEOUT=${test_timeout}"

if [ "${with_ccache}" = true ]; then
    cmake_params="${cmake_params} -DWITH_CCACHE=ON"
fi
if [ "${with_scafacos}" = true ]; then
    cmake_params="${cmake_params} -DWITH_SCAFACOS=ON"
fi

if [ "${with_fftw}" = true ]; then
    :
else
    cmake_params="-DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON ${cmake_params}"
fi

if [ "${with_coverage}" = true ]; then
    cmake_params="-DWITH_COVERAGE=ON ${cmake_params}"
fi

if [ "${with_asan}" = true ]; then
    cmake_params="-DWITH_ASAN=ON ${cmake_params}"
fi

if [ "${with_ubsan}" = true ]; then
    cmake_params="-DWITH_UBSAN=ON ${cmake_params}"
fi

if [ "${with_static_analysis}" = true ]; then
    cmake_params="-DWITH_CLANG_TIDY=ON ${cmake_params}"
fi

if [ "${with_cuda}" = true ]; then
    cmake_params="-DWITH_CUDA=ON -DWITH_CUDA_COMPILER=${with_cuda_compiler} ${cmake_params}"
else
    cmake_params="-DWITH_CUDA=OFF ${cmake_params}"
fi

command -v nvidia-smi && nvidia-smi || true
if [ "${hide_gpu}" = true ]; then
    echo "Hiding gpu from Cuda via CUDA_VISIBLE_DEVICES"
    export CUDA_VISIBLE_DEVICES=""
fi

builddir="${srcdir}/build"

outp srcdir builddir \
    make_check_unit_tests make_check_python make_check_tutorials make_check_samples make_check_benchmarks \
    cmake_params with_fftw \
    with_coverage \
    with_ubsan with_asan \
    check_odd_only \
    with_static_analysis myconfig \
    build_procs check_procs \
    with_cuda with_cuda_compiler with_ccache

echo "Creating ${builddir}..."
mkdir -p "${builddir}"
cd "${builddir}"

# load MPI module if necessary
if [ -f "/etc/os-release" ]; then
    grep -q suse /etc/os-release && . /etc/profile.d/modules.sh && module load gnu-openmpi
    grep -q 'rhel\|fedora' /etc/os-release && for f in /etc/profile.d/*module*.sh; do . "${f}"; done && module load mpi
fi

# CONFIGURE
start "CONFIGURE"

MYCONFIG_DIR="${srcdir}/maintainer/configs"
if [ "${myconfig}" = "default" ]; then
    echo "Using default myconfig."
else
    myconfig_file="${MYCONFIG_DIR}/${myconfig}.hpp"
    if [ ! -e "${myconfig_file}" ]; then
        echo "${myconfig_file} does not exist!"
        exit 1
    fi
    echo "Copying ${myconfig}.hpp to ${builddir}/myconfig.hpp..."
    cp "${myconfig_file}" "${builddir}/myconfig.hpp"
fi

cmake ${cmake_params} "${srcdir}" || exit 1
end "CONFIGURE"

# BUILD
start "BUILD"

make -k -j${build_procs} || make -k -j1 || exit ${?}

end "BUILD"

# check for exit function, which should never be called from shared library
# can't do this on CUDA though because nvcc creates a host function that just
# calls exit() for each device function, and can't do this with coverage
# because gcov 9.0 adds code that calls exit()
if [[ "${with_coverage}" == false && ( "${with_cuda}" == false || "${with_cuda_compiler}" != "nvcc" ) ]]; then
    if nm -o -C $(find . -name '*.so') | grep '[^a-z]exit@@GLIBC'; then
        echo "Found calls to exit() function in shared libraries."
        exit 1
    fi
fi

if [ "${run_checks}" = true ]; then
    start "TEST"

    # unit tests
    if [ "${make_check_unit_tests}" = true ]; then
        make -j${build_procs} check_unit_tests ${make_params} || exit 1
    fi

    # integration tests
    if [ "${make_check_python}" = true ]; then
        if [ -z "${run_tests}" ]; then
            if [ "${check_odd_only}" = true ]; then
                make -j${build_procs} check_python_parallel_odd ${make_params} || exit 1
            elif [ "${check_gpu_only}" = true ]; then
                make -j${build_procs} check_python_gpu ${make_params} || exit 1
            elif [ "${check_skip_long}" = true ]; then
                make -j${build_procs} check_python_skip_long ${make_params} || exit 1
            else
                make -j${build_procs} check_python ${make_params} || exit 1
            fi
        else
            make python_tests ${make_params}
            for t in ${run_tests}; do
                ctest --timeout 60 --output-on-failure -R "${t}" || exit 1
            done
        fi
    fi

    # tutorial tests
    if [ "${make_check_tutorials}" = true ]; then
        make -j${build_procs} check_tutorials ${make_params} || exit 1
    fi

    # sample tests
    if [ "${make_check_samples}" = true ]; then
        make -j${build_procs} check_samples ${make_params} || exit 1
    fi

    # benchmark tests
    if [ "${make_check_benchmarks}" = true ]; then
        make -j${build_procs} check_benchmarks ${make_params} || exit 1
    fi

    # installation tests
    make check_cmake_install ${make_params} || exit 1

    end "TEST"
else
    start "TEST"

    if [ "${HIP_PLATFORM}" != "hcc" ]; then
      check_proc_particle_test=${check_procs}
      if [ "${check_proc_particle_test}" -gt 4 ]; then
        check_proc_particle_test=4
      fi
      mpiexec -n ${check_proc_particle_test} ./pypresso "${srcdir}/testsuite/python/particle.py" || exit 1
    fi

    end "TEST"
fi

if [ "${with_coverage}" = true ]; then
    cd "${builddir}"
    lcov --gcov-tool "${GCOV:-gcov}" -q --directory . --ignore-errors graph --capture --output-file coverage.info # capture coverage info
    lcov --gcov-tool "${GCOV:-gcov}" -q --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
    lcov --gcov-tool "${GCOV:-gcov}" -q --remove coverage.info '*/doc/*' --output-file coverage.info # filter out docs
    # Uploading report to Codecov
    if [ -z "${CODECOV_TOKEN}" ]; then
        bash <(curl -s https://codecov.io/bash) -X gcov || echo "Codecov did not collect coverage reports"
    else
        bash <(curl -s https://codecov.io/bash) -t "${CODECOV_TOKEN}" -X gcov || echo "Codecov did not collect coverage reports"
    fi
fi

trap : 0
