#!/usr/bin/env bash
#
# Copyright (C) 2016-2022 The ESPResSo project
# Copyright (C) 2014 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
#

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

# the number of available processors depends on the CI runner
ci_procs=2
if [ "${GITLAB_CI}" = "true" ]; then
    if [[ "${OSTYPE}" == "linux-gnu"* ]]; then
        # Linux runner
        if grep -q "i7-3820" /proc/cpuinfo; then
            # communication bottleneck for more than 2 cores on Intel i7-3820
            ci_procs=2
        else
            ci_procs=4
        fi
    elif [[ "${OSTYPE}" == "darwin"* ]]; then
        # macOS runner
        ci_procs=2
    fi
elif [ "${GITHUB_ACTIONS}" = "true" ]; then
    # GitHub Actions only provide 1 core; request 2 cores to run tests
    # in parallel (OpenMPI allows oversubscription)
    ci_procs=2
else
    ci_procs=$(nproc)
fi

# handle environment variables
set_default_value srcdir "$(pwd)"
set_default_value cmake_params ""
set_default_value with_coverage false
set_default_value with_coverage_python ${with_coverage}
set_default_value with_ubsan false
set_default_value with_asan false
set_default_value with_static_analysis false
set_default_value myconfig "default"
set_default_value build_procs ${ci_procs}
set_default_value check_procs ${build_procs}
set_default_value check_odd_only false
set_default_value check_gpu_only false
set_default_value check_skip_long false
set_default_value make_check_unit_tests true
set_default_value make_check_python true
set_default_value make_check_tutorials false
set_default_value make_check_samples false
set_default_value make_check_benchmarks false
set_default_value with_fast_math false
set_default_value with_cuda false
set_default_value with_cuda_compiler "nvcc"
set_default_value with_cxx_standard 17
set_default_value build_type "RelWithAssert"
set_default_value with_ccache false
set_default_value with_hdf5 true
set_default_value with_fftw true
set_default_value with_gsl true
set_default_value with_scafacos false
set_default_value with_walberla false
set_default_value with_stokesian_dynamics false
set_default_value test_timeout 300
set_default_value hide_gpu false
set_default_value mpiexec_preflags ""

if [ "${make_check_unit_tests}" = true ] || [ "${make_check_python}" = true ] || [ "${make_check_tutorials}" = true ] || [ "${make_check_samples}" = true ] || [ "${make_check_benchmarks}" = true ]; then
    run_checks=true
else
    run_checks=false
fi

if [ "${with_coverage}" = true ]; then
    build_type="Coverage"
fi

if [ "${with_fast_math}" = true ]; then
    cmake_param_protected="-DCMAKE_CXX_FLAGS=-ffast-math -fno-finite-math-only"
fi

cmake_params="-D CMAKE_BUILD_TYPE=${build_type} -D CMAKE_CXX_STANDARD=${with_cxx_standard} -D ESPRESSO_WARNINGS_ARE_ERRORS=ON ${cmake_params}"
cmake_params="${cmake_params} -D CMAKE_INSTALL_PREFIX=/tmp/espresso-unit-tests -D ESPRESSO_INSIDE_DOCKER=ON"
cmake_params="${cmake_params} -D ESPRESSO_CTEST_ARGS=-j${check_procs} -D ESPRESSO_TEST_TIMEOUT=${test_timeout}"

if [ "${make_check_benchmarks}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_BENCHMARKS=ON"
fi

if [ "${with_ccache}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_CCACHE=ON"
fi

if [ "${with_hdf5}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_HDF5=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_HDF5=OFF"
fi

if [ "${with_fftw}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_FFTW=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_FFTW=OFF"
fi

if [ "${with_gsl}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_GSL=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_GSL=OFF"
fi

if [ "${with_scafacos}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_SCAFACOS=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_SCAFACOS=OFF"
fi

if [ "${with_stokesian_dynamics}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_STOKESIAN_DYNAMICS=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_STOKESIAN_DYNAMICS=OFF"
fi

if [ "${with_walberla}" = true ]; then
  cmake_params="${cmake_params} -D ESPRESSO_BUILD_WITH_WALBERLA_LIB=ON"
  # disable default OpenMPI CPU binding mechanism to avoid stale references to
  # waLBerla objects when multiple LB python tests run in parallel on NUMA archs
  mpiexec_preflags="${mpiexec_preflags:+$mpiexec_preflags;}--bind-to;none"
fi

if [ "${with_coverage}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_COVERAGE=ON ${cmake_params}"
fi

if [ "${with_coverage_python}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_COVERAGE_PYTHON=ON ${cmake_params}"
fi

if [ "${with_asan}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_ASAN=ON ${cmake_params}"
fi

if [ "${with_ubsan}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_UBSAN=ON ${cmake_params}"
fi

if [ "${with_static_analysis}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_CLANG_TIDY=ON ${cmake_params}"
fi

if [ "${run_checks}" = true ]; then
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_TESTS=ON"
else
    cmake_params="${cmake_params} -D ESPRESSO_BUILD_TESTS=OFF"
fi

if [ "${with_cuda}" = true ]; then
    cmake_params="-D ESPRESSO_BUILD_WITH_CUDA=ON -D ESPRESSO_CUDA_COMPILER=${with_cuda_compiler} ${cmake_params}"
else
    cmake_params="-D ESPRESSO_BUILD_WITH_CUDA=OFF ${cmake_params}"
fi

command -v nvidia-smi && nvidia-smi || true
nvidia-smi -L || true
if [ "${hide_gpu}" = true ]; then
    echo "Hiding gpu from Cuda via CUDA_VISIBLE_DEVICES"
    export CUDA_VISIBLE_DEVICES=""
fi

builddir="${srcdir}/build"

outp srcdir builddir \
    make_check_unit_tests make_check_python make_check_tutorials make_check_samples make_check_benchmarks \
    cmake_params \
    with_coverage with_coverage_python \
    with_ubsan with_asan \
    check_odd_only \
    with_static_analysis with_fast_math myconfig \
    build_procs check_procs with_walberla \
    with_cuda with_cuda_compiler with_ccache

echo "Creating ${builddir}..."
mkdir -p "${builddir}"
cd "${builddir}"

# load MPI module if necessary
if [ -f "/etc/os-release" ]; then
    grep -q suse /etc/os-release && . /etc/profile.d/modules.sh && module load gnu-openmpi
    grep -q 'rhel\|fedora' /etc/os-release && for f in /etc/profile.d/*module*.sh; do . "${f}"; done && module load mpi
    grep -q "Ubuntu 22.04" /etc/os-release && export MPIEXEC_PREFLAGS="--mca;btl_vader_single_copy_mechanism;none${mpiexec_preflags:+;$mpiexec_preflags}"
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
    if [ "${with_fast_math}" = true ]; then
        sed -i '/#define ADDITIONAL_CHECKS/d' "${builddir}/myconfig.hpp"
    fi
fi

if [ -z "${cmake_param_protected}" ]; then
  cmake "${srcdir}" ${cmake_params} || exit 1
else
  cmake "${srcdir}" ${cmake_params} "${cmake_param_protected}" || exit 1
fi
end "CONFIGURE"

# BUILD
start "BUILD"

make -k -j${build_procs} || make -k -j1 || exit ${?}

end "BUILD"

# Check for exit() function, which should never be called from a shared
# library. See details in https://github.com/espressomd/espresso/issues/2249
# Can't do this check on CUDA though because nvcc creates a host function
# that just calls exit() for each device function, and can't do this with
# coverage because gcov 9.0 adds code that calls exit(), and can't do this
# with walberla because the library calls exit() in assertions.
if [[ "${with_coverage}" == false && ( "${with_cuda}" == false || "${with_cuda_compiler}" != "nvcc" ) && "${with_walberla}" != "true" ]]; then
    if nm -o -C $(find . -name '*.so') | grep '[^a-z]exit@@GLIBC'; then
        echo "Found calls to exit() function in shared libraries."
        exit 1
    fi
fi

if [ "${run_checks}" = true ]; then
    start "TEST"

    # fail if built with CUDA but no compatible GPU was found
    if [ "${with_cuda}" = true ] && [ "${hide_gpu}" != true ]; then
        ./pypresso -c "import espressomd.cuda_init as gpu;gpu.CudaInitHandle().device = 0" || (command -v nvidia-smi && nvidia-smi || true ; exit 1)
    fi

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

    check_proc_particle_test=${check_procs}
    if [ "${check_proc_particle_test}" -gt 4 ]; then
      check_proc_particle_test=4
    fi
    mpiexec -n ${check_proc_particle_test} ./pypresso "${srcdir}/testsuite/python/particle.py" || exit 1

    end "TEST"
fi

if [ "${with_coverage}" = true ] || [ "${with_coverage_python}" = true ]; then
    start "COVERAGE"
    cd "${builddir}"

    # import codecov key
    gpg --import "${CODECOV_PUBLIC_KEY}"

    # download uploader and signatures
    curl -OSs https://uploader.codecov.io/latest/linux/codecov
    curl -OSs https://uploader.codecov.io/latest/linux/codecov.SHA256SUM
    curl -OSs https://uploader.codecov.io/latest/linux/codecov.SHA256SUM.sig

    # check uploader integrity, exit script in case of failure
    gpg --verify codecov.SHA256SUM.sig codecov.SHA256SUM
    shasum -a 256 -c codecov.SHA256SUM
    chmod +x codecov

    codecov_opts=""
    if [ "${with_coverage}" = true ]; then
        echo "Running lcov and gcov..."
        codecov_opts="${codecov_opts} --gcov"
        lcov --gcov-tool "${GCOV:-gcov}" -q --directory . --ignore-errors graph --capture --output-file coverage.info # capture coverage info
        lcov --gcov-tool "${GCOV:-gcov}" -q --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
        lcov --gcov-tool "${GCOV:-gcov}" -q --remove coverage.info '*/doc/*' --output-file coverage.info # filter out docs
    fi
    if [ "${with_coverage_python}" = true ]; then
        echo "Running python3-coverage..."
        python3 -m coverage combine testsuite/python testsuite/scripts/tutorials testsuite/scripts/samples testsuite/scripts/benchmarks
        python3 -m coverage xml
    fi
    echo "Uploading to Codecov..."
    for codecov_trial in 1 2 3; do
        codecov_errno="0"
        ./codecov ${codecov_opts} -t "${CODECOV_TOKEN}" --nonZero || codecov_errno="${?}"
        if [ "${codecov_errno}" = "0" ]; then
            break
        fi
        echo "Codecov did not upload coverage reports (return code: ${codecov_errno})" >&2
        echo "That was attempt ${codecov_trial}/3"
        echo ""
        if [ ! "${codecov_trial}" = "3" ]; then
            sleep 10s
        fi
    done

    end "COVERAGE"
fi

trap : 0
