#!/usr/bin/env bash
# Copyright (C) 2016-2018 The ESPResSo project
# Copyright (C) 2014 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

abort()
{
    echo "An error occurred. Exiting..." >&2
    echo "Command that failed: $BASH_COMMAND" >&2
    exit 1
}

trap 'abort' 0
set -e

# HELPER FUNCTIONS

# output value of env variables
function outp {
    for p in $*; do
        echo "  $p=${!p}"
    done
}

# start a block
function start {
    echo "=================================================="
    echo "START $1"
    echo "=================================================="
}

# end a block
function end {
    echo "=================================================="
    echo "END $1"
    echo "=================================================="
}

# execute and output a command
# handle environment variables
[ -z "$insource" ] && insource="false"
[ -z "$srcdir" ] && srcdir=`pwd`
[ -z "$cmake_params" ] && cmake_params=""
[ -z "$with_fftw" ] && with_fftw="true"
[ -z "$with_python_interface" ] && with_python_interface="true"
[ -z "$with_coverage" ] && with_coverage="false"
[ -z "$with_ubsan" ] && with_ubsan="false"
[ -z "$with_asan" ] && with_asan="false"
[ -z "$with_static_analysis" ] && with_static_analysis="false"
[ -z "$myconfig" ] && myconfig="default"
[ -z "$build_procs" ] && build_procs=2
[ -z "$check_procs" ] && check_procs=$build_procs
[ -z "$make_check" ] && make_check="true"
[ -z "$check_odd_only" ] && check_odd_only="false"
[ -z "$check_gpu_only" ] && check_gpu_only="false"
[ -z "$check_skip_long" ] && check_skip_long="false"
[ -z "$make_check_tutorials" ] && make_check_tutorials="false"
[ -z "$make_check_samples" ] && make_check_samples="false"
[ -z "$python_version" ] && python_version="2"
[ -z "$with_cuda" ] && with_cuda="true"
[ -z "$build_type" ] && build_type="Debug"
[ -z "$with_ccache" ] && with_ccache="false"
[ -z "$test_timeout" ] && test_timeout="300"
[ -z "$hide_gpu" ] && hide_gpu="false" 

if [ $make_check ] || [ $make_check_tutorials ] || [ $make_check_samples ]; then
  run_checks="true"
fi

# If there are no user-provided flags they
# are added according to with_coverage.
if [ -z "$cxx_flags" ]; then
    if $with_coverage; then
        cxx_flags="-Og"
        nvcc_flags="-G"
    else
        if $run_checks; then
            cxx_flags="-O3"
            nvcc_flags="-O3"
        else
            cxx_flags="-O0"
            nvcc_flags="-O0"
        fi
    fi
fi

if [ ! -z ${with_coverage+x} ]; then
  bash <(curl -s https://codecov.io/env) &> /dev/null;
fi

cmake_params="-DCMAKE_BUILD_TYPE=$build_type -DPYTHON_EXECUTABLE=$(which python$python_version) -DWARNINGS_ARE_ERRORS=ON -DTEST_NP:INT=$check_procs $cmake_params -DWITH_SCAFACOS=ON"
cmake_params="$cmake_params -DCMAKE_CXX_FLAGS=$cxx_flags -DCUDA_NVCC_FLAGS=$nvcc_flags"
cmake_params="$cmake_params -DCMAKE_INSTALL_PREFIX=/tmp/espresso-unit-tests"
cmake_params="$cmake_params -DTEST_TIMEOUT=$test_timeout"
if $with_ccache; then
  cmake_params="$cmake_params -DWITH_CCACHE=ON"
fi

command -v nvidia-smi && nvidia-smi
if [ $hide_gpu = "true" ]; then
  echo "Hiding gpu from Cuda via CUDA_VISIBLE_DEVICES"
  export CUDA_VISIBLE_DEVICES=
fi

if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

outp insource srcdir builddir \
    make_check make_check_tutorials make_check_samples \
    cmake_params with_fftw \
    with_python_interface with_coverage \
    with_ubsan with_asan \
    check_odd_only \
    with_static_analysis myconfig \
    build_procs check_procs \
    python_version with_cuda with_ccache

# check indentation of python files
pep8_command () {
    if hash pep8 2> /dev/null; then
        pep8 "$@"
    elif hash pycodestyle 2> /dev/null; then
        pycodestyle "$@"
    elif hash pycodestyle-2 2> /dev/null; then
        pycodestyle-2 "$@"
    elif hash pycodestyle-3 2> /dev/null; then
        pycodestyle-3 "$@"
    else
        echo "pep8 not found";
        exit 1
    fi
}

pep8_command --filename=*.pyx,*.pxd,*.py --select=E111 $srcdir/src/python/espressomd/
ec=$?
if [ $ec -eq 0 ]; then
    echo ""
    echo "Indentation in Python files correct..."
    echo ""
else
    echo ""
    echo "Error: Python files are not indented the right way. Please use 4 spaces per indentation level!"
    echo ""
    exit $ec
fi
# enforce style rules
pylint_command () {
    if hash pylint 2> /dev/null; then
        pylint "$@"
    elif hash pylint3 2> /dev/null; then
        pylint3 "$@"
    elif hash pylint-2 2> /dev/null; then
        pylint-2 "$@"
    elif hash pylint-3 2> /dev/null; then
        pylint-3 "$@"
    else
        echo "pylint not found";
        exit 1
    fi
}
if [ $(pylint_command --version | grep -o 'pylint.*[0-9]\.[0-9]\.[0-9]' | awk '{ print $2 }' | cut -d'.' -f2) -gt 6 ]; then
    score_option='--score=no'
else
    score_option=''
fi
pylint_command $score_option --reports=no --disable=all --enable=C1001 $(find . -name '*.py*') || { echo -e "\nOld-style classes found.\nPlease convert to new-style:\nclass C: => class C(object):\n" && exit 1; }

if ! $insource; then
    if [ ! -d $builddir ]; then
        echo "Creating $builddir..."
        mkdir -p $builddir
    fi
fi

if ! $insource; then
    cd $builddir
fi

# load MPI module if necessary
if [ -f "/etc/os-release" ]; then
    grep -q suse /etc/os-release && source /etc/profile.d/modules.sh && module load gnu-openmpi
    grep -q 'rhel\|fedora' /etc/os-release && for f in /etc/profile.d/*module*.sh; do source $f; done && module load mpi
fi

# CONFIGURE
start "CONFIGURE"

if [ $with_fftw = "true" ]; then
    :
else
    cmake_params="-DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON $cmake_params"
fi

if [ $with_python_interface = "true" ]; then
    cmake_params="-DWITH_PYTHON=ON $cmake_params"
else
    cmake_params="-DWITH_PYTHON=OFF $cmake_params"
fi

if [ $with_coverage = "true" ]; then
    cmake_params="-DWITH_COVERAGE=ON $cmake_params"
fi

if [ $with_asan = "true" ]; then
    cmake_params="-DWITH_ASAN=ON $cmake_params"
fi

if [ $with_ubsan = "true" ]; then
    cmake_params="-DWITH_UBSAN=ON $cmake_params"
fi

if [ $with_static_analysis = "true" ]; then
    cmake_params="-DWITH_CLANG_TIDY=ON $cmake_params"
fi

if [ $with_cuda = "true" ]; then
    :
else
    cmake_params="-DWITH_CUDA=OFF $cmake_params"
fi

MYCONFIG_DIR=$srcdir/maintainer/configs
if [ "$myconfig" = "default" ]; then
    echo "Using default myconfig."
else
    myconfig_file=$MYCONFIG_DIR/$myconfig.hpp
    if [ ! -e "$myconfig_file" ]; then
        echo "$myconfig_file does not exist!"
        exit 1
    fi
    echo "Copying $myconfig.hpp to $builddir/myconfig.hpp..."
    cp $myconfig_file $builddir/myconfig.hpp
fi

cmake $cmake_params $srcdir || exit 1
end "CONFIGURE"

# BUILD
start "BUILD"

make -k -j${build_procs} || make -k -j1 || exit $?

end "BUILD"

# check for exit function, which should never be called from shared library
# can't do this on CUDA though because nvcc creates a host function that just calls exit for each device funtion
if [ $with_cuda != "true" -o "$(echo $NVCC | grep -o clang)" = "clang" ]; then
    if nm -o -C $(find . -name *.so) | grep '[^a-z]exit@@GLIBC'; then
        echo "Found calls to exit() function in shared libraries."
        exit 1
    fi
fi

if $run_checks; then
    start "TEST"

    # integration and unit tests
    if $make_check; then
        if [ -z "$run_tests" ]; then
            if $check_odd_only; then
                make -j${build_procs} check_python_parallel_odd $make_params || exit 1
            elif $check_gpu_only; then
                make -j${build_procs} check_python_gpu $make_params || exit 1
            elif $check_skip_long; then
                make -j${build_procs} check_python_skip_long $make_params || exit 1
            else
                make -j${build_procs} check_python $make_params || exit 1
            fi
        else
            make python_tests $make_params
            for t in $run_tests; do
                ctest --timeout 60 --output-on-failure -R $t || exit 1
            done
        fi
        make -j${build_procs} check_unit_tests $make_params || exit 1
    fi

    # tutorial tests
    if $make_check_tutorials; then
        make -j${build_procs} check_tutorials $make_params || exit 1
    fi

    # sample tests
    if $make_check_samples; then
        make -j${build_procs} check_samples $make_params || exit 1
    fi

    # installation tests
    make check_cmake_install $make_params || exit 1

    end "TEST"
else
    start "TEST"

    if [ "$HIP_PLATFORM" != "hcc" ]; then
      mpiexec -n $check_procs ./pypresso $srcdir/testsuite/python/particle.py || exit 1
    fi

    end "TEST"
fi

if $with_coverage; then
    cd $builddir
    lcov -q --directory . --ignore-errors graph --capture --output-file coverage.info # capture coverage info
    lcov -q --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
    lcov -q --remove coverage.info '*/doc/*' --output-file coverage.info # filter out docs
    # Uploading report to CodeCov
    if [ -z "$CODECOV_TOKEN" ]; then
        bash <(curl -s https://codecov.io/bash) -X gcov || echo "Codecov did not collect coverage reports"
    else
        bash <(curl -s https://codecov.io/bash) -t "$CODECOV_TOKEN" -X gcov || echo "Codecov did not collect coverage reports"
    fi
fi

trap : 0
