#!/usr/bin/env bash
# Copyright (C) 2016 The ESPResSo project
# Copyright (C) 2014 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

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
function cmd {
    echo ">$1"
    eval $1
}

# handle environment variables
[ -z "$insource" ] && insource="false"
[ -z "$srcdir" ] && srcdir=`pwd`
[ -z "$cmake_params" ] && cmake_params=""
[ -z "$with_fftw" ] && with_fftw="true"
[ -z "$with_python_interface" ] && with_python_interface="true"
[ -z "$with_coverage" ] && with_coverage="false"
[ -z "$with_static_analysis" ] && with_static_analysis="false"
[ -z "$myconfig" ] && myconfig="default"
[ -z "$check_procs" ] && check_procs=2
[ -z "$build_procs" ] && build_procs=2
[ -z "$make_check" ] && make_check="true"

cmake_params="-DWARNINGS_ARE_ERRORS=ON -DTEST_NP:INT=$check_procs $cmake_params"

if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

outp insource srcdir builddir \
    cmake_params with_fftw \
    with_python_interface with_coverage with_static_analysis myconfig check_procs build_procs

# check indentation of python files
pep8 --filename=*.pyx,*.pxd,*.py --select=E111 $srcdir/src/python/espressomd/
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
    grep -q rhel /etc/os-release && source /etc/profile.d/modules.sh && module load mpi
fi

# CONFIGURE
start "CONFIGURE"

if [ $with_fftw = "true" ]; then
    cmake_params="$cmake_params"
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

if [ $with_static_analysis = "true" ]; then
    cmake_params="-DWITH_CLANG_TIDY=ON $cmake_params"
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

cmd "cmake $cmake_params $srcdir" || exit 1
end "CONFIGURE"

# BUILD
start "BUILD"

cmd "make -k -j${build_procs}" || cmd "make -k -j1" || exit $?

end "BUILD"

if $make_check; then
    start "TEST"

    if [ -z "$run_tests" ]; then
        cmd "make -j${build_procs} check_python $make_params" || exit 1
    else
        cmd "make python_tests $make_params"
        for t in $run_tests; do
            cmd "ctest --output-on-failure -R $t" || exit 1
        done
    fi
    cmd "make -j${build_procs} check_unit_tests $make_params" || exit 1

    end "TEST"
else
    start "TEST"

    cmd "mpiexec -n $check_procs ./pypresso $srcdir/testsuite/particle.py" || exit 1

    end "TEST"
fi

if $with_coverage; then
    cd $builddir
    lcov --directory . --capture --output-file coverage.info # capture coverage info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
    lcov --remove coverage.info '*/doc/*' --output-file coverage.info # filter out docs
    lcov --remove coverage.info '*/unit_tests/*' --output-file coverage.info # filter out unit test
    lcov --list coverage.info #debug info
    # Uploading report to CodeCov
    if [ -z "$CODECOV_TOKEN" ]; then
        bash <(curl -s https://codecov.io/bash) || echo "Codecov did not collect coverage reports"
    else
        bash <(curl -s https://codecov.io/bash) -t "$CODECOV_TOKEN" || echo "Codecov did not collect coverage reports"
    fi
fi
