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
[ -z "$srcdir" ] && srcdir=`pwd`
[ -z "$cmake_params" ] && configure_params=""
[ -z "$with_fftw" ] && with_fftw="true"
[ -z "$with_tcl" ] && with_tcl="true"
[ -z "$with_python_interface" ] && with_python_interface="true"
[ -z "$myconfig" ] && myconfig="default"
[ -z "$check_procs" ] && check_procs=2
[ -z "$make_check" ] && make_check="true"
[ -z "$make_params" ] && make_params=""

builddir=$srcdir/build

outp insource srcdir builddir \
    cmake_params with_fftw \
    with_tcl with_python_interface myconfig check_procs make_check make_params

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

# Make build dir exist.
if [ ! -d $builddir ]; then
    echo "Creating $builddir..."
    mkdir -p $builddir
fi

cd $builddir

start "CMAKE"

if $with_fftw; then
    cmake_params="$cmake_params"
else
    cmake_params="-DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON $cmake_params"
fi

if $with_tcl; then
    cmake_params="-DWITH_TCL=ON $cmake_params"
else
    cmake_params="-DWITH_TCL=OFF $cmake_params"
fi

if $with_python_interface; then
    cmake_params="-DWITH_PYTHON=ON $cmake_params"
else
    cmake_params="-DWITH_PYTHON=OFF $cmake_params"
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

cmd "cmake $cmake_params $srcdir" || exit $?
end "CMAKE"

# BUILD
start "BUILD"

cmd "make $make_params" || exit $?

end "BUILD"

if $make_check; then
    start "TEST"

    make -k check

    ec=$?

    echo "make check returned $ec"

# Did any test fail?
    if [ $ec -ne 0 ]; then
        echo "Checking for test logs."
        TCL_LOG="$builddir/testsuite/tcl/Testing/Temporary/LastTest.log"
        PYTHON_LOG="$builddir/testsuite/python/Testing/Temporary/LastTest.log"
        UNIT_LOG="$builddir/testsuite/src/core/unit_tests/Testing/Temporary/LastTest.log"

        for LOG in $TCL_LOG $PYTHON_LOG $UNIT_LOG; do
                echo $LOG         
		if [ -f $LOG ]; then
			cat $LOG;
		fi
	done         
        exit $ec
    fi

    end "TEST"
fi
