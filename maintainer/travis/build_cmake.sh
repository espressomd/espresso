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
[ -z "$insource" ] && insource="true"
[ -z "$srcdir" ] && srcdir=`pwd`
[ -z "$cmake_params" ] && cmake_params=""
[ -z "$with_fftw" ] && with_fftw="true"
[ -z "$with_tcl" ] && with_tcl="true"
[ -z "$with_python_interface" ] && with_python_interface="true"
[ -z "$myconfig" ] && myconfig="default"
[ -z "$check_procs" ] && check_procs=2
[ -z "$make_check" ] && make_check="true"

if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

outp insource srcdir builddir \
    cmake_params with_fftw \
    with_tcl with_python_interface myconfig check_procs

# check indentation of python files
#pep8 --filename=*.pyx,*.pxd,*.py --select=E111 $srcdir/src/python/espressomd/
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

if [ ! $insource ]; then
    if [ ! -d $builddir ]; then
        echo "Creating $builddir..."
        mkdir -p $builddir
    fi
fi

if [ ! $insource ]; then
    cd $builddir
fi

# CONFIGURE
start "CONFIGURE"

if [ $with_fftw = "true" ]; then
    cmake_params="$cmake_params"
else
    cmake_params="-DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON $cmake_params"
fi

if [ $with_tcl = "true" ]; then
    cmake_params="-DWITH_TCL=ON $cmake_params"
else
    cmake_params="-DWITH_TCL=OFF $cmake_params"
fi

if [ $with_python_interface = "true" ]; then
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
end "CONFIGURE"

# BUILD
start "BUILD"

cmd "make" || exit $?

end "BUILD"

if $make_check; then
    start "TEST"

    cmd "make check_python $make_params"
    ec=$?
    if [ $ec != 0 ]; then	
        cmd "cat $srcdir/testsuite/python/Testing/Temporary/LastTest.log"
        exit $ec
    fi

#    cmd "make check_tcl $make_params"
#    ec=$?
#    if [ $ec != 0 ]; then	
#        cmd "cat $srcdir/testsuite/tcl/Testing/Temporary/LastTest.log"
#        exit $ec
#    fi

    cmd "make check_unit_tests $make_params"
    ec=$?
    if [ $ec != 0 ]; then	
        cmd "cat $srcdir/src/core/unit_tests/Testing/Temporary/LastTest.log"
        exit $ec
    fi

    end "TEST"
fi
