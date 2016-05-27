#!/bin/bash --login 
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

##############################################################################################################
# configure.sh
##############################################################################################################

# in case the mpi path is not set properly run sth like:
#export PATH=/usr/lib64/mpi/gcc/openmpi/bin:$PATH

DIR=`dirname $0`
source $DIR/common.sh

# DEFAULTS
[ -z "$insource" ] && insource="true"
[ -z "$srcdir" ] && srcdir=`pwd`                                           
[ -z "$cmake_params" ] && configure_params=""  
[ -z "$myconfig" ] && myconfig="default"
# compiling on on proc is slow, lets use 4:
[ -z "$build_procs" ] && build_procs=4
# check_procs is never used...
[ -z "$check_procs" ] && check_procs=2
[ -z "$make_check" ] && make_check="true"
[ -z "$with_fftw" ] && with_fftw="true"
[ -z "$with_tcl" ] && with_tcl="true"
[ -z "$with_python_interface" ] && with_python_interface="false" #true in other build_cmake.sh
[ -z "$with_cuda" ] && with_cuda="true"
[ -z "$with_h5md" ] && with_h5md="true"

# set builddir 
if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

# create builddir if neccessary
if ! $insource; then
    if [ ! -d $builddir ]; then
        echo "Creating $builddir..."
        mkdir -p $builddir
    fi
fi

# alternatively could probably use: 'pushd $builddir' & 'popd'?
# and finally change into $builddir
# pushd $builddir
if ! $insource ; then
    cd $builddir
fi

outp insource srcdir builddir \
    configure_params with_fftw \
    with_tcl with_python_interface myconfig check_procs with_h5md with_cuda

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

# CONFIGURE
start "CONFIGURE"

if $with_fftw; then
    cmake_params="$cmake_params"
# THESE WHERE IN THE GCC SCRIPT
#     if [ -e $FFTW_HEADER ]; then
#         echo "Using FFTW => deleting mock $FFTW_HEADER..."
#         rm $FFTW_HEADER
#     fi
else
    cmake_params="-DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON $cmake_params"
#     echo "Not using FFTW => generating mock $FFTW_HEADER..."
#     echo "#error ERROR: fftw is not really present but used somewhere." \
#         > $FFTW_HEADER
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

if $with_cuda; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
    cmake_params="-DWITH_CUDA=ON $cmake_params"
    # if [ -e $CUDA_HEADER ]; then
    #     echo "Using CUDA => deleting mock $CUDA_HEADER..."
    #     rm $CUDA_HEADER
    # fi
else
    cmake_params="-DWITH_CUDA=OFF $cmake_params"
    # configure_vars="$configure_vars NVCC=/bin/false"
    # echo "Not using CUDA => generating mock $CUDA_HEADER..."
    # echo "#error ERROR: cuda is not really present but used somewhere." \
    #     > $CUDA_HEADER
fi

# The CMakeLists do not seem to have this option yet:
# if $with_h5md; then
#     cmake_params="-DWITH_HDF5=ON $cmake_params"
# else
#     cmake_params="-DWITH_HDF5=OFF $cmake_params"
# fi

MYCONFIG_DIR=$srcdir/maintainer/jenkins/configs
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

# Activate anaconda environment
cmd "source activate test"

cmd "cmake $cmake_params $srcdir" || exit $?
end "CONFIGURE"

# BUILD
start "BUILD"

cmd "make -j $build_procs" || exit $?

end "BUILD"

if $make_check; then
    start "TEST"

    cmd "make check_tcl $make_params"
    ec=$?
    if [ $ec != 0 ]; then   
        cmd "cat $srcdir/testsuite/tcl/Testing/Temporary/LastTest.log"
        exit $ec
    fi

    cmd "make check_unit_tests $make_params"
    ec=$?
    if [ $ec != 0 ]; then   
        cmd "cat $srcdir/src/core/unit_tests/Testing/Temporary/LastTest.log"
        exit $ec
    fi

    end "TEST"
fi

# popd
