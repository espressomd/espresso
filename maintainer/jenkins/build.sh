#!/bin/bash --login 
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
if ((`ps -o nice= $$` < 5)); then 
    renice -n 5 $$
fi

# DEFAULTS
[ ! -v srcdir ] && srcdir=`pwd`
[ ! -v builddir ] && builddir="$srcdir/build"
[ ! -v run_bootstrap ] && run_bootstrap="false"
[ ! -v insource ] && insource="false"
$insource && builddir=$srcdir
[ ! -v run_configure ] && run_configure="true"
[ ! -v configure_params ] && configure_params=""
[ ! -v configure_vars ] && configure_vars=""
[ ! -v with_cuda ] && with_cuda="true"
[ ! -v with_mpi ] && with_mpi="true"
[ ! -v with_fftw ] && with_fftw="true"
[ ! -v run_build ] && run_build="true"
[ ! -v myconfig ] && myconfig="default"
[ ! -v build_procs ] && build_procs=4
[ ! -v run_check ] && run_check="true"
! $with_mpi && check_procs="1"
[ ! -v check_procs ] && check_procs="4"
[ ! -v run_doc ] && run_doc="false"
[ ! -v run_dist ] && run_dist="false"
[ ! -v run_distcheck ] && run_distcheck="false"
[ ! -v distcheck_flags ] && distcheck_flags=""

function outp() {
    for p in $*; do
        echo "  $p=${!p}"
    done
}
echo "Parameters:"
outp srcdir builddir run_bootstrap insource \
    run_configure configure_params configure_vars with_cuda with_mpi with_fftw \
    run_build myconfig build_procs \
    run_check check_procs \
    run_doc run_dist run_distcheck

if [ ! -e $srcdir/configure.ac ]; then
    echo "Could not find configure.ac in $srcdir!"
    exit 1
fi

if ! $insource; then
    echo "Creating $builddir..."
    mkdir -p $builddir
fi

function start() {
  echo "START $1"
}
function end() {
    echo "END $1"
}

function cmd() {
    echo ">$1"
    eval $1
}

if $run_bootstrap; then
    start "BOOTSTRAP"
    pushd $srcdir
    ./bootstrap.sh
    popd
    end "BOOTSTRAP"
fi

# change into build dir
pushd $builddir

# CONFIGURE
if $run_configure; then
# set up configure params
    if $with_mpi; then
        configure_params="--with-mpi $configure_params"
    else
        configure_params="--without-mpi $configure_params"
    fi
    
    CUDA_HEADER=$srcdir/src/cuda.h
    if $with_cuda; then
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
        configure_params="--with-cuda $configure_params"
        configure_vars="$configure_vars NVCC=/bin/false"
        if [ -e $CUDA_HEADER ]; then
            echo "Deleting $CUDA_HEADER..."
            rm $CUDA_HEADER
        fi
    else
        configure_params="--without-cuda $configure_params"
        echo "Generating $CUDA_HEADER..."
        echo "#error ERROR: cuda is not really present but used somewhere." \
            > $CUDA_HEADER
    fi
    
    FFTW_HEADER=$srcdir/src/fftw3.h
    if $with_fftw; then
        configure_params="--with-fftw $configure_params"
        if [ -e $FFTW_HEADER ]; then
            echo "Deleting $FFTW_HEADER..."
            rm $FFTW_HEADER
        fi
    else
        configure_params="--without-fftw $configure_params"
        echo "Generating $FFTW_HEADER..."
        echo "#error ERROR: fftw is not really present but used somewhere." \
            > $FFTW_HEADER
    fi
    
    configure_vars="$configure_vars CPU_COUNT=\"$check_procs\""
    
    start "CONFIGURE"
    cmd "$srcdir/configure $configure_params"
    end "CONFIGURE"
fi

# BUILD
if $run_build; then
    start "BUILD"
    if [ "$myconfig" = "default" ]; then
        echo "Using default myconfig."
        if [ -e $builddir/myconfig.h ]; then
            echo "Deleting $builddir/myconfig.h..."
            rm $builddir/myconfig.h
        fi
    else
        myconfig_file=$srcdir/maintainer/jenkins/configs/$myconfig.h
        if [ ! -e "$myconfig_file" ]; then
            echo "$myconfig_file does not exist!"
            exit 1
        fi
        echo "Copying $myconfig.h to $builddir/myconfig.h..."
        cp $myconfig_file $builddir/myconfig.h
    fi
    
    make_params="-j $build_procs"

    cmd "make $make_params"
    end "BUILD"
fi    

# CHECK
if $run_check; then
    start "TEST"
    # something should be done after ||, otherwise Jenkins will mark
    # job as failed
    cmd "make check processors=\"$check_procs\" || CHECK_UNSTABLE=1"
    end "TEST"
fi

# DOC
if $run_doc; then
    start "DOC"
    cmd "make doc"
    end "DOC"
fi

# DIST
if $run_dist; then
    start "DIST"
    cmd "make dist dist-xz"
    end "DIST"
fi

if $run_distcheck; then
    start "DISTCHECK"
    cmd "make distcheck DISTCHECK_CONFIGURE_FLAGS=\"$distcheck_flags\""
    end "DISTCHECK"
fi

# return to calling dir
popd

echo "Finished."
exit 0
