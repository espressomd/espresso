# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
if ((`ps -o nice= $$` < 5)); then 
    renice -n 5 $$
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

function outp() {
    for p in $*; do
        echo "  $p=${!p}"
    done
}

# DIR SETTINGS
[ -z "$insource" ] && insource="true"
[ -z "$srcdir" ] && srcdir=`pwd`
if $insource; then
    builddir=$srcdir
elif [ -z "$builddir" ]; then
    builddir=$srcdir/build
fi

outp srcdir builddir insource 

if [ ! -e $srcdir/configure.ac ]; then
    echo "Could not find configure.ac in $srcdir!"
    exit 1
fi

if ! $insource; then
    if [ ! -d $builddir ]; then
        echo "Creating $builddir..."
        mkdir -p $builddir
    fi
fi

# PARALLELISM SETTINGS
if [ -n "$parallelism" ]; then
    outp parallelism

    case $parallelism in
        ("34core") 
            with_mpi=true 
            check_procs="3 4"
            ;;
        ("3core") 
            with_mpi=true 
            check_procs=3
            ;;
        ("4core")
            with_mpi=true
            check_procs=4
            ;;
        ("nompi")
            with_mpi=false
            ;;
    esac
fi

