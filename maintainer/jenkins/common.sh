#!/bin/bash

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

[ ! -v insource ] && insource="true"
[ ! -v srcdir ] && srcdir=`pwd`
if $insource; then
    builddir=$srcdir
else
    [ ! -v builddir ] && builddir="$srcdir/build"
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

