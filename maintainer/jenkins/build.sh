#!/bin/bash --login 
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
DIR=`dirname $0`
source $DIR/common.sh

# DEFAULTS
[ ! -v myconfig ] && myconfig="default"
[ ! -v build_procs ] && build_procs=4
outp myconfig build_procs

# change into build dir
pushd $builddir

# BUILD
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
cmd "make $make_params" || exit 1
end "BUILD"

popd
