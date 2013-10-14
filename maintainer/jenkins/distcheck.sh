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
[ ! -v distcheck_flags ] && distcheck_flags=""
outp distcheck_flags

# change into build dir
pushd $builddir

start "DISTCHECK"
cmd "make distcheck DISTCHECK_CONFIGURE_FLAGS=\"$distcheck_flags\"" || exit $?
end "DISTCHECK"

popd
