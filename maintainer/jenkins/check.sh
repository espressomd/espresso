#!/bin/bash --login 
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
DIR=`dirname $0`
source $DIR/common.sh

start "TEST"
[ -z "$with_mpi" ] && with_mpi="true"
! $with_mpi && check_procs=1
[ -z "$with_cuda" ] && with_cuda="true"
[ -z "$check_procs" ] && check_procs=4

outp with_mpi with_cuda check_procs

$with_cuda && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

# change into build dir
pushd $builddir

[ "$check_procs" != "1" ] &&  make_params="processors=\"$check_procs\""
# something should be done after ||, otherwise Jenkins will mark
# job as failed
cmd "make check $make_params" || CHECK_UNSTABLE=1

popd
end "TEST"
