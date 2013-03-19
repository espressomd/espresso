#!/bin/bash --login -e
BUILDDIR="$(pwd)"
SOURCEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."
source $SOURCEDIR/maintainer/jenkins/common.sh

cd $SOURCEDIR; bootstrap

cd $BUILDDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

start "CONFIGURE"
./configure --with-mpi --with-cuda CPU_COUNT="4"
end "CONFIGURE"

use_myconfig LBGPU

start "BUILD"
make -j 4
end "BUILD"

check

