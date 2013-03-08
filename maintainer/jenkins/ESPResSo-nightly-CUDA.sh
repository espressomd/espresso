#!/bin/bash --login -e
source maintainer/jenkins/common.sh

bootstrap

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

start "CONFIGURE"
./configure --with-mpi --with-cuda CPU_COUNT="4"
end "CONFIGURE"

use_myconfig LBGPU

start "BUILD"
make -j 4
end "BUILD"

check

