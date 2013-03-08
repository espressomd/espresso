#!/bin/bash --login -e
source maintainer/jenkins/common.sh
bootstrap

start "CONFIGURE"
./configure \
    CC=/home/jenkins/mpich/bin/mpicc \
    MPIEXEC=/home/jenkins/mpich/bin/mpiexec \
    CPU_COUNT=2
end "CONFIGURE"

use_myconfig maxset

start "BUILD"
make -j 2
end "BUILD"

check 
