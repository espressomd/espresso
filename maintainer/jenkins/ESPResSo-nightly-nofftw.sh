#!/bin/bash --login -e
source maintainer/jenkins/common.sh

echo "#error ERROR: fftw not really present but used somewhere." > src/fftw3.h

bootstrap

start "CONFIGURE"
./configure --without-fftw CPU_COUNT=2
end "CONFIGURE"

use_myconfig nofftw

start "BUILD"
make -j 2
end "BUILD"

check
