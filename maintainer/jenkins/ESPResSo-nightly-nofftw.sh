#!/bin/bash --login -e
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
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
