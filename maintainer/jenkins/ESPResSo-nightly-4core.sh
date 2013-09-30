#!/bin/bash --login -e
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
source maintainer/jenkins/common.sh

cd $SOURCEDIR; bootstrap

cd $BUILDDIR
start "CONFIGURE"
./configure --with-mpi CPU_COUNT="4"
end "CONFIGURE"

# copy config file
if [ "$myconfig" != default ]; then
  use_myconfig $myconfig
fi

# create mympiexec.sh
echo 'exec mpiexec --bind-to-core $@' > mympiexec.sh
chmod +x mympiexec.sh

start "BUILD"
make -j 4
end "BUILD"

check

