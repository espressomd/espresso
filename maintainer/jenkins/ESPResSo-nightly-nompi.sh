#!/bin/bash --login -e
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
source maintainer/jenkins/common.sh

bootstrap

start "CONFIGURE"
./configure --without-mpi
end "CONFIGURE"

# copy config file
if [ "$myconfig" != default ]; then
  use_myconfig $myconfig
fi

# create mympiexec.sh
echo 'exec mpiexec --bind-to-core $@' > mympiexec.sh
chmod +x mympiexec.sh

start "BUILD"
make 
end "BUILD"

check
