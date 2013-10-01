#!/bin/bash --login -e
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
source maintainer/jenkins/common.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

start "CONFIGURE"
./configure --with-cuda --with-mpi CPU_COUNT="3"
end "CONFIGURE"

# copy config file
if [ "$config" != default ]; then
  use_config $config
fi

# create mympiexec.sh
echo 'exec mpiexec --bind-to-core $@' > mympiexec.sh
chmod +x mympiexec.sh

start "BUILD"
make -j 4
end "BUILD"

check
