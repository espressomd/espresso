#!/bin/bash --login -e
BUILDDIR="$(pwd)"
SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."
source $SRCDIR/maintainer/jenkins/common.sh

bootstrap

configure --with-mpi CPU_COUNT="3"

# copy config file
if [ "$myconfig" != default ]; then
  use_myconfig $myconfig
fi

# create mympiexec.sh
echo 'exec mpiexec --bind-to-core $@' > mympiexec.sh
chmod +x mympiexec.sh

start "BUILD"
make -j 3
end "BUILD"

check
