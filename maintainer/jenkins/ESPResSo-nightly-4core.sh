#!/bin/bash --login -e
BUILDDIR="$(pwd)"
SOURCEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."
source $SOURCEDIR/maintainer/jenkins/common.sh

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

