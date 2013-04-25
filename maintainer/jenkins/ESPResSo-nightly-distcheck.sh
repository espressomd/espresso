#!/bin/bash --login -ex
source maintainer/jenkins/common.sh

# clean up dir before build
chmod -R u+w espresso-* || /bin/true
git clean -fd

bootstrap

start "CONFIGURE"
./configure
end "CONFIGURE"

start "DISTCHECK"
make distcheck DISTCHECK_CONFIGURE_FLAGS="CPU_COUNT=2"
end "DISTCHECK"
