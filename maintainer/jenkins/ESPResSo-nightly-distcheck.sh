#!/bin/bash --login -ex
# Copyright (C) 2013 Olaf Lenz
# 
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
source maintainer/jenkins/common.sh

# clean up dir before build
chmod -R u+w espresso-* || /bin/true
git clean -fd

bootstrap
export LD_LIBRARY_PATH=/usr/local/cuda/lib64

start "CONFIGURE"
./configure
end "CONFIGURE"

start "DISTCHECK"
make distcheck DISTCHECK_CONFIGURE_FLAGS="CPU_COUNT=2"
end "DISTCHECK"
