#!/bin/bash --login -e
source maintainer/jenkins/common.sh

bootstrap

start "CONFIGURE"
./configure CPU_COUNT=4
end "CONFIGURE"

use_myconfig maxset

start "BUILD"
make -j 4
end "BUILD"

check
doc
dist
