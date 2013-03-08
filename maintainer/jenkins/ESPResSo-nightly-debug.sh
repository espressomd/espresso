#!/bin/bash --login -e
source maintainer/jenkins/common.sh

bootstrap

start "CONFIGURE"
./configure
end "CONFIGURE"

use_myconfig compileonly

start "BUILD"
make -j 4
end "BUILD"
