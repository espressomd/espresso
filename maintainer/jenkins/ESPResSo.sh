#!/bin/bash --login -e
BUILDDIR="$(pwd)"
SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."
source $SRCDIR/maintainer/jenkins/common.sh

bootstrap

configure CPU_COUNT=4

use_myconfig maxset

start "BUILD"
make -j 4
end "BUILD"

check
doc
dist
