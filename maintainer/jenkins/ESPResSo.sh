#!/bin/bash --login -e
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
BUILDDIR="$(pwd)"
SOURCEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."

source maintainer/jenkins/common.sh

bootstrap

configure CPU_COUNT=4

use_myconfig maxset

start "BUILD"
make -j 4
end "BUILD"

check
doc
dist
