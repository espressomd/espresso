#!/bin/bash --login -e
BUILDDIR="$(pwd)"
SOURCEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../.."
source $SOURCEDIR/maintainer/jenkins/common.sh

cd $SOURCEDIR; bootstrap

cd $BUILDDIR


use_myconfig compileonly

start "BUILD"
make -j 4
end "BUILD"
