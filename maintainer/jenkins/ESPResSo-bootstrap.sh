#!/bin/bash --login
# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
DIR=`dirname $0`
run_bootstrap=true run_configure=false run_build=false run_check=false \
    $DIR/build.sh
