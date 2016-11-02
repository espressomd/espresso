#!/bin/bash --login 
# Copyright (C) 2015 Florian Weik
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.

DIR=`dirname $0`

source $DIR/common.sh

$DIR/../../Espresso $DIR/scripts/lj_performance.tcl

paste -d ',' lj_performance.txt > performance.txt

