#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# \
    if test $PLATFORM = OSF1; then MPIRUN=dmpirun;
# \
    else MPIRUN=mpirun; fi
# \
    exec $MPIRUN -np 8 $PLATFORM/tcl_md $0 $*

# data initialization
##################################################
puts "nproc = [setmd nproc]"
setmd npart 1000
puts "npart = [setmd npart]"
setmd box_l 25.0 35.0 12.0
puts "box =\{[setmd box]\}"
setmd procgrid 2 2 2
puts "grid = [setmd proc]"
part 1 pos 1 1 1
puts "part 1 = [part 1]"

# integration
##################################################
integrate 0
integrate 10
integrate -1

# exit
##################################################
puts "finished"
exit
