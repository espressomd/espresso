#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# \
    if test $PLATFORM = OSF1; then MPIRUN=dmpirun;
# \
    else MPIRUN=mpirun; fi
# \
    exec $MPIRUN -np 8 $PLATFORM/tcl_md $0 $*

puts "nproc = [setmd nproc]"
setmd npart 1000
puts "npart = [setmd npart]"
setmd box_l 2.5 3.5 1.2
puts "box =\{[setmd box]\}"
setmd procgrid 2 2 2
puts "grid = [setmd proc]"
part 1 pos 1 1 1
puts "part 1 = [part 1]"
puts "finished"
exit
