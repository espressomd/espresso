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
setmd box_l 10.0 8.0 6.0
puts "box =\{[setmd box]\}"
setmd procgrid 2 2 2
puts "grid = [setmd proc]"
setmd nptypes 2
puts "nptypes = [setmd nptypes]"
setmd niatypes 1
puts "niatypes = [setmd niatypes]"

# initialize and broadcast system parameters
##################################################
# done automagically

# Set up particle positions
##################################################
part 0 pos 1.0 2.0 1.0
puts "part 0 = [part 0]"
part 1 pos 1.0 1.0 1.0
puts "part 1 = [part 1]"
part 2 pos 10.0 1 1
puts "part 2 = [part 2]"
part 3 pos 6.0 1 1
puts "part 3 = [part 3]"
part 4 pos 3.0 1 1
puts "part 4 = [part 4]"
puts "npart = [setmd npart]"


# integration
##################################################
integrate init
integrate 2
integrate exit

# exit
##################################################
puts "finished"
exit
