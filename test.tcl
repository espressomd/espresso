#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# \
    if test $PLATFORM = OSF1; then MPIRUN=dmpirun;
# \
    else MPIRUN=mpirun; lamboot; fi
# \
    exec $MPIRUN -nsigs -np 8 $PLATFORM/tcl_md $0 $*

# data initialization
##################################################
puts "nproc = [setmd nproc]"
setmd box_l 10.0 8.0 16.0
puts "box =\{[setmd box]\}"
setmd procgrid 2 2 2
puts "grid = \{[setmd proc]\}"
setmd niatypes 1
puts "niatypes = [setmd niatypes]"

# setup interactions
##################################################
inter 0 1 lennard-jones 1 1 1.2 0 0
inter 1 2 lennard-jones 2 1 1.2 0 0
inter 0 2 lennard-jones 3 1 1.2 0 0

# Set up particle positions
##################################################
part 0 pos 1.0 2.0 1.0 type 0 q 1
puts "part 0 = [part 0]"
part 1 pos 1.0 1.0 1.0 type 2 q 1
puts "part 1 = [part 1]"
part 2 pos 10.0 1 1 type 1 q -2
puts "part 2 = [part 2]"
part 3 pos 6.0 1 1 type 3 q -1
puts "part 3 = [part 3]"
part 4 pos 3.0 1 1 type 0 q 1
puts "part 4 = [part 4]"


puts "nptypes = [setmd nptypes]"
puts "npart = [setmd npart]"

# pump up
##################################################

for {set i 5} { $i < 1000 } { incr i} {
    if {[expr $i % 100 == 0]} {
	puts "adding part $i"
    }
    part $i pos [expr 10*rand()] [expr 10*rand()] [expr 10*rand()] \
	q [ expr ($i % 2 == 1) ? -1 : 1 ] \
	type 0
}

# integration
##################################################
integrate init
integrate 2
integrate exit

# exit
##################################################
puts "finished"
exit
