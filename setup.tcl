#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np 8 $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs 8
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np 8 -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

set npart 100

setmd box_l 10.0 10.0 10.0

if {[setmd nproc] == 8} {
    setmd procgrid 2 2 2
}

# setup random particles
puts "setting up random particles"
expr srand([pid])
for {set i 0} { $i < $npart } { incr i} {
    part $i pos [expr 10*rand()] [expr 10*rand()] [expr 10*rand()] \
	q [ expr ($i % 2 == 1) ? -1 : 1 ] \
	type [ expr round(rand()*3) ]
}
	
# setup ramp ia
puts "setting up ramp interaction"
for {set ia1 0} { $ia1 <= 4 } { incr ia1 } {
    for {set ia2 0} { $ia2 <= 4 } { incr ia2 } {
	inter $ia1 $ia2 ramp 1 1
    }
}

set gamma 1
# write
set f [open "|gzip -c - >setup.gz" w]
for {set i 0} { $i < $npart } { incr i} {
    puts $f [part $i]
}
close $f

# relax
puts "starting ramp integration"
integrate init
for {set i 0} { $i < 100} { incr i} {
    puts "$i% done"
    integrate 1
}

# write
set f [open "|gzip -c - >config.gz" w]
writemd $f posx posy posz type q
close $f
