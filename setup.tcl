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

########### parameters
set npart 300
setmd box_l 10.0 10.0 10.0
set write finish
set write_steps 10
set mdst 1.2
set maxtime 200

setmd periodic 1 1 1
setmd bjerrum 0

source polywr.tcl

if { $write  == "yes" || $write == "finish" } {
    exec rm -f [glob -noc "configs/c*"] [glob -noc "movie/m*"]
}
 
if {[setmd n_node] == 8} {
    setmd node_grid 2 2 2
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
# no electrostatics
setmd bjerrum 0

#pairwise ramp for all particles
for {set ia1 0} { $ia1 <= 4 } { incr ia1 } {
    for {set ia2 0} { $ia2 <= 4 } { incr ia2 } {
	inter $ia1 $ia2 ramp 1.3 100
    }
}

# test for bonded interaction parameters
# friction
setmd gamma 1e4

if {"$write" == "yes" } {
    exec rm -f "configs/*" "movie/*"
}
# relax
puts "starting ramp integration"
integrate init

set cont y
for {set i 0} { $i < 100 && $cont == "y" } { incr i} {
    puts "step $i: minimum distance=[mindist] !< $mdst"
    set md [mindist]
    if {$md > $mdst || $md == "(> ramp cutoffs)"} { set cont n }
    integrate $write_steps
    if {"$write" == "yes" } {
	polywrite [format "configs/c%04d.poly" $i]
    }
}

puts "final: minimum distance=[mindist]"

# write
set f [open "|gzip -c - >config.gz" w]
writemd $f posx posy posz type q
close $f

if { $write == "finish" } {
    polywrite "configs/cfinal.poly"
}
if { $write  == "yes" || $write == "finish" } {
    catch { eval exec poly2pdb -poly [glob -noc "configs/c*"] \
	-per -psf "movie/m" >/dev/null 2>/dev/null }
}
