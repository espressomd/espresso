#!/bin/sh
# tricking...\
    PLATFORM=`uname -s`;
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np 1 $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs 1
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np 8 -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

########### parameters

set maxpart 999
set step 10
setmd box_l 10.0 10.0 10.0

set write finish
set write_steps 10
set mdst 1.2
set maxtime 200

setmd periodic 1 1 1
setmd bjerrum 0
setmd max_num_cells 512
setmd skin 0.4

source polywr.tcl

if { $write  == "yes" || $write == "finish" } {
    exec rm -f [glob -nocomplain "configs/c*"] [glob -nocomplain "movie/m*"]
}
 
if {[setmd n_node] == 8} {
    setmd node_grid 2 2 2
}
if {[setmd n_node] == 1} {
    setmd node_grid 1 1 1
}
# setup random particles
puts "setting up random particles"
for {set i 0} { $i <= $maxpart } { incr i $step} {
    part $i pos [expr 10*[tcl_rand]] [expr 10*[tcl_rand]] [expr 10*[tcl_rand]] \
	q [ expr ($i % 2 == 1) ? -1 : 1 ] \
	type [ expr round([tcl_rand]*3) ]
}
puts "total [part number] particles"

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
    catch { eval exec poly2pdb -poly [glob -nocomplain "configs/c*"] \
	-per -psf "movie/m" >/dev/null 2>/dev/null }
}
