#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=1; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else lamboot; exec /usr/lib/lam/bin/mpirun -np $NP -nsigs $PLATFORM/tcl_md $0 $*;

# \
    fi;

########### parameters
# number of particles to setup
set npart 500
# particle enumeration (1,3,5,...)
set firstpart 1
set partstep 2
# box size
setmd box_l 20.0 20.0 20.0
# number of particle types
set ntypes 2
# minimal distance of particles at finish
set mdst "guess"
# consecutive integration steps between two tests
set intsteps 100
# how many tests for minimal distance
set maxtime 200
setmd skin 1


#######################################################\

set Lx [lindex [setmd box_l] 0]
set Ly [lindex [setmd box_l] 1]
set Lz [lindex [setmd box_l] 2]

if {$mdst == "guess"} {
    set mdst [expr 0.9*pow($Lx*$Ly*$Lz/$npart,0.333333333)]
    puts "trying minimal distance of $mdst"
}

# setup random particles
puts "setting up random particles"

for {set i $firstpart; set np 0} { $np < $npart } { incr np; incr i $partstep} {
    part $i pos [expr $Lx*[tcl_rand]] [expr $Ly*[tcl_rand]] [expr $Ly*[tcl_rand]] \
	q [ expr ($i % 2 == 1) ? -1 : 1 ] \
	type [ expr round([tcl_rand]*$ntypes) ]
}

# no electrostatics
setmd bjerrum 0

#pairwise ramp for all particles
for {set ia1 0} { $ia1 <= $ntypes } { incr ia1 } {
    for {set ia2 0} { $ia2 <= $ntypes } { incr ia2 } {
	inter $ia1 $ia2 ramp $mdst 100
    }
}
# friction
setmd gamma 1e4

puts "starting ramp integration"
integrate init

set cont 1
for {set i 0} { $i < $maxtime && $cont} { incr i} {
    set md [mindist]
    puts "step $i minimum distance = $md"
    if {$md >= $mdst || $md == "(> ramp cutoffs)"} { set cont 0 }
    integrate $intsteps
}

# write
set f [open "|gzip -c - >config.gz" w]
blockfile $f write variable box_l
blockfile $f write particles "id pos type q" all
close $f

puts "wrote [part number] particles to config.gz with minimal distance [mindist]"
