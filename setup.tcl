#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $PLATFORM/Espresso $0 $*;
# \
    fi;

########### parameters
# number of particles to setup
set npart 100
# particle enumeration (1,3,5,...)
set firstpart 1
set partstep 2
# box size
setmd box_l 10.0 10.0 10.0
# number of particle types
set ntypes 1
# minimal distance of particles at finish
set mdst 2
# consecutive integration steps between two tests
set intsteps 100
# how many tests for minimal distance
set maxtime 10
# integrator
setmd skin 0.3
setmd time_step 0.0001

#######################################################\

set Lx [lindex [setmd box_l] 0]
set Ly [lindex [setmd box_l] 1]
set Lz [lindex [setmd box_l] 2]

# setup random particles
puts "setting up random particles"

for {set i $firstpart; set np 0} { $np < $npart } { incr np; incr i $partstep} {
    part $i pos [expr $Lx*[t_random]] [expr $Ly*[t_random]] [expr $Ly*[t_random]] \
	q [ expr ($i % 2 == 1) ? -1 : 1 ] \
	type [ expr round([t_random]*$ntypes) ]
}

# no electrostatics
inter coulomb 0

#pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 <= $ntypes } { incr ia1 } {
    for {set ia2 0} { $ia2 <= $ntypes } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones 3 1 1.12246 0 0
	puts [inter $ia1 $ia2]
    }
}

# cap lennard_jones
inter ljforcecap 500

# friction
setmd gamma [expr 1e4*[setmd time_step]]
setmd temp 1.

puts "starting warmup integration"


set cont 1
for {set i 0} { $i < $maxtime && $cont} { incr i} {
    set md [analyze mindist]
    puts "step $i minimum distance = $md"
    if {$md >= $mdst} { set cont 0 }

    integrate $intsteps
}

# write
set f [open "|gzip -c - >config.gz" w]
blockfile $f write variable all
blockfile $f write particles "id pos type q" all
blockfile $f write interactions
blockfile $f write bonds all
close $f
puts "wrote [setmd n_part] particles to config.gz with minimal distance [analyze mindist]"

part delete
set f [open "|gzip -cd config.gz" r]
while {[blockfile $f read auto] != "eof"} {}
close $f
puts "[setmd box_l] read [setmd n_part] particles from config.gz with minimal distance [analyze mindist]"
