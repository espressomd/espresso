#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;

#
#  This file is part of the internal section of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It should never be given outside of the institute without the explicit approval of the author.
#  It is nonetheless subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#  
set ia fene
#set ia harm
puts "Creating data for Testcase $ia.tcl, running on [setmd n_nodes] nodes:"
### set N_P     { 50   25   30   32   20   16   20    20    20    20    100   }
### set MPC     { 5    10   20   25   30   50   75    100   150   200   200   }
### set box_l   { 6.65 6.65 8.90 9.80 8.90 9.80 12.08 13.30 15.23 16.76 28.66 }
set N_p     20
set MPC     30
set bond_l  0.97
set shield  0.2
set box_l   8.90
set density 0.85

# repulsive Lennard Jones
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   0.25
set lj1_off     0.0
set lj          0

# attractive FENE / harmonic
set fene_k      30.0
set fene_r      1.5

setmd time_step 0.006
setmd skin      0.0
setmd gamma     0.0
setmd temp      0.0

set epsilon 1e-4


proc write_data {file} {
    global energy pressure
    set f [open $file "w"]
    set energy [analyze energy total]
    set pressure [analyze pressure total]
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write variable box_l
    blockfile $f write particles {id pos f}
    blockfile $f write bonds all
    close $f
}


############## $ia-specific part

setmd box_l $box_l $box_l $box_l

if {$lj==1} { inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off }
inter 0  $ia  $fene_k  $fene_r

set n_part [expr $N_p*$MPC]
puts -nonewline "Creating LJ-polymers with $ia-chains... "; flush stdout
polymer $N_p $MPC $bond_l mode SAW $shield
puts "Done."


# to ensure force recalculation
invalidate_system
integrate 0

# ensures that no other forces are on
set cureng [analyze energy $ia 0]
# tbrs
set curprs [lindex [analyze pressure $ia 0] 0]

############## end

set toteng [analyze energy total]
set totprs [analyze pressure total]

if { [expr abs($toteng - $cureng)] > $epsilon } {
    error "system has unwanted energy contributions"
}
if { [expr abs($totprs - $curprs)] > $epsilon } {
    error "system has unwanted pressure contributions"
}

set energy $cureng
set pressure $curprs
write_data "$ia\_system.data"

exit