#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Test System 3: Polyelectrolyte Solution                  #
#                 (Poor Solvent)                            #
#                                                           #
#  Created:       26.03.2003 by HL                          #
#  Last modified: 26.03.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                 pe_solution.tcl                     ="
puts "======================================================="
puts " "

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "pe_solution"
set ident "_t3"

# System parameters
#############################################################

# Number of polymers
set n_polymers 5
# length of polymers
set l_polymers 190
# distance between charged monomers
set cm_distance 3
# system density
set density 0.00002

puts [format "%d %d %d" $n_polymers $l_polymers $cm_distance]

# Interaction parameters
#############################################################

# Lennard Jones
set n_part_types  3
# repulsive LJ
set lj1_cut       1.122462048309373
set lj1_epsilon   1.0
# LJ with attractive tail
set lj2_cut       2.5
set lj2_epsilon   1.75

# FENE
set fene_cut      2.0
set fene_k        7.0

# Coulomb
set bjerrum       1.0
set accuracy      1.0e-4

# Integration parameters
#############################################################

setmd time_step 0.0125
setmd skin      0.4
setmd gamma     1.0
setmd temp      1.0

# warmup integration (with capped LJ potential)
set warm_steps   200
set warm_n_times 30
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    800
set int_n_times  1000


# More general Parameters
#############################################################

set tcl_precision 6
set vmd_output    "no"

# Initial Bond length
set bond_length   1.0
# valency of charged monomers
set cm_valency    1.0
# valency of counterions
set ci_valency    -1.0

puts "Simulate the following polyelectrolyte solution:"
puts "number of Polymers       : $n_polymers"
set n_monomers [expr ($n_polymers*$l_polymers)]
puts "number of monomers       : $n_monomers"
set polymer_charge [expr (($l_polymers+$cm_distance-1)/$cm_distance)]
puts "polymer_charge           : $polymer_charge"
set n_counterions [expr -int($polymer_charge*$n_polymers/$ci_valency)]
puts "n_counterions            : $n_counterions"
set npart [expr $n_counterions+$n_monomers]
puts "total number of particles: $npart"
set box_length [expr pow(($npart/$density),(1.0/3.0))]
puts "cubic box length         : $box_length"

# Dependent Parameters
set lj1_shift  [expr -4.0*$lj1_epsilon*(pow($lj1_cut,-12)-pow($lj1_cut,-6))]
set lj2_shift  [expr -4.0*$lj2_epsilon*(pow($lj2_cut,-12)-pow($lj2_cut,-6))]
set lj1_cutforce [expr 48*$lj1_epsilon*(0.5*pow($lj1_cut,-7)-pow($lj1_cut,-13))]
set lj2_cutforce [expr 48*$lj2_epsilon*(0.5*pow($lj2_cut,-7)-pow($lj2_cut,-13))]
#puts "LJ shifts:     1: $lj1_shift           2: $lj2_shift"
#puts "LJ Cut forces: 1: $lj1_cutforce   2: $lj2_cutforce"

# Setup
#############################################################

# simulation box
setmd box_l $box_length $box_length $box_length 

# fene
inter 0 fene $fene_k $fene_cut

# pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj1_epsilon 1.0 $lj1_cut $lj1_shift 0.0
    }
}

puts "Setup Particles (wait...)"
# polymers
polymer $n_polymers $l_polymers $bond_length mode SAW charge $cm_valency distance $cm_distance types 0 1 FENE 0
# counterions
counterions $n_counterions charge $ci_valency type 2
#puts "[part]"

set act_min_dist [analyze mindist]
puts "Placed [setmd n_part] particles with minimal distance: $act_min_dist"

#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter ljforcecap $cap 
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    set time [format "%8f" [setmd time]]
    puts -nonewline "run $i at time=$time (LJ cap=$cap) "

    integrate $warm_steps
    if { $vmd_output=="yes" } { imd positions }

    set act_min_dist [analyze mindist]
    puts -nonewline "minimal distance = $act_min_dist\r"
    flush stdout
#   Increase LJ cap
    set cap [expr $cap+10]
    inter ljforcecap $cap
    incr i
}

puts "\nWarmup done: Minimal distance [analyze mindist]"
puts "Verlet reuses: [setmd verlet_reuse]"

#############################################################
#      Integration                                          #
#############################################################

puts "Prepare integration: (tune parameters - wait...)"
inter ljforcecap 0
setmd time 0.0
# Set attractive LJ for monomer interactions (types 0 and 1
for {set ia1 0} { $ia1 < 2 } { incr ia1 } {
    for {set ia2 0} { $ia2 < 2 } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj2_epsilon 1.0 $lj2_cut $lj2_shift 0.0
    }
}
puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy mesh 16]"
puts "Interactions are now: {[inter]}"

#      Write blockfiles for restart
#############################################################
#polyBlockWrite "$name$ident.set"   {box_l time_step skin temp gamma}
#polyBlockWrite "$name$ident.start" {time} {id pos type}

# prepare observable output
set obs_file [open "$name$ident.obs" "w"]
puts $obs_file "\#$name$ident: Observables"
puts $obs_file "\#Time     R_E         R_G         R_H         E_TOT       E_KIN       E_C_k"
analyze set chains 0 $n_polymers $l_polymers
set re [analyze re]

puts "\nStart Main integration: $int_n_times times $int_steps steps"
for {set i 0} { $i < $int_n_times } { incr i} {
    set time [setmd time]
    puts -nonewline "run $i at time=$time, R_E = $re (soll 15.0, deviation [expr $re/15.0])\r"
    flush stdout

    integrate $int_steps

    set energy [analyze energy]
    set re [analyze re]
    puts $obs_file [format "%.3e %.5e %.5e %.5e %.5e %.5e %.5e" $time $re [analyze rg] [analyze rh] [lindex [lindex $energy 0] 1] [lindex [lindex $energy 1] 1] [lindex [lindex $energy [expr [llength $energy]-1]] 1] ]
    flush $obs_file
    if { $vmd_output=="yes" } { imd positions }

}

puts "\nIntegration done."

#############################################################

puts "\nFinished"