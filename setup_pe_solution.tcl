#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Setup a Polyelectrolyte solution                         #
#                                                           #
#                                                           #
#  Created:       23.01.2003 by HL                          #
#  Last modified: 23.01.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "============================================================="
puts "=              setup_pe_solution.tcl                        ="
puts "============================================================="
puts " "

#############################################################
#  Parameters                                               #
#############################################################

set vmd_output "no"

# System parameters
#############################################################

# Number of polymers
set n_polymers 5
# length of polymers
set l_polymers 238
# Bond length
set bond_length 1.0
# distance between charged monomers
set cm_distance 3
# valency of charged monomers
set cm_valency 1.0

# valency of counterions
set ci_valency -1.0

# system density
set density 0.00002

# Interaction parameters
#############################################################

set n_part_types 3

set lj1_cut       1.12246
set lj1_epsilon   1.0
set lj1_shift     1.0
set lj1_cap_start 10.0
set lj1_cap_incr  5.0

set lj2_cut       2.5
set lj2_epsilon   1.75

set fene1_cut   2.0
set fene1_k     7.0

setmd skin      0.8

# Integration parameters
#############################################################

set mdst 0.90
setmd time_step 0.01
set intsteps 100
set maxtime_w  100
set maxtime_i1 100
set maxtime_i2 100000

setmd gamma 1.0
setmd temp 1.0

# Other parameters
#############################################################
set tcl_precision 5
set mypi 3.141592653589793

#############################################################
#  Setup System                                             #
#############################################################

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
setmd box_l $box_length $box_length $box_length 

# Interaction setup
#############################################################

# fene
inter 1 fene $fene1_k $fene1_cut

# pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj1_epsilon 1 $lj1_cut $lj1_shift 0
    }
}
set lj1_cap $lj1_cap_start
setmd lj_force_cap $lj1_cap
setmd bjerrum 0.0

# Particle setup
#############################################################

# check existing warmup configuration
if { [file exists "pe_warmup.gz"]} {
    set f [open "|gzip -cd pe_warmup.gz" r]
    while {[blockfile $f read auto] != "eof" } {}
    close $f
    puts "Found file pe_warmup.gz with:"
    puts "n_part = [part number]"
    puts "grid   = \{[setmd node_grid]\}"
    puts "box    = \{[setmd box]\}" 
    if { $npart != [part number] } {
	puts "Found configuration differs from this setup. Create new Startconfiguration"
	set warmup 0
    } else {
	set warmup 1
    }
} else {
    set warmup 0
}

# create new start configuration if necessary
if { $warmup==0 } {
    puts "Create new Start Configuration..."
    # Polymer setup:
    set part_id 0
    for {set i 0} { $i < $n_polymers } {incr i} {
	set posx [expr $box_length*[tcl_rand]]
	set posy [expr $box_length*[tcl_rand]]
	set posz [expr $box_length*[tcl_rand]]
	part $part_id pos $posx $posy $posz \
	    q $cm_valency type 1
	incr part_id 
	for {set n 1} { $n < $l_polymers } {incr n} {
	    set phi [expr 2.0*$mypi*[tcl_rand]]
	    set theta [expr $mypi*[tcl_rand]]
	    set vecx [expr $bond_length*sin($phi)*sin($theta)]
	    set vecy [expr $bond_length*cos($phi)*sin($theta)]
	    set vecz [expr $bond_length*cos($theta)]
	    set posx [expr $posx+$vecx]
	    set posy [expr $posy+$vecy]
	    set posz [expr $posz+$vecz]
	    part $part_id pos $posx $posy $posz \
		q [expr ($n%$cm_distance == 0) ? $cm_valency : 0] \
		type  [expr ($n%$cm_distance == 0) ? 1 : 0] \
		bond 1 [expr $part_id-1]
	    incr part_id 
	}
    }
    # Counterion setup
    for {set i 0} { $i < $n_counterions } {incr i} {
	set posx [expr $box_length*[tcl_rand]]
	set posy [expr $box_length*[tcl_rand]]
	set posz [expr $box_length*[tcl_rand]]
	part $part_id pos $posx $posy $posz \
	    q $ci_valency type 2
	#puts "[part $part_id]"
	incr part_id 
    }
} else {
# setup the bonds...
    set part_id 0
    for {set i 0} { $i < $n_polymers } {incr i} {
	incr part_id 
	for {set n 1} { $n < $l_polymers } {incr n} {
	    part $part_id bond 1 [expr $part_id-1]
	    incr part_id 
	}
    }
}

#############################################################
#  Prepare IMD Connection                                   #
#############################################################

if { $vmd_output=="yes" } {
    writepsf pe_sol.psf
    writepdb pe_sol.pdb
}

#############################################################
#  Warmup Integration                                       #
#############################################################
if { $warmup == 0 } {
    puts "start warmup integration"
    set cont 1
    for {set i 0} { $i < $maxtime_w && $cont} { incr i} {
	set md [analyze mindist]
	puts "step $i minimum distance = $md, force cap = $lj1_cap"
	if {$md >= $mdst} { set cont 0 }
	setmd lj_force_cap $lj1_cap
	integrate $intsteps
	puts "[analyze energy lj 1 1]"
	puts "[analyze energy]"
	set lj1_cap [expr $lj1_cap + $lj1_cap_incr]
    }
    # write
    set f [open "|gzip -c - >pe_warmup.gz" w]
    blockfile $f write variable box_l
    blockfile $f write particles "id pos type q" all
    close $f
    puts "wrote [part number] particles to pe_warmup.gz with minimal distance [analyze mindist]"
}

#############################################################
#      Integration                                          #
#############################################################

# turn off potentil caps
setmd lj_force_cap 0
# electrostatics
inter coulomb 1.5 p3m 70.0 16 3 0.026609

puts "[inter]"

if { $vmd_output=="yes" } {
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    puts "opened port $port"
}

#puts "[analyze]"

#puts [part 1 print bonds]


#analyze
set f [open "pe_obs.dat" w]
puts $f "\#    RE            RG"


puts "start integration 1 (good solvent)"
for {set i 0} { $i < $maxtime_i1 } { incr i} {
    puts "step $i at time=[setmd time]"
    integrate $intsteps
    puts "[analyze energy]"
    if { $vmd_output=="yes" } { imd positions }

}

puts "Set attractive LJ between monomers"
inter 0 0 lennard-jones $lj2_epsilon 1 $lj2_cut 0 0
inter 0 1 lennard-jones $lj2_epsilon 1 $lj2_cut 0 0
#inter 1 0 lennard-jones $lj2_epsilon 1 $lj2_cut 0 0
inter 1 1 lennard-jones $lj2_epsilon 1 $lj2_cut 0 0

puts "[inter]"

puts "start integration 2 (poor solvent)"
for {set i 0} { $i < $maxtime_i2 } { incr i} {
    puts "step $i at time=[setmd time]"
    integrate $intsteps
}

close $f

if { $vmd_output=="yes" } { puts "imd: [imd disconnect]" }
