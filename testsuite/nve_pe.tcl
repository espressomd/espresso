# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
# 
#############################################################
#                                                           #
#  Test System: Single PE chain fixed at end                #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "LENNARD_JONES"
require_feature "ELECTROSTATICS"
require_feature "FFTW"
require_feature "BOND_ANGLE_COSINE"
require_feature "EXTERNAL_FORCES"
require_feature "ADRESS" off

puts "----------------------------------------------"
puts "- Testcase nve_pe.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

proc hexconvert { vec shift_vec d_space} {
    set hvec {{1. 0.5 0.5} {0.0 0.8660254 0.28867513} {0. 0. 0.81649658}}
    set rvec {0 0 0}
    set dim 3
    for {set j 0} { $j < $dim } {incr j} {
        for {set i 0} { $i < $dim } {incr i} {
	    lset rvec $j [expr [lindex $rvec $j] + [lindex $vec $i] * [lindex [lindex $hvec $j] $i]]
        }
    }
    lsqr $rvec
    for {set j 0} { $j < $dim } {incr j} {
        lset rvec $j [expr ([lindex $d_space $j] * [lindex $rvec $j] + [lindex $shift_vec $j])]
    }
    return $rvec
}

proc create_chain { part_id rvec p_length b_length} {
    set posx [lindex $rvec 0]
    set posy [lindex $rvec 1]
    set posz [lindex $rvec 2]
    for {set j 0} { $j < $p_length } {incr j} {
        part $part_id pos $posx $posy $posz
        incr part_id
        set posz [expr $posz + $b_length]
    }
    return $part_id
}

if { [setmd n_nodes] == 3 || [setmd n_nodes] == 6 } {
    puts "Testcase nve_pe.tcl does not run on 3 or 6 nodes"
    exit 0
}


# System parameters
#############################################################

set n_poly 1
set p_length 10
set b_length 1.00
set d_max 2.24
set d_space [list [expr $d_max] [expr $d_max] 1]
set density 0.00001

# Interaction parameters
#############################################################

set ljr_cut       1.12246204831
set ljr_eps       1.0

set fene_r        2.0
set fene_k        7.0

set bend_k        10.0
set accuracy      1.0e-6

set bjerrum 2.0

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5
set int_steps    1000

# Other parameters
#############################################################
set tcl_precision 10
set mypi          3.141592653589793
set ener_tolerance 0.1

#############################################################
#  Setup System                                             #
#############################################################

set n_part [expr $n_poly * $p_length * (1.) ]
set volume [expr $n_part/$density]
set sphere_rad [expr pow((3.0*$volume)/(4.0*$mypi),1.0/3.0)]
set  box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set shift_vec [list [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0]]

setmd box_l     $box_l $box_l $box_l
setmd periodic  1 1 1
setmd time_step $time_step
setmd skin      $skin
thermostat off

# Interaction setup
#############################################################

# repulsive LJ for all
set ljr_shift  0.25
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0

# FENE 
inter 0 fene $fene_k $fene_r

# Stiffness
inter 1 angle $bend_k

#############################################################
#  Create a bundle of n_poly chains                         #
#############################################################

set part_id 0
set vec {0 0 0}
set rvec [hexconvert $vec $shift_vec $d_space]
set part_id [create_chain $part_id $rvec $p_length $b_length]

# setup bonds etc.
for {set i 0} { $i < [expr $n_poly * $p_length ]} {incr i} {
    part $i type 0 q 1
    if { $i < [expr $p_length / 2.] } {part $i q -1}
    if { [expr $i % $p_length ]  != 0 } {
        part [expr $i-1] bond 0 $i
        if { [expr $i % $p_length]  != 1 } {
            part [expr $i-1] bond 1 [expr $i-2] [expr $i]
        }
    }
}

part 0 fix 1 1 1 
part 1 v 1 1 1

# particle numbers
set n_mono [expr $n_poly * $p_length]
set n_part [expr $n_mono ]

# Coulomb tuning 
#puts -nonewline "P3M Parameter Tuning ... please wait\r"
#flush stdout
#puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy]"

#Use pretuned p3m parameters:
inter coulomb 2 p3m 125.07 8 6 0.0195788 9.47835e-07

# Status report
puts "$n_poly PE chain of length $p_length and charge distance $b_length"
#puts "Constraints:\n[constraint]"
puts "Interactions:\n[inter]\n"
puts "Particles:\n[part]\n"


#############################################################
#      Integration                                          #
#############################################################

puts "Initial Energy: [analyze energy]\n"
set ini_energy [lindex [lindex [analyze energy] 0] 1 ]
integrate $int_steps
set fin_energy [lindex [lindex [analyze energy] 0] 1 ]
puts "Final Energy: [analyze energy]\n"
set error [expr abs( $fin_energy - $ini_energy) / $ini_energy * 100.]
puts "Energy deviation in NVE simulation: $error %"

if { $error > $ener_tolerance } {
    error_exit "energy deviation greater than $ener_tolerance % "
} else {
    puts "Alles in Ordnung :) "
}

exit 0
