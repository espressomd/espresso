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
#  Test System: NaCl crystal - Madelung constant            #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "FFTW"
require_max_nodes_per_side 2

puts "----------------------------------------------"
puts "- Testcase madelung.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

#############################################################
#  Parameters                                               #
#############################################################

# The Madelung constant
set madelung_nacl 1.747564594633182190636212035544397403481
set force_epsilon 1e-14

# System identification: 
set name  "madelung"
set ident "_t3"

# System parameters
#############################################################

# grid constant for NaCl crystal
set nacl_const 2.0
# replications of the unit cell in each direcrtion
set nacl_repl  4

set charge 1.0

# for further checks shift the grid in the simulation box
set x0 0.0
set y0 0.0
set z0 0.0

# Interaction parameters coulomb p3m)
#############################################################

set bjerrum   1.0
set accuracy  1.0e-6
# if you do not like the automatic tuning...
set p3m_alpha 0.408905
set p3m_rcut  7.95 
set p3m_mesh  16
set p3m_cao   7

# Integration parameters
#############################################################
setmd time_step 0.01
setmd skin      0.03
thermostat off

set int_steps   1


# Other parameters
#############################################################
set tcl_precision 12

#############################################################
#  Setup System                                             #
#############################################################

set box_l [expr 2.0*$nacl_const*$nacl_repl]
setmd box_l $box_l $box_l $box_l

# Particle setup
#############################################################

set part_id 0
set nc $nacl_const
for {set i 0} { $i < $nacl_repl } {incr i} {
    set xoff [expr $x0+($i * 2.0*$nacl_const)]
    for {set j 0} { $j < $nacl_repl } {incr j} {    
	set yoff [expr $y0+($j * 2.0*$nacl_const)]
	for {set k 0} { $k < $nacl_repl } {incr k} {
	    set zoff [expr $z0+($k * 2.0*$nacl_const)]
	    
	    part $part_id pos [expr $xoff] [expr $yoff] [expr $zoff] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff] [expr $zoff] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff+$nc] [expr $zoff] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff+$nc] [expr $zoff] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff] [expr $zoff+$nc] type 1 q -$charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff] [expr $zoff+$nc] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff] [expr $yoff+$nc] [expr $zoff+$nc] type 0 q $charge
	    incr part_id
	    part $part_id pos [expr $xoff+$nc] [expr $yoff+$nc] [expr $zoff+$nc] type 1 q -$charge
	    incr part_id

	}
    }
}

#puts -nonewline "P3M Parameter Tuning ... please wait\r"
#flush stdout
#puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy mesh 32]"
# Use pretuned p3m parameters:
inter coulomb $bjerrum p3m 7.97000e+00 32 6 4.48505e-01

# print system specification
set n_part [setmd n_part]
puts "NaCl crystal with $n_part particle (grid $nc, replic $nacl_repl)"
#puts "in cubic box of size $box_l. Crystal offset: $x0 $y0 $z0"
#puts "Interactions: \n{ [inter] }"

#writepdb "$name$ident.pdb"

#puts "\nIntegrate $int_steps step(s)"
integrate $int_steps

# calculate Madelung constant
set energy "[analyze energy]"
puts "\nEnergy: $energy"
set calc_mad_const [expr -2.0*$nacl_const*[lindex [lindex $energy 0] 1] / $n_part]
puts "\nExact value of the Madelung constant: $madelung_nacl"
puts "Madelung constant:                    $calc_mad_const"
puts "relative error:                       [expr ($calc_mad_const-$madelung_nacl)/$madelung_nacl]   (<$accuracy)"

# check forces
set max_force 0.0
for {set i 0} { $i < $n_part } {incr i} {
    set force [part $i print f]
    set fx [lindex $force 0]
    if { $fx > $max_force } { set max_force $fx }
    set fy [lindex $force 1]
    if { $fy > $max_force } { set max_force $fy }
    set fz [lindex $force 2]
    if { $fz > $max_force } { set max_force $fz }
}
puts "Maximal force component:              $max_force"

if { [expr ($calc_mad_const-$madelung_nacl)/$madelung_nacl] > $accuracy } {
    error_exit "Madelung constant failed to reach accuracy goal"
}
if { $max_force > $force_epsilon } {
    error_exit "Forces larger than $force_epsilon occured"
}

set cur_pres [analyze pressure]
puts "Total Pressure:           $cur_pres"
set cur_pres [lindex [lindex $cur_pres 0] 1]
set der_pres [expr [lindex [lindex $energy 0] 1]/(3*pow([lindex [setmd box_l] 0],3))]
set err_pres [expr ($der_pres - $cur_pres)/$cur_pres]
puts "Derived Virial Pressure U/3V:        $der_pres   ( relative error: $err_pres )"
if { $err_pres > $accuracy } {
    error_exit "pressure derivation failed to reach accuracy goal"
}

exit 0
