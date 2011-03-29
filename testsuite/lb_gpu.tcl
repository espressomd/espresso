# Copyright (C) 2010,2011 The ESPResSo project
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

#############################################################
#                                                           #
# Basic tests of the Lattice Boltzmann implementation       #
#                                                           #
# 1) check conservation of fluid mass                       #
# 2) check conservation of total momentum                   #
# 3) measure temperature of colloid (no auto-check so far)  #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "LB_GPU"
require_feature "LENNARD_JONES"

puts "----------------------------------------"
puts "- Testcase lb_gpu.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

#############################################################
# Procedures                                                #
#############################################################

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} {blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write interactions
    blockfile $f write particles {id type q pos v f}
    blockfile $f write bonds
    close $f
}

#############################################################
# Parameters                                                #
#############################################################
set errf [lindex $argv 1]


require_feature "LB_GPU"
#require_feature "LB"

# Integration parameters
#############################################################
set int_steps     	1
puts $int_steps
set int_times		1000
puts $int_times

set time_step     0.01
set tau           0.01

set agrid         1.0

set box_l         64
set box_l_x       128.0

set dens          1.0
set viscosity     3.0
#set bulk_viscosity 0.1
set friction      5.0
#equ velo 0.025
set temp          0.0

set skin          0.5

#set mom_prec      1.e-3
#set mass_prec     1.e-3
#set temp_prec     5.e-3

# Other parameters
#############################################################

if { [ catch {
#############################################################
# System Setup                                              #
#############################################################

setmd time_step $time_step
setmd skin $skin
puts [setmd time_step]
# Simulation boxlo
#############################################################
setmd box_l $box_l $box_l $box_l
puts [setmd box_l]
setmd periodic 1 1 1
#cellsystem domain_decomposition -no_verlet_list
#cellsystem nsquare

#inter 0 0 lennard-jones 1 1 12.5 0 0
#puts [inter]

puts [setmd cell_grid]
#part 0 pos 10 10 10
# Fluid
#############################################################
lbfluid gpu dens $dens visc $viscosity agrid $agrid tau $tau friction $friction ext_force 0.5 0 0
#bulk_viscosity $bulk_viscosity ext_force 0.03 0 0
thermostat lb $temp

#lbnode_exf_gpu 0 0 0 set force 10 0 0
#	for { set i 0 } { $i < 32 } { incr i } {
#		for { set j 0 } { $j < 32 } { incr j } {
#		lbnode_exf_gpu 0 $i $j set force 1 0 0
#		}
#	}
lbboundary wall normal 0 0 1 dist 0
lbboundary wall normal 0 0 -1 dist [expr -$box_l+1]

lbboundary wall normal 0 1 0 dist 0
lbboundary wall normal 0 -1 0 dist [expr -$box_l+1]
# Particles
#############################################################
# load colloid from file
#read_data "~/espresso/testsuite/lb_system.data"

#part 0 pos 10 10 10
#part 0 v 0. 0. -10.5
#set k	0
#set m	5
#for { set m 0 } { $m < 1 } { incr m } {
#	for { set i 0 } { $i < 64 } { incr i } {
#		for { set j 0 } { $j < 64 } { incr j } {
#			
		 	#if {sqrt(pow($j-15., 2) + pow($i-15.,2)) < 4.} { continue }
#			part $k pos $m. $i. $j.
#			incr k
#		}
#	}
#}
#part deleteall
#part 0 pos 0 5 15
#puts [setmd n_part]
# here you can create the necessary snapshot
#write_data "lb_system.data"

# give the colloid a kick
#for { set i 0 } { $i < [setmd n_part] } { incr i } { 
#	part $i fix 
#        part $i v 0. 0. 0.
#puts "fixing part [ part 0 print fix ]"
 #   set vx [lindex [part $i print v] 0]
  #  set vy [lindex [part $i print v] 1]
   # set vz [lindex [part $i print v] 2]
    #part $i v [expr 1.0+$vx] $vy $vz
#}

# determine initial fluid mass and total momentum (fluid is at rest)
set fluidmass [expr $dens*pow($box_l,3)]
set tot_mom { 0.0 0.0 0.0 }
#for { set i 0 } { $i < [setmd n_part] } { incr i } {
#    lset tot_mom 0 [expr [lindex $tot_mom 0]+[lindex [part $i print v] 0]]
#    lset tot_mom 1 [expr [lindex $tot_mom 1]+[lindex [part $i print v] 1]]
#    lset tot_mom 2 [expr [lindex $tot_mom 2]+[lindex [part $i print v] 2]]
#}

set max_dmass 0.0
set max_dmx   0.0
set max_dmy   0.0
set max_dmz   0.0

#############################################################
# Integration                                               #
#############################################################
puts "Running at temperature T=[setmd temp] [thermostat]"

set max_dmass 0.0
set max_dmx   0.0
set max_dmy   0.0
set max_dmz   0.0

set avg_temp  0.0
set var_temp  0.0

#puts "fluid: [analyze fluid momentum]"
#puts "parts: [analyze momentum particles]"
#prepare_vmd_connection "vmd_data"
#after 3000
puts [time {
for { set i 1 } { $i <= $int_times } { incr i } {

    #puts -nonewline "Loop $i of $int_times starting at time [format %f [setmd time]]\n"; flush stdout
#\n
#puts [part 0 print pos v force]
	#for { set k 0 } { $k < $box_l } { incr k } {
	#	for { set j 0 } { $j < $box_l } { incr j } {
	#		lbnode_exf_gpu 5 $k $j set force 0.01 0 0
	#	}
	#}
	#lbprint field lb_field/datei$i.vtk
    integrate $int_steps
    #puts [part 0 print pos]
    #puts -nonewline [analyze energy]

#imd positions
#after 50
    # check fluid mass conservation
    #set dmass [expr abs([analyze fluid mass]-$fluidmass)]
    #if { $dmass > $mass_prec } {
#	error "mass deviation too large $dmass"
 #   }
  #  if { $dmass > $max_dmass } { set max_dmass $dmass }

    # check total momentum conservation
    #set mom [analyze momentum]
    #puts "fluid: [analyze fluid momentum]"
    #puts "parts: [analyze momentum particles]"
   # set dmx [expr abs([lindex $mom 0]-[lindex $tot_mom 0])]
    #set dmy [expr abs([lindex $mom 1]-[lindex $tot_mom 1])]
    #set dmz [expr abs([lindex $mom 2]-[lindex $tot_mom 2])]
    #if { $dmx > $mom_prec || $dmy > $mom_prec || $dmz > $mom_prec } {
#	error "momentum deviation too large $mom $tot_mom $dmx $dmy $dmz"
 #   }
  #  if { $dmx > $max_dmx } { set max_dmx $dmx }
   # if { $dmy > $max_dmy } { set max_dmy $dmy }
    #if { $dmz > $max_dmz } { set max_dmz $dmz }

    # temperature of the colloid
    #set temp [expr 2.0/3.0*[analyze energy kinetic]/[setmd n_part]]
    #set avg_temp [expr $avg_temp+$temp]
    #set var_temp [expr $var_temp+$temp*$temp]

    # temperature of the fluid
    #set fluid_temp [analyze fluid temp]

}    
puts "\n"
}]
#############################################################
# Analysis and Verification                                 #
#############################################################
#set avg_temp [expr $avg_temp/$int_times]
#set var_temp [expr $var_temp/$int_times - $avg_temp*$avg_temp]
#set rel_temp_error [expr abs(($avg_temp-[setmd temp])/[setmd temp])]
lbprint v vtk field.vtk
puts "\n"
#puts "Maximal mass deviation $max_dmass"
#puts "Maximal momentum deviation in x $max_dmx, in y $max_dmy, in z $max_dmz"

#puts "\nAverage temperature $avg_temp (relative deviation $rel_temp_error)\n"
#if { $rel_temp_error > $temp_prec } {
#    error "relative temperature deviation too large"
#}

} res ] } {
    error_exit $res
}

exec rm -f $errf
#after 10000
exit 0

#############################################################
