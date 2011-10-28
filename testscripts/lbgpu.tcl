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
#
#############################################################
#                                                           #
# Basic tests of the Lattice Boltzmann implementation       #
#                                                           #
# 1) check conservation of fluid mass                       #
# 2) check conservation of total momentum                   #
# 3) measure temperature of colloid (no auto-check so far)  #
#                                                           #
#############################################################
require_feature "LB_GPU"
require_feature "LENNARD_JONES"

puts "----------------------------------------"
puts "- Testcase lb.tcl running on [format %02d [setmd n_nodes]] nodes  -"
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
# Integration parameters
#############################################################
set int_steps     10
set int_times     1

set time_step     0.02
set tau           0.02

set agrid         1.0

set box_l         30.0

set dens          0.85
set viscosity     3.0
set friction      1.0

set temp          1.0

set skin          0.5

set mom_prec      1.e-2
set mass_prec     1.e-8
set temp_prec     0.7

# Other parameters
#############################################################

#############################################################
# System Setup                                              #
#############################################################
puts "cudadevices: [cuda list]"
setmd time_step $time_step
setmd skin $skin

# Simulation box
#############################################################
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1
cellsystem domain_decomposition -no_verlet_list

# Fluid
#############################################################
lbfluid gpu dens $dens visc $viscosity agrid $agrid tau $tau
lbfluid friction $friction

thermostat lb $temp

# Particles
#############################################################
# load colloid from file
#read_data "lb_system.data"
# here you can create the necessary snapshot
#write_data "lb_system.data"
#part 0 pos 5 5 5 
# give the colloid a kick
#for { set i 0 } { $i < [setmd n_part] } { incr i } { 
#    set vx [lindex [part $i print v] 0]
#    set vy [lindex [part $i print v] 1]
#    set vz [lindex [part $i print v] 2]
#    part $i v $vx [expr 1.0+$vy] $vz
#}

# determine initial fluid mass and total momentum (fluid is at rest)
set fluidmass [expr $dens*pow($box_l,3)]
puts "#particles: [setmd n_part]\n"
#puts "tot: $tot_mom"
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

for { set i 1 } { $i <= $int_times } { incr i } {

    puts -nonewline "Loop $i of $int_times starting at time [format %f [setmd time]]\n"
    integrate $int_steps


    # temperature of the fluid
    set fluid_temp [analyze fluid temp]

}    

#############################################################
# Analysis and Verification                                 #
#############################################################

puts "fluid temperature [analyze fluid temp] (relative deviation [expr $fluid_temp-1.0])\n"
#if { $rel_temp_error > $temp_prec } {
#    error "relative temperature deviation too large"
#}


exit 0

#############################################################
