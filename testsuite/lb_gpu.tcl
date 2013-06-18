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
require_feature "ADRESS" off

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

# Integration parameters
#############################################################
set int_steps     100
set int_times     20

set time_step     0.005
set tau           0.02

set agrid         1.0

set box_l         30.0

set dens          0.85
set viscosity     30.0
set friction      2.0

set temp          1.0

set skin          0.5

set mom_prec      1.e-2
set mass_prec     1.e-8
set temp_confidence 10

# Other parameters
#############################################################

if { [ catch {
#############################################################
# System Setup                                              #
#############################################################

# make real random draw
set cmd "t_random seed"
for {set i 0} {$i < [setmd n_nodes]} { incr i } {
    lappend cmd [expr [pid] + $i] }
eval $cmd

setmd time_step $time_step
setmd skin $skin

# Simulation box
#############################################################
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1

# Particles
#############################################################
# load colloid from file
read_data "lb_system.data"
thermostat langevin 1. 1.
integrate 1000
kill_particle_motion
thermostat off
#part 0 pos 10 10 10
# here you can create the necessary snapshot
#write_data "lb_system.data"

# Fluid
#############################################################

set components [setmd lb_components]
if {$components==1} {
  lbfluid gpu agrid $agrid dens $dens visc $viscosity agrid $agrid tau $tau
  lbfluid friction $friction
}
if {$components==2} {
  lbfluid gpu agrid $agrid dens $dens $dens visc $viscosity $viscosity agrid $agrid tau $tau sc_coupling 0.0 0.0 0.0 
  lbfluid friction $friction $friction
  exit 0
}



thermostat lb $temp

# give the colloid a kick
for { set i 0 } { $i < [setmd n_part] } { incr i } { 
    set vx [lindex [part $i print v] 0]
    set vy [lindex [part $i print v] 1]
    set vz [lindex [part $i print v] 2]
    part $i v [expr .1+$vx] $vy $vz
}

# determine initial fluid mass and total momentum (fluid is at rest)
set fluidmass [expr $dens*pow($box_l,3)]
set tot_mom { 0.0 0.0 0.0 }
for { set i 0 } { $i < [setmd n_part] } { incr i } {
    lset tot_mom 0 [expr [lindex $tot_mom 0]+[lindex [part $i print v] 0]]
    lset tot_mom 1 [expr [lindex $tot_mom 1]+[lindex [part $i print v] 1]]
    lset tot_mom 2 [expr [lindex $tot_mom 2]+[lindex [part $i print v] 2]]
}

## warm up particle and fluid
integrate 1000

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
set avg_fluid_temp  0.0
set var_fluid_temp  0.0

#puts "fluid: [analyze fluid momentum]"
#puts "parts: [analyze momentum particles]"

for { set i 1 } { $i <= $int_times } { incr i } {

    puts -nonewline "Loop $i of $int_times starting at time [format %f [setmd time]]\r"; flush stdout
    integrate $int_steps

    # check fluid mass conservation
    #set dmass [expr abs([analyze fluid mass]-$fluidmass)]
    #if { $dmass > $mass_prec } {
	#error "mass deviation too large $dmass"
    #}
    #if { $dmass > $max_dmass } { set max_dmass $dmass }

    # check total momentum conservation
    set p_mom [analyze momentum particles]
    set f_mom [analyze fluid momentum]
 
    set momx [expr [lindex $p_mom 0]+[lindex $f_mom 0]]
    set momy [expr [lindex $p_mom 1]+[lindex $f_mom 1]]
    set momz [expr [lindex $p_mom 2]+[lindex $f_mom 2]]

    #puts "fluid: [analyze fluid momentum]"
    #puts "parts: [analyze momentum particles]"
    set dmx [expr abs( $momx-[lindex $tot_mom 0])]
    set dmy [expr abs( $momy-[lindex $tot_mom 1])]
    set dmz [expr abs( $momz-[lindex $tot_mom 2])]
    if { $dmx > $mom_prec || $dmy > $mom_prec || $dmz > $mom_prec } {
	error "momentum deviation too large $p_mom $f_mom $dmx $dmy $dmz"
    }
    if { $dmx > $max_dmx } { set max_dmx $dmx }
    if { $dmy > $max_dmy } { set max_dmy $dmy }
    if { $dmz > $max_dmz } { set max_dmz $dmz }

    # temperature of the colloid
    set temp [expr 2.0/[degrees_of_freedom]*[analyze energy kinetic]/[setmd n_part]]
    set avg_temp [expr $avg_temp+$temp]
    set var_temp [expr $var_temp+$temp*$temp]

    # temperature of the fluid
    set fluid_temp [analyze fluid temp]

    set avg_fluid_temp [expr $avg_fluid_temp+$fluid_temp]
    set var_fluid_temp [expr $var_fluid_temp+$fluid_temp*$fluid_temp]
}

#############################################################
# Analysis and Verification                                 #
#############################################################
set tcl_precision 6

puts "NOTE: this is a statistical test, which can fail,"
puts "even if everything works correctly. However, the"
puts "chance is really SMALL, so if you see an error"
puts "here, consider a bug in the thermostat.\n"

set avg_temp [expr $avg_temp/$int_times]
set avg_fluid_temp [expr $avg_fluid_temp/$int_times]

set temp_dev [expr sqrt(2.0/([setmd n_part]*[degrees_of_freedom]))]
set temp_prec [expr $temp_confidence*$temp_dev/sqrt($int_times)]

puts "\n"
puts "maximal mass deviation $max_dmass"
puts "maximal momentum deviation in x $max_dmx, in y $max_dmy, in z $max_dmz\n"

puts "average temperature         $avg_temp"
puts "fluid temperature           $avg_fluid_temp"
puts "maximally accepted deviations are $temp_prec"

if {[expr abs($avg_temp - [setmd temp])] > $temp_prec} {
    error "relative particle temperature deviation too large"
}
if { [expr abs($fluid_temp - [setmd temp])] > $temp_prec} {
    error "relative fluid temperature deviation too large"
}

} res ] } {
    error_exit $res
}

exit 0

#############################################################
