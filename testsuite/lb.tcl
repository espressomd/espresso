#!/bin/sh
#
# tricking... the line after this comment is interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
#
# $Id$
#
# This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
# It is therefore subject to the ESPResSo license agreement which you
# accepted upon receiving the distribution and by which you are
# legally bound while utilizing this file in any form or way. 
# There is NO WARRANTY, not even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# You should have received a copy of that license along with this
# program; if not, refer to http://www.espresso.mpg.de/license.html
# where its current version can be found, or write to
# Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 
# 55021 Mainz, Germany.
# Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#############################################################
#                                                           #
# Basic tests of the Lattice Boltzmann implementation       #
#                                                           #
# 1) check conservation of fluid mass                       #
# 2) check conservation of total momentum                   #
# 3) measure temperature of colloid (no auto-check so far)  #
#                                                           #
#############################################################

puts "----------------------------------------"
puts "- Testcase lb.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

#############################################################
# Procedures                                                #
#############################################################
proc error_exit {error} {
    global errf
    #set f [open $errf "w"]
    puts "Error occured: $error"
    #close $f
    exit -666
}

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

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	set f [open $errf "w"]
	puts $f "not compiled in: $feature"
	close $f
	exit -42
    }
}

#############################################################
# Parameters                                                #
#############################################################
set errf [lindex $argv 1]

require_feature "LB"

# Integration parameters
#############################################################
set int_steps     10
set int_times     1000

set time_step     0.005
set tau           0.02

set agrid         1.0

set box_l         30.0

set dens          0.85
set viscosity     3.0
set friction      20.0

set temp          1.0

set skin          0.5

set mom_prec      1.e-8
set mass_prec     1.e-12
set temp_prec     5.e-3

# Other parameters
#############################################################

if { [ catch {
#############################################################
# System Setup                                              #
#############################################################

setmd time_step $time_step
setmd skin $skin

# Simulation box
#############################################################
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1
cellsystem domain_decomposition -no_verlet_list

# Fluid
#############################################################
lbfluid dens $dens visc $viscosity agrid $agrid tau $tau
lbfluid friction $friction

setmd temperature $temp

# Particles
#############################################################
# load colloid from file
read_data "lb_system.data"

# here you can create the necessary snapshot
#write_data "lb_system.data"

# give the colloid a kick
for { set i 0 } { $i < [setmd n_part] } { incr i } { 
    set vx [lindex [part $i print v] 0]
    set vy [lindex [part $i print v] 1]
    set vz [lindex [part $i print v] 2]
    part $i v [expr 1.0+$vx] $vy $vz
}

# determine initial fluid mass and total momentum (fluid is at rest)
set fluidmass [expr $dens*pow($box_l,3)]
set tot_mom { 0.0 0.0 0.0 }
for { set i 0 } { $i < [setmd n_part] } { incr i } {
    lset tot_mom 0 [expr [lindex $tot_mom 0]+[lindex [part $i print v] 0]]
    lset tot_mom 1 [expr [lindex $tot_mom 1]+[lindex [part $i print v] 1]]
    lset tot_mom 2 [expr [lindex $tot_mom 2]+[lindex [part $i print v] 2]]
}

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

    puts -nonewline "Loop $i of $int_times starting at time [format %f [setmd time]]\r"; flush stdout
    integrate $int_steps

    # check fluid mass conservation
    set dmass [expr abs([analyze fluid mass]-$fluidmass)]
    if { $dmass > $mass_prec } {
	error "mass deviation too large $dmass"
    }
    if { $dmass > $max_dmass } { set max_dmass $dmass }

    # check total momentum conservation
    set mom [analyze momentum]
    #puts "fluid: [analyze fluid momentum]"
    #puts "parts: [analyze momentum particles]"
    set dmx [expr abs([lindex $mom 0]-[lindex $tot_mom 0])]
    set dmy [expr abs([lindex $mom 1]-[lindex $tot_mom 1])]
    set dmz [expr abs([lindex $mom 2]-[lindex $tot_mom 2])]
    if { $dmx > $mom_prec || $dmy > $mom_prec || $dmz > $mom_prec } {
	error "momentum deviation too large $mom $tot_mom $dmx $dmy $dmz"
    }
    if { $dmx > $max_dmx } { set max_dmx $dmx }
    if { $dmy > $max_dmy } { set max_dmy $dmy }
    if { $dmz > $max_dmz } { set max_dmz $dmz }

    # temperature of the colloid
    set temp [expr 2.0/3.0*[analyze energy kinetic]/[setmd n_part]]
    set avg_temp [expr $avg_temp+$temp]
    set var_temp [expr $var_temp+$temp*$temp]

    # temperature of the fluid
    set fluid_temp [analyze fluid temp]

}    

#############################################################
# Analysis and Verification                                 #
#############################################################
set avg_temp [expr $avg_temp/$int_times]
set var_temp [expr $var_temp/$int_times - $avg_temp*$avg_temp]
set rel_temp_error [expr abs(($avg_temp-[setmd temp])/[setmd temp])]

puts "\n"
puts "Maximal mass deviation $max_dmass"
puts "Maximal momentum deviation in x $max_dmx, in y $max_dmy, in z $max_dmz"

puts "\nAverage temperature $avg_temp (relative deviation $rel_temp_error)\n"
#if { $rel_temp_error > $temp_prec } {
#    error "relative temperature deviation too large"
#}

} res ] } {
    error_exit $res
}

exec rm -f $errf

exit 0

#############################################################
