#!/bin/sh
#
# tricking... the line after this comment is interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
#
# $Id: lb.tcl,v 2.3 2007-11-26 15:03:50 stuehn Exp $
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
# Copyright (C) 2016 The ESPResSo project
# Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

source "tests_common.tcl"

require_feature "SD"
require_max_nodes_per_side 1
require_cudadevice
require_feature "SD_NOT_PERIODIC" "off"
#require_feature "BD"


#############################################################
# Parameters                                                #
#############################################################
set int_steps     1
set int_times  	  13
set warmup        100

set time_step     0.01

set skin          0.5
set charge        1
set F		  0.0
set yukawa_energy_cut     0.1

set rc 1.5

# Other parameters
set gamma 500
set start 1
set thermo "sd"
set kappa_t 3.0
set temp 1
set pdens 0.1
set box_l 20.0
set seed [pid]
set radius 0.3
set viscosity 1
set box_l_y       $box_l
set box_l_z       $box_l
set box_l_x       $box_l
set N  [expr int($pdens*$box_l_x*$box_l_y*$box_l_z)]
set wz [expr pow(3./(4*3.1415926535897931*$pdens),(1./3))]
set lambda_t [expr $wz/$kappa_t]
if {$temp > "0"} {set bjerrum [expr $gamma*$wz]
   set cutoff [expr -log($yukawa_energy_cut/$bjerrum)*$lambda_t]
}
set temp_steps 0
puts "bjerrum: $bjerrum\n"
puts "lambda: $lambda_t\n"
puts "cutoff: $cutoff\n"
puts "number of parts: $N\n"
puts "temperatur: $temp\n"
puts "number of temp steps: $temp_steps\n"
#############################################################
#                                                           #
#############################################################

#############################################################
# System Setup                                              #
#############################################################

setmd time_step [expr $time_step/100]
setmd skin $skin
t_random seed $seed
# Simulation boxlo
#############################################################
setmd box_l $box_l_x $box_l_y $box_l_z
setmd periodic 1 1 1

#############################################################
#                                                           #
#############################################################
# Fluid
#############################################################

#cuda setdevice $dev


setmd sd_radius $radius

setmd sd_visc $viscosity

#############################################################
# Particles
#############################################################


for { set j 0 } { $j < $N } { incr j } {
	part $j pos [expr [t_random]*$box_l_x] [expr [t_random]*$box_l_y] [expr [t_random]*$box_l_z] type 0 q $charge
}



#############################################################
# interaction                                               #
#############################################################
if { 1 } { 
	if { $start == 1} {
		set lambda 0.1
		
		while { $lambda < $lambda_t } {
		    inter coulomb $bjerrum dh [expr 1.0/$lambda] $cutoff
		    integrate 7
		    set lambda [expr $lambda*1.01]
		}
	    }
	    setmd min_global_cut [expr $rc+0.25]
	    inter coulomb $bjerrum dh [expr 1.0/$lambda_t] $cutoff
    }

setmd skin 0.08
#setmd cell_grid 4 4 4
#puts "Tuning cells .."
#tune_cells
#puts [setmd skin]
#puts [setmd cell_grid]

################warmup######################################
setmd time_step 0.01
thermostat bd $temp


set logname "sd_thermalization.log"
set log [open $logname w]

for {set i 0} {$i < $warmup} {incr i} {
   integrate_sd [expr $int_steps*100]
   set tmp [analyze energy coulomb]
   puts $log $tmp
}


#############################################################
# start measurement for bd                                  #
#############################################################


set energy_bd 0
for {set i 0} {$i < $int_times} {incr i} {
   integrate_sd $int_steps
   set energy_bd [expr $energy_bd + $tmp ]
}
flush $log
puts "Brownian Dynamics part finished, starting Stokesian Dynamics ..."
thermostat off
thermostat sd $temp
for {set i 0} {$i < $int_times} {incr i} {
   integrate_sd $int_steps
   set tmp [analyze energy coulomb]
   puts $log $tmp
}
flush $log
set energy_sd 0
for {set i 0} {$i < $int_times} {incr i} {
   integrate_sd $int_steps
   set tmp [analyze energy coulomb]
   puts $log $tmp
   set energy_sd [expr $energy_sd + $tmp ]
}

set energy_bd [expr $energy_bd / $int_times]
set energy_sd [expr $energy_sd / $int_times]
#puts $energy_sd
#puts $energy_bd

set error [expr abs(($energy_sd - $energy_bd) / ($energy_sd + $energy_bd )*2)]
set prec 0.01
if {$error > $prec} {
    puts "Error was $error. This is bigger than $prec, which is the limit. It can happen sometimes that this fails as this is a stochastic test.
If it happens often, there is probably something wrong."
    puts "You can also increase the \$int_times variable to increase the accuracy of the test, or check the log file $logname for trends."
    error_exit
}

puts "Success. The error was $error."
flush $log
close $log


ok_exit

#############################################################
