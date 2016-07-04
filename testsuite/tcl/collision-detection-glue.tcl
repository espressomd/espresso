# Copyright (C) 2011,2012,2016 The ESPResSo project
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
#  Test collision detection in the "glue to surface>" mode
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "COLLISION_DETECTION"
require_feature "ADRESS" off
require_max_nodes_per_side {1 1 1}

puts "---------------------------------------------------------------"
puts "- Testcase collision-detection-glue.tcl running on 1 nodes"
puts "---------------------------------------------------------------"

# Setup
setmd box_l 19 19 19

thermostat off
setmd time_step 0.01
inter 2 harmonic 1 1
inter 3 harmonic 1 0.0001
setmd skin 0
part 0 pos 0 0 0 type 2 
part 1 pos 5.5 0 0 type 3
part 2 pos 3 0 0 type 0
inter 2 3 lennard-jones 0.001 5.6 5.7 auto
# analyze the bonding structure for pair bonds
proc analyze_topology {bond_type {check_others 0}} {
    set bonded ""
    for {set i 0} {$i <= [setmd max_part]} {incr i} {
	foreach bond [lindex [part $i pr bond] 0] {
	    if {[lindex $bond 0] == $bond_type} {
		set j [lindex $bond 1]
		if {$i < $j} {
		    lappend bonded "$i $j"
		} {
		    lappend bonded "$j $i"
		}
	    } {
		if {$check_others} {
		    error_exit "bond $bond at particle $i of unexpected type found"
		}
	    }
	}
    }
    return [lsort $bonded]
}

# Test default setting
if { "[on_collision]" != "off" } {
    error_exit "Collision detection should be off by default."
}

# Test switching it off
on_collision off
if { "[on_collision]" != "off" } {
    error_exit "Disabling collision_detection does not work"
}

# Make sure, it doesn't do anything when turned off
integrate 0

set bonds [analyze_topology "" 1]
if {$bonds != ""} {
    error_exit "Bonds were created when collision detection was off." 
}

# Check setting of parameters
setmd min_global_cut 5.0
on_collision glue_to_surface 5.5 2 3 1 2 3 4 1.0

set res [on_collision]

if { ! ($res == "glue_to_surface 5.500000 2 3 1 2 3 4 1.000000") } {
    error_exit "Setting collision_detection parameters does not work: $res"
}

# Check the actual collision detection
integrate 0 recalc_forces

# Check bonds between colliding particles
set bonds [analyze_topology 2]
if {$bonds != "{0 1}"} {
    error_exit "bond between colliding particles not created as it should: bonds are $bonds"
}

# Check bond between the glued particle and the virtual site
set bonds [analyze_topology 3]
if {$bonds != "{0 3}"} {
    error_exit "bond between glued particle and vs not created as it should: bonds are $bonds"
}

# Check if the virtual site has the correct settings
if { [part 3 print virtual] != 1 } {
 error_exit "The supposed virtual particle doesn't have the virtual flag set."
}
# Is the vs attached correctly
set vs_info [lrange [part 3 print vs_relative] 0 1]
if {$vs_info != "1 4.500000" } {
 error_exit "the vs_relative params are wrong: $vs_info"
}

# Integrate again and make sure, no extra bonds are added
# enforce force recalculation
integrate 0 recalc_forces

# Check, whether the bonds are still correct, not doubled
set bonds [analyze_topology 2]
if {$bonds != "{0 1}"} {
    error_exit "After 2nd run: bond between colliding particles not created as it should: bonds are $bonds"
}

# Check bond between the glued particle and the virtual site
set bonds [analyze_topology 3]
if {$bonds != "{0 3}"} {
    error_exit "After 2nd run: bond between glued particle and vs not created as it should: bonds are $bonds"
}



# Check whether the number of particles is correct (3 normal +1 vs =4)
if {[setmd n_part] != 4} {
    error_exit "Incorrect number of particles [setmd n_part] in the simulation. Too many or too few virtual sites were created."
}

# Check the particle type of the virtual sites
if { ! ([part 3 print type]==1) } {
    error_exit "type of virtual sites is incorrect."
}
if { ! ([part 0 print type]==4) } {
    error_exit "type of glued particle is incorrect: [part 0 print type]."
}

# test exception, generating another collision
part 2 pos 11 0 0 type 2
on_collision exception glue_to_surface 5.5 2 3 1 2 3 4 1.0
if {![catch {integrate 1} err]} {
    error_exit "no exception was thrown at collision, although requested"
}

set bonds ""
foreach exception [lrange $err 1 end] {
    if {[lrange $exception 1 3] != "collision between particles"} {
	error_exit "unexpected exception $exception"
    }
    lappend bonds "[lindex $exception 4] [lindex $exception 6]"
}
set bonds [lsort $bonds]

if {$bonds != "{1 2}"} {
    error_exit "exception bonds $bonds wrong"
}

# Check, whether the bonds are also correct
# Between centers of colliding particles
set bonds [analyze_topology 2]
if {$bonds != "{0 1} {1 2}"} {
    error_exit "bonds not correctly created: bonds are $bonds"
}
# Between glued particle and vs
set bonds [analyze_topology 3]
if {$bonds != "{0 3} {2 4}"} {
    error_exit "bonds not correctly created: bonds are $bonds"
}

exit 0
