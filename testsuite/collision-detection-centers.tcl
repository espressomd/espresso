# Copyright (C) 2011,2012,2013 The ESPResSo project
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
#  Test collision detection with binding of centers of colliding particles
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "COLLISION_DETECTION"

puts "---------------------------------------------------------------"
puts "- Testcase collision-detection-centers.tcl running on [setmd n_nodes] nodes"
puts "---------------------------------------------------------------"

# Setup
setmd box_l 10 10 10

thermostat off
setmd time_step 0.01
inter 3 harmonic 2 2
inter 0 0 lennard-jones 0.0001 2 2.1 auto
inter 7 harmonic 1 1
setmd skin 0
part 0 pos 9   0 0 
# place close to boundary to check pbc and processor boundaries
part 1 pos 0.5 0 0
part 2 pos 3.0 0 0

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
setmd min_global_cut 1.0
on_collision bind_centers 2.0 7

set res [on_collision]
if { ! ( ([lindex $res 0] == "bind_centers") && (abs([lindex $res 1]-2) <1E-5) && ([lindex $res 2] == 7)) } {
    error_exit "Setting collision_detection parameters for bind_centers does not work"
}

# Check the actual collision detection
integrate 0

# Check, whether the bonds are correct
set bonds [analyze_topology 7 1]
if {$bonds != "{0 1}"} {
    error_exit "bond not created as it should: bonds are $bonds"
}

# Integrate again and make sure, no extra bonds are added
integrate 0 recalc_forces

# Check, whether the bonds are still correct, not doubled
set bonds [analyze_topology 7 1]
if {$bonds != "{0 1}"} {
    error_exit "bond doubled: bonds are $bonds"
}

# Exchange particles 1 and 0 in positions and make sure, no extra bonds are added
part 1 pos 9   0 0 
part 0 pos 0.5 0 0

integrate 0

# Check, whether the bonds are still correct, not doubled
set bonds [analyze_topology 7 1]
if {$bonds != "{0 1}"} {
    error_exit "bond double on exchange: bonds are $bonds"
}


# test exception, generating another collision
part 2 pos 7 0 0
on_collision exception bind_centers 2.0 3

if {![catch {integrate 0} err]} {
    error_exit "no exception was thrown at collision, although requested"
}

set bonds ""
foreach exception [lrange $err 2 end] {
    if {[string is integer $exception]} { continue }
    if {[lrange $exception 0 2] != "collision between particles"} {
	error_exit "unexpected exception $exception"
    }
    lappend bonds "[lindex $exception 3] [lindex $exception 5]"
}
set bonds [lsort $bonds]

if {$bonds != "{0 1} {1 2}"} {
    error_exit "exception bonds $bonds wrong"
}

# Check, whether the bonds are correct
# old bonds should not have changed
set bonds [analyze_topology 7]
if {$bonds != "{0 1}"} {
    error_exit "bonds of other type have changed unexpectedly: bonds are $bonds"
}

set bonds [analyze_topology 3]
if {$bonds != "{0 1} {1 2}"} {
    error_exit "bonds not correctly created: bonds are $bonds"
}

exit 0
