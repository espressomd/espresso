# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
#  
# this file is part of espresso.
#  
# espresso is free software: you can redistribute it and/or modify
# it under the terms of the gnu general public license as published by
# the free software foundation, either version 3 of the license, or
# (at your option) any later version.
#  
# espresso is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  see the
# gnu general public license for more details.
#  
# you should have received a copy of the gnu general public license
# along with this program.  if not, see <http://www.gnu.org/licenses/>. 

# 
#############################################################
#                                                           #
#  test collision detection with binding of centers of colliding particles
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "COLLISION_DETECTION"
if {[has_feature "LEES_EDWARDS"]} {
    require_max_nodes_per_side 2
}

puts "---------------------------------------------------------------"
puts "- Testcase collision-detection-centers.tcl running on [setmd n_nodes] nodes"
puts "---------------------------------------------------------------"

# setup
setmd box_l 10 10 10

thermostat off
setmd time_step 0.01
inter 3 harmonic 2 2
inter 0 0 lennard-jones 0.0001 2 2.1 auto
inter 7 harmonic 2 1
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

# test default setting
if { "[on_collision]" != "off" } {
    error_exit "collision detection should be off by default."
}

# test switching it off
on_collision off
if { "[on_collision]" != "off" } {
    error_exit "disabling collision_detection does not work"
}

# make sure, it doesn't do anything when turned off
integrate 0

set bonds [analyze_topology "" 1]
if {$bonds != ""} {
    error_exit "bonds were created when collision detection was off." 
}

# check setting of parameters
setmd min_global_cut 1.0
on_collision bind_centers 2.0  7

set res [on_collision]
if { ! ( ([lindex $res 0] == "bind_centers") && (abs([lindex $res 1]-2) <1e-5) && ([lindex $res 2] == 7)) } {
    error_exit "setting collision_detection parameters for bind_centers does not work"
}

# check the actual collision detection
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

# check, whether the bonds are still correct, not doubled
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
foreach exception [lrange $err 1 end] {
    if {[regexp { ERROR: collision between particles (\d+) and (\d+)} $exception -> id1 id2]} {
        lappend bonds "$id1 $id2"
    } else {
	error_exit "unexpected exception $exception"
    }
    
}
set bonds [lsort $bonds]

if {$bonds != "{0 1} {1 2}"} {
    error_exit "exception bonds $bonds wrong, expected {0 1} {1 2}"
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
