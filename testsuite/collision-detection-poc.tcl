# Copyright (C) 2010,2011 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
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

require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "COLLISION_DETECTION"
require_max_nodes_per_side {1 1 1}


# Setup
setmd box_l 10 10 10

thermostat off
setmd time_step 0.01
inter 0 0 lennard-jones 0.1 1.01 1.01 auto
inter 0 harmonic 1 1
inter 1 harmonic 1 0.0001
setmd skin 0
part 0 pos 0 0 0 
part 1 pos 1 0 0

part 2 pos 3 0 0
integrate 0


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
#integrate 1
#if {! ([part 0 print bond]== "{ } " && [part 1 print bond] == "{ } ") } {
# error_exit "Bonds were created when collision detection was off." 
#}

# Check setting of parameters
on_collision bind_at_point_of_collision 1.0 0 1 1
set res [on_collision]
if { ! ( ([lindex $res 0] == "bind_at_point_of_collision") && (abs([lindex $res 1]-1) <1E-5) && ([lindex $res 2] == 0) && ([lindex $res 3] == 1) && ([lindex $res 4] == 1)) } {
  error_exit "Setting collision_detection parameters for bind_centers does not work"
}

# Check the actual collision detection
integrate 1 

# Get bond info
set bond1 [part 0 print bond]
set bond2 [part 1 print bond]
set bond3 [part 2 print bond]


# Check that the single particle did not get a bond
if {$bond3 != "{ } "} {
 error_exit "3rd particle should not get a bond."
}

# Check, whether the bonds are correct
if {!((($bond1=="{ {0 1} } ") && ($bond2=="{ } ")) || (($bond2=="{ {0 0} } ") && ($bond1=="{ } "))) } {
 error_exit "Bond between first 2 particles incorrect: $bond1, $bond2."
}

# Integrate again and make sure, no extra bonds are added
integrate 1


# Check, whether the bonds are correct
if {!((($bond1=="{ {0 1} } ") && ($bond2=="{ } ")) || (($bond2=="{ {0 0} } ") && ($bond1=="{ } "))) } {
 error_exit "Bond between first 2 particles incorrect."
}

# Check whether two virtual sites have been created
if {[setmd n_part] !=5} {
 error_exit "Incorrect number of particles in the simulation. Too many or too few virtual sites were created."
}

# Check the bonds between virtual sites
set bond1 [part 3 print bonds]
set bond2 [part 4 print bonds]

if {!((($bond1=="{ {1 4} } ") && ($bond2=="{ } ")) || (($bond2=="{ {1 3} } ") && ($bond1=="{ } "))) } { 
 error_exit "Bonds between the virtual sites are incorrect."
}

# Check the particle type of the virtual sites
if { (! ([part 3 print type]==1 && [part 4 print type]==1))} {
 error_exit "
}
