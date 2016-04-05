# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
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
require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "BOND_ANGLE"

puts "---------------------------------------------------------------"
puts "- Testcase collision-detection-angular.tcl running on [setmd n_nodes] nodes"
puts "---------------------------------------------------------------"

# Setup
setmd box_l 10 10 10

thermostat off
setmd time_step 0.01
setmd skin 0

# pair bond
inter 2 harmonic 1 1
# angular bonds
for {set a 0} {$a <= 180} {incr a} {
    inter [expr 3 + $a] angle_harmonic 200.0 [expr $a * [PI] / 180]
}

# triangle around x=0
part 0 pos -0.366 2 0 
# place close to boundary to check pbc and processor boundaries
part 1 pos 0.5 2.5 0
part 2 pos 0.5 1.5 0

# reverse triangle
part 3 pos 0.366 2 5 
# place close to boundary to check pbc and processor boundaries
part 4 pos -0.5 2.5 5
part 5 pos -0.5 1.5 5

# line
part 6 pos 7 7 7
part 7 pos 8 7 7
part 8 pos 9 7 7

# Check setting of parameters
setmd min_global_cut 1.0
on_collision bind_three_particles 1.0 2 3 180

set res [on_collision]
if { ! ( ([lindex $res 0] == "bind_three_particles") && (abs([lindex $res 1]-1) <1E-5)
         && ([lindex $res 2] == 2) && ([lindex $res 3] == 3) && ([lindex $res 4] == 180)) } {
    error_exit "Setting collision_detection parameters for bind_centers does not work"
}

# Check the actual collision detection
integrate 0

for {set p 0} {$p < 9} {incr p} {
    puts $p-->[part $p pr bonds]
}

puts "-------------------------------------------------------------------"
# Integrate again and make sure, no extra bonds are added
integrate 0 recalc_forces

for {set p 0} {$p < 9} {incr p} {
    puts $p-->[part $p pr bonds]
}

exit 0
