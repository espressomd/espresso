# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
source "tests_common.tcl"

require_feature ELECTROSTATICS 
require_feature MASS
require_feature EXTERNAL_FORCES

proc veccompare { a b } {
#  puts "$a $b"
  if { [ llength $a ]  != [ llength $b ] } {
    return 0
  }
  for { set i 0 } { $i < [ llength $a ] } { incr i } {
    set A [ lindex $a $i ] 
    set B [ lindex $b $i ] 
    if { !( $A == $B || [ expr abs($A - $B) < 1e-6 ]) } {
      return 0
    }
  }
  return 1
}


puts "---------------------------------------------------------------"
puts "- Testcase observable.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 8 8 8
setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 1 0 0 v 2 0 0 type 0 q 1 mass 1 ext_force 1 0 0 
part 1 pos 2 0 0 v 4 0 0 type 0 q -1 ext_force -2 0 0
part 2 pos 1 1 0 v 2 2 0 type 1
part 3 pos 2 1 0 v 4 2 0 type 1
part 4 pos 3 0 0 v 6 0 0 type 0 mass 2
integrate 0


############### Observable particle_position ##########################
## Here we are checking if particle specifications are working
# position particle 0
set p0 [ observable new particle_positions id 0 ]
if { ![ veccompare [ observable $p0 print ] { 1 0 0 } ] }  {
  error "P0 is not working"
}
# position particle 0, 3
set p1 [ observable new particle_positions id { 2 3 } ]
if { ![ veccompare [ observable $p1 print ] { 1 1 0 2 1 0} ] }  {
  error "P1 is not working"
}
# positions particle 0 1
set p2 [ observable new particle_positions type { 0 } ]
if { ![ veccompare [ observable $p2 print ] { 1 0 0 2 0 0 3 0 0 } ] }  {
  error "P2 is not working"
}
# position particle 0-3
set p3 [ observable new particle_positions type { 0 1 } ]
if { ![ veccompare [ observable $p3 print ] { 1 0 0 2 0 0 1 1 0 2 1 0  3 0 0 } ] }  {
  error "P3 is not working"
}


############## Observable particle velocities #######################

set v0 [ observable new particle_velocities id 0 ]
if { ![ veccompare [ observable $v0 print ] { 2 0 0 } ] }  {
  error "V0 is not working"
}
# position particle 0, 3
set v1 [ observable new particle_velocities id { 2 3 } ]
if { ![ veccompare [ observable $v1 print ] { 2 2 0 4 2 0} ] }  {
  error "V1 is not working"
}
# positions particle 0 1
set v2 [ observable new particle_velocities type { 0 } ]
if { ![ veccompare [ observable $v2 print ] { 2 0 0 4 0 0 6 0 0 } ] }  {
  error "V2 is not working"
}
# position particle 0-3
set v3 [ observable new particle_velocities type { 0 1 } ]
if { ![ veccompare [ observable $v3 print ] { 2 0 0 4 0 0 2 2 0 4 2 0  6 0 0 } ] }  {
  error "V3 is not working"
}
# position all particles
set v4 [ observable new particle_velocities all ]
if { ![ veccompare [ observable $v4 print ] { 2 0 0 4 0 0 2 2 0 4 2 0  6 0 0 } ] }  {
  error "V3 is not working"
}

############# Observable com_position ##########

set com_pos1 [ observable new com_position id { 0 1 } ]
if { ![ veccompare [ observable $com_pos1 print ] { 1.5 0 0 } ] }  {
  error "com_pos1 is not working"
}
set com_pos2 [ observable new com_position type 0 ]
if { ![ veccompare [ observable $com_pos2 print ] { 2.25 0 0 } ] }  {
  error "com_pos2 is not working"
}


############# Observable com_position ##########
set com_vel1 [ observable new com_velocity id { 0 1 } ]
if { ![ veccompare [ observable $com_vel1 print ] { 3 0 0 } ] }  {
  error "com_vel1 is not working"
}
set com_vel2 [ observable new com_velocity type 0 ]
if { ![ veccompare [ observable $com_vel2 print ] { 4.5 0 0 } ] }  {
  error "com_vel2 is not working"
}


############# Observable dipole_moment #####################
set dipm [ observable new dipole_moment id { 0 1 } ]
if { ![ veccompare [ observable $dipm print ] { -1 0 0 } ] }  {
  error "dipm is not working"
}

############# Observable current       #####################
set current [ observable new current id { 0 1 } ]
if { ![ veccompare [ observable $current print ] { -2 0 0 } ] }  {
  error "current is not working"
}

############# Observable particle_currents #####################
set particle_currents [ observable new particle_currents id { 0 1 } ]
if { ![ veccompare [ observable $particle_currents print ] { 2 0 0 -4 0 0 } ] }  {
  error "particle_currents is not working"
}

############# Observable stress tensor ##########
set stress_tensor [ observable new stress_tensor ] 
if { ![ veccompare [ observable $stress_tensor print ] [ lreplace [ lindex [ analyze stress ] 0 ] 0 0 ] ] } {
  error "stress_tensor is not working"
}

############# Observable com force and particle force ###### 
set particle_forces [ observable new particle_forces id { 0 1 } ] 
if { ![ veccompare [ observable $particle_forces print ] { 1 0 0 -2 0 0 } ] } {
  error "particle_forces is not working"
}

set com_force [ observable new com_force id { 0 1 } ] 
if { ![ veccompare [ observable $com_force print ] { -1 0 0 } ] } {
  error "particle_forces is not working"
}

############# Observable tcl command #######################
proc p1 {} {
  return { 1 2 3 }
}
set tclcommand [ observable new tclcommand 3 "p1" ] 
puts [ p1 ]
puts "This is the output [ observable $tclcommand print ]"
if { ![ veccompare [ observable $tclcommand print ] { 1 2 3 } ] } {
  error "tclcommand is not working"
}

exit 0
