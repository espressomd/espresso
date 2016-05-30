# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#  Test virutal sites                                       #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "EXTERNAL_FORCES"
require_feature "VIRTUAL_SITES_NO_VELOCITY" off
require_feature "THERMOSTAT_IGNORE_NON_VIRTUAL" off
require_feature "ELECTROSTATICS"
require_feature "LENNARD_JONES"
require_feature "ROTATIONAL_INERTIA"

puts "---------------------------------------------------------------"
puts "- Testcase virtual-sites.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "---------------------------------------------------------------"

setmd box_l 10 10 10
cellsystem domain_decomposition 
setmd time_step 0.01 
setmd skin 1
thermostat off

# check that the min global_cutoff works
setmd min_global_cut 1
if {[setmd max_range] < [setmd min_global_cut] + [setmd skin] - 0.001} {
    error_exit "max cut is too small ( [setmd max_cut_nonbonded] vs. [setmd min_global_cut] + [setmd skin]) "
}
puts "OK: max cut is [setmd max_range], should not be smaller than [setmd min_global_cut] + [setmd skin]"

part 0 pos 5 5 5 quat 1.0 0.0 0.0 0.0 omega_lab 1 2 3
part 1 pos 5 5 6 virtual 1 vs_auto_relate_to 0 
part 2 pos 5 5 4 virtual 1 vs_auto_relate_to 0 
puts [part 1 print pos ]
puts [part 2 print pos ]


integrate 0
if { ([part 1 print pos] != "5.0 5.0 6.0") || ([part 2 print pos] != "5.0 5.0 4.0")} {
 error_exit  "Error: Virutal sites incorrectly positioned: [part 1 print pos] [part 2 print pos]"
} else {
 puts "OK: Position of virtual sites"
}

# Verify that velocity of virtual sites is omega_(central particle) times r
set r [vecsub [part 1 print pos] [part 0 print pos]]
set omega [part 0 print omega_lab]
set v [veccross_product3d $omega $r]
if {[veclen [vecsub $v [part 1 print v]]] >0.0001 } {
 error_exit "Error: Particle 1 velocity incorrect."
} else {
 puts "OK: Velocity particle 1" 
}

set r [vecsub [part 2 print pos] [part 0 print pos]]
set omega [part 0 print omega_lab]
set v [veccross_product3d $omega $r]
if {[veclen [vecsub $v [part 2 print v]]] >0.0001 } {
 error_exit "Error: Particle 2 velocity incorrect. Is [part 2 print v], should be $v."
} else {
 puts "OK: Velocity particle 2" 
}

# Make sure, virtual particles actual motion is parallel to their velocities
set pos1 [part 1 print pos]
set pos2 [part 2 print pos]
set v1 [part 1 print v]
set v2 [part 2 print v]

integrate 1
set diff1 [vecscale 100 [vecsub [part 1 print pos] $pos1]]
set diff2 [vecscale 100 [vecsub [part 2 print pos] $pos2]]
if { [vecdot_product $diff1 $v1] <0.9 } {
 error_exit "Error: Particle 1 moved in wrong direction." 
} else {
 puts "OK: Actual motion of particle 1 is in the right direction. ( $diff1 $v1 ) "
}
if { [vecdot_product $diff2 $v2] <0.9 } {
 error_exit "Error: Particle 1 moved in wrong direction." 
} else {
 puts "OK: Actual motion of particle 2 is in the right direction. ( $diff2 $v2 ) "
}

# Test transformation of forces accumulating on virtual sites
set f1 "3 4 5"
set f2 "4 5 6" 
part 1 ext 3 4 5
part 2 ext 4 5 6
integrate 0
set t [part 0 print torque_lab]
set f [part 0 print f]
if { [veclen [vecsub [vecadd $f1 $f2] $f]] >1E-4 } {
 error_exit "Error: Force on central particle should be [vecadd $f1 $f2] but is $f" 
} else {
 puts "OK: Force on central particle"
}

set r1 [vecsub [part 1 print pos] [part 0 print pos]]
set r2 [vecsub [part 2 print pos] [part 0 print pos]]
set t1 [veccross_product3d $r1 $f1]
set t2 [veccross_product3d $r2 $f2]

if {[veclen [vecsub [vecadd $t1 $t2] $t]] >1E-5 } {
 error_exit "Error: Torque on central particle should be [vecadd $t1 $t2] but is $t" 
} else {
 puts "OK: Torque on central particle"
}



# Testing wether cellsystem and periodic boundaries are handled correctly
part delete

setmd time_step 0.005
cellsystem domain_decomposition 

setmd skin 0.5
setmd box_l 20 20 20 
setmd min_num_cells 27
setmd max_num_cells 100000
setmd periodic 1 1 1
thermostat langevin 1 1 


inter coulomb 1 dh 0.45 0.45
inter 1 1 lennard-jones 10 0.41 0.42

# check that reseting the max cut works
setmd min_global_cut 0
if {[setmd max_range] > 0.45 + [setmd skin] + 0.001} {
    error_exit "max cut ([setmd max_cut_nonbonded]) was not reduced correctly"
}
puts "OK: max cut is [setmd max_range], should not be bigger than 0.45 + [setmd skin]"

setmd min_global_cut 0.21
part 0 pos 2 2 2 v -1 0 0 type 0
part 1 pos 1.8 2 2 virtual 1 vs_auto_relate_to 0 q 1 type 1
part 2 pos 2.2 2 2 virtual 1 vs_auto_relate_to 0 q 1 type 1

set error 0
puts "Checking "
for {set i 0} {$i<10000} {incr i } {
 integrate 1 
 if { ([analyze energy coulomb] <1.5) || ([analyze energy nonbonded 1 1 ] < 12) } {
  puts "$i [analyze energy coulomb] [analyze energy nonbonded 1 1]"
  puts "Error: Lost interaction between particles"
  set error 1
 }
 set v [part 0 print v]
 set v1 [part 1 print v]
 set v2 [part 2 print v]
 if { [veclen [vecsub [vecadd $v1 $v2] [vecscale 2 $v]]] >1E-5 } {
  error_exit "Error: Velocities of outer particles do not add up to twice the velocity of center of mass."
  puts "[vecadd $v1 $v2] $v"
  set error 1
 }
  set r [vecsub [part 1 print pos] [part 0 print pos]]
  set omega [part 0 print omega_lab]
  
  if {abs([veclen $r]-0.2)  >1E-5} {
   error_exit "Distance between particle 0 and 1 incorrect."
  }
  
  set r [vecsub [part 2 print pos] [part 0 print pos]]
  set omega [part 0 print omega_lab]
  
  if {abs([veclen $r]-0.2)  >1E-5} {
   error_exit "Distance between particle 0 and 2 incorrect."
  }
}
puts "OK: Handling of periodic boundaries"
puts "OK: Velocities of outer particles add up to velocity of center of mass"


