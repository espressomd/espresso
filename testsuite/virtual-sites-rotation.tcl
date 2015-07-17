
# Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
setmd box_l 20 20 20 
setmd min_global_cut 5

require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "THERMOSTAT_IGNORE_NON_VIRTUAL" off
require_feature "ROTATIONAL_INERTIA"

puts "---------------------------------------------------------------"
puts "- Testcase virtual-sites-rotation.tcl running on 1 nodes"
puts "---------------------------------------------------------------"
cellsystem nsquare -no_verlet_list
set mass 200
set j1 60.
set j2 80. 
set j3 100.

part 0 pos 0 0 0
part 1 pos 5 0 0 vs_auto_relate_to 0 virtual 1
part 1 rinertia $j1 $j2 $j3 mass $mass


set kT 1.5
set halfkT 0.75
thermostat langevin $kT 1

# no need to rebuild Verlet lists, avoid it
setmd skin 1.0
setmd time_step 0.01

set n 1.


set ox2 0.
set oy2 0.
set oz2 0.



set ox2 0.
set oy2 0.
set oz2 0.


set loops 18000 
puts "Thermalizing..."
integrate 10000
puts "Measuring..."

for {set i 0} {$i <$loops} {incr i} {
 integrate 100
 # Get kinetic energy in each degree of freedom for all particles
  set p 1
  set o [part $p print omega_lab]
  set ox2 [expr $ox2 +pow([lindex $o 0],2)]
  set oy2 [expr $oy2 +pow([lindex $o 1],2)]
  set oz2 [expr $oz2 +pow([lindex $o 2],2)]
}

set tolerance 0.1

set Eox [expr 0.5 * $j1 *$ox2/$n/$loops]
set Eoy [expr 0.5 * $j2 *$oy2/$n/$loops]
set Eoz [expr 0.5 * $j3 *$oz2/$n/$loops]

set do [expr 1./3. *($Eox +$Eoy +$Eoz)/$halfkT-1.]

puts "1/2 kT = $halfkT"
puts "rotation: $Eox $Eoy $Eoz"

puts "Deviation in rotational energy: $do"

if { abs($do) > $tolerance } {
 error "Relative deviation in rotational energy too large: $do"
}

exit 0
