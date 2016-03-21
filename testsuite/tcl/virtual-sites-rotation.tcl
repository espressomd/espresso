
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
set box 20
setmd box_l $box $box $box
setmd min_global_cut 1.1

require_feature "VIRTUAL_SITES_RELATIVE"
require_feature "THERMOSTAT_IGNORE_NON_VIRTUAL" off
require_feature "ROTATIONAL_INERTIA"
require_feature "MASS"

puts "---------------------------------------------------------------"
puts "- Testcase virtual-sites-rotation.tcl running on  [format %02d [setmd n_nodes]] nodes  -"
puts "---------------------------------------------------------------"
cellsystem nsquare -no_verlet_list

set mass [expr rand() *20]
set j1 [expr rand() * 20]
set j2 [expr rand() * 20]
set j3 [expr rand() * 20]


set kT 1.5
set halfkT 0.75
thermostat langevin $kT 1

# no need to rebuild Verlet lists, avoid it
setmd skin 0.0
setmd time_step 0.01

set n 4
for {set p 0 } { $p < [expr $n*2] } { incr p 2} { 
  part $p pos [expr $box/$n] 0 0
  part [expr $p+1] pos [expr $box/$n+1] 0 0 vs_auto_relate_to 0 virtual 1
  part [expr $p+1] rinertia $j1 $j2 $j3 mass $mass
}

set ox2 0.
set oy2 0.
set oz2 0.


set loops 4000 
puts "Thermalizing..."
integrate 1000
puts "Measuring..."

if { [catch {
  for {set i 0} {$i <$loops} {incr i} {
   integrate 50
   # Get kinetic energy in each degree of freedom for all particles
    for {set p 1 } { $p < [setmd n_part] } { incr p 2} { 
       set o [part $p print omega_body]
       set ox [lindex $o 0]
       set oy [lindex $o 1]
       set oz [lindex $o 2]
       set ox2 [expr $ox2 +$ox*$ox]
       set oy2 [expr $oy2 +$oy*$oy]
       set oz2 [expr $oz2 +$oz*$oz]
    }
  }
  
  set tolerance 0.1
  
  set Eox [expr 0.5 * $j1 *$ox2/$n/$loops]
  set Eoy [expr 0.5 * $j2 *$oy2/$n/$loops]
  set Eoz [expr 0.5 * $j3 *$oz2/$n/$loops]
  
  set do [expr 1./3. *($Eox +$Eoy +$Eoz)/$halfkT-1.]
  
  puts "1/2 kT = $halfkT"
  puts "rotation: $Eox $Eoy $Eoz"
  
  if { abs($do) > $tolerance } {
   error "Relative deviation in rotational energy too large: [expr abs($do)]"
  } else { 
   puts "Deviation in rotational energy: [expr abs($do)], this is OK (less than [expr 100*$tolerance]%)"
  }
  
} res ] } { 
  error_exit $res
}

exit 0
