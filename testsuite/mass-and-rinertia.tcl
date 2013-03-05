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
source "tests_common.tcl"

require_feature "MASS"
require_feature "ROTATIONAL_INERTIA"

# Cchecks if every degree of freedom has 1/2 kT of energy, even when
# mass and inertia tensor are active

set kT 1.5
set halfkT 0.75
thermostat langevin $kT 1

# no need to rebuild Verlet lists, avoid it
setmd skin 1.0
setmd time_step 0.01

set n 500
set mass 200
set j1 300
set j2 400
set j3 500
for {set i 0} {$i<$n} {incr i} {
  part $i pos 0 0 0 rinertia $j1 $j2 $j3 mass $mass
}

set vx2 0.
set vy2 0.
set vz2 0.
set ox2 0.
set oy2 0.
set oz2 0.


set loops 3
puts "Thermalizing..."
integrate 300
puts "Measuring..."

for {set i 0} {$i <$loops} {incr i} {
 integrate 100
 # Get kinetic energy in each degree of freedom for all particles
 for {set p 0} {$p <$n} {incr p} {
  set v [part $p print v]
  set o [part $p print omega]
  set ox2 [expr $ox2 +pow([lindex $o 0],2)]
  set oy2 [expr $oy2 +pow([lindex $o 1],2)]
  set oz2 [expr $oz2 +pow([lindex $o 2],2)]
  set vx2 [expr $vx2 +pow([lindex $v 0],2)]
  set vy2 [expr $vy2 +pow([lindex $v 1],2)]
  set vz2 [expr $vz2 +pow([lindex $v 2],2)]
 }
}

set tolerance 0.1
set Evx [expr 0.5 * $mass *$vx2/$n/$loops]
set Evy [expr 0.5 * $mass *$vy2/$n/$loops]
set Evz [expr 0.5 * $mass *$vz2/$n/$loops]

set Eox [expr 0.5 * $j1 *$ox2/$n/$loops]
set Eoy [expr 0.5 * $j2 *$oy2/$n/$loops]
set Eoz [expr 0.5 * $j3 *$oz2/$n/$loops]

set dv [expr 1./3. *($Evx +$Evy +$Evz)/$halfkT-1.]
set do [expr 1./3. *($Eox +$Eoy +$Eoz)/$halfkT-1.]

puts "1/2 kT = $halfkT"
puts "translation: $Evx $Evy $Evz rotation: $Eox $Eoy $Eoz"

puts "Deviation in translational energy: $dv"
puts "Deviation in rotational energy: $do"

if { abs($dv) > $tolerance } {
 error "Relative deviation in translational energy too large: $dv"
}
if { abs($do) > $tolerance } {
 error "Relative deviation in translational energy too large: $dv"
}

exit 0