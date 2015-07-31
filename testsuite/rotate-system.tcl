
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

source "tests_common.tcl"

set tcl_precision 14

proc vectorsTheSame {a b} {
 set tol 1E-4
 set diff [vecsub $a $b]
 if { [veclen $diff] > $tol } {
  return 0
 }
 return 1
}

require_feature "PARTIAL_PERIODIC"
require_max_nodes_per_side 1

puts "----------------------------------------------"
puts "- Testcase rotate-system.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

set l 9
setmd box_l $l $l $l
setmd periodic 0 0 0

set n 100

set oldPos ""
for {set i 0} {$i <$n} {incr i} {
 part $i pos [expr $l *[t_random]] [expr $l *[t_random]] [expr $l *[t_random]]
 set oldPos "$oldPos {[part $i print pos]}"
}

set theta [expr acos(2*[t_random]-1.)]
set phi [expr [t_random] *2. *3.142]
set axis "[expr sin($theta) *cos($phi)] [expr sin($theta)*sin($phi)] [expr cos($theta)]"
set alpha [expr [t_random] *2. *3.142]

puts "Rotating: phi=$phi, theta=$theta, alpha=$alpha, Axis=$axis"

set oldCom [analyze centermass 0]
puts "Center of mass berfore rotation: $oldCom"

rotate_system $phi $theta $alpha

set newCom [analyze centermass 0]
puts "Center of mass after rotation: $newCom"

if { ! [vectorsTheSame $oldCom $newCom] } {
 error "Center of mass changed during rotation"
}


set newCom [analyze centermass 0]
for {set i 0} {$i < $n} {incr i} {
 set newPos [part $i print pos]

 # Calculate expected position from position before rotation
 
 # Shift center of mass to origin
 set expPos [vecsub [lindex $oldPos $i] $oldCom]
 
 # Rotate using Rodrigues formula
 set a [vecscale [expr cos($alpha)] $expPos]
 set b [vecscale [expr sin($alpha)] [veccross_product3d $axis $expPos]]
 set c [vecscale [expr (1. -cos($alpha)) *[vecdot_product $axis $expPos]] $axis]
 set expPos [vecadd $a [vecadd $b $c]]

 # Re-add center of mass
 set expPos [vecadd $oldCom $expPos] 

 # Compare
 if { ! [vectorsTheSame $newPos $expPos]} {
  puts "Particle $i, new position is $newPos, but $expPos expected. Old pos: [lindex $oldPos $i]"
  puts "  Diff: [vecsub $expPos $newPos], angle=[expr acos([vecdot_product $newPos $expPos] / [veclen $newPos] /[veclen $expPos])]"
  error "Expected and measured positions don't match"
 }
}

