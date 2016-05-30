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

source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "EXTERNAL_FORCES"
require_feature "FFTW"
# we can't reach our precision with more than 2 CPUS per side
require_max_nodes_per_side {2 2 2}

puts "---------------------------------------------------------------"
puts "- Testcase iccp3m.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 10 10 10
setmd time_step 0.01
setmd skin 0.5
thermostat off

part 000 pos 5 5 1 q 1
set counter 1

set normals [ list ] 
set areas [ list ] 
set normals [ list ] 
set epsilons [ list ] 
for { set i 0 } { $i < 10 } { incr i } {
  for { set j 0 } { $j < 10 } { incr j } {
    part $counter pos $i $j 0 fix 1 1 1 q 1.
    lappend normals { 0 0 1. }
    lappend areas 1 
    lappend epsilons .5 
    incr counter 
  }
}
for { set i 0 } { $i < 10 } { incr i } {
  for { set j 0 } { $j < 10 } { incr j } {
    part $counter pos $i $j 5 fix 1 1 1 q 1.
    lappend normals { 0 0 -1. }
    lappend areas 1 
    lappend epsilons .1 
    incr counter 
  }
}

puts "[ inter coulomb 1. p3m tunev2 accuracy 1e-3 mesh 32 cao 4 ]"

iccp3m 200 eps_out 1 max_iterations 60 convergence 1e-1 relax 0.7 areas $areas normals $normals epsilons $epsilons first_id 1
if {[catch {
    integrate 0
    integrate 100
} res]} {
    error_exit "iccp3m: caught error $res"
}

set refpos 1.0429519727877221 
set refforce 0.07816245324943571
set pos [ lindex [ part 0 print pos ] 2 ] 
set force [ lindex [ part 0 print force ] 2 ]
if { abs($refpos-$pos) > 1e-4 } {
    error_exit "iccp3m: position is wrong"
}
if { abs($refforce - $force) > 1e-4 } {
    error_exit "iccp3m: force is wrong"
}

exit 0
