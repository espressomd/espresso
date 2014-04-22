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

puts "---------------------------------------------------------------"
puts "- Testcase correlation.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"
setmd box_l 10 10 10
setmd time_step 0.01
setmd skin 1
thermostat langevin 1 1 
thermostat off
part 0 pos 0 0 0 v 1 2 3

observable new particle_positions all

correlation new obs1 0 dt 0.1 tau_max 10 corr_operation square_distance_componentwise compress1 linear
correlation new obs1 0 dt 0.01 tau_max 10 corr_operation square_distance_componentwise tau_lin 10 
integrate 1000
correlation 0 autoupdate start
correlation 1 autoupdate start
integrate 20000

#set corr [ correlation ]
#puts "correlations: \n$corr";
set corr [ correlation 1 print ]
#puts "obs: [observable 0 print]"
#puts "corr: \n$corr";
for { set i 0 } { $i < [ llength $corr ]} { incr i } {
  set t [ lindex [ lindex $corr $i ] 0 ] 
  if { abs([ lindex [ lindex $corr $i ] 2 ] -  $t*$t ) > 0.0001 } {
    error "test failed1"
  }
  if { abs([ lindex [ lindex $corr $i ] 3 ] -  2*2*$t*$t ) > 0.0001 } {
    error "test failed2, "
  }
  if { abs([ lindex [ lindex $corr $i ] 4 ] -  3*3*$t*$t ) > 0.0001 } {
    error "test failed3"
  }
}

exit 0
