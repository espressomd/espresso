# Copyright (C) 2012,2014,2015,2016 The ESPResSo project
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
puts "- Testcase external_potential.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

require_max_nodes_per_side 1

setmd box_l 6. 6. 6.
thermostat off
setmd skin 1.0
setmd time_step 0.05
cellsystem domain_decomposition -no_verlet_list
#cellsystem nsquare
#cellsystem layered
external_potential tabulated file "harmonic_potential.dat" scale "1." 
part 0 pos 4.100  3.000  3.000
set initial_energy [ analyze energy total ] 
for { set i 0 } { $i < 1000 } { incr i } {
  integrate 100
  set current_energy [ analyze energy total  ] 
  puts "$current_energy [ expr $current_energy - $initial_energy ]"
  if { abs($current_energy - $initial_energy ) > 0.02}  {
    puts "Error in energy conservation"
    puts "The deviations are larger than 0.02"
    error_exit
  }
}

puts "Test successful"
ok_exit
