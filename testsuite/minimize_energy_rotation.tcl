#############################################################
#                                                           #
#  Lennard Jones Liquid                                     #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
setmd min_global_cut 0.1

require_feature "ROTATION"
require_feature "ROTATION_PER_PARTICLE"
require_feature "DIPOLES"
setmd skin 0
setmd time_step 0.01
thermostat off


part 0 pos 0 0 0 dip 2 3 4
set l [veclen "2 3 4"]
constraint ext_magn_field -1 1 1
minimize_energy 1E-8 10000 0.01 0.01
set target "-1 1 1"
set target [vecscale [expr $l /[veclen $target]] $target]
puts "Expected: $target, got [part 0 print dip]"

if { [veclen [vecsub [part 0 print dip] $target]] > 0.0001 } {
  error "Dipole moment did not align to external field for unconstraint rotation."
}

# Prevent rotation of the particle
part 0 pos 0 0 0 dip 2 3 4 rotation 0
set l [veclen "2 3 4"]
constraint ext_magn_field -1 1 1
minimize_energy 1E-8 10000 0.01 0.01
set target "2 3 4"
set target [vecscale [expr $l /[veclen $target]] $target]
puts "Expected: $target, got [part 0 print dip]"

if { [veclen [vecsub [part 0 print dip] $target]] > 0.0001 } {
  error "Dipole moment did not align to external field for constraint rotation (rotation 0)."
}


# Freeze only rotation around x and y axes. This should still not rotate at all,
# but tests a different part of the code
part 0 pos 0 0 0 dip 2 3 4 rotation 8
set l [veclen "2 3 4"]
constraint ext_magn_field -1 1 1
minimize_energy 1E-8 10000 0.01 0.01
set target "2 3 4"
set target [vecscale [expr $l /[veclen $target]] $target]
puts "Expected: $target, got [part 0 print dip]"

if { [veclen [vecsub [part 0 print dip] $target]] > 0.0001 } {
  error "Dipole moment did not align to external field for partially constraint rotation (rotation 8)."
}


