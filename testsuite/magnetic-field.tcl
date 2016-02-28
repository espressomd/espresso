# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
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

source "tests_common.tcl"

require_feature "DIPOLES"
require_feature "ROTATION"
require_feature "CONSTRAINTS"

# Test for the correctness of energy of a dipole in a magnetic field and torque on 
# a dipole in a magnetic field

# Setup
thermostat off
setmd skin 0
setmd time_step 0.01

part 0 pos 0 0 0 dip 1 2 3
constraint ext_magn_field 4 5 6

# Check energy
set E [analyze energy total]
if { abs($E+32.) > 1E-11 } {
 error "Incorrect energy for dipole in magnetic field: $E should be -32"
}

# Check torque
integrate 0 
set T [part 0 print torque_lab]
set expected [veccross_product3d {1 2 3} {4 5 6}]
if { [veclen [vecsub $expected $T]] >1E-6 } {
 error  "Torque for dipole in magnetic field incorrect. Should be $expected but is $T"
}

puts "Energy and torque for a dipole in a homogeneous magnetic field are correct."

exit 0
