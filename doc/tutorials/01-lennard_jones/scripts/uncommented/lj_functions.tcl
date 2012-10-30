# Copyright (C) 2010,2011,2012 The ESPResSo project
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
#  Functions (TCL procedures used in lj_tutorial )          #
#                                                           # 
#############################################################

proc rescale_velocities { target_temperature particle_number } {

  set energies [analyze energy]
  set kinetic [lindex $energies 1 1]
  set factor [expr sqrt(1.5*$target_temperature*$particle_number/$kinetic)]

  for {set i 0} {$i<$particle_number} {incr i} {
    set vel [part $i print v]
    part $i v [expr [lindex $vel 0]*$factor] [expr [lindex $vel 1]*$factor] [expr [lindex $vel 2]*$factor]
  }
}

proc save_sim {cfile parinfo range} {
  blockfile $cfile write variable all
  blockfile $cfile write tclvariable all
  blockfile $cfile write particles $parinfo $range
  blockfile $cfile write interactions
  blockfile $cfile write bonds
  blockfile $cfile write random
  blockfile $cfile write seed
  blockfile $cfile write bitrandom
  blockfile $cfile write bitseed
}
