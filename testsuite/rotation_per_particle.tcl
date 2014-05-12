# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

require_feature "ROTATION"
require_feature "ROTATION_PER_PARTICLE"

puts "----------------------------------------------"
puts "- Testcase rotation_per_particle.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

set epsilon 5e-4
thermostat langevin 1 1

setmd time_step 0.01
setmd skin 0


part 0 pos 0 0 0

# Check that rotation is on by default and that the rotation al degrees are thermalized
set rot [part 0 print rotation]
if {$rot != "1"} {
 error "Rotation not on by default!"
}

integrate 100
set quat [part 0 print quat]
set zero "1. 0. 0. 0."
if {[veclen [vecsub $quat $zero]] < $epsilon} {
 error "Rotational degrees not thermalized when rotation is on"
}

# Check that the rotational properties don't change any more, when rotation is off
set quat [part 0 print quat]
set omega [part 0 print omega_lab]
part 0 rotation 0
integrate 1000
set quatN [part 0 print quat]
set omegaN [part 0 print omega_lab]

if {[veclen [vecsub $quat $quatN]] >$epsilon} {
 error "Quaternions changed even when rotation is off"
}

if {[veclen [vecsub $omega $omegaN]] > $epsilon} {
 error "omegaernions changed even when rotation is off"
}

exit 0
