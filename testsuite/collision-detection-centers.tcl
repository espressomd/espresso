# Copyright (C) 2010,2011 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
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
#  Test collision detection with binding of centers of colliding particles
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "COLLISION_DETECTION"

setmd box_l 10 10 10

thermostat off
sdetmd time_step 0.01
setmd skin 0
part 0 pos 0 0 0 
part 1 pos 1 0 0

paaart 2 pos 3 0 0


if { "[on_collision]" != "off" } {
  error_exit "Collision detection should be off by default."
}

on_collision off
if { "[on_collision]" != "off" } {
  error_exit "Disabling collision_detection does not work"
}

on_collision bind_centers 1.0 1
if { "[on_collision]" != "bind_centers 1.0 1" } {
  error_exit "Setting collision_detection parameters for bind_centers does not work"
}

integrate 0 

part 0 print bond
part 1 print bond
part 2 print bond

