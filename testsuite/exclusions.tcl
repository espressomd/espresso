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
#############################################################
#                                                           #
#  Test exclusions
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "EXCLUSIONS"

setmd skin 0.1
setmd time_step 0.01
thermostat off
part 0 pos 0 0 0
part 1 pos 0.01 0.01 0.1
inter 0 0 lennard-jones 1 0.1 0.112 auto
integrate 0
set E [analyze energy nonbonded 0 0]
set f0 [part 0 print f]
set f1 [part 1 print f]
if { [veclen $f0 ] <1 } {
 error_exit "Force not there while exclusion off."
}
if { [veclen $f1 ] <1 } {
 error_exit "Force not there while exclusion off."
}
if { $E <0.1 } {
 error_exit "Energy not there while exclusion off."
}

part 0 exclude 1

integrate 0
set E [analyze energy nonbonded 0 0]
set f0 [part 0 print f]
set f1 [part 1 print f]
if { [veclen $f0 ] >0.001 } {
 error_exit "Force found for a pair with an exclusion."
}
if { [veclen $f1 ] >0.001 } {
 error_exit "Force found for a pair with an exclusion."
}
if { $E >0.001 } {
 error_exit "Energy found for a pair with an exclusion."
}

exit 0
