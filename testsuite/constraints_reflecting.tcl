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
# Basic tests of the Lattice Boltzmann implementation       #
#                                                           #
# 1) check conservation of fluid mass                       #
# 2) check conservation of total momentum                   #
# 3) measure temperature of colloid (no auto-check so far)  #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "CONSTRAINTS"
require_feature "LENNARD_JONES"

puts "----------------------------------------"
puts "- Testcase constraint_reflecting.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

#############################################################
# Procedures                                                #
#############################################################

# Integration parameters
#############################################################

set time_step     0.01
set box_l         10.0
set skin          0.5

set prec     5.e-5

# Other parameters
#############################################################

if { [ catch {
#############################################################
# System Setup                                              #
#############################################################

setmd time_step $time_step

# Simulation box
#############################################################
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1
setmd skin $skin

thermostat off

# A Wall constraint
constraint wall normal 1. 0. 0. dist 0.5 reflecting 1 type 0

# Fake interaction 
inter 0 0 lennard-jones 0 1 1 0 0

# After the 1000 integration steps, the particle should be at x=6. y=10.
set refx 6.
set refy 10.

# Just one particle
#############################################################

part 0 pos 5 0 0 v -1 1 0

for { set i 0 } { $i < 100 } { incr i } {
  integrate 10
}

set pos [ part 0 print pos ]
puts "position after the reflection process $pos"
if { abs([ lindex $pos 0 ] - $refx ) > $prec } {
  error "The x position of the particle is wrong"
}
if { abs([ lindex $pos 1 ] - $refy ) > $prec } {
  error "The y position of the particle is wrong"
}
puts "The reflection test was passed!"

# To check bounce back we reset the system
constraint delete 0
constraint wall normal 1. 0. 0. dist 0.5 reflecting 2 type 0
part 0 pos 5 0 0 v -1 1 0

# After the 1000 integration steps, the particle should be at x=6. y=-0.98
set refx 6.
set refy -0.98

for { set i 0 } { $i < 100 } { incr i } {
  integrate 10
}

set pos [ part 0 print pos ]
puts "position after the bounce back $pos"
if { abs([ lindex $pos 0 ] - $refx ) > $prec } {
  error "The x position of the particle is wrong"
}
if { abs([ lindex $pos 1 ] - $refy ) > $prec } {
  error "The y position of the particle is wrong"
}
puts "The bounce back test was passed!"



} res ] } {
    error_exit $res
}

exit 0

#############################################################
