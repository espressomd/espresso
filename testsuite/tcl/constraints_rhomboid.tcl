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

require_feature "LENNARD_JONES"
require_feature "CONSTRAINTS"

puts "---------------------------------------------------------"
puts "- Testcase constraint_rhomboid.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "---------------------------------------------------------"

if { [ catch {
#############################################################
# System Setup                                              #
#############################################################

  setmd time_step 0.05
  setmd skin 0.1
  thermostat off
  
# Simulation box
#############################################################
  setmd box_l 15 15 15
  constraint rhomboid corner 5 5 5 b 5 0 0 a 0 5 0 c 0 0 5 direction outside type 1

# Particles
#############################################################
  part 0 pos 7.5 7.5 0 v 0 0 1 type 0
  part 1 pos 0 7.5 7.5 v 1 0 0 type 0
  part 2 pos 7.5 0 7.5 v 0 1 0 type 0
  part 3 pos 7.5 7.5 15 v 0 0 -1 type 0
  part 4 pos 15 7.5 7.5 v -1 0 0 type 0
  part 5 pos 7.5 15 7.5 v 0 -1 0 type 0

  part 6 pos 0 0 0 v 1 1 1 type 0
  part 7 pos 0 15 0 v 1 -1 1 type 0
  part 8 pos 0 15 15 v 1 -1 -1 type 0
  part 9 pos 0 0 15 v 1 1 -1 type 0
  part 10 pos 15 0 0 v -1 1 1 type 0
  part 11 pos 15 15 0 v -1 -1 1 type 0
  part 12 pos 15 15 15 v -1 -1 -1 type 0
  part 13 pos 15 0 15 v -1 1 -1 type 0

  part 14 pos 7.5 0 0 v 0 1 1 type 0
  part 15 pos 0 7.5 0 v 1 0 1 type 0
  part 16 pos 0 0 7.5 v 1 1 0 type 0
  part 17 pos 7.5 15 15 v 0 -1 -1 type 0
  part 18 pos 15 7.5 15 v -1 0 -1 type 0
  part 19 pos 15 15 7.5 v -1 -1 0 type 0
  part 20 pos 7.5 0 15 v 0 1 -1 type 0
  part 21 pos 0 7.5 15 v 1 0 -1 type 0
  part 22 pos 15 0 7.5 v -1 1 0 type 0
  part 23 pos 0 15 7.5 v 1 -1 0 type 0
  part 24 pos 7.5 15 0 v 0 -1 1 type 0
  part 25 pos 15 7.5 0 v -1 0 1 type 0

  inter 0 1 lennard-jones 1. 1. 1.1225 0.25 0

#############################################################
# Integration                                               #
#############################################################

  integrate 150

#############################################################
# Analysis and Verification                                 #
#############################################################

  set pos { {7.5 7.5 0.583545561232}
            {0.583545561232 7.5 7.5}
            {7.5 0.583545561232 7.5}
            {7.5 7.5 14.4164544388}
            {14.4164544388 7.5 7.5}
            {7.5 14.4164544388 7.5}
            {1.41248232369 1.41248232369 1.41248232369}
            {1.41248232369 13.5875176763 1.41248232369}
            {1.41248232369 13.5875176763 13.5875176763}
            {1.41248232369 1.41248232369 13.5875176763}
            {13.5875176763 1.41248232369 1.41248232369}
            {13.5875176763 13.5875176763 1.41248232369}
            {13.5875176763 13.5875176763 13.5875176763}
            {13.5875176763 1.41248232369 13.5875176763}
            {7.5 1.15757330101 1.15757330101}
            {1.15757330101 7.5 1.15757330101}
            {1.15757330101 1.15757330101 7.5}
            {7.5 13.842426699 13.842426699}
            {13.842426699 7.5 13.842426699}
            {13.842426699 13.842426699 7.5}
            {7.5 1.15757330101 13.842426699}
            {1.15757330101 7.5 13.842426699}
            {13.842426699 1.15757330101 7.5}
            {1.15757330101 13.842426699 7.5}
            {7.5 13.842426699 1.15757330101}
            {13.842426699 7.5 1.15757330101}
          }

  set sum 0.0
  
  for {set i 0} {$i <= 25} {incr i} {
    set part_pos [part $i print pos]
    
    for {set k 0} {$k < 3} {incr k} {
      set sum [expr $sum + [expr abs([lindex $part_pos $k]-[lindex [lindex $pos $i] $k])]]
    }
  }
  
  if {$sum > 1.e-8} {
    error "Combined particle trajectory deviation of [format %1.3E $sum] exceeds limit of 1.000E-08 after reflection on rhomboid."
  }

} res ] } {
   error_exit $res
}

exit 0
