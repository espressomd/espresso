# place particles to generate a certain dihedral angle and plot the
# dihedral energy against this angle.
#
# The script creates three setups, with the line connecting the center
# parallel to one of the three major axis. All three columns should give
# the same output.
#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2012 The ESPResSo project
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

# box size
set box_l 10.0
setmd box_l $box_l $box_l $box_l

# The dihedral potential including parameters.
inter 1 dihedral 2 3 0

# minimal distance of the particles, to avoid areas with large forces
# where the convergence is bad
set min_dist 0

# particles are always in a sphere of this diameter
# choose bigger than the cutoff to test the shifting
set max_dist 1.4

# number of tests
set tries 10000

# below you typically don't need to change anything
##################################################

# create the particles so that we can create the bond
for {set i 0} {$i < 4} {incr i} {
    part $i pos 0 0 0
}
part 1 bond 1 0 2 3

setmd periodic 0 0 0
thermostat off
setmd time_step 1
setmd skin 0

puts "# angle E1 E2"
for {set angle 0} {$angle < 2*[PI]} {set angle [expr $angle + 0.01]} {
    part 1 pos 0 0 0
    part 2 pos 1 0 0
    part 3 pos 1 1 0
    part 0 pos 0 [expr cos($angle)] [expr sin($angle)]
    set e1 [analyze energy total]

    part 1 pos 0 0 0
    part 2 pos 0 1 0
    part 3 pos 0 1 1
    part 0 pos [expr sin($angle)] 0 [expr cos($angle)]
    set e2 [analyze energy total]

    part 1 pos 0 0 0
    part 2 pos 0 0 1
    part 3 pos 1 0 1
    part 0 pos [expr cos($angle)] [expr sin($angle)] 0
    set e3 [analyze energy total]

    puts "$angle $e1 $e2 $e3"

}

exit