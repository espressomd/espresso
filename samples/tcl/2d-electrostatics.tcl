# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#    Max-Planck-Institute for Polymer Research, Theory Group
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

#######################################################################
#
# minimal example of using MMM2D / P3M+ELC to simulate a salt melt
# between walls. No data is recorded or printed, just to demonstrate
# the setup of the electrostatics methods.
#
#######################################################################

# choose the electrostatics algorithm, mmm2d or elc
set method "elc"
thermostat langevin 1.0 1.0
setmd time_step 0.01
setmd skin 0.3
set n_pairs 100
set box_l 10
# z needs to be a bit smaller for ELC, not a problem with MMM2D
set box_l_z 9

# always a cubic box, we make z shorter by a wall
setmd box_l $box_l $box_l $box_l

#    create particles in a slab
#
# Keeps a distance of 1 to both walls, to
# avoid collisions with the wall during warmup.
# Due to the LJ potential, the particle centers will
# maintain this typical distance also during the
# actual run.
####################################################

for {set p 0} {$p < $n_pairs} {incr p} {
    part $p pos [expr $box_l*rand()] [expr $box_l*rand()] [expr 1 + ($box_l_z - 2)*rand()] q +1
    part [expr $p + $n_pairs] pos [expr $box_l*rand()] [expr $box_l*rand()] [expr 1 + ($box_l_z - 2)*rand()] q -1
}

#    confining wall
#########################################

constraint wall normal  0  0  1 dist 0 type 1
constraint wall normal  0  0 -1 dist [expr -$box_l_z] type 1

#    wall excluded volume interaction
#########################################

inter 0 1 lennard-jones 1.0 1.0 1.12246 auto 0.0

#    warmup
#
# In order to avoid overlaps, gradually make the
# potential steeper. But only particle-particle, not
# particle-wall.
####################################################

set capradius 1.123

inter forcecap individual

set cnt 0
while {$capradius > 0.9} {
    if {$cnt % 10 == 0} {
	puts "capradius = $capradius"
    }
    incr cnt
    set capradius [expr $capradius*0.99]
    inter 0 0 lennard-jones 1.0 1.0 1.12246 auto 0.0 $capradius
    
    integrate 10
}

# turn force capping off
inter forcecap 0
inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

#    Coulomb interactions
#########################################

if {$method == "mmm2d"} {
    #########################################
    #        setup using MMM2D
    #########################################

    setmd periodic 1 1 0
    # ideally, this should be tuned to best performance. Just try different values 
    cellsystem layered [expr 6/[setmd n_nodes]]
    inter coulomb 1.0 mmm2d 1e-4
} {
    #########################################
    #      setup using P3M + ELC 
    #########################################

    cellsystem domain_decomposition

    setmd periodic 1 1 1
    puts [inter coulomb 1.0 p3m tunev2 mesh 32 accuracy 1e-4]
    inter coulomb elc 1e-4 [expr $box_l - $box_l_z]
}

# just simulate
#########################################

while {1} {
    puts "time: [setmd time]"
    integrate 100
}
exit 0
