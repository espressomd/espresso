#
# print the maximal error in the numerical derivation of a potential
# against the displacement dr. This is a central derivate, so we
# expect at least a dr^3 scaling, if the forces are continuous.
#
# This can be used to check whether the forces and energy functions of
# a potential fit.
#
#######################################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
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

# box size, should be larger than the potentials cutoff
set box_l 10.0
setmd box_l $box_l $box_l $box_l

# number of interaction partners, 2 for nonbonded
set n_part 4

# create the particles so that you can create bonds
for {set i 0} {$i < $n_part} {incr i} {
    part $i pos 0 0 0
}

# The potential to test, including parameters.
#
# - can be bonded or nonbonded
# - all particles have type 0
inter 1 dihedral 1 1 0
part 1 bond 1 0 2 3

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

setmd periodic 0 0 0
thermostat off
setmd time_step 1
setmd skin 0

puts "# dr maximal deviation"
puts "# should scale like dr^3 for sufficiently small dr"
for {set dr 0.01} {$dr > 0.0001} {set dr [expr $dr*0.5]} {
    set maxDev 0
    for {set r 0} {$r < $tries} {incr r} {
	# 1. position
	while {1} {
	    # particles in a sphere around $bp
	    set bp "[expr $box_l*rand()] [expr $box_l*rand()] [expr $box_l*rand()]"
	    for {set i 0} {$i < $n_part} {incr i} {
		while {1} {
		    set vec ""
		    for {set j 0} {$j<3} {incr j} {
			lappend vec [expr (2.0*rand() - 1.0)*$max_dist/2] }
		    if { [veclen $vec] <= $max_dist/2 } { break }
		}
		set p1($i) [vecadd $bp $vec]
		eval part $i pos $p1($i)
	    }
	    if {[analyze mindist] > $min_dist} { break }
	}

	set e1 [analyze energy total]
	integrate 0
	for {set i 0} {$i < $n_part} {incr i} {
	    set f1($i) [part $i print force]
	}

	# 2., slightly displaced position

	for {set i 0} {$i < $n_part} {incr i} {
	    set dx($i) [vec_random $dr]	    
	    set p2($i) [vecadd $p1($i) $dx($i)]
	    eval part $i pos $p2($i)
	}

	set e2 [analyze energy total]
	integrate 0
	for {set i 0} {$i < $n_part} {incr i} {
	    set f2($i) [part $i print force]
	}

	# Check forces
	set dE [expr $e2 - $e1]
	set fdotdx 0
	for {set i 0} {$i < $n_part} {incr i} {
	    # average force
	    set force [vecscale 0.5 [vecadd $f1($i) $f2($i)]]
	    set fdotdx [expr $fdotdx + [vecdot_product $force $dx($i)]]
	}

	set dev [expr abs($dE + $fdotdx)]
	if {$dev > $maxDev} { set maxDev $dev }
    }
    puts "$dr $maxDev"
}

exit