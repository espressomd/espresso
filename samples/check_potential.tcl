# a test for energy/force correspondence for simple radial potentials
# 
# takes the numerical derivate of the potential and allows to compare
# with the force at randomly chosen points.
#
# Copyright (C) 2010,2012 The ESPResSo project
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

# the potential to test, including parameters
set potential "hertzian 1 1"
# displacement for discrete derivative
set dr 0.0001
# system box, should be big enough for the cutoff
set box_l 10
# number of checks
set tries 10000

setmd box_l $box_l $box_l $box_l
eval inter 0 0 $potential
thermostat off
setmd time_step 1
setmd skin 0

part 0 pos 0 0 0
for {set r 0} {$r < $tries} {incr r} {
    set p "[expr $box_l*rand()] [expr $box_l*rand()] [expr $box_l*rand()]"
    set distvec [bond_vec_min "0 0 0" $p]
    set dir     [vecnorm $distvec]
    
    eval part 1 pos $p
    set e1 [analyze energy total]
    integrate 0
    set f1 [part 1 print force]
    
    eval part 1 pos [vecadd $p [vecscale $dr $dir]]
    set e2 [analyze energy total]
    integrate 0
    set f2 [part 1 print force]
    
    set force [vecscale 0.5 [vecadd $f1 $f2]]

    # expected value from numerical derivative
    set rad_force_expect [expr -($e2 - $e1)/$dr]
    # measured value along axis
    set rad_force [vecdot_product $force $dir]
    set rad_force_diff [expr $rad_force_expect - $rad_force]
    # non-radial part of the force, should be zero
    set force_expect [vecscale $rad_force_expect $dir]
    set force_excess [veclen [vecsub $force_expect $force]]

    set dist [veclen $distvec]

    puts "$dist $rad_force_diff $force_excess"
}
