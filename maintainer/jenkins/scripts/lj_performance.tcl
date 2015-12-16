#############################################################
#                                                           #
#  Lennard Jones Liquid                                     #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013,2014 The ESPResSo project
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
#  Parameters                                               #
#############################################################

# System parameters
#############################################################

# 10 000  Particles
set box_l   10.7437
set density 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246

# Integration parameters
#############################################################

setmd time_step 0.01
setmd skin      0.4
thermostat langevin 1.0 1.0

# warmup integration (with capped LJ potential)
set warm_steps   100
set warm_n_times 30
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    1000
set int_n_times  10

# Other parameters
#############################################################
set tcl_precision 6

#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

setmd box_l $box_l $box_l $box_l

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut auto

# Particle setup
#############################################################

set volume [expr $box_l*$box_l*$box_l]
set n_part [expr floor($volume*$density)]

for {set i 0} { $i < $n_part } {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
 
    part $i pos $posx $posy $posz type 0
}

set act_min_dist [analyze mindist]
setmd max_num_cells 2744

# set LJ cap
set cap 20
inter forcecap $cap

# Warmup Integration Loop
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    integrate $warm_steps

    # Warmup criterion
    set act_min_dist [analyze mindist]

#   Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}

inter forcecap 0

set sum_e 0
set sum_f 0

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time] "

    set fit [lindex [time { integrate $int_steps }] 0]
    set eit [lindex [time { analyze energy }] 0]

    set sum_e [expr $sum_e + $eit]
    set sum_f [expr $sum_f + $fit]
    
}

set fd [open "lj_performance.txt" "w"]
puts $fd "LJ_Force,LJ_Energy"
puts $fd "[expr $sum_f /$int_n_times / 1000.], [expr $sum_e / 1000.]"
close $fd

# terminate program
exit
