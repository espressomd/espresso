#############################################################
#                                                           #
#  Lennard Jones Liquid                                     #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

puts " "
puts "======================================================="
puts "=       h5_lj.tcl                                 ="
puts "======================================================="
puts " "

#############################################################
#  Parameters                                               #
#############################################################
set name  "lj_liquid"

set n_part 500
set box_l   10.7437
setmd box_l $box_l $box_l $box_l


set int_n_times  500
set int_steps    10
setmd time_step 0.01

set warm_steps   100
set warm_n_times 30

setmd skin      0.4
set min_dist     0.9

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut auto
thermostat langevin 1.0 1.0


#############################################################
#  Start positions                                          #
#############################################################
for {set i 0} { $i < $n_part } {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
 
    part $i pos $posx $posy $posz type 0
}

#############################################################
#  Warmup Integration                                       #
#############################################################
set act_min_dist [analyze mindist]
#setmd max_num_cells 2744
puts "\nStart warmup integration:"
set cap 20
inter forcecap $cap

# Warmup Integration Loop
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    integrate $warm_steps
    # Warmup criterion
    set act_min_dist [analyze mindist]
    puts "run $i at time=[setmd time] (LJ cap=$cap) min dist = $act_min_dist\r"
    flush stdout
	# Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}


#############################################################
#  H5MD file and datasets                                   #
#############################################################
if { [file exists "h5_lj.h5"] != 1 } {
    # Intitialize H5MD dataset. If it already exists it will be opened and extended
	h5md_init "h5_lj.h5" 
    # Initialize user defined H5MD 1D-observable e.g. energy
    h5md_observable1D_init "energy"
} else {
    h5md_init "h5_lj.h5" 
}

#############################################################
#      Integration                                          #
#############################################################
puts "\nStart integration: run $int_n_times times $int_steps steps"
inter forcecap 0

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time]\r"
    integrate $int_steps
    set E [analyze energy kinetic]
	# write to user defined H5MD observable
    h5md_observable1D_write "energy" $E
	# write all particle positions to H5MD dataset
    h5md_write_positions
}
# Terminate program
puts "\n\nFinished"
# Close H5MD-dataset
h5md_close
exit
