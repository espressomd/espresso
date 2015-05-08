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
setmd max_num_cells 2744
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
#   Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}


#############################################################
#  H5MD file and datasets                                   #
#############################################################
#File
h5mdfile H5Fcreate "h5_lj.h5" 
#Groups
h5mdfile H5Gcreate2 "particles"
h5mdfile H5Gcreate2 "observables"



#Dataset positions
h5mdfile H5Screate_simple type double dims 1 $n_part 3
h5mdfile H5Pset_chunk dims 5 $n_part 3
h5mdfile H5Dcreate2 "/particles/pos"
h5mdfile H5Dopen2 "/particles/pos"
for { set k 0 } { $k < $n_part } { incr k } {
    set pos [part $k pr p]
    h5mdfile H5_write_value value [lindex $pos 0] index 0 $k 0
    h5mdfile H5_write_value value [lindex $pos 1] index 0 $k 1
    h5mdfile H5_write_value value [lindex $pos 2] index 0 $k 2
}
h5mdfile H5Dwrite
#Dataset energy
h5mdfile H5Screate_simple type double dims 1 1
h5mdfile H5Pset_chunk dims 5 1
h5mdfile H5Dcreate2 "/observables/energy"
h5mdfile H5Dopen2 "/observables/energy"
set E [analyze energy kinetic]
h5mdfile H5_write_value value $E index 0 0
h5mdfile H5Dwrite

#############################################################
#      Integration                                          #
#############################################################
puts "\nStart integration: run $int_n_times times $int_steps steps"
inter forcecap 0

for {set i 0} { $i < $int_n_times } { incr i} {
    #puts -nonewline "run $i at time=[setmd time] "
    integrate $int_steps
    set E [analyze energy kinetic]
    
    #WRITE H5MD
    #Positions
    h5mdfile H5Dopen2 "/particles/pos"
    h5mdfile H5Dextend dims [expr $i+2] $n_part 3
    h5mdfile H5Sselect_hyperslab offset [expr $i+1] 0 0
    h5mdfile H5Screate_simple type double dims 1 $n_part 3
    for { set k 0 } { $k <= [setmd max_part] } { incr k } {
	set pos [part $k pr p]
	h5mdfile H5_write_value value [lindex $pos 0] index 0 $k 0
	h5mdfile H5_write_value value [lindex $pos 1] index 0 $k 1
	h5mdfile H5_write_value value [lindex $pos 2] index 0 $k 2
	
    }
    h5mdfile H5Dwrite
    #Energy
    h5mdfile H5Dopen2 "/observables/energy"
    h5mdfile H5Dextend dims [expr $i+2] 1
    h5mdfile H5Sselect_hyperslab offset [expr $i+1] 0
    h5mdfile H5Screate_simple type double dims 1 1
    h5mdfile H5_write_value value $E index 0 0
    h5mdfile H5Dwrite

}

h5mdfile H5Pclose
h5mdfile H5Dclose
h5mdfile H5Sclose
h5mdfile H5Gclose
h5mdfile H5Fclose
h5mdfile H5_free_memory
# terminate program
puts "\n\nFinished"
exit
