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
puts "=       lj_liquid.tcl                                 ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "lj_liquid"
set ident "_s1"

# 
set vmd_output "yes"

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
set int_n_times  100

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

puts "Simulate $n_part particles in a cubic simulation box "
puts "[setmd box_l] at density $density"
puts "Interactions:\n[inter]"
set act_min_dist [analyze mindist]
puts "[part 0]"
puts "[part 1]"
puts "Start with minimal distance $act_min_dist"

setmd max_num_cells 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# prepare vmd connection
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident$"
}

#open Observable file
set obs_file [open "$name$ident.obs" "w"]
puts $obs_file "\# System: $name$ident"
puts $obs_file "\# Time\tE_tot\tE_kin\tE_LJ"
flush $obs_file

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter forcecap $cap

# Warmup Integration Loop
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {

    integrate $warm_steps

    # Visualization
    if { $vmd_output=="yes" } { imd positions }

    # Warmup criterion
    set act_min_dist [analyze mindist]
    puts "run $i at time=[setmd time] (LJ cap=$cap) min dist = $act_min_dist\r"
    flush stdout

#   write observables
    puts $obs_file "[setmd time] [analyze energy total] [analyze energy kinetic] [analyze energy nonbonded 0 0]"

#   Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}

# Just to see what else we may get from the c code
puts "\nro variables:"
puts "cell_grid     [setmd cell_grid]" 
puts "cell_size     [setmd cell_size]" 
puts "local_box_l   [setmd local_box_l]" 
puts "max_cut       [setmd max_cut]" 
puts "max_part      [setmd max_part]" 
puts "max_range     [setmd max_range]" 
puts "max_skin      [setmd max_skin]" 
puts "n_nodes       [setmd n_nodes]" 
puts "n_part        [setmd n_part]" 
puts "n_part_types  [setmd n_part_types]" 
puts "periodicity   [setmd periodicity]" 
puts "transfer_rate [setmd transfer_rate]" 
puts "verlet_reuse  [setmd verlet_reuse]" 

# write configurations file
set trajectory [open "$name$ident.config" "w"]

blockfile $trajectory write variable {box_l time_step skin}
blockfile $trajectory write interactions
blockfile $trajectory write integrate
blockfile $trajectory write thermostat
flush $trajectory

#############################################################
#      Integration                                          #
#############################################################
puts "\nStart integration: run $int_n_times times $int_steps steps"

inter forcecap 0

puts [analyze energy]

set j 0
for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time] "

    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }
#   write observables
    set energies [analyze energy]
    puts $obs_file "[setmd time] [analyze energy total] [analyze energy kinetic] [analyze energy nonbonded 0 0]"
    flush $obs_file
    puts -nonewline "temp = [expr [lindex $energies 1 1]/(([degrees_of_freedom]/2.0)*[setmd n_part])]\r"
    flush stdout

#   write intermediate configuration
    if { $i%10==0 } {
	blockfile $trajectory write particles
	flush $trajectory
	incr j
    }
}

# write end configuration
blockfile $trajectory write particles
blockfile $trajectory write variable { time }
close $trajectory

close $obs_file

# terminate program
puts "\n\nFinished"
exit
