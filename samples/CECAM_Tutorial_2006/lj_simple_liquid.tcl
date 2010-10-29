#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
# 
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
# 
#############################################################
#                                                           #
#  Sample System 1: Lennard Jones Liquid                    #
#                                                           #
#                                                           #
#  Created:       17.03.2003 by HL                          #
#  Last modified: 18.03.2003 by HL                          #
#  Last modified: 26.11.2006 by TS                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=       Sample Script 1: lj_liquid_simple.tcl                ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "lj_liquid_simple"
set ident "_s1"

# 
set vmd_output "no"

# System parameters
#############################################################

set box_l   10.0
set density 0.5

# Interaction parameters (repulsive Lennard Jones)
#############################################################

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   $lj1_eps

# Integration parameters
#############################################################

setmd time_step 0.01
setmd skin      0.4
set temperature 1.0
set gamma 1.0

thermostat set langevin $temperature $gamma

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

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0

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
puts $obs_file "\# Time\tE_tot\tE_kin\t..."

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter ljforcecap $cap

# Warmup Integration Loop
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {

    integrate $warm_steps

    # Visualization
    if { $vmd_output=="yes" } { imd positions }

    # Warmup criterion
    set act_min_dist [analyze mindist]
    puts -nonewline "run $i at time=[setmd time] (LJ cap=$cap) min dist = $act_min_dist\r"
    flush stdout

#   write observables
    puts $obs_file "{ time [setmd time] } [analyze energy]"

#   Increase LJ cap
    set cap [expr $cap+10]
    inter ljforcecap $cap
    incr i
}

# write parameter file
# polyBlockWrite "$name$ident.set" {box_l time_step skin temp gamma } "" 

#############################################################
#      Integration                                          #
#############################################################
puts "\nStart integration: run $int_n_times times $int_steps steps"

inter ljforcecap 0

puts [analyze energy]

set j 0
for {set i 0} { $i < $int_n_times } { incr i} {


    puts -nonewline "run $i at time=[setmd time] "

    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }
#   write observables
    set energies [analyze energy]
    puts $obs_file "{ time [setmd time] } $energies"
    puts -nonewline "temp = [expr [lindex $energies 1 1]/(1.5*[setmd n_part])]\r"
    
    if {$i > [expr $int_n_times/2]} {analyze append}
    flush stdout

#   write intermediate configuration
#    if { $i%10==0 } {
#	polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
#	incr j
#    }

}

# write end configuration
polyBlockWrite "$name$ident.end" {time box_l} {id pos type}

close $obs_file

# calculate radial distribution function and write result to a file
set part_type_1 0
set part_type_2 0
set rdf_from 0.5
set rdf_to 3.5
set rdf_n_bins 100
set rdf_file [open "$name$ident.rdf" "w"]
puts $rdf_file "[analyze <rdf> $part_type_1 $part_type_2 $rdf_from $rdf_to $rdf_n_bins]"
close $rdf_file

# write standard .pdb-file, that can be read be many graphics programs (e.g. molekel)
writepdb "$name$ident.pdb"

puts "\n\nFinished"

# terminate program
exit
