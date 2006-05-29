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
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 
#############################################################
#                                                           #
#  Sample System 6: Charged Ideal Gas                       #
#                                                           #
#                                                           #
#  Created:       27.02.2004 by HL                          #
#  Last modified: 29.02.2004 by HL                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=       Sample Script 6: pressure.tcl                 ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"


#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "pressure"
set ident "_s6"

# Visualization?  File I/O?  Retune Electrostatics?
set vmd_output  "no"
set write_stuff "yes"
set retune      "no"

# System parameters
#############################################################

# 10 000  Particles
set n_part      10000
set density     0.001
set volume      [expr $n_part/$density]
set box_l       [expr pow($volume,1./3.)]

# Interaction parameters (electrostatics)
#############################################################

set bjerrum     1
set accuracy    1.0e-04
set kappa       [expr sqrt(4*[PI]*$bjerrum*1.0*$density)]

# Integration parameters
#############################################################

setmd time_step 0.01
setmd skin      0.5
setmd gamma     1.0
setmd temp      1.0

# integration
set int_steps    100
set int_n_times  1000
set int_avg      100

# Other parameters
#############################################################
set   tcl_precision 6
setmd max_num_cells 2744


#############################################################
#  Setup System                                             #
#############################################################

# System setup
#############################################################

puts "System Setup:"
setmd box_l $box_l $box_l $box_l

# Particle setup
#############################################################

puts -nonewline "    Setting up $n_part particles at random positions... "; flush stdout
for {set i 0} { $i < $n_part } {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
    part $i pos $posx $posy $posz type 0
    if { $i % 2 == 0 } { part $i q +1.0 } else { part $i q -1.0 }
}
puts "Done ([setmd n_part] particles now)."

puts -nonewline "    Checking that total charge is zero... "; flush stdout
set total_charge 0
for {set i 0} { $i < $n_part } {incr i} { 
    set total_charge [expr $total_charge + [part $i print q]] }
if { abs($total_charge) > $accuracy } {
    puts "Failed.\nERROR: System has non-zero total charge! Aborting... "; exit }
puts "Done (found $total_charge as overall charge)."

# Interaction setup
#############################################################

puts "    Creating electrostatic interactions... "
if { $retune=="yes" } {
    puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy]"
} else {
    # inter coulomb 2.0 p3m 3.29062 64 6 0.973715 8.15277e-05
    # inter coulomb 1.0 p3m 45.9159 32 3 0.0420705 8.64147e-05
    inter coulomb 1 p3m 26.4304 32 4 0.0859363 9.39165e-05
}


#############################################################
#  Integration                                              #
#############################################################

puts "\nSimulate $n_part particles in a cubic simulation box of length $box_l at density $density"
puts "Interactions:\n[inter]"
set  act_min_dist [analyze mindist]
puts "Starting with a minimal distance $act_min_dist and a screening-length of kappa=$kappa "

# prepare vmd connection
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident$"
}

#open Observable file
if { $write_stuff=="yes" } { 
    set obs_file [open "$name$ident.obs" "w"]
    puts $obs_file "\# System: $name$ident"
    puts $obs_file "\# Time\tE_tot\tE_kin\t..."
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

# write parameter file
if { $write_stuff=="yes" } { polyBlockWrite "$name$ident.set" {box_l time_step skin temp gamma } "" }


puts "\nStart integration: run $int_n_times times $int_steps steps, averaging the pressure over the past $int_avg measurements"
set j 0; set avg_i 0; set pc_avg 0.0; set pDH_avg 0.0
for {set i 0} { $i < $int_avg } {incr i} { lappend avg_pc 0.0; lappend avg_pDH 0.0 }
for {set i 0} { $i < $int_n_times } { incr i} {

    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }

#   derive observables
    puts -nonewline "    Run $i at time=[setmd time] "; flush stdout
    set ener [analyze energy]
    set pres [analyze pressure]
    set temp [expr [lindex [lindex $ener 1] 1]/(1.5*[setmd n_part])]
    set ptot [expr [lindex [lindex $pres 0] 1]]
    set p_c  [expr [lindex [lindex $pres 2] 1]]
    set p_DH [expr -$temp*pow($kappa,3)/(24*[PI])]
    puts -nonewline "T=$temp; p=$ptot, p_c=$p_c, p_DH=$p_DH "; flush stdout

#   average and check derivations
    set ind_avg [expr $avg_i % $int_avg]
    set pc_avg  [expr $pc_avg  - [lindex $avg_pc  $ind_avg] + $p_c ]
    set pDH_avg [expr $pDH_avg - [lindex $avg_pDH $ind_avg] + $p_DH]
    set avg_pc  [lreplace $avg_pc  $ind_avg $ind_avg $p_c ]
    set avg_pDH [lreplace $avg_pDH $ind_avg $ind_avg $p_DH]
    incr avg_i; if { $avg_i >= $int_avg } { set avg_n $int_avg } else { set avg_n $avg_i }
    set tmp_pc [expr $pc_avg/$avg_n]; set tmp_pDH [expr $pDH_avg/$avg_n]; set tmp_err [expr ($pc_avg-$pDH_avg)/(0.01*$pDH_avg)]
    puts -nonewline "(<p_c>=$tmp_pc, <p_DH>=$tmp_pDH => $tmp_err\% error) \r"; flush stdout

#   write observables
    if { $write_stuff=="yes" } { puts $obs_file "{ time [setmd time] } { pressure $ptot $p_c $p_DH $avg_n $tmp_pc $tmp_pDH $tmp_err } $ener $pres" }

#   write intermediate configuration
    if { $i%10==0 } { if { $write_stuff=="yes" } { polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type} }; incr j }
}

# write end configuration
if { $write_stuff=="yes" } { polyBlockWrite "$name$ident.end" {time box_l} {id pos type} }

close $obs_file

# terminate program
puts "\n\nFinished"
exit
