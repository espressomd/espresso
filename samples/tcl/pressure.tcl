#############################################################
#                                                           #
#  Charged Ideal Gas                                        #
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
puts "=       pressure.tcl                                  ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "pressure"
set ident "_s6"

# Visualization?  File I/O?
set vmd_output  "no"
set write_stuff "yes"

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
thermostat langevin 1.0 1.0

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
puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy]"


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

# write configurations file
if { $write_stuff=="yes" } {
    set trajectory [open "$name$ident.config" "w"]

    blockfile $trajectory write variable {box_l time_step skin}
    blockfile $trajectory write interactions
    blockfile $trajectory write integrate
    blockfile $trajectory write thermostat
    flush $trajectory
}


puts "\nStart integration: run $int_n_times times $int_steps steps, averaging the pressure over the past $int_avg measurements"
set j 0; set avg_i 0; set pc_avg 0.0; set pDH_avg 0.0
for {set i 0} { $i < $int_avg } {incr i} { lappend avg_pc 0.0; lappend avg_pDH 0.0 }
for {set i 0} { $i < $int_n_times } { incr i} {

    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }

#   derive observables
    puts -nonewline "    Run $i at time=[setmd time] "; flush stdout
    set temp [expr [analyze energy kinetic]/(([degrees_of_freedom]/2.0)*[setmd n_part])]
    set ptot [analyze pressure total]
    set p_c  [analyze pressure coulomb]
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
    if { $write_stuff=="yes" } { puts $obs_file "{ time [setmd time] } { pressure $ptot $p_c $p_DH $avg_n $tmp_pc $tmp_pDH $tmp_err }" }

#   write intermediate configuration
    if { $i%10==0 } { if { $write_stuff=="yes" } {
	blockfile $trajectory write particles
	flush $trajectory
    }
    }
}

# write end configuration
if { $write_stuff=="yes" } {
    blockfile $trajectory write particles
    blockfile $trajectory write variable { time }
    close $trajectory
}

close $obs_file

# terminate program
puts "\n\nFinished"
exit
