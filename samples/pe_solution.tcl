#############################################################
#                                                           #
#  Polyelectrolyte Solution                                 #
#  (Poor Solvent)                                           #
#                                                           #
#############################################################
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
#  

puts " "
puts "======================================================="
puts "=       Sample script 3: pe_solution.tcl              ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "pe_solution"
set ident "_s3"

# System parameters
#############################################################

# Number of polymers
set n_polymers 5
# length of polymers
set l_polymers 190
# distance between charged monomers
set cm_distance 3
# system density
set density 0.00002

puts [format "%d %d %d" $n_polymers $l_polymers $cm_distance]

# Interaction parameters
#############################################################

# Lennard Jones
set n_part_types  3
# repulsive LJ
set lj1_cut       1.122462048309373
set lj1_epsilon   1.0
# LJ with attractive tail
set lj2_cut       2.5
set lj2_epsilon   1.75

# FENE
set fene_cut      2.0
set fene_k        7.0

# Coulomb
set bjerrum       1.5
set accuracy      1.0e-4

# Integration parameters
#############################################################

setmd time_step 0.0125
setmd skin      0.4
thermostat langevin 1.0 1.0

# warmup integration (with capped LJ potential)
set warm_steps   200
set warm_n_times 30
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    200
set int_n_times  5000


# More general Parameters
#############################################################

set tcl_precision 6
set vmd_output    "no"

# Initial Bond length
set bond_length   1.0
# valency of charged monomers
set cm_valency    1.0
# valency of counterions
set ci_valency    -1.0

puts "Simulate the following polyelectrolyte solution:"
puts "number of polymers       : $n_polymers"
set n_monomers [expr ($n_polymers*$l_polymers)]
puts "number of monomers       : $n_monomers"
set polymer_charge [expr (($l_polymers+$cm_distance-1)/$cm_distance)]
puts "polymer_charge           : $polymer_charge"
set n_counterions [expr -int($polymer_charge*$n_polymers/$ci_valency)]
puts "n_counterions            : $n_counterions"
set npart [expr $n_counterions+$n_monomers]
puts "total number of particles: $npart"
set box_length [expr pow(($npart/$density),(1.0/3.0))]
puts "cubic box length         : $box_length"

# Dependent Parameters
set lj1_shift  [calc_lj_shift 1 $lj1_cut]
set lj2_shift  [calc_lj_shift 1 $lj2_cut]

# Setup
#############################################################

# simulation box
setmd box_l $box_length $box_length $box_length 

# fene
inter 0 fene $fene_k $fene_cut

# pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj1_epsilon 1.0 $lj1_cut $lj1_shift 0.0
    }
}

puts "Setup Particles (wait...)"
# polymers
polymer $n_polymers $l_polymers $bond_length mode PSAW charge $cm_valency distance $cm_distance types 0 1 FENE 0
# counterions
counterions $n_counterions charge $ci_valency type 2
#puts "[part]"
set n_part [setmd n_part]
for { set i 0 } { $i < $n_part } { incr i } {
    puts "[part $i]"
}

set act_min_dist [analyze mindist]
puts "Placed [setmd n_part] particles with minimal distance: $act_min_dist"

#############################################################
#  prepare vmd connection                                       #
#############################################################
    if { $vmd_output=="yes" } {
	puts -nonewline "\nWrite psf and pdb for VMD connection... "; flush stdout
	writepsf "$name$ident.psf" ; writepdb "$name$ident.pdb"
	puts -nonewline "Output created, establishing link... "; flush stdout
	for {set port 10000} { $port < 65000 } { incr port } {
	    catch {imd connect $port} res
	    if {$res == ""} break
	}
	if { $port==65000 } { puts "Failed." } else { puts "Done (now listening at port $port)." 
	   # puts "    What you have to do now for a VMD connection:"
	   # puts "    (1) Start vmd in current directory (best before running the script)."
	   # puts "    (2) Enter on vmd command line: 'mol load psf $name$ident.psf pdb $name$ident.pdb'"
	    set HOSTNAME [exec hostname]
	   # puts "    (3) Enter on vmd command line: 'imd connect $HOSTNAME $port'"
	   # puts "    (4) To have the chains coloured individually, set 'Coloring-Method' to 'ResName' in the 'Graphics'-menu"
	    imd listen 0
	    set vmdout_file [open "vmdoutput.script" "w"]
	    puts $vmdout_file "mol load psf $name$ident.psf pdb $name$ident.pdb"
	    puts $vmdout_file "rotate stop"
	    puts $vmdout_file "imd connect $HOSTNAME $port"
	    close $vmdout_file
	    exec vmd -e vmdoutput.script &
	}
    }

#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter forcecap $cap 
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    set time [format "%8f" [setmd time]]
    puts -nonewline "run $i at time=$time (LJ cap=$cap)\r "
    flush stdout

    puts [analyze energy]
    integrate $warm_steps

    if { $vmd_output=="yes" } { imd positions }

    set act_min_dist [analyze mindist]
    puts -nonewline "minimal distance = $act_min_dist\r"
    flush stdout
#   Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}

puts "\nWarmup done: Minimal distance [analyze mindist]"
puts "Verlet reuses: [setmd verlet_reuse]"

#############################################################
#      Integration                                          #
#############################################################

puts "Prepare integration: (tune parameters - wait...)"
inter forcecap 0
setmd time 0.0
# Set attractive LJ for monomer interactions (types 0 and 1
for {set ia1 0} { $ia1 < 2 } { incr ia1 } {
    for {set ia2 0} { $ia2 < 2 } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj2_epsilon 1.0 $lj2_cut $lj2_shift 0.0
    }
}
#puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy mesh 16]"
puts "[inter coulomb $bjerrum scafacos_p3m]"
puts "Interactions are now: {[inter]}"

#      Write blockfiles for restart
#############################################################
set trajectory [open "$name$ident.config" "w"]

blockfile $trajectory write variable {box_l time_step skin}
blockfile $trajectory write interactions
blockfile $trajectory write integrate
blockfile $trajectory write thermostat
flush $trajectory

# prepare observable output
set obs_file [open "$name$ident.obs" "w"]
puts $obs_file "\#$name$ident: Observables"
puts $obs_file "\#Time     R_E         R_G         R_H         E_TOT       E_KIN       E_C_k"
analyze set chains 0 $n_polymers $l_polymers
set re_full [analyze re]
set re [lindex $re_full 0]
set j 0
puts "\nStart Main integration: $int_n_times times $int_steps steps"
for {set i 0} { $i <= $int_n_times } { incr i} {
    set time [setmd time]
    puts -nonewline "run $i at time=$time, R_E = $re_full (soll 15.0, deviation [expr $re/15.0])\r"
    flush stdout

    integrate $int_steps

    set f [open "$name$ident\_dist.dat" w]
    puts -nonewline $f "\#"
    puts $f "[analyze distribution 2 {0 1} 1.0 300.0 30 1 1]"
    flush $f
    close $f

    set energy [analyze energy]
    set re_full [analyze re]
    set re [lindex $re_full 0]
    set rg [lindex [analyze rg] 0]
    set rh [lindex [analyze rh] 0]

    puts $obs_file [format "%.3e %.5e %.5e %.5e %.5e %.5e %.5e" $time $re $rg $rh [lindex [lindex $energy 0] 1] [lindex [lindex $energy 1] 1] [lindex [lindex $energy [expr [llength $energy]-1]] 1] ]
    flush $obs_file
    if { $vmd_output=="yes" } { imd positions }
    #   write intermediate configuration
    if { $i%50==0 } {
	blockfile $trajectory write particles
	flush $trajectory
	incr j
    }
}

puts "\nIntegration done."

#############################################################

puts "\nFinished"
exit
