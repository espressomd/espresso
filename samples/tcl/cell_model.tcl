#############################################################
#                                                           #
#  Simulation in a spherical cell                           #
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

source ../scripts/bundle.tcl

puts " "
puts "======================================================="
puts "=                 cell_model.tcl                      ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "bundle"
set ident "_t5"

# System parameters
#############################################################

# number of Polyelectrolytes
set n_poly        3

# Polyelectrolyte backbone
set l_back        20
set c_pat         {0 1 1}

# Hairs
set h_dist        3
set h_length      3
set h_dens        0.5

# Counterions
set ci_val       -1

# Spherical simulation cell
set density       0.0005

# Interaction parameters
#############################################################

set ljr_cut       1.12246204831
set ljr_eps       1.0

set lja_cut       2.5
set lja_eps       1.7

set fene_r        2.0
set fene_k        7.0

set bend_k        30.0

# The real value for PPP would be around bjerrum=1.7
set bjerrum       3.0
set accuracy      1.0e-4

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5

set warm_steps   100
set warm_n_times 30
set min_dist     0.9

set ci_steps     1000
set ci_n_times   10

set int_steps    1000
set int_n_times  1000

# Other parameters
#############################################################
set tcl_precision 6
set mypi          3.141592653589793
set vmd_output    "no"
set vmd_wait      3

#############################################################
#  Setup System                                             #
#############################################################

# Particle numbers
# Backbone:      id: 0 - (n_mono-1)
# Counterions:   id: n_mono - (hair_start-1)
# Hairs:         id: hair_start - (n_part-1)  
#############################################################
set n_mono [expr $n_poly * $l_back]
set n_hair [expr $n_poly * $h_length*int(($l_back/(1.0*$h_dist))+0.999)]
set n_ci 0
set n_charges 0
for {set i 0} { $i < $l_back } {incr i} {
    set pat_num [expr $i%[llength $c_pat]]
    if { [lindex $c_pat $pat_num] != 0 } { 
	set n_ci [ expr $n_ci + [lindex $c_pat $pat_num] ] 
	incr n_charges
    }
}
set n_ci [expr $n_ci*$n_poly]
set n_charges [expr $n_charges*$n_poly]
set back_charge $n_ci
set n_ci [expr -$n_ci/$ci_val]
set n_charges [expr $n_charges+$n_ci]
if { $n_ci < 0 } { puts "Impossible ci valency" exit }
set hair_start [expr $n_mono+$n_ci] 
set n_part [expr $n_mono + $n_ci + $n_hair]

# Simulation box
set volume [expr $n_part/$density]
set sphere_rad [expr pow((3.0*$volume)/(4.0*$mypi),1.0/3.0)]
set  box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set center [list [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0]]
constraint sphere center [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0] radius $sphere_rad type 4 

setmd box_l     $box_l $box_l $box_l
setmd periodic  0 0 0
setmd time_step $time_step
setmd skin      $skin
thermostat langevin 1.0 1.0
cellsystem NSQUARE

# Interaction setup
#############################################################

# type 0 uncharged backbone, type 1 counterion, type 2 charged backbone, type 3 hair
set max_part_type 3

# repulsive LJ for all
set ljr_shift  $ljr_eps
for {set i 0} { $i <= $max_part_type } {incr i} {
    for {set j $i} { $j <= $max_part_type } {incr j} {
    inter $i $j lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
    }
}
# attractive LJ for hairs
set lja_shift [expr -4.0*$lja_eps*(pow($lja_cut,-12)-pow($lja_cut,-6))]
inter 3 3 lennard-jones $lja_eps 1.0 $lja_cut $lja_shift 0
# FENE 
inter 0 fene $fene_k $fene_r
# Stiffness
inter 1 angle $bend_k

#sphere constraint interactions
for {set i 0} { $i <= $max_part_type } {incr i} {
    inter $i 4 lennard-jones 1.0 1.0 1.12246 1.0 0.
}

# Particle setup
#############################################################

# Parameters for the counterion two state setuup
set f_ppp 0
for { set i 0 } {$i < [llength $c_pat] } { incr i } { set f_ppp [expr $f_ppp + [lindex $c_pat $i] ] }
set f_ppp [expr $f_ppp/double([llength $c_pat])]
set xi_bundle [expr $n_poly*$f_ppp*$bjerrum]
# the factor 0.9 is an empiric correction due too the finite length of the bundle
set man_frac [expr 0.9*(1.0-1.0/$xi_bundle)]
set man_rad  [expr sqrt($sphere_rad)]

# setup bundle cylinder
set cyl_rad [expr 1.0+sqrt(($n_hair)/($mypi*$l_back*$h_dens))]

bundle_backbone_setup   $n_poly $l_back $cyl_rad $c_pat $h_dist $center {0 2} 0
bundle_counterion_two_state_setup $n_ci $sphere_rad $ci_val $center 1 $n_mono $man_rad $l_back $man_frac
bundle_hair_setup       $n_poly $l_back $cyl_rad $h_dist $h_length $center 3 $hair_start 0

# Control Loop
for {set i 0} { $i < $n_part } {incr i} {
    puts "[part $i]"
}

# Checks
if { $l_back > [expr 1.8*$sphere_rad] } { 
    puts "Sphere too small!" 
    exit 
}

# Status report                                             #
#############################################################
puts "Simulate $n_part particles in a spherical cell with radius $sphere_rad"
puts "    Micelle:         $n_poly PPPs of length $l_back. Cylinder radius: $cyl_rad"
puts "    Charge pattern:  { $c_pat }, total backbone charge:  $back_charge"
puts "    Hairs:           Length $h_length, distanvce on backbone: $h_dist monomer."
puts "    Counterions:     $n_ci counterions with valency $ci_val."
puts "    CI Setup:        Fraction $man_frac in cylinder with radius $man_rad around the bundle"
puts "    Particles:       Total: $n_part, Charged: $n_charges, Uncharged [expr $n_part-$n_charges]"
puts "\nConstraints:\n[constraint]"

# write configurations file
set trajectory [open "$name$ident.config" "w"]

blockfile $trajectory write variable {box_l time_step skin}
blockfile $trajectory write interactions
blockfile $trajectory write thermostat
flush $trajectory

#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident" $vmd_wait 1
}

#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration: max $warm_n_times times $warm_steps steps"
integrate 0
set act_min_dist [analyze mindist]
set cap 20  
inter forcecap $cap
# fix backbone monomers
for {set i 0} { $i < $warm_n_times  && $act_min_dist < $min_dist } { incr i} {
    puts -nonewline "Run $i at time=[setmd time] min_dist=$act_min_dist\r"
    flush stdout
    integrate $warm_steps

    blockfile $trajectory write particles
    flush $trajectory

    if { $vmd_output=="yes" } { imd positions }
    set act_min_dist [analyze mindist]
    set cap [expr $cap+10]
    inter forcecap $cap  
}

puts "Warmup done: Minimal distance [analyze mindist]"

#############################################################
#  Counterion Equilibration                                #
#############################################################
for {set i 0} { $i < $n_mono } {incr i} { part $i fix }

puts "\nStart Counterion Equilibration: $ci_n_times times $ci_steps steps"
# remove LJ cap
inter forcecap 0
# Start Coulomb interaction
inter coulomb $bjerrum dh 0.0 [expr 2.0*$sphere_rad]

puts "Interactions:\n [inter]"

for {set i 0} { $i < $ci_n_times  } { incr i} {
    puts -nonewline "Run $i at time=[setmd time]\r"
    flush stdout
    integrate $ci_steps

    blockfile $trajectory write particles
    flush $trajectory
   
    if { $vmd_output=="yes" } { imd positions }
}

puts "Counterion Equilibration done."

#############################################################
#      Integration                                          #
#############################################################

# free backbone monomers
for {set i 0} { $i < $n_mono } {incr i} { part $i unfix }

puts "\nStart integration: $int_n_times times $int_steps steps"

analyze set chains 0 $n_poly $l_back
for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "Run $i at time=[setmd time], re = [analyze re]\r"
    flush stdout
    integrate $int_steps

    blockfile $trajectory write particles
    flush $trajectory

    if { $vmd_output=="yes" } { imd positions }
}

puts "\nIntegration done."

#############################################################
puts "\nFinished"

exit
