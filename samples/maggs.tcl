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

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	set f [open $errf "w"]
	puts $f "not compiled in: $feature"
	close $f
	exit -42
    }
}

puts "----------------------------------------"
puts "- Testcase maggs.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"
require_feature "ELECTROSTATICS"

set use_warmup "yes"

# System size
set density  0.07
set temp     1.0
set num_part 500

# Tuning parameters
set time_step   0.01
set skin        0.3
set max_cells   16


# Interaction parameters
#############################################################

# Lennard-Jones
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]

# Coulomb
set bjerrum       5.0
set f_mass        0.01
set mesh          24

# Integration parameters
#############################################################

set gamma     1.5

# warmup integration (with capped LJ potential)
set warm_steps         [expr int(1./$time_step)]
set warm_n_times       40
set warm_n_min         10
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    2000
set int_n_times  50

set tcl_precision 7

#############################################################
#  Setup System                                             #
#############################################################

# setup new configuration
set box_l [expr pow( ($num_part/$density), (1.0/3.0) )]
setmd box_l $box_l $box_l $box_l

for {set i 0} { $i < $num_part } {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
    
    part $i pos $posx $posy $posz
    if {[expr $i % 2 == 0]} { part $i q +1.0 type 1} else { part $i q -1.0 type 0}
    
}

# Settings 
inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
inter 1 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
inter 1 1 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0

set max_cells3d [expr $max_cells*$max_cells*$max_cells]
setmd max_num_cells $max_cells3d

setmd time_step $time_step
setmd skin      $skin      
thermostat langevin $temp $gamma
integrate 0

#############################################################
#  Warmup Integration                                       #
#############################################################
set act_min_dist [analyze mindist]

set cap 20
inter forcecap $cap
set i 0
while { $i < $warm_n_times && ( $act_min_dist < $min_dist || $i < $warm_n_min ) } {
    integrate $warm_steps
    set act_min_dist [analyze mindist]
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}
inter forcecap 0

# maggs requires domain decompostion with no verlet lists
cellsystem domain_decomposition -no_verlet_list
inter coulomb $bjerrum maggs $f_mass $mesh 0.

set act_min_dist [analyze mindist]

set i 0
while { $i < $warm_n_times && ( $act_min_dist < $min_dist || $i < $warm_n_min ) } {
    integrate $warm_steps
    set act_min_dist [analyze mindist]
    incr i
}

#############################################################
#      Integration                                          #
#############################################################

set accepted_error 0.1

set r_bins [expr $num_part/10]

#open file with correct RDF
set f_correct  [open "maggs_correct_rdf.dat" r]
for {set i 0} { $i < $r_bins } {incr i} {
    gets $f_correct in_string
    set correct_rdf($i) [lindex $in_string 1]
}
close $f_correct

for {set i 0} { $i < $int_n_times } { incr i} {
    integrate $int_steps  
    set rdf_pp [analyze rdf 1 1 0.5 7.0 $r_bins]
    # Add up rdf values to result
    if { $i == 0 } {
	for { set j 0 } {$j < $r_bins} {incr j } {
	    set bin($j)  "[lindex [lindex [lindex $rdf_pp 1] $j] 0]"
	    set result_pp($j) [expr [lindex [lindex [lindex $rdf_pp 1] $j] 1]]
	}	
    } else {
	for { set j 0 } {$j < $r_bins} {incr j } {
	    set result_pp($j) [expr $result_pp($j) + [expr [lindex [lindex [lindex $rdf_pp 1] $j] 1]]]
	}
    }
    # compare test results with correct data
    set max_error $accepted_error
    for { set j 0 } {$j < $r_bins} {incr j } {
	set temp_val [expr $result_pp($j)/[expr $i + 1]]
	set temp_error [expr abs($correct_rdf($j)-$temp_val)]
	if { $temp_error > $max_error } { 
	    set max_error $temp_error 
	}
    } 
    
    puts -nonewline "done [expr $i+1] at time=[setmd time] with error=$max_error \r"
    flush stdout
    
}
if { $max_error > $accepted_error } {
    puts "\nRDF error is too large"
} else {
    puts "\nMEMD seems to have no errors"
} 


