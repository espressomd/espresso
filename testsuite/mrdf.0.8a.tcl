#!/bin/sh

#############################################################
#                                                           #
#  Test 1 :       Pair distribution function (Maggs)        #
#                                                           #
#                                                           #
#  Created:       06.05.2004 by Beemaster                   #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                maggs.rdf.tcl                        ="
puts "======================================================="
puts " "

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "mrdf"
set ident ".0.8a"

set use_warmup "yes"

# System size
set density  0.07
set temp     1.0
set num_part 500

# Tuning parameters
set time_step   0.01
set skin        0.3
set max_cells   16


set vmd_output    "no"
set vmd_wait      3

# Interaction parameters
#############################################################

# Lennard-Jones
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [expr 0.25*$lj1_eps]

# Coulomb
set bjerrum       5.0
set f_mass        0.01
set mesh          24

# Integration parameters
#############################################################

set gamma     1.5

# warmup integration (with capped LJ potential)
set warm_steps         [expr int(1./$time_step)]
set field_warm_steps   [expr 5*$warm_steps]
set warm_n_times       100
set warm_n_min         40
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9

# integration
set int_steps    [expr int(10/$time_step)]
set int_n_times  2000

set tcl_precision 7

set system_name "$name$ident"
puts "Simulating $system_name "

#############################################################
#  Setup System                                             #
#############################################################

# setup new configuration
set setup_time [lindex [time {
    set box_l [expr pow( ($num_part/$density), (1.0/3.0) )]
    setmd box_l $box_l $box_l $box_l
    
    for {set i 0} { $i < $num_part } {incr i} {
	set posx [expr $box_l*[t_random]]
	set posy [expr $box_l*[t_random]]
	set posz [expr $box_l*[t_random]]

	part $i pos $posx $posy $posz
	if {[expr $i % 2 == 0]} { part $i q +1.0 type 1} else { part $i q -1.0 type 0}
	
    }
} ] 0]
puts "Particle setup done ($setup_time microseconds)"

# Settings 
set init_time [lindex [time {

    inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
    inter 1 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
    inter 1 1 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0

    set max_cells3d [expr $max_cells*$max_cells*$max_cells]
    setmd max_num_cells $max_cells3d
    
    setmd time_step $time_step
    setmd skin      $skin      
    setmd gamma     $gamma    
    setmd temp      $temp    
    integrate 0
} ] 0]
puts "Initialization done ($init_time  microseconds; dt=[setmd time_step], skin=[setmd skin], gamma=[setmd gamma], temp=[setmd temp])"


#############################################################
#  Warmup Integration                                       #
#############################################################
if { ![file exists "$system_name.warmup"] && $use_warmup == "yes" } {
    set warmup_time [lindex [time {
	set act_min_dist [analyze mindist]
	puts "\nStart warmup integration I:"
	puts "At maximum $warm_n_times times $warm_steps steps (minimum: $warm_n_min)"
	puts "Stop if minimal distance is larger than $min_dist"
	
	set cap 20
	inter ljforcecap $cap
	
	set i 0
	while { $i < $warm_n_times && ( $act_min_dist < $min_dist || $i < $warm_n_min ) } {
	    puts -nonewline "run [expr $i+1] at time=[setmd time] (LJ cap=$cap) "
	    
	    integrate $warm_steps

	    set act_min_dist [analyze mindist]
	    puts -nonewline "minimal distance = $act_min_dist\r"
	    flush stdout
	    set cap [expr $cap+10]
	    inter ljforcecap $cap
	    incr i
	}
	inter ljforcecap 0
    } ] 0]
    set warm_totsteps [expr $i*$warm_steps]

puts "\nPrepare integration: (tune parameters - wait...)"
inter ljforcecap 0
# maggs requires domain decompostion with no verlet lists
cellsystem domain_decomposition -no_verlet_list
puts "[inter coulomb $bjerrum maggs $f_mass $mesh 0.]"
puts "Interactions are now: {[inter]}"
    set warmup_time [expr $warmup_time + [lindex [time {
	set act_min_dist [analyze mindist]
	puts "\nStart warmup integration II:"
	puts "At maximum $warm_n_times times $field_warm_steps steps (minimum: $warm_n_min)"
	puts "Stop if minimal distance is larger than $min_dist"
	
	set i 0
	while { $i < $warm_n_times && ( $act_min_dist < $min_dist || $i < $warm_n_min ) } {
	    puts -nonewline "run [expr $i+1] at time=[setmd time] (full LJ, full ES) "
	    
	    integrate $field_warm_steps
	    set act_min_dist [analyze mindist]
	    puts -nonewline "minimal distance = $act_min_dist\r"
	    flush stdout
	    incr i
	}
    } ] 0]]
    set warm_totsteps [expr $warm_totsteps + $i*$warm_steps]
}



# Variable check point
puts "\nbox_l         [setmd box_l]" 
puts "cell_grid     [setmd cell_grid]" 
puts "cell_size     [setmd cell_size]" 
puts "max_cut       [setmd max_cut]" 
puts "max_range     [setmd max_range]" 
puts "max_skin      [setmd max_skin]" 
puts "n_part        [setmd n_part]" 
puts "n_part_types  [setmd n_part_types]" 
puts "verlet_reuse  [setmd verlet_reuse]" 


#############################################################
#      Integration                                          #
#############################################################

puts "\nStart integration: run $int_n_times times $int_steps steps"
set r_bins [expr $num_part/10]

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run [expr $i+1] at time=[setmd time] \r"
    flush stdout
    
    integrate $int_steps  

    set rdf_pp [analyze rdf 1 1 0.5 7.0 $r_bins]
    set rdf_pn [analyze rdf 1 0 0.5 7.0 $r_bins]
# Add up rdf values to result
    if { $i == 0 } {
	for { set j 0 } {$j < $r_bins} {incr j } {
	    set bin($j)  "[lindex [lindex [lindex $rdf_pp 1] $j] 0]"
	    set result_pp($j) [expr [lindex [lindex [lindex $rdf_pp 1] $j] 1]]
	    set result_pn($j) [expr [lindex [lindex [lindex $rdf_pn 1] $j] 1]]
	}	
    } else {
	for { set j 0 } {$j < $r_bins} {incr j } {
	    set result_pp($j) [expr $result_pp($j) + [expr [lindex [lindex [lindex $rdf_pp 1] $j] 1]]]
	    set result_pn($j) [expr $result_pn($j) + [expr [lindex [lindex [lindex $rdf_pn 1] $j] 1]]]
	}
    }
    # write averaged results to file
    set f [open "$system_name.pp.rdf" w]
    for {set j 0} { $j < $r_bins } { incr j} {
	set temp_val [expr $result_pp($j)/[expr $i + 1]]
	puts $f [format "%.3e %.6e" [expr $bin($j)] $temp_val]
    }
    flush $f
    close $f
    set f [open "$system_name.pn.rdf" w]
    for {set j 0} { $j < $r_bins } { incr j} {
	set temp_val [expr $result_pn($j)/[expr $i + 1]]
	puts $f [format "%.3e %.6e" [expr $bin($j)] $temp_val]
    }
    flush $f
    close $f
}

puts "\nFinished"
