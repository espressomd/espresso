#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Simulation in a spherical cell                           #
#  5 PPPs N=61 lb-elj scan                                  #
#                                                           #
#  Created:       09.04.2003 by HL                          #
#  Last modified: 25.04.2003 by HL                          #
#                                                           #
#############################################################

source /people/pckr41/limbach/Espresso/scripts/bundle.tcl

puts " "
puts "======================================================="
puts "=          Bundle Simulations: PPP 2 N 61             ="
puts "=                lb0.50elj2.00.tcl                    ="
puts "======================================================="
puts " "

#puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################


# System parameters
#############################################################

# number of Polyelectrolytes
set n_poly        2

# Polyelectrolyte backbone
set l_back        61
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
set lja_eps       2.0

set fene_r        2.0
set fene_k        7.0

set bend_k        30.0

# The real value for PPP would be around bjerrum=1.7
set bjerrum       0.5

# System identification: 
set name  "lb[format %.2f $bjerrum]elj[format %.2f $lja_eps]"
puts "$name"
set ident ""

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5

set warm_steps   100
set warm_n_times 30
set min_dist     0.9

set ci_steps     1000
set ci_n_times   200

set int_steps    1000
set int_n_times  1000

set cfg_write_steps 10000

# Absolute times:
set time     0
set ci_time  [expr $ci_steps*$ci_n_times*$time_step]
set cfg_cnt 0
# If you wish to restart, use maxtime...
set restart 1
set max_time 2000

# Other parameters
#############################################################
set tcl_precision 6
set mypi          3.141592653589793
set vmd_output    "no"
set vmd_wait      1

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

# check for existing configurations
puts "Searching for existing configs ..."
set i 0
while { $restart == 1 } {
    if { [file exists "$name$ident.[format %04d $i]" ] } {
	incr i
    } else { set restart 0 }
}
if { $i > 1 } {
    set i [expr $i-1]
    set file [open "$name$ident.[format %04d $i]" r]
    while { [blockfile $file read auto] != "eof" } {}
    close $file
    set time [setmd time]
    puts "Restart Simulation with $name$ident.[format %04d $i] at t=$time" 
    set step       [expr $i*$cfg_write_steps]
    set last_write [expr $i*$cfg_write_steps]
    set cfg_cnt [expr $i+1]
    set restart 1	  
} else { 
    set restart 0     
    set step 0
    set last_write 0
}

# Simulation box
set volume [expr $n_part/$density]
set sphere_rad [expr pow((3.0*$volume)/(4.0*$mypi),1.0/3.0)]
set  box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set center [list [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0]]
constraint sphere [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0] $sphere_rad 1.0 1.0 1.12246 1.0

if { $restart == 0 } { setmd box_l     $box_l $box_l $box_l }

setmd periodic  0 0 0
setmd time_step $time_step
setmd skin      $skin
setmd gamma     1.0
setmd temp      1.0

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
set lja_shift [expr -(pow($lja_cut,-12)-pow($lja_cut,-6))]
inter 3 3 lennard-jones $lja_eps 1.0 $lja_cut $lja_shift 0
# FENE 
inter 0 fene $fene_k $fene_r
# Stiffness
inter 1 angle $bend_k

# Particle setup
#############################################################

# Parameters for the counterion two state setuup
set f_ppp 0
for { set i 0 } {$i < [llength $c_pat] } { incr i } { set f_ppp [expr $f_ppp + [lindex $c_pat $i] ] }
set f_ppp [expr $f_ppp/double([llength $c_pat])]
set xi_bundle [expr $n_poly*$f_ppp*$bjerrum]
# the factor 0.9 is an empiric correction due too the finite length of the bundle
if { $xi_bundle > 1.0 } {
    set man_frac [expr 0.9*(1.0-1.0/$xi_bundle)]
} else { set man_frac 0 }
set man_rad  [expr sqrt($sphere_rad)]

# setup bundle cylinder
set cyl_rad [expr 1.0+sqrt(($n_hair)/($mypi*$l_back*$h_dens))]

if { $restart == 0 } { 
    bundle_backbone_setup   $n_poly $l_back $cyl_rad $c_pat $h_dist $center {0 2} 0
    bundle_counterion_two_state_setup $n_ci $sphere_rad $ci_val $center 1 $n_mono $man_rad $l_back $man_frac
    bundle_hair_setup       $n_poly $l_back $cyl_rad $h_dist $h_length $center 3 $hair_start 0
}

# Control Loop
for {set i 0} { $i < $n_part } {incr i} {
#    puts "[part $i]"
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

if { $n_part != [setmd n_part] } {
    puts "Particle number mismatch $n_part != [setmd n_part]"
}

#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident" $vmd_wait 1
}

#############################################################
#  Warmup Integration                                       #
#############################################################

if { $restart == 0 } {
    puts "\nStart warmup integration: max $warm_n_times times $warm_steps steps"
    set act_min_dist [analyze mindist]
    set cap 20  
    inter ljforcecap $cap
    # fix backbone monomers
    for {set i 0} { $i < $warm_n_times  && $act_min_dist < $min_dist } { incr i} {
	puts -nonewline "Run $i at time=[setmd time] min_dist=$act_min_dist\r"
	flush stdout
	integrate $warm_steps
	if { $vmd_output=="yes" } { imd positions }
	set act_min_dist [analyze mindist]
	set cap [expr $cap+10]
	inter ljforcecap $cap  
    }
    puts "Warmup done: Minimal distance [analyze mindist]"

} else { puts "Skip warmup (Minimal distance [analyze mindist])" }


#############################################################
#  Counterion Equilibration                                #
#############################################################
if { $time < $ci_time } {
    if { $restart == 1 } {
	set ci_n_times [expr int(($ci_time-$time)/($ci_steps*$time_step))+1]
    } 
    for {set i 0} { $i < $n_mono } {incr i} { part $i fix }
    
    puts "\nStart Counterion Equilibration: $ci_n_times times $ci_steps steps"
    # remove LJ cap
    inter ljforcecap 0
    # Start Coulomb interaction
    inter coulomb $bjerrum dh 0.0 [expr 2.0*$sphere_rad]
    
    puts "Interactions:\n [inter]"
    
    if { $restart == 0 } {
	polyBlockWrite "$name$ident.[format %04d $cfg_cnt]" {time box_l} {id pos type}
	incr cfg_cnt
    }
    for {set i 0} { $i < $ci_n_times  } { incr i} {
	puts -nonewline "Run $i at time=[setmd time]\r"
	flush stdout
	integrate $ci_steps
	if { $vmd_output=="yes" } { imd positions }
	
	set step [expr $step+$ci_steps]
	if { $step >= [expr $last_write+$cfg_write_steps] } {
	    polyBlockWrite "$name$ident.[format %04d $cfg_cnt]" {time box_l} {id pos type}
	    set last_write [expr $cfg_cnt*$cfg_write_steps]
	    incr cfg_cnt
	}
    }
    puts "Counterion Equilibration done."
} else { puts "Skip Counterion Equilibration" }


#############################################################
#      Integration                                          #
#############################################################

# free backbone monomers
for {set i 0} { $i < $n_mono } {incr i} { part $i unfix }

if { $restart == 0 } {
    puts "\nStart integration: $int_n_times times $int_steps steps"
} else {
    set int_n_times [expr int(($max_time-$time)/($int_steps*$time_step))+1]
    puts "\nStart integration: $int_n_times times $int_steps steps"
}

analyze set chains 0 $n_poly $l_back
for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "Run $i at time=[setmd time], re = [analyze re]\r"
    flush stdout
    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }

    set step [expr $step+$int_steps]
    if { $step >= [expr $last_write+$cfg_write_steps] } {
	polyBlockWrite "$name$ident.[format %04d $cfg_cnt]" {time box_l} {id pos type}
	set last_write [expr $cfg_cnt*$cfg_write_steps]
	incr cfg_cnt
    }
}

puts "\nIntegration done."

#############################################################
puts "\nFinished"