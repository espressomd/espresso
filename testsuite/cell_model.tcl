#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Template for Bundle Simulations                          #
#                                                           #
#                                                           #
#  Created:       23.01.2003 by HL                          #
#  Last modified: 23.01.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                bundle_temp.tcl                      ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "bundle"
set ident "temp"

# System parameters
#############################################################

set n_poly        4

set l_back        40
set c_dist        3
set c_num         2

set h_dist        3
set h_length      3
set h_dens        0.5

set sphere_rad    25.0
#set density       0.001

# Interaction parameters
#############################################################

set ljr_cut       1.12246204831
set ljr_eps       1.0

set lja_cut       2.5
set lja_eps       1.7

set fene_r        2.0
set fene_k        7.0

set bend_k        10.0

set bjerrum       3.0
set accuracy      1.0e-4

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5

set warm_steps   100
set warm_n_times 30
set min_dist     0.9

set ci_steps     100
set ci_n_times   100

set int_steps    100
set int_n_times  100

# Other parameters
#############################################################
set tcl_precision 6
set mypi          3.141592653589793
set vmd_output    "yes"
set vmd_wait      10

#############################################################
#  Setup System                                             #
#############################################################

set box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set center      [expr $box_l/2.0]
constraint sphere $center $center $center $sphere_rad 1.0 1.0 1.12246 1.0

setmd box_l     $box_l $box_l $box_l
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
set lja_shift [expr -4.0*$lja_eps*(pow($lja_cut,-12)-pow($lja_cut,-6))]
inter 3 3 lennard-jones $lja_eps 1.0 $lja_cut $lja_shift 0
# FENE 
inter 0 fene $fene_k $fene_r
# Stiffness
inter 1 angle $bend_k

# Particle setup
#############################################################

# particle numbers
set n_mono [expr $n_poly * $l_back]
set n_hair [expr $n_poly * $h_length*int(($l_back/(1.0*$h_dist))+0.999)]
set n_ci   [expr $n_poly * $c_num * ($l_back/$c_dist)]
set tmp [expr ($l_back%$c_dist)-$c_dist+$c_num]
if { $tmp > 0 } { set n_ci [expr $n_ci + $n_poly*$tmp] }
set n_part [expr $n_mono + $n_ci + $n_hair]

# setup bundle cylinder
set cyl_rad [expr 1.0+sqrt(($n_hair)/($mypi*$l_back*$h_dens))]

set part_id 0
set s_rad2 [expr $sphere_rad*$sphere_rad]

# Backbone setup
set off1  [expr -($h_dist/2.0-0.5)]
set alpha [expr -(2.0*$mypi/$n_poly)]
for {set i 0} { $i < $n_poly } {incr i} {
    set posx [expr $center - $l_back/2.0 + $off1 + ($i%$h_dist)]
    set posy [expr $center + $cyl_rad*sin($i*$alpha)]
    set posz [expr $center + $cyl_rad*cos($i*$alpha)]
    for {set j 0} { $j < $l_back } {incr j} {
	part $part_id pos $posx $posy $posz type 0 fix
	if { ($j%$c_dist) >= ($c_dist-$c_num) } {
	    part $part_id type 2 q 1.0
	}
	if { $j > 0 } {
	    part [expr $part_id-1] bond 0 $part_id
	}
	if { $j > 1 } {
	    part [expr $part_id-1] bond 1 [expr $part_id-2] [expr $part_id]
	}
	set posx [expr $posx+1.0]
	incr part_id
    }
}

# Counterions 
for {set i 0} { $i < $n_ci } {incr i} {
    set dist [expr $s_rad2 + 10.0]
    while { $dist > [expr $s_rad2-(2.0*$sphere_rad-1.0)] } {
        set posx [expr 2*$sphere_rad*[t_random]-$sphere_rad]
        set posy [expr 2*$sphere_rad*[t_random]-$sphere_rad]
        set posz [expr 2*$sphere_rad*[t_random]-$sphere_rad]
        set dist [expr $posx*$posx + $posy*$posy + $posz*$posz]
    }
    part $part_id pos [expr $posx+$center] [expr $posy+$center] [expr $posz+$center] type 1 q -1.0
    incr part_id
}

# Hair setup
set back_id 0
for {set i 0} { $i < $n_poly } {incr i} {
    for {set j 0} { $j < $l_back } {incr j} {
	if { $j%$h_dist == 0 } {
	    set back_pos [part $back_id print pos]
	    set posx [lindex $back_pos 0]
	    for {set k 1} { $k <= $h_length } {incr k} {
		set posy [expr $center + ($cyl_rad-$k)*sin($i*$alpha)]
		set posz [expr $center + ($cyl_rad-$k)*cos($i*$alpha)]
		part $part_id pos $posx $posy $posz type 3
		if { $k == 1 } {
		    part $part_id bond 0 $back_id
		} else {
		    part $part_id bond 0 [expr $part_id-1]
		}
		incr part_id
	    }
	}
	incr back_id
    }
}

for {set i 0} { $i < $part_id } {incr i} {
   # puts "[part $i]"
}

# Status report
puts "Simulate $n_part particles in a spherical cell with radius $sphere_rad"
puts "$n_poly PE backbones of length $l_back and charge distance $c_dist"
puts "with hairs of length $h_length every $h_dist monomer"
puts "neutralized by $n_ci counterions."
puts "Backbone monomers: $n_mono ($n_ci charged), Hair monomers: $n_hair."
puts "Setup in centered cylinder with radius $cyl_rad."
puts "Constraints:\n[constraint]"

#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    writepsf "$name$ident.psf"
    writepdb "$name$ident.pdb"
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    set HOSTNAME [exec hostname]
    set vmdout_file [open "vmdoutput.script" "w"]
    puts $vmdout_file "mol load psf $name$ident.psf pdb $name$ident.pdb"
    puts $vmdout_file "rotate stop"
    puts $vmdout_file "imd connect $HOSTNAME $port"
    close $vmdout_file
    puts "PSF and PDB written. IMD connection on port $port"
    puts "Start VMD in the same directory on the machine you like :"
    puts "vmd -e vmdoutput.script &"
    imd listen $vmd_wait
}

#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration: max $warm_n_times times $warm_steps steps"
set act_min_dist [analyze mindist]
set cap 20  
inter ljforcecap $cap

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

#############################################################
#  Counterion Equilibration                                #
#############################################################

puts "\nStart Counterion Equilibration: $ci_n_times times $ci_steps steps"
# remove LJ cap
inter ljforcecap 0
# Start Coulomb interaction
inter coulomb $bjerrum dh 0.0 [expr 2.0*$sphere_rad]

puts "Interactions:\n [inter]"

for {set i 0} { $i < $ci_n_times  } { incr i} {
    puts -nonewline "Run $i at time=[setmd time]\r"
    flush stdout
    integrate $ci_steps
    if { $vmd_output=="yes" } { imd positions }
}

puts "Counterion Equilibration done."

#############################################################
#      Integration                                          #
#############################################################

# free backbone monomers
for {set i 0} { $i < $n_mono } {incr i} {
    part $i v 0 0 0 unfix 
}

#for {set i 0} { $i < $part_id } {incr i} {
 #   puts "[part $i]"
#}
analyze set chains 0 $n_poly $l_back
puts "\nStart integration: $int_n_times times $int_steps steps"

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "Run $i at time=[setmd time], re = [analyze re]\r"
    flush stdout
    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }
}

puts "\nIntegration done."

#############################################################
puts "\nFinished"