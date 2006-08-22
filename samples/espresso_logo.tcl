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
#  System: ESPRESSO LOGO                                    #
#  Even though this system is just for fun                  #
#  it is actually a real polyelectrolyte solution           #
#                                                           #
#  Created:       17.03.2003 by HL                          #
#  Last modified: 18.03.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=    Sample script 5: espresso_logo.tcl               ="
puts "======================================================="
puts " "

require_feature "LENNARD_JONES" "EXTERNAL_FORCE"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "espresso"
set ident "_s5"

# 
set vmd_output "yes"

# System parameters
#############################################################

set box_l 20.0

set cup_top_circ 21
set cup_bot_circ 15
set cup_height   6

set saucer_circ  30

set n_steam 5
set l_steam 6

# Interaction parameters (repulsive Lennard Jones)
#############################################################

# Lennard Jones
set n_part_types  3
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   $lj1_eps

# FENE
set fene_cut      2.0
set fene_k        7.0

# Coulomb
set bjerrum       2.5
set accuracy      1.0e-4


# Integration parameters
#############################################################

setmd time_step 0.01
setmd skin      1.0
thermostat langevin 1.0 1.0

# warmup integration (with capped LJ potential)
set fixed_steps   10
set fixed_n_times 100

# integration
set int_steps    10
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

# fene
inter 0 fene $fene_k $fene_cut
# lj
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
        inter $ia1 $ia2 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0.0
    }
}
# coulomb: after particle setup

# Particle setup
#############################################################

set mypi   3.141592653589793

# espresso cup (polyelectrolytes)
set pid 0
for {set i 0} { $i < $cup_height } {incr i} {
    set circ [expr $cup_bot_circ + $i*($cup_top_circ-$cup_bot_circ)/double($cup_height-1)]
    set rad  [expr $circ/(2.0*$mypi)]
    set alpha [expr 2.0*$mypi/int($circ)]
    set posy [expr 2.0+$i] 
    for {set j 0} { $j < [expr int($circ)] } {incr j} {
        set posx [expr $box_l/2.0 + $rad*sin($j*$alpha+($mypi/2.0))]
        set posz [expr $box_l/2.0 + $rad*cos($j*$alpha+($mypi/2.0))]
	part $pid pos $posx $posy $posz type 0 q 1.0
	if { $j > 0 } { part $pid bond 0 [expr $pid-1] }
	incr pid
    }   
}

set rad [expr $cup_bot_circ/(2.0*$mypi)]
set posy 2.0
puts "cup: bottom_rad $rad"
while { $rad > 1.0 } {
    set rad   [expr $rad-0.9]
    set circ  [expr 2.0*$mypi*$rad]
    set alpha [expr 2.0*$mypi/int($circ)]
    for {set j 0} { $j < [expr int($circ)] } {incr j} {
        set posx [expr $box_l/2.0 + $rad*sin($j*$alpha+($mypi/2.0))]
        set posz [expr $box_l/2.0 + $rad*cos($j*$alpha+($mypi/2.0))]
	part $pid pos $posx $posy $posz type 0 q 1.0
	if { $j > 0 } { part $pid bond 0 [expr $pid-1] }
	incr pid
    }   
}

# cup handle (one more polyelectrolyte)
set hand_rad  [expr ($cup_height-4.0)/sqrt(2.0)]
set hand_circ [expr (1.5*$mypi*$hand_rad)]
set hand_xoff [expr ($cup_bot_circ+$cup_top_circ)/(4.0*$mypi)+1.2]
set hand_yoff [expr 2.0+$cup_height/2.0 -0.2]
set alpha [expr 2.0*$mypi/int(4.0*$hand_circ/3.0)]
set beta  [expr sin(($cup_top_circ-$cup_bot_circ)/(2.0*$mypi*$cup_height-1))]
set beta [expr $beta - $mypi/8.0]
puts "handle rad=$hand_rad circ=$hand_circ beta $beta"
set posz [expr ($box_l/2.0)+0.5]
for {set i 0} { $i < [expr int($hand_circ)]} {incr i} {
    set posx [expr $hand_xoff+$box_l/2.0 + $hand_rad*sin($i*$alpha+$beta)]
    set posy [expr $hand_yoff + $hand_rad*cos($i*$alpha+$beta)]
    part $pid pos $posx $posy $posz type 0 q 1.0
    if { $i > 0 } { part $pid bond 0 [expr $pid-1] }
    incr pid
}

# saucer (counterions)
set s_rad_o [expr $saucer_circ/(2.0*$mypi)]
set s_rad_i [expr $cup_bot_circ/(2.0*$mypi)]
set n_saucer [expr int($s_rad_o-$s_rad_i) + 1]
puts "saucer: r_out=$s_rad_o, r_in=$s_rad_i, rings=$n_saucer"
set n_ci 0
for { set i 0 } { $i < $n_saucer } { incr i } {
    set n_ci [expr $n_ci + int($saucer_circ-($i*2.0*$mypi))]
}
set ci_val [expr -$pid/double($n_ci)]
puts "Parts: total [expr $pid+$n_ci], cup $pid, saucer $n_ci (valency $ci_val)"
for { set i 0 } { $i < $n_saucer } { incr i } {
    set rad  [expr $s_rad_o-$i]
    set alpha [expr 2.0*$mypi/int($saucer_circ-($i*2.0*$mypi))]
    set posy [expr 2.3 - 0.5*$i]
    for {set j 0} { $j < int($saucer_circ-($i*2.0*$mypi)) } {incr j} {
        set posx [expr $box_l/2.0 + $rad*sin($j*$alpha)]
        set posz [expr $box_l/2.0 + $rad*cos($j*$alpha)]
	part $pid pos $posx $posy $posz type 1 q $ci_val
	incr pid
    }
}

for {set i 0} { $i < $pid} {incr i} { part $i fix }

# steam (some neutral polymers)
for {set i 0} { $i < $n_steam } {incr i} {
    set posy [expr $cup_height] 
    set rad  [expr ($cup_top_circ-12.0)/(2.0*$mypi)]
    set alpha [expr 2.0*$mypi/int($n_steam)]
    set posx [expr $box_l/2.0 + $rad*sin($i*$alpha)]
    set posz [expr $box_l/2.0 + $rad*cos($i*$alpha)]
    for {set j 0} { $j < $l_steam } {incr j} {
	part $pid pos $posx $posy $posz type 2 
	if { $j > 0 } { 
	    part $pid bond 0 [expr $pid-1]
	    part $pid ext_force 0.0 5.0 0.0
	}
	if { $j == 0 } { part $pid fix }
	incr pid
	set posy [expr $posy+1.0]
    }
}

#rerun only when you have changed some parameters
#puts "[inter coulomb $bjerrum p3m tune accuracy 0.0005]"
inter coulomb $bjerrum p3m 9.0 16 4 2.92863e-01

puts "Interactions:\n[inter]"
set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

#############################################################
#  Fixed Integration                                       #
#############################################################


# vmd connection
writepsf "espresso.psf"
writepdb "espresso.pdb"
for {set port 10000} { $port < 65000 } { incr port } {
    catch {imd connect $port} res
    if {$res == ""} break
}
set HOSTNAME [exec hostname]
set vmdout_file [open "vmd_start.script" "w"]
puts $vmdout_file "axes location off"
puts $vmdout_file "mol load psf espresso.psf pdb espresso.pdb"
puts $vmdout_file "rock y by -2"
puts $vmdout_file "mol modstyle 0 0 CPK 1.500000 0.600000 8.000000 6.000000"
puts $vmdout_file "mol modcolor 0 0 SegName"
puts $vmdout_file "imd connect $HOSTNAME $port"
close $vmdout_file

exec vmd -e vmd_start.script &

#############################################################
#      Integration                                          #
#############################################################

puts "\nStart integration with fixed cup:"
puts "Run $fixed_n_times times $fixed_steps steps"

set j 0
for {set i 0} { $i < $fixed_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time]\r"
    integrate $fixed_steps
    polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
    incr j
    if { $vmd_output=="yes" } { imd positions }
    flush stdout
    incr i
}

puts "\nStart integration with unfixed particles:"
puts "Run $int_n_times times $int_steps steps"
for {set i 0} { $i < $pid} {incr i} { part $i unfix }

for {set i 0} { $i < $int_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time] \r"
    flush stdout

    integrate $int_steps
    polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
    incr j

    if { $vmd_output=="yes" } { imd positions }
}

exec rm vmd_start.script
puts "\n\nFinished"
exit
