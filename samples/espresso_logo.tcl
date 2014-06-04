#############################################################
#                                                           #
#  System: ESPRESSO LOGO                                    #
#  Even though this system is just for fun                  #
#  it is actually a real polyelectrolyte solution           #
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
puts "=              espresso_logo.tcl                      ="
puts "======================================================="
puts " "

require_feature "LENNARD_JONES" "EXTERNAL_FORCES" "ELECTROSTATICS"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "espresso"
set ident "_s5"

# Visualisation: online (imd), or offline (vtf-file)
set vmd_online "no"
set vmd_offline "yes"

# name of the vmd start script
set vmd_script "espresso_logo_vmd.tcl"

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
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]

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

set vtf_bonds ""

# espresso cup (polyelectrolytes)
append vtf_bonds "\#cup\n"
set pid 0
set start_pid 0
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
    set end_pid [expr $pid-1]
    append vtf_bonds "bond " $start_pid "::" $end_pid "," $start_pid ":" $end_pid "\n"
    set start_pid $pid
}

append vtf_bonds "\#cup bottom\n"
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
    set end_pid [expr $pid-1]
    append vtf_bonds "bond " $start_pid "::" $end_pid "," $start_pid ":" $end_pid "\n"
    set start_pid $pid
}

# cup handle (one more polyelectrolyte)
append vtf_bonds "\#handle\n"
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
set end_pid [expr $pid-1]
append vtf_bonds "bond " $start_pid "::" $end_pid "\n"
set start_pid $pid

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
append vtf_bonds "\#steam\n"
set start_pid $pid
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
    set end_pid [expr $pid-1]
    append vtf_bonds "bond " $start_pid "::" $end_pid "\n"
    set start_pid $pid
}


#rerun only when you have changed some parameters
#puts "[inter coulomb $bjerrum p3m tune accuracy 0.0005]"
inter coulomb $bjerrum p3m tune accuracy 2.92863e-01 mesh 36

puts "Interactions:\n[inter]"
set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

#############################################################
#  Fixed Integration                                       #
#############################################################

# find free port for IMD
if { $vmd_online == "yes" } then {
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    if {$port == 65000} {
	error "could not open a port for IMD, please try offline"
    }
}

# create VMD start script
if { $vmd_offline == "yes" || $vmd_online == "yes" } then {
    puts "Creating $vmd_script..."
    set vmdout_file [open $vmd_script "w"]
    if { $vmd_online == "yes" } then {
	puts $vmdout_file "mol load vtf espresso_logo_online.vtf"
	puts $vmdout_file "imd connect localhost $port"
    } elseif { $vmd_offline == "yes" } then {
	puts $vmdout_file "mol load vtf espresso_logo.vtf"
    }

    puts $vmdout_file {
display resize 600 600
render options POV3 vmd_povray +W600 +H600 +ua -I%s +FN &
axes location off
color Display Background 8
# cup blue
color Name O 0
# saucer red
color Name N 1 
# steam silver
color Name S 6

animate goto 44
scale to 0.12
translate to 0 0.4 0.5

mol addrep 0

# cup and saucer
mol modselect 0 0 "not name S"
mol modstyle 0 0 CPK 3 0.3 8 6
mol modmaterial 0 0 Glossy

# steam
mol modselect 1 0 "name S"
mol modstyle 1 0 CPK 3 0.3 8 6
mol modmaterial 1 0 Glass2
}
    
    close $vmdout_file
    puts "$vmd_script finished."
}

# start vmd
if { $vmd_online == "yes" } then {
    set vtf_file [open "espresso_logo_online.vtf" w]
    writevsf $vtf_file
    puts $vtf_file $vtf_bonds
    writevcf $vtf_file
    close $vtf_file

    # start VMD
    exec vmd -e $vmd_script &
    # wait for VMD to connect
    imd listen 100000
}


# create vtf-file
if { $vmd_offline == "yes" } then {
    set vtf_file [open "espresso_logo.vtf" w]
    writevsf $vtf_file
    # connect cup
    puts $vtf_file $vtf_bonds
    writevcf $vtf_file
}



#############################################################
#      Integration                                          #
#############################################################

puts "\nStart integration with fixed cup:"
puts "Run $fixed_n_times times $fixed_steps steps"

for {set i 0} { $i < $fixed_n_times } { incr i} {
    puts -nonewline "run $i at time=[setmd time]\r"
    integrate $fixed_steps
    if { $vmd_online=="yes" } { imd positions }
    if { $vmd_offline=="yes" } then { writevcf $vtf_file }
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

    if { $vmd_online=="yes" } { imd positions }
    if { $vmd_offline=="yes" } then { writevcf $vtf_file }
}

if { $vmd_online=="yes" } { exec rm vmd_start.script }
if { $vmd_offline=="yes" } then { close $vtf_file }

puts "\n\nFinished"
exit
