# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
source "tests_common.tcl"

require_feature "MASS"
require_feature "LENNARD_JONES"

puts "-------------------------------------------"
puts "- Testcase mass.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------"

####################### PROCEDURES

## read and write data

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos type mass v}
    close $f
}

## creates the data for the mass testcase

proc cr_v {arg1 arg2} {
    set delta [expr $arg2 - $arg1]
    set v [expr [t_random]*$delta + $arg1] 
    return $v 
}

proc cr_t {arg1} {
    set arg [expr int([t_random]*($arg1-1))]
    return $arg
}

############################## SYSTEM PARAMETERS

setmd time_step 0.001
setmd skin 0.5
set eps_mom 1e-6
set eps_en 0.2
set step 1000
set step_count 10
set rcut_min 1.0
set rcut_max 2.0
set nt 10

set watch 0
set count_hits 0
set write_system 0

################################################

proc system_setup {nt rcut_min rcut_max step} {
    
    ###################System setup
    set v_min 1.0
    set v_max 20.0
    set m_min 1.0
    set m_max 10.0
    set np 1000
    set sik 5.5
    set density 0.74
    ###############################

    if {$np/2.0 != int($np/2.0)} {error "Number of particles must be even"}
    #calculate boxlength
 
    set sum 0
    set drcut [expr ($rcut_max - $rcut_min) / ($nt - 1.0)]
    for {set i 0} {$i < $nt} {incr i} {
	 set rcut_i [expr pow(($rcut_min + $i*$drcut),3)]
	 set sum [expr $sum+$rcut_i]
    }
    set V [expr (3.0/$density) * ($np/2.0) * ($sum/$nt) ]
    set bl [expr pow($sik*$V,1.0/3.0)]
    #create particles with "counterions"
    setmd box_l [expr $bl/3.0] $bl $bl 
    counterions [expr int($np*0.5)] start 0 mode SAW $rcut_max
    set np1 [setmd max_part]
    set list [part]
    part deleteall
    counterions [expr int($np*0.5)] start 0 mode SAW $rcut_max
    setmd box_l $bl $bl $bl
    foreach i $list {
	part [expr [setmd max_part]+1] pos [expr [lindex $i 2]+[expr $bl*2.0/3.0]] [lindex $i 3] [lindex $i 4]
    }
    unset list 
    # set mass, type and velocity
    set dm [expr ($m_max - $m_min) / ($nt - 1.0)]
    for {set i 0} {$i <= $np1} {incr i} {
	#set type and mass
	set t [cr_t $nt]
	set m [expr $m_min + $t*$dm]
	part $i type $t mass $m
	part [expr $i+($np/2)] type $t mass $m
	#set velocities
	set veloc [cr_v $v_min $v_max]
	part $i v $veloc 0.0 0.0
	part [expr $i+($np/2)] v [expr -1*$veloc] 0.0 0.0
    }
}

##  calculating the momentum
proc momentum {arg1 arg2} {
    set mom "0 0 0"
    for {set i $arg1} {$i <= $arg2} {incr i} {
	set m_tmp [part $i print mass]
	set vel_tmp [part $i print v]
	set mom_tmp [vecscale $m_tmp $vel_tmp]
	set mom [vecadd $mom $mom_tmp]
    }
    return $mom
}

##counting hits
proc count { list } {
    set res 0
    for {set i 0} {$i <= [setmd max_part]} {incr i} {
	set v_new [part $i print v]
	lappend v_old [lindex $list $i 12] [lindex $list $i 13] [lindex $list $i 14] 
	set diff [vecsub $v_new $v_old] 
	foreach e $diff { if {$e > 0.5} {incr res}; break}
	unset v_old
    }
    return $res
}

############################## START SIMULATION
##### Lennard-Jones, no thermostat

if { [catch {  
    thermostat off

    #here you can create the necessary snapshot
    if {$write_system == 1} {
	system_setup $nt $rcut_min $rcut_max $step
	write_data "mass_system.data"
    } else {read_data "mass_system.data"}
    set drcut [expr ($rcut_max - $rcut_min) / ($nt - 1.0)]
    set lj_eps 1.0
    set lj_off 0.0
    for {set i 0} {$i < $nt} {incr i } {
	for {set j $i} {$j < $nt} {incr j} {
	    if {$i == $j} {
		set lj_sig [expr $rcut_min + $i*$drcut]
	    } else { 
		set lj_sig [expr ($rcut_min + $i*$drcut + $rcut_min + $j*$drcut) / 2.0]
	    }
	    set lj_cut [expr pow(2.0,1.0/6.0) * $lj_sig]
	    set lj_shf [expr -(pow($lj_sig/($lj_cut-$lj_off),12.0) - pow($lj_sig/($lj_cut-$lj_off),6.0))]
	    inter $i $j lennard-jones $lj_eps $lj_sig $lj_cut $lj_shf $lj_off
	}
    }
    set mom [momentum 0 [setmd max_part]]
    set en_old [analyze energy total]
    puts "The momentum is $mom."
    puts "The total energy is $en_old."
    puts "The kinetic energy is [analyze energy kinetic]."
    
    #here you can watch the simuluation and count the hits 
    if {$count_hits == 1} {
	    set nh 0
	    if {$watch == 1} {prepare_vmd_connection vmd 3000}
	    for {set i 0} {$i <= $step} {incr i $step_count} {
		set list [part]
		integrate $step_count
		if {$watch == 1} {imd position}
		set nh [expr $nh + [count $list]] 
	    }
	    set av_nh [expr $nh*$step_count/$step]
     } elseif {$count_hits == 0} {
	    if {$watch == 1} {
		prepare_vmd_connection vmd 3000
		for {set i 0} {$i <= $step} {incr i $step_count} {integrate $step_count; imd position}
	    } else {integrate $step}
     }
#### Analyze energy and momentum
    set mom [momentum 0 [setmd max_part]]
    set en_new [analyze energy total]
    set en_err [expr abs(($en_new-$en_old)/$en_old) * 100]
    puts "The deviation of the momentum is [veclen $mom]"
    puts "The relative error of the energy is $en_err %."
    if { [veclen $mom] > $eps_mom } {
        puts "deviation of the momentum in direction x is [lindex $mom 0]"
	puts "deviation of the momentum in direction y is [lindex $mom 1]"
        puts "deviation of the momentum in direction z is [lindex $mom 2]"
        error "deviation of momentum too large"
    }
    if { $en_err > $eps_en } { error "The error of the relative energy is too large"}
    puts "The total energy now is $en_new."
    puts "The kinetic energy now is [analyze energy kinetic]."
    if {$count_hits == 1} { puts "The average number of hits was $av_nh"}
########################################################
} res ] } {
    error_exit $res
}

exit 0
