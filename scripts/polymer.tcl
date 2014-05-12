#############################################################
#                                                           #
# polymer.tcl                                               #
# ===========                                               #
#                                                           #
# Script to setup polymer chains & networks.                #
#                                                           #
#############################################################
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
proc collision_tcl { i posx posy posz shield box_length } {
# Checks whether particle <i> at coordinates (<posx>, <posy>, <posz>) collides
# with any other particle j < i due to a minimum image distance smaller than <shield>.
    set min_dist [mindist $posx $posy $posz]
    if { [string is double $min_dist] } { 
	if { $min_dist > $shield } { return 0 } else { return 1 }
    } else {
	set min_dist 30000
	for {set j 0} {$j<$i} {incr j} {
	    set posj [part $j print pos]; set posjx [lindex $posj 0]; set posjy [lindex $posj 1]; set posjz [lindex $posj 2]
	    set dx [expr ($posx-$posjx) - round(($posx-$posjx)/$box_length)*$box_length]
	    set dy [expr ($posy-$posjy) - round(($posy-$posjy)/$box_length)*$box_length]
	    set dz [expr ($posz-$posjz) - round(($posz-$posjz)/$box_length)*$box_length]
	    set min_dist [min $min_dist [expr $dx*$dx + $dy*$dy + $dz*$dz]]
	}
	if { $min_dist > [sqr $shield] } { return 0 } else { return 1 }
    }
}



proc polymer_tcl { MPC bond_length box_length {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} {p7 NA} {p8 NA} {p9 NA} {p10 NA} {p11 NA} {p12 NA} {p13 NA} {p14 NA} {p15 NA} } {
# polymer_tcl <MPC> <bond_length> <box_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charge <val_cM>] [distance <cM_dist>] [types <type_nM> [<type_cM>]] [FENE <type_FENE>]
# Creates one polymer chain within the simulation box,
# and returns how often the attempt to place a monomer failed in the worst case.
# Parameters:  <MPC>         = monomers per chain
#              <bond_length> = length of the bonds between two monomers
#              <box_length>  = length of the simulation box
#              <part_id>     = particle number of the start monomer (defaults to '0')
#              <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
#              <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
#              <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
#              <val_cM>      = valency of charged monomers (defaults to '1')
#              <cM_dist>     = distance between two charged monomers' indices (defaults to '1')
#              <type_{n|c}P> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1')
#              <type_FENE>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0')
    set param [list $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 $p14 $p15]
    set part_id 0; set mode "SAW"; set shield 0; set max_try 30000; set val_cM 1; set cM_dist 1; set type_nM 0; set type_cM 1; set type_FENE 0
    for {set i 0} {$i < 15} {incr i} {
	switch [lindex $param $i] {
	    "start" { incr i; set part_id [lindex $param $i] 
		      if { ![string is integer $part_id] } { puts "Index of polymer chain's first monomer must be integer (got: $part_id)!"
			  puts "Aborting...\n"; exit } }
	    "mode"  { incr i; set mode [lindex $param $i]
		      if { $mode=="SAW" } { 
			  set shield [lindex $param [expr $i+1]]
			  if { [string is double $shield] } {
			      if { $shield < 0 } { puts "The SAW-shield must be >0 (got: $shield)!\nAborting...\n"; exit 
			      } else { incr i 
				  set max_try [lindex $param [expr $i+1]]
				  if { [string is integer $max_try] } { incr i } else { set max_try 30000 } }
			  } else { set shield 0 }
		      } elseif { $mode!="RW" } { puts "The mode for the chain's setup must be either 'SAW' or 'RW' (got: $mode)!\nAborting...\n"; exit } }
	    "charge" {incr i; set val_cM [lindex $param $i]
		      if { ![string is integer $val_cM] } { puts "The charge of the chain's monomers must be integer (got: $val_cM)!"
			  puts "Aborting...\n"; exit } }
	    "distance" { incr i; set cM_dist [lindex $param $i]
		      if { ![string is integer $cM_dist] } { puts "The distance between two charged monomers' indices must be integer (got: $cM_dist)!"
			  puts "Aborting...\n"; exit } }
	    "types" { incr i; set type_nM [lindex $param $i]
		      if { ![string is integer $type_nM] } { puts "The type-\# of neutral monomers must be integer (got: $type_nM)!"
			  puts "Aborting...\n"; exit }
		      set type_cM [lindex $param [expr $i+1]]
		      if { ![string is integer $type_cM] } { set type_cM $type_nM } else { incr i } }
	    "FENE"  { incr i; set type_FENE [lindex $param $i]
		      if { ![string is integer $type_FENE] } { puts "The type-\# of the FENE-interaction must be integer (got: $type_FENE)!"
			  puts "Aborting...\n"; exit } }
	    default { if { [lindex $param $i]!="NA" } {
		          puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	            }
	}
    }
    set max_cnt 0
    for {set cnt2 0} {$cnt2<$max_try} {incr cnt2} {
	for {set cnt 0} {$cnt<$max_try} {incr cnt} {
	    set posx [expr $box_length*[tcl_rand]]
	    set posy [expr $box_length*[tcl_rand]]
	    set posz [expr $box_length*[tcl_rand]]
	    if { $mode!="SAW" || [expr [collision_tcl $part_id $posx $posy $posz $shield $box_length]==0] } { break }
	}
	if {$cnt >= $max_try} { puts "Failed to find a suitable place for the start-monomer for $cnt times!\nAborting...\n"; exit }
	part $part_id pos $posx $posy $posz q $val_cM type $type_cM
	incr part_id; set max_cnt [max $cnt $max_cnt]
	for {set n 1} { $n < $MPC } {incr n} {
	    for {set cnt 0} {$cnt<$max_try} {incr cnt} {
		set theta [expr     [PI]*[tcl_rand]]
		set phi   [expr 2.0*[PI]*[tcl_rand]]
		set posx  [expr $posx + $bond_length*sin($theta)*cos($phi)]
		set posy  [expr $posy + $bond_length*sin($theta)*sin($phi)]
		set posz  [expr $posz + $bond_length*cos($theta)]
		if { $mode!="SAW" || [expr [collision_tcl $part_id $posx $posy $posz $shield $box_length]==0] } { break }
	    }
	    if { $cnt >= $max_try } { 
		puts -nonewline "Warning! Attempt \# [expr $cnt2+1] to build a polymer failed after $cnt unsuccessful trials to place monomer $n! " 
		puts "Retrying by re-setting the start-monomer..."
		set part_id [expr $part_id - $n]; set n 0; break }
	    part $part_id pos $posx $posy $posz q [expr ($n%$cM_dist == 0) ? $val_cM : 0] \
		type  [expr ($n%$cM_dist == 0) ? $type_cM : $type_nM ] \
		bond $type_FENE [expr $part_id-1]
	    incr part_id; set max_cnt [max $cnt $max_cnt]
	}
	if { $mode!="SAW" || $n>0 } { break}
    }
    if {$cnt2 >= $max_try} { puts "Failed to place current polymer chain in the simulation box for $cnt2 times!\nAborting...\n"; exit
    } else { return [max $max_cnt $cnt2] }
}



proc polymers_tcl { N_P MPC bond_length box_length {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} {p7 NA} {p8 NA} {p9 NA} {p10 NA} {p11 NA} {p12 NA} {p13 NA} {p14 NA} {p15 NA} } {
# polymers_tcl <N_P> <MPC> <bond_length> <box_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charge <val_cM>] [distance <cM_dist>] [types <type_nM> [<type_cM>]] [FENE <type_FENE>]
# Creates <N_P> polymer chains within the simulation box using the 'polymer_tcl' command,
# and returns how often the attempt to place a polymer failed in the worst case.
# Parameters:  <N_P>         = number of polymer chains to create
#              <MPC>         = monomers per chain
#              <bond_length> = length of the bonds between two monomers
#              <box_length>  = length of the simulation box
#              <part_id>     = particle number of the start monomer (defaults to '0')
#              <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
#              <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
#              <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
#              <val_cM>      = valency of charged monomers (defaults to '1')
#              <cM_dist>     = distance between two charged monomers' indices (defaults to '1')
#              <type_{n|c}P> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1')
#              <type_FENE>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0')
    set param [list $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 $p14 $p15]
    set part_id 0
    for {set i 0} {$i < 15} {incr i} {
	switch [lindex $param $i] {
	    "start" { incr i; set part_id [lindex $param $i] 
		      if { ![string is integer $part_id] } { puts "Index of the first monomer must be integer (got: $part_id)!\nAborting...\n"; exit } }
	    default { }
	}
    }
    for { set i 0 } { $i < [min 15 [llength $param]] } {incr i} {
	if { [lindex $param $i]=="start" } { incr i } else { eval set p[expr $i+1] [lindex $param $i] }
    }
    if { [llength $param]==15 } { set param [lreplace $param end-1 end] }
    set max_cnt 0
    for { set i 0 } { $i < $N_P } { incr i } {
	set max_cnt [max $max_cnt [polymer_tcl $MPC $bond_length $box_length start [expr $part_id+$i*$MPC] \
				      $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 ] ]
    }
    return $max_cnt
}



proc counterions_tcl { N_CI box_length {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} {p7 NA} {p8 NA} {p9 NA} {p10 NA} } {
# counterions_tcl <N_CI> <box_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charge <val_CI>] [type <type_CI>]
# Creates <N_CI> counterions of charge <val_CI> within the simulation box of cubic size <box_length>,
# and returns how often the attempt to place a particle failed in the worst case.
# Parameters:  <N_CI>        = number of counterions to create
#              <box_length>  = length of the simulation box
#              <part_id>     = particle number of the first counterion (defaults to '[setmd n_part]')
#              <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
#              <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
#              <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
#              <val_CI>      = valency of the counterions (defaults to '-1')
#              <type_CI>     = type number of the counterions to be used with "part" (default to '2')
    set param [list $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10]
    set part_id [setmd n_part]; set mode "SAW"; set shield 0; set max_try 30000; set val_CI -1; set type_CI 2
    for {set i 0} {$i < 10} {incr i} {
	switch [lindex $param $i] {
	    "start" { incr i; set part_id [lindex $param $i] 
		      if { ![string is integer $part_id] } { puts "Index of the first counterion must be integer (got: $part_id)!\nAborting...\n"; exit } }
	    "mode"  { incr i; set mode [lindex $param $i]
		      if { $mode=="SAW" } { 
			  set shield [lindex $param [expr $i+1]]
			  if { [string is double $shield] } {
			      if { $shield < 0 } { puts "The SAW-shield must be >0 (got: $shield)!\nAborting...\n"; exit 
			      } else { incr i 
				  set max_try [lindex $param [expr $i+1]]
				  if { [string is integer $max_try] } { incr i } else { set max_try 30000 } }
			  } else { set shield 0 }
		      } elseif { $mode!="RW" } { puts "The mode for the chain's setup must be either 'SAW' or 'RW' (got: $mode)!\nAborting...\n"; exit } }
	    "charge" {incr i; set val_CI [lindex $param $i]
		      if { ![string is integer $val_CI] } { puts "The charge of the counterions must be integer (got: $val_CI)!"
			  puts "Aborting...\n"; exit } }
	    "type"  { incr i; set type_CI [lindex $param $i]
		      if { ![string is integer $type_CI] } { puts "The type-\# of the counterions must be integer (got: $type_CI)!"
			  puts "Aborting...\n"; exit } }
	    default { if { [lindex $param $i]!="NA" } {
		          puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	            }
	}
    }
    set max_cnt 0; set part_0 $part_id
    for {set i 0} { $i < $N_CI } {incr i} {
	for {set cnt 0} {$cnt<$max_try} {incr cnt} {
	    set posx [expr $box_length*[tcl_rand]]
	    set posy [expr $box_length*[tcl_rand]]
	    set posz [expr $box_length*[tcl_rand]]
	    if { $mode!="SAW" || [expr [collision_tcl $part_id $posx $posy $posz $shield $box_length]==0] } { break }
	}
	if {$cnt >= $max_try} { puts "Failed to find a suitable place for counterion [expr $part_id-$part_0] after $cnt trials!\nAborting...\n"; exit }
	part $part_id pos $posx $posy $posz q $val_CI type $type_CI
	incr part_id; set max_cnt [max $cnt $max_cnt]
    }
    return $max_cnt
}



proc salt_tcl { N_pS N_nS box_length {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} {p7 NA} {p8 NA} {p9 NA} {p10 NA} {p11 NA} {p12 NA} } {
# salt_tcl <N_pS> <N_nS> <box_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charges <val_pS> [<val_nS>]] [types <type_pS> [<type_nS>]]
# Creates <N_pS> positively and <N_nS> negatively charged salt ions of charge <val_pS> and <val_nS> within the simulation box of cubic size <box_length>,
# and returns how often the attempt to place a particle failed in the worst case.
# Parameters:  <N_pS>/<N_nS> = number of salt ions to create
#              <box_length>  = length of the simulation box
#              <part_id>     = particle number of the first salt ion (defaults to '[setmd n_part]')
#              <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
#              <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
#              <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
#              <val_{p|n}S>  = valencies of the salt ions (default to '1' and '-1', respectively); if <val_nS> is not given, <val_nS> = -1*<val_pS>
#              <type_{p|n}S> = type numbers to be used with "part" (default to '3' and '4'); if <type_nS> is not given, <type_nS> = <type_pS> is assumed
    set param [list $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12]
    set part_id [setmd n_part]; set mode "SAW"; set shield 0; set max_try 30000; set val_pS 1; set val_nS -1; set type_pS 3; set type_nS 4
    for {set i 0} {$i < 12} {incr i} {
	switch [lindex $param $i] {
	    "start" { incr i; set part_id [lindex $param $i] 
		      if { ![string is integer $part_id] } { puts "Index of the first counterion must be integer (got: $part_id)!\nAborting...\n"; exit } }
	    "mode"  { incr i; set mode [lindex $param $i]
		      if { $mode=="SAW" } { 
			  set shield [lindex $param [expr $i+1]]
			  if { [string is double $shield] } {
			      if { $shield < 0 } { puts "The SAW-shield must be >0 (got: $shield)!\nAborting...\n"; exit 
			      } else { incr i 
				  set max_try [lindex $param [expr $i+1]]
				  if { [string is integer $max_try] } { incr i } else { set max_try 30000 } }
			  } else { set shield 0 }
		      } elseif { $mode!="RW" } { puts "The mode for the chain's setup must be either 'SAW' or 'RW' (got: $mode)!\nAborting...\n"; exit } }
	    "charges" {incr i; set val_pS [lindex $param $i]
		      if { ![string is integer $val_pS] } { puts "The charge of the positive salt ions must be integer (got: $val_pS)!"
			  puts "Aborting...\n"; exit }
		      set val_nS [lindex $param [expr $i+1]]
		      if { ![string is integer $val_nS] } { set val_nS [expr (-1)*$val_pS] } else { incr i } }
	    "types" { incr i; set type_pS [lindex $param $i]
		      if { ![string is integer $type_pS] } { puts "The type-\# of positive salt ions must be integer (got: $type_pS)!"
			  puts "Aborting...\n"; exit }
		      set type_nS [lindex $param [expr $i+1]]
		      if { ![string is integer $type_nS] } { set type_nS $type_pS } else { incr i } }
	    default { if { [lindex $param $i]!="NA" } {
		          puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	            }
	}
    }
    set max_cnt 0; set part_0 $part_id
    for {set i 0} { $i < $N_pS } {incr i} {
	for {set cnt 0} {$cnt<$max_try} {incr cnt} {
	    set posx [expr $box_length*[tcl_rand]]
	    set posy [expr $box_length*[tcl_rand]]
	    set posz [expr $box_length*[tcl_rand]]
	    if { $mode!="SAW" || [expr [collision_tcl $part_id $posx $posy $posz $shield $box_length]==0] } { break }
	}
	if {$cnt >= $max_try} { 
	    puts "Failed to find a suitable place for positive salt ion [expr $part_id-$part_0] after $cnt trials!\nAborting...\n"; exit }
	part $part_id pos $posx $posy $posz q $val_pS type $type_pS
	incr part_id; set max_cnt [max $cnt $max_cnt]
    }
    for {set i 0} { $i < $N_nS } {incr i} {
	for {set cnt 0} {$cnt<$max_try} {incr cnt} {
	    set posx [expr $box_length*[tcl_rand]]
	    set posy [expr $box_length*[tcl_rand]]
	    set posz [expr $box_length*[tcl_rand]]
	    if { $mode!="SAW" || [expr [collision_tcl $part_id $posx $posy $posz $shield $box_length]==0] } { break }
	}
	if {$cnt >= $max_try} { 
	    puts "Failed to find a suitable place for negative salt ion [expr $part_id-$part_0] after $cnt trials!\nAborting...\n"; exit }
	part $part_id pos $posx $posy $posz q $val_nS type $type_nS
	incr part_id; set max_cnt [max $cnt $max_cnt]
    }
    return $max_cnt
}



proc velocities_tcl { N_T v_max {p1 NA} {p2 NA} } {
# velocities_tcl <N_T> <v_max> [start <part_id>]
# Sets the velocities of <N_T> particles to a random value [-vmax,vmax],
# and returns the averaged velocity when done.
# Parameters:  <N_T>         = number of particles of which the velocities should be set
#              <v_max>       = maximum velocity to be used
#              <part_id>     = particle number of the first of the <N_T> particles (defaults to '0')
    set param [list $p1 $p2]
    set part_id 0
    for {set i 0} {$i < 2} {incr i} {
	switch [lindex $param $i] {
	    "start" { incr i; set part_id [lindex $param $i] 
		      if { ![string is integer $part_id] } { puts "Index of the first particle must be integer (got: $part_id)!\nAborting...\n"; exit } }
	    default { if { [lindex $param $i]!="NA" } {
		          puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	            }
	}
    }
    set v_av 0
    for {set i $part_id} { $i < [expr $N_T+$part_id] } {incr i} {
	set r     [expr   $v_max*[tcl_rand]]
	set theta [expr     [PI]*[tcl_rand]]
	set phi   [expr 2.0*[PI]*[tcl_rand]]
	set vx  [expr $r*sin($theta)*cos($phi)]
	set vy  [expr $r*sin($theta)*sin($phi)]
	set vz  [expr $r*cos($theta)]
	part $i v $vx $vy $vz
	set v_av [expr $v_av + sqrt($vx*$vx + $vy*$vy + $vz*$vz)]
    }
    return [expr $v_av/$N_T]
}

