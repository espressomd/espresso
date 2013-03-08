#############################################################
#                                                           #
#  Diamond Hydrogel Networks                                #
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

puts " "
puts "==================================================="
puts "=                  diamond.tcl                    ="
puts "==================================================="
puts " "

puts "Program Information: \n[code_info]\n"



#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "diamond"
set ident "_t6"

# On 'yes' connects to 'vmd' visualizing current configuration
set vmd_output "no"


# System parameters
#############################################################

set N_node  8
set N_P     16 
set MPC     20
set N_CI    328
set N_T     [expr $N_node + $N_P*$MPC + $N_CI]
set packing { 0.45 0.4 0.1 0.05 0.025 0.01 0.0075 0.005 0.0025 0.001 0.00075  }


# Interaction parameters
#############################################################

# repulsive Lennard Jones
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.122462048309373
set lj1_shift   [expr -4.0*$lj1_eps*(pow($lj1_cut,-12)-pow($lj1_cut,-6))]

# attractive FENE
set fene_k      15.551
set fene_r      1.25

# electrostatics
set bjerrum     1.75241
set accuracy    1.0e-1



# Integration parameters
#############################################################

setmd time_step 0.012
setmd skin      0.4
thermostat langevin 1.0 1.0

# warmup integration (with capped LJ potential) until particles are at least $min_dist apart (but at least $min_loop loops)
set warm_step   200
set warm_loop   300
set warm_cap1   10
set warm_incr   25
set min_dist    0.9
set min_loop    200

# integration (with full LJ potential) for $int_time
set int_step    1000
set int_time    [expr 20000*[setmd time_step]]


# Other parameters
#############################################################

set tcl_precision  6
set random_seeds  { }



#############################################################
#  Setup System                                             #
#############################################################

if { [file exists "$name$ident.DHN"] } {
    set DHN_file [open "$name$ident.DHN" "a"]
} else {
    set DHN_file [open "$name$ident.DHN" "w"]
    puts -nonewline $DHN_file "ID int_time pack_i density a_cube re rg rh re2/rg2 Temp mindist "
    if { $bjerrum > 0.0 } { 
	puts $DHN_file "p_total D(p_total) p_FENE D(p_FENE) p_lj D(p_lj) p_c D(p_c) p_ideal p_osmotic"
    } else { puts $DHN_file "p_total D(p_total) p_FENE D(p_FENE) p_lj D(p_lj) p_ideal p_osmotic" }
    flush $DHN_file
}


# Random number generator setup
#############################################################

if { [llength $random_seeds] > 0 } { eval t_random seed $random_seeds }


# Particle & interaction setup
#############################################################

set i 1
foreach pack_i $packing {
  set n_p_i $N_P;
  set mpc_i $MPC;
  set int_time_i $int_time; 
  set int_step_i $int_step
  set n_part $N_T
  set name_i "$name[format %02d $i]"

  set density [expr $pack_i*(6/[PI])]
  set a_cube  [expr pow($N_T/$density,1./3.)]
  setmd box_l $a_cube $a_cube $a_cube
  set bond_l  [expr sqrt(3*[sqr [expr $a_cube/4.]])/($mpc_i+1)]

  puts "\n====================\n=== System [format %2d $i]/[llength $packing] ===\n===================="
  puts "\nSimulate a Diamond Hydrogel Network with $N_node tetra-functional nodes connected by $n_p_i polymer chains with $mpc_i monomers each"
  puts "    in a cubic simulation box of length [setmd box_l] at a density $density which corresponds to a network packing fraction $pack_i."

#  if { ![file exists "$name_i$ident.wrm" ] && ![file exists "$name_i$ident.crl" ] && ![file exists "$name_i$ident.end" ] } {

  inter 0 HARMONIC  $fene_k  $fene_r
#	inter 0 FENE      $fene_k  $fene_r
  for {set ia1 0} {$ia1<3} {incr ia1} { 
     for {set ia2 0} {$ia2<3} {incr ia2} { 
      inter $ia1 $ia2 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
    }
  }
  puts "Interaction Bonds: [inter 0]"
 	puts -nonewline "Creating diamond polymers of (initial) bond-length $bond_l... "; flush stdout
  puts "Done ([diamond $a_cube $bond_l $mpc_i counterions $N_CI charges 1. 1. -1.])."
	
 	if { $vmd_output=="yes" } {
    puts -nonewline "Write psf and pdb for VMD connection... "; flush stdout
    prepare_vmd_connection "$name_i$ident"; puts "Done."
  } 
# } 


#############################################################
#  Warm-up Integration (with capped LJ-potential)           #
#############################################################

	puts -nonewline "\nStart warm-up integration (capped LJ-interactions) for maximal [expr $warm_step*$warm_loop] timesteps in $warm_loop loops; "
	puts "stop if minimal distance is larger than $min_dist after at least [expr $warm_step*$min_loop] timesteps."
        setmd time 0;
        set tmp_cap $warm_cap1;
        inter ljforcecap $tmp_cap;
        inter coulomb 0.0
	set obs_file [open "$name_i$ident.obs1" "w"]
	puts $obs_file "t mindist re rg rh Temp"
	puts $obs_file "[setmd time] [analyze mindist] [analyze re 8 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp]"
	puts -nonewline "    Analysis at t=[setmd time]: mindist=[analyze mindist], "
	puts "re=[lindex [analyze re] 0], rg=[lindex [analyze rg] 0], rh=[lindex [analyze rh] 0], T=[setmd temp]."
	for { set j 0 } { $j < $warm_loop } { incr j } {
	    integrate $warm_step;
      set tmp_dist [analyze mindist]
	    if { $vmd_output=="yes" } { imd positions }
	    puts -nonewline "    \[$i\] Step [expr ($j+1)*$warm_step]/[expr $warm_step*$warm_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "LJ's cap = $tmp_cap, Temp = $tmp_Temp"; flush stdout
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp"
	    puts -nonewline ", mindist=[analyze mindist], re=[lindex [analyze re] 0], rg=[lindex [analyze rg] 0], rh=[lindex [analyze rh] 0]...\r"
	    flush stdout
	    if { ($tmp_dist >= $min_dist) && ($j > $min_loop) } { break }
	    inter forcecap $tmp_cap; set tmp_cap [expr $tmp_cap + $warm_incr]
	}
	flush $obs_file; close $obs_file

#############################################################
#      Integration                                          #
#############################################################

	setmd time 0; set int_loop [expr int($int_time_i/([setmd time_step]*$int_step_i)+0.56)]; set tmp_step 0
	puts "\nStart integration (full interactions) with timestep [setmd time_step] until time t>=$int_time_i (-> $int_loop loops). "
	puts -nonewline "    Activating electrostatics... "; flush stdout
	puts "Done.\n[inter coulomb $bjerrum p3m tune accuracy $accuracy]"
	puts -nonewline "Remove capping of LJ-interactions... "; flush stdout; inter forcecap 0; puts "Done."
	set sfx "[expr int(ceil(log10($int_loop*$int_step_i)))+1]d"
  set tmp_start 0; set obs_file [open "$name_i$ident.obs2" "w"]
  set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 0]
  puts $obs_file "t mindist re dre re2 dre2 rg drg rg2 drg2 rh drh Temp p p2 ideal pid FENE pf pf2 lj plj plj2 coulomb pc pc2"
  puts $obs_file "[setmd time] [analyze mindist] [analyze re 8 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp] $ptot"
  puts -nonewline "    Analysis at t=[setmd time]: mindist=[analyze mindist], "
  puts "re=[lindex [analyze re] 0], rg=[lindex [analyze rg] 0], rh=[lindex [analyze rh] 0], T=[setmd temp], p=$p1."
  analyze append
  #checkpoint_set "$name_i$ident.[eval format %0$sfx 0]" "all" "tmp_step"
  for { set j $tmp_start } { $j < $int_loop } { incr j } {
    integrate $int_step_i; set tmp_dist [analyze mindist]
    if { $vmd_output=="yes" } { imd positions }
    set tmp_step [expr ($j+1)*$int_step_i]
    puts -nonewline "    \[$i\] Step $tmp_step/[expr $int_step_i*$int_loop] (t=[setmd time]): "; flush stdout
    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "Temp = $tmp_Temp"; flush stdout
    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 0]
    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp $ptot"
    set tmp_conf [analyze append]; flush $obs_file
    puts -nonewline ", mindist=[analyze mindist], re=[lindex [analyze re] 0], "
    puts -nonewline "rg=[lindex [analyze rg] 0], rh=[lindex [analyze rh] 0], p=$p1...\r"; flush stdout 
	}
	puts -nonewline "\n    Integration complete; saving checkpoint to '$name_i$ident.end'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.end" "-" "-"; puts "Done."; close $obs_file

	puts -nonewline "\nFinished with current system; "

  puts -nonewline "Cleaning up for next system... "; flush stdout; 
#  part deleteall;
  analyze remove;
  setmd time 0;
  inter coulomb 0.0; 
  incr i; 
  puts "Done.\n"
}
}



# Wrapping up
puts -nonewline "Closing files... "; flush stdout
close $DHN_file
puts -nonewline "Compressing outputs to save disk-space... "; flush stdout
eval exec gzip *.obs2 *.ps *00
puts "Done."
puts "\nThe Diamond Hydrogel Networks Testcase is now complete.\nThanks for watching, and Good Night!\n"
exit
