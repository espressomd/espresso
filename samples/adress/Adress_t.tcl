######################################################################
#                                                                    #
#  System: AdResS simulation of a Liquid of tetrahedral molecules    #
#                                                                    #
######################################################################
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
puts "=              Script:  Adress_t.tcl                  ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################
source tetrahedral.tcl

set n_mol 2520
# System identification: 
set name  "mono"
set ident "_adress_$n_mol"
set sfx "8d"

# vmd options
set vmd_output "no"

# basic parameters
set box_x 36
set box_y 20
set box_z 20

# particle density
set tetra_density [expr double($n_mol)/double($box_x*$box_y*$box_z)]
setmd box_l $box_x $box_y $box_z

# adress setup
set hyb_width 12.0
set ex_width [expr ($box_x-2.0*$hyb_width)/2.0]
set topo_kind 2
set ex_center [expr 0.5*$box_x]
set wf 0
#adress set topo 1 width 1.0
adress set topo $topo_kind width [expr $ex_width*0.5] $hyb_width center x $ex_center wf $wf 

#thermodynamic force
puts "Setting up thermodynamic forces"
set s_pref 0.0
thermodynamic_force 1 "thermo_force.tab" $s_pref
puts "Done"

#############################################################
# Integration parameters                                    #
#############################################################
           
setmd time_step 0.005
setmd skin      0.4

# langevin thermostat
set system_temp 1.0
set gamma 0.5
thermostat langevin $system_temp $gamma

puts "Thermostat has been set"
puts "[thermostat]"

# warmup integration (with capped LJ potential)
set warm_steps   1000
set warm_n_times 20
# do the warmup until the particles have at least the distance min_dist
set min_dist   0.9

# integration and equilibration
set int_steps    1000
set int_n_times  1000
set equi_steps   1000
set equi_n_times 100

#############################################################
# Other parameters                                          #
#############################################################
set tcl_precision 6
set checkpoint  $int_steps
      
# Initial Bond length
set bond_length   1.0
           
# parameters for density profile and rdf calculation
set nbins_dpr 1000
set nbins_rdf 100
set rmin 0
set rmax 10
set dir 0
set type_tetra 1

#####################           
# Interaction setup #
#####################
           
#bonds
inter 0 tabulated bond bond_tetra.tab

# atomistic : repulsive Lennard Jones potential
set lj_eps     1.0
set lj_sig     1.0
set lj_cut     [expr 1.122462048309373*$lj_sig]
set lj_shift   0.25

inter 0 0 lennard-jones $lj_eps $lj_sig $lj_cut $lj_shift 0

#coarse-grained : tabulated potential
inter 1 1 tabulated cg_tetra.tab

##################
# Particle setup #
##################
           
set n_part_per_mol 4
set real_n_part_per_mol 5
           
set volume [expr $box_x*$box_y*$box_z]
set n_ats [expr int($n_mol*$n_part_per_mol)]
set n_tetra [expr $n_ats+$n_mol]

#particle mass
set mass_ex 1.0
set mass_cg 4.0   

# creation of molecules
           
for {set i 0} { $i < $n_mol } {incr i} {
    set posx [expr $box_x*[t_random]]
    set posy [expr $box_y*[t_random]]
    set posz [expr $box_z*[t_random]]
    
    create_pyramide $i $posx $posy $posz $bond_length $n_part_per_mol $mass_ex
    
    #creates com particles
    set cm_id  [expr $i*($n_part_per_mol+1) + $n_part_per_mol]
    part $cm_id pos $posx $posy $posz mass $mass_cg type 1 virtual 1
    
    set mol 0

    set j [expr ($n_part_per_mol+1)*$i]
    
    lappend mol $cm_id
   
    for {set k $j} { $k < [expr $j+4] } {incr k} {
	lappend mol $k
    }   
    lappend topo $mol
}

puts "Syncronizing topology..."
eval analyze set $topo 
analyze set topo_part_sync

# print topology
set topo_file [open "$name$ident.top" w]
puts $topo_file "$topo"
close $topo_file

puts "AdResS simulation of a liquid of tetrahedral particles."
puts "$n_mol molecules in a simulation box of dimensions [setmd box_l] at density $tetra_density"
puts "Interactions:\n[inter]"

setmd max_num_cells [expr [setmd n_part]/[setmd n_nodes]]

set systemTime [clock seconds] 
#############################################################
#  Warm-up Integration (with capped LJ and TAB potentials)  #
#############################################################
set n_part_exp [expr int($n_tetra) ]

if { [file exists "$name$ident.wrm" ] } {
    set n_part_exp [expr int($n_tetra ) ]
    set inp [open "$name$ident.wrm" r]
    puts -nonewline "\nSkipping warm-up integration: Existing checkpoint found (currently reading it... "; flush stdout
    while { [set type [blockfile $inp read auto]] != "eof" } { 
	puts "read $type" 
    }
    close $inp; puts "done) with [setmd n_part] particles ($n_part_exp expected)."
    if { $n_part_exp != [setmd n_part] } {
	puts "WARNING: Configuration does not correspond to current case $i!"; exit      
    }
    puts "Syncronizing topology..."
    eval analyze set $topo 
    analyze set topo_part_sync
} else {
    puts "\nStart warmup integration (capped LJ and MORSE interactions):"
    puts "At maximum $warm_n_times times $warm_steps steps"
    puts "Stop if minimal distance is larger than $min_dist"
    setmd time 0.0
    set act_min_dist [analyze mindist 0 0]
    puts "Start with minimal distance $act_min_dist"
    
    #open Observable files
    
    set obs_file [open "$name$ident.wrm.obs" "w"]
    puts $obs_file "t E_kin mindist bond_ave"
    
    # set LJ cap
    set cap 10
    inter forcecap $cap
    puts "Interactions:\n[inter]"
    # Warmup Integration Loop
    set i 0
    while { $i < $warm_n_times && $act_min_dist < $min_dist } {
	
	integrate $warm_steps
		
	puts "Average bond length: [bond_length_ave $n_mol $n_part_per_mol]" 
	
	# Warmup criterion
	set act_min_dist [analyze mindist 0 0]
		
	puts -nonewline "run $i at time=[setmd time] (LJ cap=$cap) min dist = $act_min_dist\r"
	flush stdout
	
	puts $obs_file "[setmd time] [analyze energy_kinetic 1] $act_min_dist [bond_length_ave $n_mol $n_part_per_mol]"
	
	#   Increase LJ cap
	set cap [expr $cap+10]
	inter forcecap $cap
	incr i
    }
    
    # write everything to disk (set checkpoint)
    puts -nonewline "\n    Warm-up complete; saving checkpoint to '$name$ident.wrm'... ";flush stdout
    polyBlockWrite "$name$ident.wrm" {all} {id pos v f type}; puts "Done."
    flush $obs_file; close $obs_file
}

puts ""
puts -nonewline "Remove capping of LJ interactions... "; flush stdout; inter forcecap 0; puts "Done."

puts -nonewline "Switch on capping of LJ interactions... "; flush stdout; inter forcecap 1000000000; puts "Done."

# Just to see what else we may get from the c code
puts "\nro variables:"
puts "cell_grid     [setmd cell_grid]" 
puts "cell_size     [setmd cell_size]" 
puts "local_box_l   [setmd local_box_l]" 
puts "max_cut       [setmd max_cut]" 
puts "max_part      [setmd max_part]" 
puts "max_range     [setmd max_range]" 
puts "max_skin      [setmd max_skin]" 
puts "n_nodes       [setmd n_nodes]" 
puts "n_part        [setmd n_part]" 
puts "n_part_types  [setmd n_part_types]" 
puts "periodicity   [setmd periodicity]" 
puts "transfer_rate [setmd transfer_rate]" 
puts "verlet_reuse  [setmd verlet_reuse]" 
puts "time_step     [setmd time_step]"
puts "adress_vars   [setmd adress_vars]"


#############################################################
#      Equilibration                                        #
#############################################################

if { [file exists "$name$ident.equi.gz" ] } {
    set n_part_exp [expr int($n_tetra)]
    set inp [open "|gzip -cd $name$ident.equi.gz" r]
    puts -nonewline "\nSkipping equilibration: Existing checkpoint found (currently reading it... "; flush stdout
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp; puts "done) with [setmd n_part] particles ($n_part_exp expected)."
    if { $n_part_exp != [setmd n_part] } { puts "WARNING: Configuration does not correspond to current case!"; exit      }
    puts "Syncronizing topology..."
    eval analyze set $topo 
    analyze set topo_part_sync
} else {
    setmd time 0; set tmp_step 0
    
    puts "\nStart equilibration: run $equi_n_times times $equi_steps steps"
    puts "Interactions:\n[inter]"
    
    if { [file exists "$name$ident.equi.chk" ] } {
	puts -nonewline "    Checkpoint found (currently reading it..."; flush stdout
	checkpoint_read "$name$ident.equi"
	set tmp_start [expr int([setmd time]/[setmd time_step]/$equi_steps)]
	if { [expr $tmp_step/$equi_steps] != $tmp_start } { 
	    puts "failed: Checkpoint corrupt, time_step is wrong! Expected: $tmp_start, got: [expr $tmp_step/$equi_steps])"; exit 
	}
	puts "done) at time [setmd time]: Skipping ahead to timestep [expr int($tmp_step+1)] in loop $tmp_start!"
	#open Observable files
	set obs_file [open "$name$ident.equi.obs" "a"]
	
	puts "Syncronizing topology..."
	eval analyze set $topo 
	analyze set topo_part_sync
	
    } else {
	set tmp_start 0
	
	# open Observable files
	set obs_file [open "$name$ident.equi.obs" "w"]
	
	# write observables
	puts $obs_file "t T nEXs nCGs nHYBs"
	checkpoint_set "$name$ident.equi.[eval format %0$sfx 0].gz" "all" "-"
    }
    
    for { set j $tmp_start } { $j < $equi_n_times} { incr j } {
	integrate $equi_steps; set tmp_step [expr ($j+1)*$equi_steps]
	if { $vmd_output=="yes" } { imd positions }
	puts "    Step $tmp_step/[expr $equi_steps*$equi_n_times] (t=[setmd time]): "; flush stdout
	
	set currenttime [clock seconds]
	set timedif [expr $currenttime - $systemTime]
	if {$timedif >=118800} {
	    set retval 7
	    puts ".... Run finished with exit code $retval" ;
	    set out_restart [open "restart_flag.dat" "w"]
	    puts $out_restart "$retval" ;
	    exit $retval
	}   
	
	set temperature [expr (2.0/3.0)*([analyze energy_kinetic 1]+[analyze energy_kinetic 3])/double($n_mol+$n_mon)] 
	
	set n_tetra_reg [n_regime 1 $n_tetra $box_x $hyb_width 0]
	
	# write and print observables
	
	puts $obs_file "[setmd time] $temperature $n_tetra_reg"
	puts "Analysis at t=[setmd time]: T=$temperature"
	puts "EX  : [lindex $n_tetra_reg 0]"
	puts "CG  : [lindex $n_tetra_reg 1]"
	puts "HYB : [lindex $n_tetra_reg 2]"
	
	puts " " 
	
	flush stdout
	flush $obs_file
	
	# set partial checkpoint (will have previous 'configs' by [analyze append] => averages will be correct)
	if { [expr $tmp_step % $checkpoint]==0 } {
	    puts "\r    Step $tmp_step: Checkpoint at time [setmd time]... "; flush stdout; flush $obs_file
	    checkpoint_set "$name$ident.equi.[eval format %0$sfx $tmp_step].gz" 0 "tmp_step" "-"
	}
    }
        
    close $obs_file
    
    # derive ensemble averages
    set avg [calcObsAv $name$ident.equi.obs {1 2 3 4}]
    
    set temp_ave [lindex $avg 2 0]
    set ntetraex_ave [lindex $avg 2 1]
    set ntetracg_ave [lindex $avg 2 2]
    set ntetrahy_ave [lindex $avg 2 3]
    
    set temp_err [lindex $avg 3 0]
    set ntetraex_err [lindex $avg 3 1]
    set ntetracg_err [lindex $avg 3 2]
    set ntetrahy_err [lindex $avg 3 3]
    
    puts "" 
    puts "Quantities averaged over [lindex $avg 0] configurations:"
    puts ""
    puts "Average temperature: temp_ave +- $temp_err"
    puts "Average number of explicit molecules       : $ntetraex_ave +- $ntetraex_err"
    puts "Average number of coarse-grained molecules : $ntetracg_ave +- $ntetracg_err"
    puts "Average number of hybrid molecules         : $ntetrahy_ave +- $ntetrahy_err"
    puts ""
    
    set obs_file_2 [open "$name$ident.equi.ave" "w"]
    puts $obs_file_2 "adress_vars   [setmd adress_vars]"
    puts $obs_file_2 "$n_mol tetrahedral molecules; box size : [setmd box_l]"
    puts $obs_file_2 ""
    puts $obs_file_2 "Average temperature: temp_ave +- $temp_err"
    puts $obs_file_2 "Average number of explicit molecules       : $ntetraex_ave +- $ntetraex_err"
    puts $obs_file_2 "Average number of coarse-grained molecules : $ntetracg_ave +- $ntetracg_err"
    puts $obs_file_2 "Average number of hybrid molecules         : $ntetrahy_ave +- $ntetrahy_err"
    puts $obs_file_2 ""  
    close $obs_file_2

    # write everything to disk (set checkpoint)
    # (the whole configs-array is not included here for space constraints (it may exceed 1700MB),
    #  it is however stored fractionally in the partial checkpoints, so use 'checkpoint_read' to restore it)
    puts -nonewline "\n    Equilibration complete; saving checkpoint to '$name$ident.equi.gz'... ";flush stdout
    polyBlockWrite "$name$ident.equi.gz" {all} {id pos v f type}; puts "Done."
}


#############################################################
#      Integration                                          #
#############################################################

if { [file exists "$name$ident.end.gz" ] } {
    set n_part_exp [expr int($n_tetra) ]
    set inp [open "|gzip -cd $name$ident.end.gz" r ]
    puts -nonewline "\nSkipping integration: Existing checkpoint found (currently reading it... "; flush stdout
    while { [blockfile $inp read auto] != "eof" } {}
    close $inp; puts "done) with [setmd n_part] particles ($n_part_exp expected)."
    if { $n_part_exp != [setmd n_part] } { puts "WARNING: Configuration does not correspond to current case!"; exit
    }
    puts "Syncronizing topology..."
    eval analyze set $topo 
    analyze set topo_part_sync
} else {
    setmd time 0; set tmp_step 0
    
    puts "\nStart integration: run $int_n_times times $int_steps steps"
    puts "Interactions:\n[inter]"
    
    if { [file exists "$name$ident.chk" ] } {
	puts -nonewline "    Checkpoint found (currently reading it..."; flush stdout
	checkpoint_read "$name$ident"
	set tmp_start [expr int([setmd time]/[setmd time_step]/$int_steps)]
	if { [expr $tmp_step/$int_steps] != $tmp_start } { 
	    puts "failed: Checkpoint corrupt, time_step is wrong! Expected: $tmp_start, got: [expr $tmp_step/$int_steps])"; exit 
	}
	puts "done) at time [setmd time]: Skipping ahead to timestep [expr int($tmp_step+1)] in loop $tmp_start!"
	set obs_file [open "$name$ident.obs" "a"]
	puts "Syncronizing topology..."
	eval analyze set $topo 
	analyze set topo_part_sync
    }  else {
	set tmp_start 0
	
	set obs_file [open "$name$ident.obs" "w"]
	set temperature [expr (2.0/3.0)*([analyze energy_kinetic 1]+[analyze energy_kinetic 3])/double($n_mol+$n_mon)] 
	
	set n_tetra_reg [n_regime 1 [expr $n_tetra] $box_x $hyb_width 0]
	
	#write observables
	puts $obs_file "t T nEXs nCGs nHYs"
	puts $obs_file "[setmd time] $temperature $n_tetra_reg"
	puts "Analysis at t=[setmd time]: T=$temperature"
	puts "EX  : [lindex $n_tetra_reg 0]"
	puts "CG  : [lindex $n_tetra_reg 1]"
	puts "HYB : [lindex $n_tetra_reg 2]"
	puts " " 
	
	analyze append; checkpoint_set "$name$ident.[eval format %0$sfx 0].gz" 0 "tmp_step"
    }
    
    for { set j $tmp_start } { $j < $int_n_times} { incr j } {
	integrate $int_steps; set tmp_step [expr ($j+1)*$int_steps]
	if { $vmd_output=="yes" } { imd positions }
	puts "    Step $tmp_step/[expr $int_steps*$int_n_times] (t=[setmd time]): "; flush stdout
	
	set currenttime [clock seconds]
	set timedif [expr $currenttime - $systemTime]
	if {$timedif >=118800} {
	    set retval 7
	    puts ".... Run finished with exit code $retval" ;
	    set out_restart [open "restart_flag.dat" "w"]
	    puts $out_restart "$retval" ;
	    exit $retval
	}   
	
	set temperature [expr (2.0/3.0)*([analyze energy_kinetic 1])/double($n_mol)] 
	set n_sol_reg [n_regime 1 [expr $n_solv+$n_solu] $box_x $hyb_width 0]
	
	#write observables
	puts $obs_file "[setmd time] $temperature $n_tetra_reg"
	puts "Analysis at t=[setmd time]: T=$temperature"
	puts "EX  : [lindex $n_tetra_reg 0]"
	puts "CG  : [lindex $n_tetra_reg 1]"
	puts "HYB : [lindex $n_tetra_reg 2]"
	puts " " 
	
	flush stdout;analyze append; flush $obs_file
	
	# set partial checkpoint (will have previous 'configs' by [analyze append] => averages will be correct)
	if { [expr $tmp_step % $checkpoint]==0 } {
	    puts "\r    Step $tmp_step: Checkpoint at time [setmd time]... "; flush stdout
	    checkpoint_set "$name$ident.[eval format %0$sfx $tmp_step].gz" 0 "tmp_step" "-"
	}
    }
    close $obs_file
    
    puts "Calculating ensemble averages..."
    # derive ensemble averages
    set avg [calcObsAv $name$ident.obs {1 2 3 4}]
    
    set temp_ave [lindex $avg 2 1]
    set ntetraex_ave [lindex $avg 2 2]
    set ntetracg_ave [lindex $avg 2 3]
    set ntetrahy_ave [lindex $avg 2 4]
    
    set temp_err [lindex $avg 3 1]
    set ntetraex_err [lindex $avg 3 2]
    set ntetracg_err [lindex $avg 3 3]
    set ntetrahy_err [lindex $avg 3 4]
    
    puts "" 
    puts "Quantities averaged over [lindex $avg 0] configurations:"
    puts "Average temperature: $temp_ave +- $temp_err"
    puts "Average number of explicit molecules       : $ntetraex_ave +- $ntetraex_err"
    puts "Average number of coarse-grained molecules : $ntetracg_ave +- $ntetracg_err"
    puts "Average number of hybrid molecules         : $ntetrahy_ave +- $ntetrahy_err"
    puts ""
    
    set obs_file_2 [open "$name$ident.ave" "w"]
    puts $obs_file_2 "adress_vars   [setmd adress_vars]"
    puts $obs_file_2 "$n_mol tetrahedral molecules; box size : [setmd box_l]"
    puts $obs_file_2 ""
    puts $obs_file_2 "Average temperature: temp_ave +- $temp_err"
    puts $obs_file_2 "Average number of explicit molecules       : $ntetraex_ave +- $ntetraex_err"
    puts $obs_file_2 "Average number of coarse-grained molecules : $ntetracg_ave +- $ntetracg_err"
    puts $obs_file_2 "Average number of hybrid molecules         : $ntetrahy_ave +- $ntetrahy_err"
    puts $obs_file_2 ""  
    close $obs_file_2
    
    #radial distribution function
    set mol_rdf [analyze <rdf> 1 1 $rmin $rmax $nbins_rdf $int_n_times]
    set obs_rdf [open "$name$ident.rdf" "w"]
    for {set gg 0 } {$gg < $nbins_rdf } {incr gg } {
	puts $obs_rdf "[lindex [lindex [lindex $mol_rdf 1 ] $gg] 0] [lindex [ lindex [lindex $mol_rdf 1 ] $gg] 1 ]"
    }
    close $obs_rdf
    
    #density profiles
    set x_dens_prof [analyze <density_profile> $nbins_dpr $tetra_density  $dir $int_n_times $type_tetra ] 
    set obs_dprof [open "$name$ident.dpr" "w"]
    for {set gg 0 } {$gg < $nbins_dpr } {incr gg } {
	puts $obs_dprof "[lindex [lindex [lindex $x_dens_prof 0] $gg ] 0] [lindex [lindex [lindex $x_dens_prof 0]  $gg ] 1]"
    }
    close $obs_dprof
        
    puts -nonewline "\n    Integration complete; saving checkpoint to '$name$ident.end.gz'... ";flush stdout
    polyBlockWrite "$name$ident.end.gz" {all} {id pos v f type}; puts "Done."
}

# terminate program
puts "\n\nFinished"
set retval 0
puts ".... Run finished with exit code $retval" ;
set out_restart [open "restart_flag.dat" "w"]
puts $out_restart "$retval" ;
exit $retval
