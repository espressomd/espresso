#############################################################
#                                                           #
# Kremer-Grest's Linear Polymer Melts                       #
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
puts "==================================================="
puts "=              kremerGrest.tcl                    ="
puts "==================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "kremerGrest"
set ident "_s2"

# On 'yes' connects to 'vmd' visualizing current configuration
set vmd_output "no"


# System parameters
#############################################################

set N_P     { 50   25   30   32   20   16   20    20    20    20    100   }
set MPC     { 5    10   20   25   30   50   75    100   150   200   200   }
set bond_l  0.97
set shield  0.2
set box_l   { 6.65 6.65 8.90 9.80 8.90 9.80 12.08 13.30 15.23 16.76 28.66 }
set density 0.85


# Interaction parameters
#############################################################

# repulsive Lennard Jones
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]
set lj1_off     0.0

# attractive FENE
set fene_k      30.0
set fene_r      1.5


# Integration parameters
#############################################################

setmd time_step 0.006
setmd skin      0.4
thermostat langevin 1.0 0.5

# warmup integration (with capped LJ potential) until particles are at least $min_dist apart
set warm_step   200
set warm_loop   300
set warm_cap1   10
set warm_incr   25
set min_dist    0.90

# integration (with full LJ potential) for $int_time
set int_step  {  500  500   500  1000   500  2000  5000  10000 20000 10000 25000 }
set int_time  { 3000 7200 12000 60000 14400 42000 48000 120000 90000 60000 30000 }


# Other parameters
#############################################################

set tcl_precision  6
set random_seeds  { }
set obs           { mindist re rg rh g123 }
set re2           { 5.2  13.1 29.7 37.8 46.7 82.7 118.3 163.8 263.8 250.7 300.3 }
set rg2           { 0.92  2.2  5.0  6.3  7.7 13.3  20.1  27.5  42.5  46.1  53.6 }
set pKG           { 5.55 5.20 5.04 4.97 4.99 4.93  4.90  4.90  4.87  4.88  4.84 }
set checkpoint    100000



#############################################################
#  Setup System                                             #
#############################################################

if { [file exists "$name$ident.KKG"] } {
    set KKG_file [open "$name$ident.KKG" "a"]
    set VIR_file [open "$name$ident.VIR" "a"]
} else {
    set KKG_file [open "$name$ident.KKG" "w"]
    puts $KKG_file "ID N_P MPC box_l int_time re re-reKG reKG rg rg-rgKG rgKG re2/rg2 Temp mindist p D(p) p_osm pKG p-pKG"; flush $KKG_file
    set VIR_file [open "$name$ident.VIR" "w"]
    puts $VIR_file "ID N_P MPC box_l int_time p_total D(p_total) p_FENE D(p_FENE) p_lj D(p_lj) p_ideal p_osmotic pKG p-pKG"; flush $VIR_file
}

# Random number generator setup
#############################################################

if { [llength $random_seeds] > 0 } { eval t_random seed $random_seeds }


# Interaction & Particle setup
#############################################################

set i 1
foreach n_p_i $N_P  mpc_i $MPC  box_l_i $box_l  int_time_i $int_time  rg2_i $rg2  int_step_i $int_step  pKG_i $pKG {
    setmd box_l $box_l_i $box_l_i $box_l_i

    inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off
    inter 0   FENE          $fene_k  $fene_r

    puts "\n====================\n=== System [format %2d $i]/[llength $N_P] ===\n===================="
    puts "\nSimulate $n_p_i polymer chains with $mpc_i monomers of (initial) bond-length $bond_l each" 
    puts "    in a cubic simulation box of length [setmd box_l] at density $density with gamma [setmd gamma] and temp [setmd temp]."
    puts "Interactions:\n   [inter]"
    set n_part [expr $n_p_i*$mpc_i]; set name_i "$name[format %02d $i]"
    if { ![file exists "$name_i$ident.wrm" ] && ![file exists "$name_i$ident.end" ] } {
	puts -nonewline "Creating LJ-polymers... "; flush stdout
	polymer $n_p_i $mpc_i $bond_l mode SAW $shield
	puts "Done."
    }


# Prepare connection to 'vmd'
#############################################################

    if { $vmd_output=="yes" } {
	puts -nonewline "\nWrite psf and pdb for VMD connection... "; flush stdout
	writepsf "$name_i$ident.psf" $n_p_i $mpc_i; writepdb "$name_i$ident.pdb"
	puts -nonewline "Output created, establishing link... "; flush stdout
	for {set port 10000} { $port < 65000 } { incr port } {
	    catch {imd connect $port} res
	    if {$res == ""} break
	}
	if { $port==65000 } { puts "Failed." } else { puts "Done (now listening at port $port)." 
	    puts "    What you have to do now for a VMD connection:"
	    puts "    (1) Start vmd in current directory (best before running the script)."
	    puts "    (2) Enter on vmd command line: 'mol load psf $name_i$ident.psf pdb $name_i$ident.pdb'"
	    set HOSTNAME [exec hostname]
	    puts "    (3) Enter on vmd command line: 'imd connect $HOSTNAME $port'"
	    puts "    (4) To have the chains coloured individually, set 'Coloring-Method' to 'ResName' in the 'Graphics'-menu"
	    imd listen 0
	}
    }


#############################################################
#  Warm-up Integration (with capped LJ-potential)           #
#############################################################

    if { [file exists "$name_i$ident.wrm" ] } {
	set inp [open "$name_i$ident.wrm" r]
	puts -nonewline "\nSkipping warm-up integration: Existing checkpoint found (currently reading it... "; flush stdout
	while { [blockfile $inp read auto] != "eof" } {}
	close $inp; puts "done) with [setmd n_part] particles ([expr $n_p_i*$mpc_i] expected)."
	if { [expr $n_p_i*$mpc_i] != [setmd n_part] } { puts "WARNING: Configuration does not correspond to current case $i!"; exit }
    } else {
	puts -nonewline "\nStart warm-up integration (capped LJ-interactions) for maximal [expr $warm_step*$warm_loop] timesteps in $warm_loop loops; "
	puts "stop if minimal distance is larger than $min_dist."
	setmd time 0; set tmp_cap $warm_cap1; inter forcecap $tmp_cap
	set obs_file [open "$name_i$ident.obs1" "w"]
	puts $obs_file "t mindist re rg rh Temp"
	puts $obs_file "[setmd time] [analyze mindist] [analyze re 0 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp]"
	puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp]."
	for { set j 0 } { $j < $warm_loop } { incr j } {
	    integrate $warm_step; set tmp_dist [analyze mindist]
	    if { $vmd_output=="yes" } { imd positions }
	    puts -nonewline "    \[$i\] Step [expr ($j+1)*$warm_step]/[expr $warm_step*$warm_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/([degrees_of_freedom]/2.0)]; puts -nonewline "LJ's cap = $tmp_cap, Temp = $tmp_Temp"; flush stdout
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp"
	    puts -nonewline ", mindist=[analyze mindist], re=[lindex [analyze re] 0], rg=[lindex [analyze rg] 0], rh=[analyze rh]...\r"; flush stdout
	    if { $tmp_dist >= $min_dist } { break }
	    inter forcecap $tmp_cap; set tmp_cap [expr $tmp_cap + $warm_incr]
	}
	# write everything to disk (set checkpoint)
	puts -nonewline "\n    Warm-up complete; saving checkpoint to '$name_i$ident.wrm'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.wrm" "-"; puts "Done."
	flush $obs_file; close $obs_file
    }
 

#############################################################
#      Integration                                          #
#############################################################

    if { [file exists "$name_i$ident.end" ] } {
	set inp [open "$name_i$ident.end" r]
	puts -nonewline "Skipping integration: Existing checkpoint found (currently reading it... "; flush stdout
	while { [blockfile $inp read auto] != "eof" } {}
	close $inp; puts "done) with [setmd n_part] particles ([expr $n_p_i*$mpc_i] expected)."
	if { [expr $n_p_i*$mpc_i] != [setmd n_part] } { puts "WARNING: Configuration does not correspond to current case $i!"; exit }
    } else {
	setmd time 0; set int_loop [expr int($int_time_i/([setmd time_step]*$int_step_i)+0.56)]; set tmp_step 0
	puts -nonewline "\nStart integration (full interactions) with timestep [setmd time_step] until time t>=$int_time_i (-> $int_loop loops); "
	puts "aiming for re^2 = [lindex $re2 [expr $i-1]], rg^2 = [lindex $rg2 [expr $i-1]], and p = $pKG_i."
	puts -nonewline "    Remove capping of LJ-interactions... "; flush stdout; inter forcecap 0; puts "Done."
	set sfx "[expr int(ceil(log10($int_loop*$int_step_i)))+1]d"
	if { [file exists "$name_i$ident.chk" ] } {
	    puts -nonewline "    Checkpoint found (currently reading it... "; flush stdout
	    checkpoint_read "$name_i$ident"
	    set tmp_start [expr int([setmd time]/[setmd time_step]/$int_step_i)]
	    if { [expr $tmp_step/$int_step_i] != $tmp_start } { 
		puts "failed: Checkpoint corrupt, time_step is wrong! Expected: $tmp_start, got: [expr $tmp_step/$int_step_i])"; exit 
	    }
	    puts "done) at time [setmd time]: Skipping ahead to timestep [expr int($tmp_step+1)] in loop $tmp_start!"
	    set obs_file [open "$name_i$ident.obs2" "a"]
	    set tmp_dist [analyze mindist]; set tmp_re [analyze re 0 $n_p_i $mpc_i]; set tmp_rg [analyze rg]; set tmp_rh [analyze rh]
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
	    puts -nonewline "    Analysis at t=[setmd time]: mindist=$tmp_dist, T=[setmd temp], "
	    puts "re^2=[lindex $tmp_re 2], rg^2=[lindex $tmp_rg 2], rh=[lindex $tmp_rh 0], p=$p1."
	} else {
	    set tmp_start 0; set obs_file [open "$name_i$ident.obs2" "w"]
	    set tmp_dist [analyze mindist]; set tmp_re [analyze re 0 $n_p_i $mpc_i]; set tmp_rg [analyze rg]; set tmp_rh [analyze rh]
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
	    puts $obs_file "t mindist re dre re2 dre2 rg drg rg2 drg2 rh drh Temp pressure p ideal pid F# FENE pf lj# lj# lj plj "
	    puts $obs_file "[setmd time] $tmp_dist $tmp_re $tmp_rg $tmp_rh [setmd temp] $ptot"
	    puts -nonewline "    Analysis at t=[setmd time]: mindist=$tmp_dist, T=[setmd temp], "
	    puts "re^2=[lindex $tmp_re 2], rg^2=[lindex $tmp_rg 2], rh=[lindex $tmp_rh 0], p=$p1."
	    analyze append; checkpoint_set "$name_i$ident.[eval format %0$sfx 0]" "all" "tmp_step"
	}
	for { set j $tmp_start } { $j < $int_loop } { incr j } {
	    integrate $int_step_i; set tmp_step [expr ($j+1)*$int_step_i]
	    if { $vmd_output=="yes" } { imd positions }
	    puts -nonewline "    \[$i\] Step $tmp_step/[expr $int_step_i*$int_loop] (t=[setmd time]): "; flush stdout
	    set tmp_dist [analyze mindist]; set tmp_re [analyze re]; set tmp_rg [analyze rg]; set tmp_rh [analyze rh]
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 1]
	    set tmp_Temp [expr [analyze energy kin]/$n_part/([degrees_of_freedom]/2.0)]
	    puts $obs_file "[setmd time] $tmp_dist $tmp_re $tmp_rg $tmp_rh $tmp_Temp $ptot"
	    puts -nonewline "mindist=$tmp_dist, T=$tmp_Temp"; flush stdout
	    set tmp_conf [analyze append]; flush $obs_file
	    # set partial checkpoint (will have previous 'configs' by [analyze append] => averages will be correct)
	    if { [expr $tmp_step % $checkpoint]==0 } {
		puts -nonewline "\r    \[$i\] Step $tmp_step: Checkpoint at time [setmd time]... "; flush stdout; flush $obs_file
		checkpoint_set "$name_i$ident.[eval format %0$sfx $tmp_step]" [expr int($checkpoint/$int_step_i)] "tmp_step" "-"
		puts -nonewline "set (with <re^2>=[lindex [analyze <re>] 2], <rg^2>=[lindex [analyze <rg>] 2] averaged over $tmp_conf configurations"
		puts ", <p>=[lindex [nameObsAv $name_i$ident.obs2 p] 1])."
	    } else { 
		puts -nonewline ", re^2=[lindex $tmp_re 2], rg^2=[lindex $tmp_rg 2], rh=[lindex $tmp_rh 0], p=$p1...\r"; flush stdout 
	    }
	}
	# write everything to disk (set checkpoint)
	# (the whole configs-array is not included here for space constraints (it may exceed 1700MB),
	#  it is however stored fractionally in the partial checkpoints, so use 'checkpoint_read' to restore it)
	puts -nonewline "\n    Integration complete; saving checkpoint to '$name_i$ident.end'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.end" "-" "-"; puts "Done."; close $obs_file

	puts -nonewline "\nFinished with current system: \n    "; flush stdout
	# derive ensemble averages
	set avg [nameObsAv $name_i$ident.obs2 { Temp mindist p pid pf plj }]
	set tmp_Temp [lindex $avg 1]; set tmp_min [lindex $avg 2]
	set p1 [lindex $avg 3]; set pid [lindex $avg 4]; set p_os [expr $p1/$pid]
	set pf1 [lindex $avg 5]; set plj1 [lindex $avg 6] 
	set d_p12 [lindex $avg 7]; set d_pf12 [lindex $avg 9]; set d_plj12 [lindex $avg 10]
	set d_pKG [expr ($p1-$pKG_i)/$pKG_i]
	set tmp_re [lindex [analyze <re>] 2]; set tmp_rg [lindex [analyze <rg>] 2]
	set tmp_reKG [lindex $re2 [expr $i-1]]; set tmp_rgKG [lindex $rg2 [expr $i-1]]
	set tmp_divE [expr ($tmp_re-$tmp_reKG)/$tmp_reKG]; set tmp_divG [expr ($tmp_rg-$tmp_rgKG)/$tmp_rgKG]
	set tmp_rat2 [expr $tmp_re/$tmp_rg]
	puts -nonewline "<re^2> = $tmp_re ([expr 100*$tmp_divE]% -> $tmp_reKG), "
	puts -nonewline "<rg^2> = $tmp_rg ([expr 100*$tmp_divG]% -> $tmp_rgKG), "
	puts "<re^2>/<rg^2> = $tmp_rat2 (RW=6);"
	puts "    <Temp> = $tmp_Temp, <p> = $p_os*p_id = $p1+-$d_p12 (=[expr 100*$d_p12/$p1]% error / [expr 100*$d_pKG]% -> $pKG_i)."
	# append ensemble averages to .KKG-file
	puts -nonewline $KKG_file "$i $n_p_i $mpc_i $box_l_i $int_time_i "
	puts -nonewline $KKG_file "$tmp_re $tmp_divE $tmp_reKG $tmp_rg $tmp_divG $tmp_rgKG "
	puts $KKG_file "$tmp_rat2 $tmp_Temp $tmp_min $p1 $d_p12 $p_os $pKG_i $d_pKG"; flush $KKG_file
	# sort <g1>, <g2>, and <g3> into .g123-file
	set outG [open "$name_i$ident.g123" "w"]
	for {set gx 1} {$gx<=3} {incr gx} { eval set tmp_g$gx [list [analyze <g$gx>]] }
	for {set gt 0} {$gt<[llength $tmp_g1]} {incr gt} { 
	    puts $outG "[expr $gt*[setmd time_step]] [lindex $tmp_g1 $gt] [lindex $tmp_g2 $gt] [lindex $tmp_g3 $gt]"
	}
	close $outG
	# look at pressure and internal distances
	puts -nonewline $VIR_file "$i $n_p_i $mpc_i $box_l_i $int_time_i "
	puts $VIR_file "$p1 $d_p12 $pf1 $d_pf12 $plj1 $d_plj12 $pid $p_os $pKG_i $d_pKG"; flush $VIR_file
	puts -nonewline "    Analyzing internal distances... "; flush stdout
	set outI [open "$name_i$ident.idf" "w"]; set tmp_idf [analyze <internal_dist>]
	for {set gt 0} {$gt<[llength $tmp_idf]} {incr gt} { puts $outI "$gt [lindex $tmp_idf $gt]" }
	close $outI; puts "Done."
	# create gnuplots
	puts -nonewline "Creating a gnuplot from current results... "; flush stdout
	plotObs $name_i$ident.obs2 {1:13 1:5 1:9 1:11 1:2 1:15} titles {Temp re^2 rg^2 rh mindist p} labels [concat "time (tau)" "$name_i$ident.obs2"]
	plotObs $name_i$ident.g123 {1:2 1:3 1:4} titles {<g1> <g2> <g3>} labels [concat "time (tau)" "$name_i$ident.g123"] scale "logscale xy"
	plotObs $name_i$ident.idf {1:2} titles {<internal_dist>} labels [concat "|i-j|" "$name_i$ident.idf"] scale "logscale xy"
	lappend plotted "$name_i$ident.obs2"; lappend plotted "$name_i$ident.g123"; lappend plotted "$name_i$ident.idf"
	puts "Done."
    }
    puts -nonewline "Cleaning up for next system... "; flush stdout; 
    part deleteall; analyze remove; setmd time 0; incr i; puts "Done.\n"
}
# Final gnuplots
puts -nonewline "Creating a gnuplot of the averaged quantities... "; flush stdout
plotObs $name$ident.KKG {3:6 3:8 3:9 3:11} titles {"<re^2>" "re^2_KG" "<rg^2>" "rg^2_KG"} labels { "monomers per chain" } scale "logscale xy"
plotObs $name$ident.VIR {3:6 3:8 3:10 3:12 3:13 3:14} titles {"<p>" "<p_FENE>" "<p_lj>" "<p_ideal>" "<p_osmotic>" "<pKG>"} labels { "monomers per chain" } scale "logscale x"
lappend plotted "$name$ident.KKG"; lappend plotted "$name$ident.VIR"; puts "Done."
# puts -nonewline "Combining all plots into '$name_i$ident.final.ps'... "; flush stdout
# plotJoin $plotted "$name_i$ident.final.ps"; puts "Done."
# Wrapping up
puts -nonewline "Closing files... "; close $KKG_file; close $VIR_file; puts "Done."
puts "\nThe Kremer-Grest-Testcase is now complete.\nThanks for watching, and Good Night!\n"

exit
