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
#
#  This file is part of the internal section of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It should never be given outside of the institute without the explicit approval of the author.
#  It is nonetheless subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#  


puts "Program Information: \n[code_info]\n"


#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "analysis_system"
set ident ".data"


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
set lj1_shift   0.25
set lj1_off     0.0

# attractive FENE
set fene_k      30.0
set fene_r      1.5


# Integration parameters
#############################################################

setmd time_step 0.006
setmd skin      0.4
setmd gamma     0.5
setmd temp      1.0

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

set random_seeds  { }
set obs           { mindist re rg rh g123 }
set re2           { 5.2  13.1 29.7 37.8 46.7 82.7 118.3 163.8 263.8 250.7 300.3 }
set rg2           { 0.92  2.2  5.0  6.3  7.7 13.3  20.1  27.5  42.5  46.1  53.6 }
set pKG           { 5.55 5.20 5.04 4.97 4.99 4.93  4.90  4.90  4.87  4.88  4.84 }
set checkpoint    100000



#############################################################
#  Setup System                                             #
#############################################################


# Random number generator setup
#############################################################

if { [llength $random_seeds] > 0 } { eval t_random seed $random_seeds }


# Interaction & Particle setup
#############################################################

set i 1
foreach n_p_i $N_P  mpc_i $MPC  box_l_i $box_l  int_time_i $int_time  rg2_i $rg2  int_step_i $int_step  pKG_i $pKG {
    if {$n_p_i != 20 || $mpc_i != 30} { continue }

    setmd box_l $box_l_i $box_l_i $box_l_i

    inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off
    inter 0   FENE          $fene_k  $fene_r

    puts "\n====================\n=== System [format %2d $i]/[llength $N_P] ===\n===================="
    puts "\nSimulate $n_p_i polymer chains with $mpc_i monomers of (initial) bond-length $bond_l each" 
    puts "    in a cubic simulation box of length [setmd box_l] at density $density with gamma [setmd gamma] and temp [setmd temp]."
    puts "Interactions:\n   [inter]"
    set n_part [expr $n_p_i*$mpc_i]; set name_i "$name\_BB"
    # set name_i "$name[format %02d $i]"
    if { ![file exists "$name_i$ident.wrm" ] && ![file exists "$name_i$ident.end" ] } {
	puts -nonewline "Creating LJ-polymers... "; flush stdout
	polymer $n_p_i $mpc_i $bond_l mode SAW $shield
	puts "Done."
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
	setmd time 0; set tmp_cap $warm_cap1; inter ljforcecap $tmp_cap
	set obs_file [open "$name_i$ident.obs1" "w"]
	puts $obs_file "t mindist re rg rh Temp"
	puts $obs_file "[setmd time] [analyze mindist] [analyze re 0 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp]"
	puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp]."
	for { set j 0 } { $j < $warm_loop } { incr j } {
	    integrate $warm_step; set tmp_dist [analyze mindist]
	    puts -nonewline "    \[$i\] Step [expr ($j+1)*$warm_step]/[expr $warm_step*$warm_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "LJ's cap = $tmp_cap, Temp = $tmp_Temp"; flush stdout
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp"
	    puts -nonewline ", mindist=[analyze mindist], re=[lindex [analyze re] 0], rg=[lindex [analyze rg] 0], rh=[lindex [analyze rh] 0]...\r"; flush stdout
	    if { $tmp_dist >= $min_dist } { break }
	    inter ljforcecap $tmp_cap; set tmp_cap [expr $tmp_cap + $warm_incr]
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
#if {$int_step_i == 500 && $int_time_i == 14400} { set int_loop 2000 }
if {$int_step_i == 500 && $int_time_i == 14400} { set int_loop 10; set int_step_i 100000 }
setmd temp 0.0; setmd gamma 0.0
	puts -nonewline "\nStart integration (full interactions) with timestep [setmd time_step] until time t>=$int_time_i (-> $int_loop loops); "
	puts "aiming for re = [expr sqrt([lindex $re2 [expr $i-1]])], rg = [expr sqrt([lindex $rg2 [expr $i-1]])], and p = $pKG_i."
	puts -nonewline "    Remove capping of LJ-interactions... "; flush stdout; inter ljforcecap 0; puts "Done."
	set sfx "[expr int(ceil(log10($int_loop*$int_step_i)))+1]d"
	if { [file exists "$name_i$ident.chk" ] } {
	    puts -nonewline "    Checkpoint found (currently reading it... "; flush stdout
	    checkpoint_read "$name_i$ident"
	    set tmp_start [expr int([setmd time]/[setmd time_step]/$int_step_i)]
	    if { [expr $tmp_step/$int_step_i] != $tmp_start } { 
		puts "failed: Checkpoint corrupt, time_step is wrong! Expected: $tmp_start, got: [expr $tmp_step/$int_step_i])"; exit 
	    }
	    puts "done) at time [setmd time]: Skipping ahead to timestep [expr int($tmp_step+1)] in loop $tmp_start!"
	    set obs_file [open "$name_i$ident.obs2" "a"]; analyze set chains 0 $n_p_i $mpc_i
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 0]
	    puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp], p=$p1."
	} else {
	    set tmp_start 0; set obs_file [open "$name_i$ident.obs2" "w"]
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 0]
	    puts $obs_file "t mindist re dre re2 dre2 rg drg rg2 drg2 rh drh Temp p p2 ideal pid FENE pf pf2 lj plj plj2"
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re 0 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp] $ptot"
	    puts "    Analysis at t=[setmd time]: rg^2=[lindex [analyze rg] 2], re^2=[lindex [analyze re] 2], rh=[lindex [analyze rh] 0], p=$p1, mindist=[analyze mindist], T=[setmd temp]."
	    set observables ""; set listables ""; set sys_obs [list box_l gamma periodicity skin temperature time time_step]
	    analyze append; checkpoint_set "$name_i$ident.[eval format %0$sfx 0].gz" "all" "observables listables" "all" "$sys_obs" "-"
	}
	for { set j $tmp_start } { $j < $int_loop } { incr j } {
	    integrate $int_step_i; set tmp_dist [analyze mindist]
	    set tmp_step [expr ($j+1)*$int_step_i]
	    puts -nonewline "    \[$i\] Step $tmp_step/[expr $int_step_i*$int_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "Temp = $tmp_Temp"; flush stdout
	    set ptot [eval concat [eval concat [analyze pressure]]]; set p1 [lindex $ptot 0]
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp $ptot"
	    set tmp_conf [analyze append]; flush $obs_file
	    # set partial checkpoint (will have previous 'configs' by [analyze append] => averages will be correct)
	    if { [expr $tmp_step % $checkpoint]==0 } {
		puts -nonewline "\r    \[$i\] Step $tmp_step: Checkpoint at time [setmd time]... "; flush stdout; flush $obs_file
		set observables ""; set listables ""; set sys_obs [list box_l gamma periodicity skin temperature time time_step]
		set mindist1 [analyze mindist]; set mindist2 [analyze mindist 0 0]; lappend observables $mindist1 $mindist2
		set nbhood [analyze nbhood 13 2.5]; lappend listables $nbhood
		set distto [analyze distto 13]; lappend observables $distto
		set energy [analyze energy total]; lappend observables $energy
		set pressure [analyze pressure total]; lappend observables $pressure
		set re [analyze re]; set re_av [analyze <re>]; lappend listables $re $re_av
		set rg [analyze rg]; set rg_av [analyze <rg>]; lappend listables $rg $rg_av
		set rh [analyze rh]; set rh_av [analyze <rh>]; lappend listables $rh $rh_av
		set idf [analyze internal_dist]; set idf_av [analyze <internal_dist>]; lappend listables $idf $idf_av
		set bdf [analyze bond_dist index 13]; set bdf_av [analyze <bond_dist> index 13]; lappend listables $bdf $bdf_av
		set bondl [analyze bond_l]; set bondl_av [analyze <bond_l>]; lappend listables $bondl $bondl_av
		set gff [analyze formfactor 1 10 10]; set gff_av [analyze <formfactor> 1 10 10]; lappend listables $gff $gff_av
		set g1v [analyze <g1>]; set g2v [analyze <g2>]; set g3v [analyze <g3>]; lappend listables $g1v $g2v $g3v
		checkpoint_set "$name_i$ident.[eval format %0$sfx $tmp_step].gz" [expr int($checkpoint/$int_step_i)] "observables listables" "all" "$sys_obs" "-"
		puts -nonewline "set (with <re>=[lindex [analyze <re>] 0], <rg>=[lindex [analyze <rg>] 0] averaged over $tmp_conf configurations"
		puts ", <p>=[lindex [nameObsAv $name_i$ident.obs2 p] 1])."
	    } else { puts -nonewline ",  rg^2=[lindex [analyze rg] 2], re^2=[lindex [analyze re] 2], rh=[lindex [analyze rh] 0], p=$p1, mindist=[analyze mindist], T=[setmd temp]...\r"; 
		flush stdout }
	}
	# write everything to disk (set checkpoint)
	# (the whole configs-array is not included here for space constraints (it may exceed 1700MB),
	#  it is however stored fractionally in the partial checkpoints, so use 'checkpoint_read' to restore it)
	puts -nonewline "\n    Integration complete; saving checkpoint to '$name_i$ident.end'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.end" "-" "-"; puts "Done."; close $obs_file

	puts -nonewline "\nFinished with current system: "
# 	# derive ensemble averages
# 	set avg [nameObsAv $name_i$ident.obs2 { Temp mindist p p2 pid pf pf2 plj plj2 }]
# 	set tmp_Temp [lindex $avg 1]; set tmp_min [lindex $avg 2]
# 	set p1 [lindex $avg 3]; set p2 [lindex $avg 4]; set pid [lindex $avg 5]; set p_os [expr $p1/$pid]
# 	set pf1 [lindex $avg 6]; set pf2 [lindex $avg 7]; set plj1 [lindex $avg 8]; set plj2 [lindex $avg 9]
# 	set d_p12 [expr sqrt(abs($p2 - $p1*$p1)/([lindex $avg 0]-1))]
# 	set d_pf12 [expr sqrt(abs($pf2 - $pf1*$pf1)/([lindex $avg 0]-1))]
# 	set d_plj12 [expr sqrt(abs($plj2 - $plj1*$plj1)/([lindex $avg 0]-1))]
# 	set d_pKG [expr ($p1-$pKG_i)/$pKG_i]
# 	set tmp_re [lindex [analyze <re>] 0]; set tmp_rg [lindex [analyze <rg>] 0]
# 	set tmp_reKG [expr sqrt([lindex $re2 [expr $i-1]])]; set tmp_rgKG [expr sqrt([lindex $rg2 [expr $i-1]])]
# 	set tmp_divE [expr ($tmp_re-$tmp_reKG)/$tmp_reKG]; set tmp_divG [expr ($tmp_rg-$tmp_rgKG)/$tmp_rgKG]
# 	set tmp_rat2 [expr $tmp_re*$tmp_re/($tmp_rg*$tmp_rg)]
# 	puts -nonewline "<re> = $tmp_re ([expr 100*$tmp_divE]% -> $tmp_reKG), "
# 	puts -nonewline "<rg> = $tmp_rg ([expr 100*$tmp_divG]% -> $tmp_rgKG), "
# 	puts "<re2>/<rg2> = $tmp_rat2 (RW=6);"
# 	puts "    <Temp> = $tmp_Temp, <p> = $p_os*p_id = $p1+-$d_p12 (=[expr 100*$d_p12/$p1]% error / [expr 100*$d_pKG]% -> $pKG_i)."
## 	# append ensemble averages to .KKG-file
## 	puts -nonewline $KKG_file "$i $n_p_i $mpc_i $box_l_i $int_time_i "
## 	puts -nonewline $KKG_file "$tmp_re $tmp_divE $tmp_reKG $tmp_rg $tmp_divG $tmp_rgKG "
## 	puts $KKG_file "$tmp_rat2 $tmp_Temp $tmp_min $p1 $d_p12 $p_os $pKG_i $d_pKG"; flush $KKG_file
# 	# sort <g1>, <g2>, and <g3> into .g123-file
# 	set outG [open "$name_i$ident.g123" "w"]
# 	for {set gx 1} {$gx<=3} {incr gx} { eval set tmp_g$gx [list [analyze <g$gx>]] }
# 	for {set gt 0} {$gt<[llength $tmp_g1]} {incr gt} { 
# 	    puts $outG "[expr $gt*[setmd time_step]] [lindex $tmp_g1 $gt] [lindex $tmp_g2 $gt] [lindex $tmp_g3 $gt]"
# 	}
# 	close $outG
# 	# look at pressure and internal distances
## 	puts -nonewline $VIR_file "$i $n_p_i $mpc_i $box_l_i $int_time_i "
## 	puts $VIR_file "$p1 $d_p12 $pf1 $d_pf12 $plj1 $d_plj12 $pid $p_os $pKG_i $d_pKG"; flush $VIR_file
# 	set outI [open "$name_i$ident.idf" "w"]; set tmp_idf [analyze <internal_dist>]
# 	for {set gt 0} {$gt<[llength $tmp_idf]} {incr gt} { puts $outI "$gt [lindex $tmp_idf $gt]" }
# 	close $outI
# 	# create gnuplots
# 	puts -nonewline "Creating a gnuplot from current results... "; flush stdout
# 	plotObs $name_i$ident.obs2 {1:6 1:3 1:4 1:5 1:2 1:7} titles {Temp re rg rh mindist p} labels [concat "time (tau)" "$name_i$ident.obs2"]
# 	plotObs $name_i$ident.g123 {1:2 1:3 1:4} titles {<g1> <g2> <g3>} labels [concat "time (tau)" "$name_i$ident.g123"] scale "logscale xy"
# 	plotObs $name_i$ident.idf {1:2} titles {<internal_dist>} labels [concat "|i-j|" "$name_i$ident.idf"] scale "logscale xy"
# 	lappend plotted "$name_i$ident.obs2"; lappend plotted "$name_i$ident.g123"; lappend plotted "$name_i$ident.idf"
	puts "Done."
    }
puts [analyze mindist]
    puts -nonewline "Cleaning up for next system... "; flush stdout; 
    part deleteall; analyze remove; setmd time 0; incr i; puts "Done.\n"
}
# Final gnuplots
# puts -nonewline "Creating a gnuplot of the averaged quantities... "; flush stdout
## plotObs $name$ident.KKG {3:6 3:8 3:9 3:11} titles {"<re>" "reKG" "<rg>" "rgKG"} labels { "monomers per chain" } scale "logscale xy"
## plotObs $name$ident.VIR {3:6 3:8 3:10 3:12 3:13 3:14} titles {"<p>" "<p_FENE>" "<p_lj>" "<p_ideal>" "<p_osmotic>" "<pKG>"} labels { "monomers per chain" } scale "logscale x"
## lappend plotted "$name$ident.KKG"; lappend plotted "$name$ident.VIR"; puts "Done."
# # puts -nonewline "Combining all plots into '$name_i$ident.final.ps'... "; flush stdout
# # plotJoin $plotted "$name_i$ident.final.ps"; puts "Done."
# Wrapping up
## puts -nonewline "Closing files... "; close $KKG_file; close $VIR_file; puts "Done."
puts "\nThe Kremer-Grest-Testcase is now complete.\nThanks for watching, and Good Night!\n"

exit