#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Test System 2: Kremer-Grest's Linear Polymer Melts       #
#                                                           #
#                                                           #
#  Created:       21.03.2003 by BAM                         #
#                                                           #
#############################################################

puts " "
puts "==================================================="
puts "=                kremerGrest.tcl                  ="
puts "==================================================="
puts " "

puts "Program Information: \n[code_info]\n"



#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "kremerGrest"
set ident "_t2"

# On 'yes' connects to 'vmd' visualizing current configuration
set vmd_output "yes"


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
set lj1_shift   $lj1_eps

# attractive FENE
set fene_k      7.0
set fene_r      2.0


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
set min_dist    0.92

# integration (with full LJ potential) for $int_time
set int_step    500
set int_time  { 3000 7200 12000 60000 14400 42000 48000 120000 90000 60000 30000 }


# Other parameters
#############################################################

set tcl_precision  6
set random_seeds  { }
set obs           { mindist re rg rh g123 }
set re2           { 5.2  13.1 29.7 37.8 46.7 82.7 118.3 163.8 263.8 250.7 300.3 }
set rg2           { 0.92  2.2  5.0  6.3  7.7 13.3  20.1  27.5  42.5  46.1  53.6 }
set checkpoint    100000



#############################################################
#  Setup System                                             #
#############################################################

if { [file exists "$name$ident.KKG"] } {
    set KKG_file [open "$name$ident.KKG" "a"]
} else {
    set KKG_file [open "$name$ident.KKG" "w"]
    puts $KKG_file "ID N_P MPC box_l int_time re re-reKG reKG rg rg-rgKG rgKG re2/rg2 Temp mindist"; flush $KKG_file
}

# Random number generator setup
#############################################################

if { [llength $random_seeds] > 0 } { eval t_random seed $random_seeds }


# Interaction & Particle setup
#############################################################

set i 1
foreach n_p_i $N_P  mpc_i $MPC  box_l_i $box_l  int_time_i $int_time  rg2_i $rg2  {
    setmd box_l $box_l_i $box_l_i $box_l_i

    inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
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
	setmd time 0; set tmp_cap $warm_cap1; inter ljforcecap $tmp_cap
	set obs_file [open "$name_i$ident.obs1" "w"]
	puts $obs_file "t mindist re rg rh Temp"
	puts $obs_file "[setmd time] [analyze mindist] [analyze re 0 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp]"
	puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp]."
	for { set j 0 } { $j < $warm_loop } { incr j } {
	    integrate $warm_step; set tmp_dist [analyze mindist]
	    if { $vmd_output=="yes" } { imd positions }
	    puts -nonewline "    \[$i\] Step [expr ($j+1)*$warm_step]/[expr $warm_step*$warm_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "LJ's cap = $tmp_cap, Temp = $tmp_Temp"; flush stdout
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp"
	    puts -nonewline ", mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh]...\r"; flush stdout
	    if { $tmp_dist >= $min_dist } { break }
	    inter ljforcecap $tmp_cap; set tmp_cap [expr $tmp_cap + $warm_incr]
	}
	# write everything to disk (set checkpoint)
	puts -nonewline "\n    Warm-up complete; saving checkpoint to '$name_i$ident.wrm'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.wrm" "-"; puts "Done."
	close $obs_file
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
	setmd time 0; set int_loop [expr int($int_time_i/([setmd time_step]*$int_step)+0.56)]
	puts -nonewline "\nStart integration (full interactions) with timestep [setmd time_step] until time t>=$int_time_i (-> $int_loop loops); "
	puts "aiming for re = [expr sqrt([lindex $re2 [expr $i-1]])] and rg = [expr sqrt([lindex $rg2 [expr $i-1]])]."
	puts -nonewline "    Remove capping of LJ-interactions... "; flush stdout; inter ljforcecap 0; puts "Done."
	if { [file exists "$name_i$ident.chk" ] } {
	    set inp [open "$name_i$ident.chk" r]
	    puts -nonewline "    Checkpoint found (currently reading it... "; flush stdout
	    while { [blockfile $inp read auto] != "eof" } {}
	    set tmp_start [expr int([setmd time]/[setmd time_step]/$int_step)]
	    if { [expr $j+1] != $tmp_start } { puts "failed: Checkpoint corrupt, time_step is wrong! Expected: $tmp_start, got: [expr $j+1])"; exit }
	    close $inp; puts "done) at time [setmd time]: Skipping ahead to timestep [expr int([setmd time]/[setmd time_step]+1)] in loop $tmp_start!"
	    set obs_file [open "$name_i$ident.obs2" "a"]; analyze set chains 0 $n_p_i $mpc_i
	    puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp]."
	} else {
	    set tmp_start 0; set obs_file [open "$name_i$ident.obs2" "w"]
	    puts $obs_file "t mindist re rg rh Temp"
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re 0 $n_p_i $mpc_i] [analyze rg] [analyze rh] [setmd temp]"
	    puts "    Analysis at t=[setmd time]: mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh], T=[setmd temp]."
	}
	analyze append
	for { set j $tmp_start } { $j < $int_loop } { incr j } {
	    integrate $int_step; set tmp_dist [analyze mindist]
	    if { $vmd_output=="yes" } { imd positions }
	    puts -nonewline "    \[$i\] Step [expr ($j+1)*$int_step]/[expr $int_step*$int_loop] (t=[setmd time]): "; flush stdout
	    set tmp_Temp [expr [analyze energy kin]/$n_part/1.5]; puts -nonewline "Temp = $tmp_Temp"; flush stdout
	    puts $obs_file "[setmd time] [analyze mindist] [analyze re] [analyze rg] [analyze rh] $tmp_Temp"
	    set tmp_conf [analyze append]
	    # set partial checkpoint (will have previous 'configs' by [analyze append] => averages will be correct)
	    if { [expr (($j+1)*$int_step) % $checkpoint]==0 } {
		puts -nonewline "\r    \[$i\] Step [expr ($j+1)*$int_step]: Checkpoint at time [setmd time]... "; flush stdout
		polyBlockWriteAll "$name_i$ident.chk" "j"
		puts "set (with <re>=[analyze <re>], <rg>=[analyze <rg>] averaged over $tmp_conf configurations)."
	    } else { puts -nonewline ", mindist=[analyze mindist], re=[analyze re], rg=[analyze rg], rh=[analyze rh]...\r"; flush stdout }
	}
	# write everything to disk (set checkpoint)
	puts -nonewline "\n    Integration complete; saving checkpoint to '$name_i$ident.end'... ";flush stdout
	polyBlockWriteAll "$name_i$ident.end" "-"; puts "Done."; close $obs_file

	puts -nonewline "\nFinished with current system; "
	# derive ensemble averages
	lappend what "[calcObsAv $name_i$ident.obs2 { 1 2 3 4 5 } ]"
	set avg [findObsAv { Temp mindist } [lindex $what end]]; set tmp_Temp [lindex $avg 1]; set tmp_min [lindex $avg 2]
	set tmp_re [analyze <re>]; set tmp_rg [analyze <rg>]
	set tmp_reKG [expr sqrt([lindex $re2 [expr $i-1]])]; set tmp_rgKG [expr sqrt([lindex $rg2 [expr $i-1]])]
	set tmp_divE [expr ($tmp_re-$tmp_reKG)/$tmp_reKG]; set tmp_divG [expr ($tmp_rg-$tmp_rgKG)/$tmp_rgKG]
	set tmp_rat2 [expr $tmp_re*$tmp_re/($tmp_rg*$tmp_rg)]
	puts -nonewline "<re> = $tmp_re ([expr 100*$tmp_divE]% -> $tmp_reKG), "
	puts -nonewline "<rg> = $tmp_rg ([expr 100*$tmp_divG]% -> $tmp_rgKG), "
	puts "<re2>/<rg2> = $tmp_rat2 (RW=6), <Temp> = $tmp_Temp."
	# append ensemble averages to .KKG-file
	puts -nonewline $KKG_file "$i $n_p_i $mpc_i $box_l_i $int_time_i "
	puts -nonewline $KKG_file "$tmp_re $tmp_divE $tmp_reKG $tmp_rg $tmp_divG $tmp_rgKG "
	puts $KKG_file "$tmp_rat2 $tmp_Temp $tmp_min"; flush $KKG_file
	# sort <g1>, <g2>, and <g3> into .g123-file
	set outG [open "$name_i$ident.g123" "w"]
	for {set gx 1} {$gx<=3} {incr gx} { eval set tmp_g$gx [list [analyze <g$gx>]] }
	for {set gt 0} {$gt<[llength $tmp_g1]} {incr gt} { 
	    puts $outG "[expr $gt*[setmd time_step]] [lindex $tmp_g1 $gt] [lindex $tmp_g2 $gt] [lindex $tmp_g3 $gt]"
	}
	close $outG 
	# create gnuplots
	puts -nonewline "Creating a gnuplot from current results... "; flush stdout
	plotObs $name_i$ident.obs2 {1:6 1:3 1:4 1:5 1:2} titles {Temp re rg rh mindist} labels [concat "time (tau)" "$name_i$ident.obs2"]
	plotObs $name_i$ident.g123 {1:2 1:3 1:4} titles {<g1> <g2> <g3>} labels [concat "time (tau)" "$name_i$ident.g123"] scale "logscale xy"
	lappend plotted "$name_i$ident.obs2"; lappend plotted "$name_i$ident.g123"
	puts "Done."
    }
    puts -nonewline "Cleaning up for next system... "; flush stdout; 
    part deleteall; analyze remove; setmd time 0; incr i; puts "Done.\n"
}
# Final gnuplots
puts -nonewline "Creating a gnuplot of the averaged quantities... "; flush stdout
plotObs $name$ident.DKG {3:6 3:8 3:9 3:11} titles {"<re>" "reKG" "<rg>" "rgKG"} labels { "monomers per chain" } scale "logscale xy"
lappend plotted "$name$ident.DKG"; puts "Done."
# puts -nonewline "Combining all plots into '$name_i$ident.final.ps'... "; flush stdout
# plotJoin $plotted "$name_i$ident.final.ps"; puts "Done."
# Wrapping up
puts -nonewline "Closing files... "; close $KKG_file; puts "Done."
puts "\nThe Kremer-Grest-Testcase is now complete.\nThanks for watching, and Good Night!\n"
