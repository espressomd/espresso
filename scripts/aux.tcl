#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#  
#############################################################
#                                                           #
# aux.tcl                                                   #
# =======                                                   #
#                                                           #
# Several additional auxiliary functions for Espresso.      #
#                                                           #
# Created:       01.10.2002 by BAM                          #
#                                                           #
#############################################################




#
# timeStamp
# ---------
# 
# Squeezes a prefix into a given filename and adds the 
# current date as postfix before the suffix.
#
# Input:
# - complete path 'destination'
# - 'prefix' to be squeezed between path- and file-name
# - 'postfix' to be added before 'suffix'; 
#   if it is '-1', the current date will be used
# 
# Created:       01.10.2002 by BAM
# 
#############################################################

proc timeStamp { destination prefix postfix suffix } {

    # Customize 'destination' to include informations about the current time and the Deserno-wishlist used
    puts -nonewline "    Building a customized destination file name... "
    flush stdout

    # Squeeze a timestamp just before the 'suffix'
    if { [string compare $postfix -1]==0 } {
	set tmp_var [clock format [clock seconds] -format %y%m%d]
    } else {
	set tmp_var $postfix
    }
    set destination "[string range $destination 0 [string last .$suffix $destination]]$tmp_var.$suffix"

    # Squeeze the 'prefix' in between path and filename
    set tmp_var [expr [string first [lindex [split $destination \"/\"] end] $destination]-1]
    set destination "[string range $destination 0 $tmp_var]$prefix\_[lindex [split $destination \"/\"] end]"
    puts "Done."

    return $destination
}




#
# polyBlockWrite
# --------------
# 
# Writes current polymer configuration to disk,
# including all bonds and interactions there are,
# using Axa's blockfile-format.
#
# Input:
# - complete path 'destination';
#   if the filename ends with '.gz', the file will be compressed
# - a list of 'Espresso' parameters to be saved (out of node_grid|box_l|niatypes|time_step|skin|gamma|bjerrum|...
#   ...p3m_alpha|p3m_r_cut|p3m_mesh|p3m_cao|p3m_epsilon|p3m_mesh_offset|max_num_cells|periodicity);
#   if an empty list '{}' is supplied, no parameters and no interactions are written
#   Default value: All the above mentioned parameters.
# - a string containing which informations (out of pos|type|q|v|f) on the particles should be saved to disk;
#   if an empty string ' "" 'is provided, no particles, and no bonds are written
#   Default vaule: All the above mentioned informations.
# 
# Created:       08.11.2002 by BAM
# 
#############################################################

proc polyBlockWrite { destination {write_param "all"} {write_part "id pos type q v f"} } {

    # Open output-file - compressed, if desired
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" w]
    } else {
	set f [open "$destination" "w"]
    }
    
    # Write parameters and interactions, if desired
    if { "$write_param" != "{}" } {
	foreach j $write_param {
	    blockfile $f write variable $j
	}
	blockfile $f write interactions
	blockfile $f write integrate
	blockfile $f write thermostat
    }
    
    # Write particles and bonds, if desired
    if { "$write_part" != "" } {
	blockfile $f write particles $write_part all
	blockfile $f write bonds all
    }

    # Close file
    flush $f; close $f
}


proc polyBlockWriteAll { destination {tclvar "all"} {cfg "1"} {rdm "random"} } {
    polyBlockWrite $destination
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" a] } else { set f [open "$destination" "a"] }
    # Write tcl-variables, if desired
    if { "$tclvar" != "-" } { foreach j $tclvar { blockfile $f write tclvariable $j } }
    
    # Write seed/status of random number generator, if desired
    if { "$rdm" != "-" } { foreach j $rdm { blockfile $f write $j } }
    flush $f
    # Write stored analysis-configurations, if desired
    if { "$cfg" != "-" }  { blockfile $f write configs }

    # Close file
    flush $f; close $f
}

proc polyBlockWriteTclvar { destination {tclvar "all"} } {
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" a] } else { set f [open "$destination" "a"] }
    # Write tcl-variables, if desired
    if { "$tclvar" != "-" } { foreach j $tclvar { blockfile $f write tclvariable $j } }
    # Close file
    flush $f; close $f
}

# reads the file named '$source' into Espresso
proc polyBlockRead { source } {
    set inp [open "$source" r]; while { [eof $inp] != 1 } { blockfile $inp read auto }; close $inp
}


proc checkpoint_set { destination { cnt "all" } { tclvar "all" } { ia "all" } { var "all" } { ran "all" } } {
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" w]
	set chk [open "[join [lrange [split $destination .] 0 [expr [llength [split $destination .]]-3]] .].chk" a]
    } else { 
	set f [open "$destination" "w"]
	set chk [open "[join [lrange [split $destination .] 0 [expr [llength [split $destination .]]-2]] .].chk" "a"]
    }
    if { "$var" != "-" } { blockfile $f write variable $var }
    if { "$tclvar" != "-" } { foreach j $tclvar { blockfile $f write tclvariable $j } }
    if { "$ia" != "-" } { blockfile $f write interactions; blockfile $f write integrate; blockfile $f write thermostat }
    set part_write "id pos type q v f "
    if { [regexp "ROTATION" [code_info]]} { lappend part_write "quat omega torque " }
    if { [regexp "CONSTRAINTS" [code_info]]} { lappend part_write "fix" }
    if { [regexp "EXTERNAL_FORCES" [code_info]]} { lappend part_write "ext_force" }
    blockfile $f write particles "$part_write"
    blockfile $f write bonds
    if { "$ran" != "-" } { blockfile $f write random }
    if { "$ran" != "-" } { blockfile $f write bitrandom }
    flush $f
    if { "$cnt" == "all" } { blockfile $f write configs } else { blockfile $f write configs $cnt }
    flush $f; close $f
    puts $chk "$destination"; flush $chk; close $chk
    invalidate_system
}


# this command is intended to read all the checkpoints listed in a .chk-file
# parameters are 'read_all_chks' == 1 (if all checkpoints should be read), != 1 (if only the last one should be read),
# 'write_pdb' == 0 (do nothing, ignore 'pdb_sfx'), != 0 (if all [e.g. 459] checkpoints should be converted to $write_pdb0000.pdb ... $write_pdb0459.pdb),
# and 'pdb_sfx' (giving the number of digits to be used in enumbering the .pdb-files)
proc checkpoint_read { origin { read_all_chks 1 } { write_pdb 0 } { pdb_sfx 5 }} {
    if { [file exists "$origin.chk"] } { 
	set chk [open "$origin.chk" "r"] 
    } elseif { 
	[file exists "$origin"] } { set chk [open "$origin" "r"] 
    } else { 
	puts "ERROR: Could not find checkpoint-list $origin!\nAborting..."; exit 
    }
    if { $write_pdb !=0 } { set pdb_ind 0; set pdb_sfx [join [list "%0" $pdb_sfx "d"] ""] }
    if { $read_all_chks } {
	while { [eof $chk]==0 } { 
	    if { [gets $chk source] > 0 } {
		if { [string compare [lindex [split $source "."] end] "gz"]==0 } { 
		    set f [open "|gzip -cd $source" r] 
		} else { 
		    set f [open "$source" "r"] 
		}
		part deleteall
		while { [blockfile $f read auto] != "eof" } {}
		puts -nonewline "."; flush stdout; # puts "read $source"
		close $f
		if { $write_pdb !=0 } { if {$pdb_ind==0} {writepsf "$write_pdb.psf"}; writepdb "$write_pdb[format $pdb_sfx $pdb_ind].pdb"; incr pdb_ind }
	    }
	}
    } else {
	puts " Reading just one chkpnt " 
	set tmp_chk "NA"
	while { [eof $chk]==0 } { 
	    if { [gets $chk source] > 0 } { 
		set tmp_chk $source 
	    }   
	}
	if { $tmp_chk == "NA" } { 
	    puts "ERROR: Didn't find any checkpoints! Aborting..."; exit 
	} else { 
	    set source $tmp_chk 
	}
	if { [string compare [lindex [split $source "."] end] "gz"]==0 } { 
	    set f [open "|gzip -cd $source" r]
	} else { 
	    set f [open "$source" "r"] 
	}
	part deleteall
	while { [blockfile $f read auto] != "eof" } {}
	puts -nonewline "."; flush stdout; # puts "read $source"
	close $f
	if { $write_pdb !=0 } { if {$pdb_ind==0} {writepsf "$write_pdb.psf"}; writepdb "$write_pdb[format $pdb_sfx $pdb_ind].pdb"; incr pdb_ind }
    }
    close $chk
}


proc polyConfMovWriteInit { write prfx polyConfAux } {
    if { $write=="yes" || $write=="movie" } {
	set param [list output_path output_prfx write_param write_part movie_path movie_prfx N_P MPC N_CI N_pS N_nS]
	for {set i 0} {$i < [llength $polyConfAux]} {incr i} { eval set [lindex $param $i] [lindex $polyConfAux $i] }
	set tmp_file "$prfx.gz"; set tmp_file "$output_prfx$tmp_file"
	puts -nonewline "    Saving configuration to $tmp_file... "; flush stdout
	polyBlockWrite "$output_path$tmp_file" "$write_param" "$write_part"
	if { $write=="movie" } {
	    set tmp_file $prfx; set tmp_file "$movie_prfx$tmp_file"
	    puts -nonewline "create movie $tmp_file... "; flush stdout
	    writepsf "$movie_path$tmp_file.psf" $N_P $MPC $N_CI $N_pS $N_nS
	    if {[llength $polyConfAux] > 6} {
		writepdb "$movie_path$tmp_file.pdb"
	    } else { writepdb "$movie_path$tmp_file.pdb" }
	}
	puts "Done."
    }
}

proc polyConfMovWrite { write prfx digit step polyConfAux } {
    if { ($write == "yes" || $write=="movie") } {
	set param [list output_path output_prfx write_param write_part movie_path movie_prfx]
	for {set i 0} {$i < [llength $polyConfAux]} {incr i} { eval set [lindex $param $i] [lindex $polyConfAux $i] }
	set tmp_format "d"; set tmp_format "$prfx%0$digit$tmp_format"
	set tmp_file [eval format "$output_prfx$tmp_format.gz" $step]
	puts -nonewline " (saving $tmp_file - "; flush stdout
	polyBlockWrite "$output_path$tmp_file" "$write_param" "$write_part"
	if { $write=="movie" } {
	    set tmp_file [eval format "$movie_prfx$tmp_format" $step]
	    puts -nonewline "$tmp_file - "; flush stdout
	    if {[llength $polyConfAux] > 6} {
		writepdb "$movie_path$tmp_file.pdb"
	    } else { writepdb "$movie_path$tmp_file.pdb" }
	}
	puts -nonewline "done)"; flush stdout
    }
}


proc analysisInit { stat stat_out N_P MPC simtime { noted "na" } { notedD "na" } } {
    if {[llength $stat]>0} {
	puts -nonewline "    Analyzing at t=$simtime: "; flush stdout
	puts -nonewline "Init all ([analyze set chains 0 $N_P $MPC]): "; flush stdout
	puts -nonewline $stat_out "t "; set tmp_stat "$simtime "
	foreach i $stat {
	    if { $i == "g123" } { puts -nonewline "init g123 ([analyze g123 -init 0 $N_P $MPC ]), "; flush stdout }
	    set tmp_var [analyze $i]
	    puts -nonewline $stat_out "$i "; set tmp_stat "$tmp_stat $tmp_var"
	    puts -nonewline "$i=$tmp_var; "; flush stdout
	}
	if { $noted != "na" } { puts -nonewline $stat_out "$noted "; set tmp_stat "$tmp_stat $notedD" }
	puts $stat_out "\n$tmp_stat"; flush $stat_out
	puts "Done."
    }
}

proc analysis { stat stat_out N_P MPC simtime { noted "na" } } {
    puts -nonewline $stat_out "$simtime "
    foreach i $stat {
	set tmp_stat [analyze $i]
	puts -nonewline ", $i=$tmp_stat"; flush stdout
	puts -nonewline $stat_out "$tmp_stat "
    }
    if { $noted != "na" } { puts -nonewline $stat_out "$noted " }
    puts $stat_out " "; flush $stat_out
}



#
# tune_cells
# ----------
# 
# Tunes the skin and cell-grid to balance verlet_reuse vs. 
# organisational overhead.
#
#############################################################

proc tune_cells { { int_steps 1000 } { min_cells 1 } { max_cells 30 } { tol_cells 0.099 } } {
    set msg_string "Tuning current system with [setmd n_part] particles and initial skin=[setmd skin] and cell_grid=[setmd cell_grid] to:"
    array set int_times [list "0 0 0" 0.0]
    proc eval_cells { i num_cells } {
	upvar int_steps int_steps  int_times int_times  msg_string msg_string
	setmd max_num_cells [expr int(pow($num_cells,3))]
	setmd skin 0.0
	set msg_string "$msg_string\n    run [expr $i+1] (t=[setmd time]; max_cells=[format %.3f $num_cells],"
	set msg_string "$msg_string skin=[format %.3g [setmd skin [expr 0.999*[setmd max_skin]]]],"
	set msg_string "$msg_string grid=[setmd cell_grid], cell_size=[format %.3f [lindex [set tmp_cell [setmd cell_size]] 0]]"
	set msg_string "$msg_string [format %.3f [lindex $tmp_cell 1]] [format %.3f [lindex $tmp_cell 2]], max_cut=[format %.3f [setmd max_cut]],"
	set msg_string "$msg_string max_range=[format %.3f [setmd max_range]])"
	if { [llength [set int_time [array get int_times "[setmd cell_grid]"]]]==2 } {
	    set integ_time [lindex $int_time 1]
	} else {
	    set integ_time [lindex [time {
		integrate $int_steps
	    } ] 0]
	    array set int_times [list "[setmd cell_grid]" $integ_time]
	}
	set msg_string "$msg_string => [expr $integ_time/1.0e6]s (vr=[format %.3f [setmd verlet_reuse]])"
	return $integ_time
    }
    set i 0; set ax $min_cells; set bx [expr 0.5*($min_cells+$max_cells)]; set cx $max_cells
    set GR 0.61803399; set GC [expr 1.0-$GR]; set x0 $ax; set x3 $cx
    if { abs($cx-$bx) > abs($bx-$ax) } { set x1 $bx; set x2 [expr $bx + $GC*($cx-$bx)] } else { set x2 $bx; set x1 [expr $bx - $GC*($bx-$ax)] }
    set f1 [eval_cells $i $x1]; incr i 
    set f2 [eval_cells $i $x2]; incr i
    while { abs($x3-$x0) > $tol_cells*(abs($x1) + abs($x2))} {
	if { $f2 < $f1 } {
	    set x0 $x1; set x1 $x2; set x2 [expr $GR*$x1 + $GC*$x3]
	    set f1 $f2; set f2 [eval_cells $i $x2]; incr i
	} else {
	    set x3 $x2; set x2 $x1; set x1 [expr $GR*$x2 + $GC*$x0]
	    set f2 $f1; set f1 [eval_cells $i $x1]; incr i
	}
    }
    set msg_string "$msg_string\nFinal set: "
    if { $f1 < $f2 } { set fs [eval_cells $i $x1] } else { set fs [eval_cells $i $x2] }
    return $msg_string
}



#
# stopParticles
# -------------
# 
# Sets the velocities and forces of all particles currently
# stored in Espresso to zero.
#
#############################################################

proc stop_particles { } { 
    for { set i 0} { $i < [setmd n_part] } {incr i} {
	part $i v 0 0 0
	part $i f 0 0 0
    }
}
proc stopParticles { } {
    puts -nonewline "        Setting all particles' velocities and forces to zero... "; flush stdout
    set old_i 0; set old_e [setmd n_part]; set old [expr 10*$old_e]
    for { set i 0} { $i < $old_e } {incr i} {
	if { [expr $old_i*100]>=$old } { set old [expr $old + 10*$old_e]; puts -nonewline ".,."; flush stdout }; incr old_i
	part $i v 0 0 0
	part $i f 0 0 0
    }
    puts ".,. Done (stopped $old_e particles and as many forces)."
}

#
# center of mass motion
# ---------------------
#
# Calculate the center of mass motion of the system stored in Espresso
#
#############################################################
proc system_com { } {
   set npart [setmd n_part]
    
    # calculate center of mass motion
    set com_vel { 0 0 0 }
    set part_cnt 0
    set i 0
    while { $part_cnt < $npart } {
	if { [part $i] != "na" } {
	    set part_vel [part $i print v] 
	    set com_vel  [vecadd $com_vel $part_vel]
	    incr part_cnt
	}
	incr i
    }
    return [vecscale [expr 1.0/$npart] $com_vel]
}

#
# galileiTransformParticles 
# -------------------------
# 
# Perform a Galilei Transformation on all Particles stored 
# in Espresso such that the center of mass motion of the 
# system is zero. Returns old center of mass motion.
# Assumes particles of equal mass!
#
#############################################################

proc galileiTransformParticles { } {
    set com_vel [system_com]
    
    # subtract center of mass motion from particle velocities
    set npart [setmd n_part]; set part_cnt 0; set i 0
    while { $part_cnt < $npart } {
	if { [part $i] != "na" } {
	    set part_vel [vecsub [part $i print v] $com_vel]
	    part $i v [lindex $part_vel 0] [lindex $part_vel 1] [lindex $part_vel 2]
	    incr part_cnt
	}
	incr i
    }
    return $com_vel
}

#
# Prepare Connection to VMD
# -------------------------
# 
#
#############################################################

proc prepare_vmd_connection { {filename "vmd"} {wait "0"} {start "1" } } {
    writepsf "$filename.psf"
    writepdb "$filename.pdb"
    for {set port 10000} { $port < 65000 } { incr port } {
	catch {imd connect $port} res
	if {$res == ""} break
    }
    set HOSTNAME [exec hostname]
    set vmdout_file [open "vmd_start.script" "w"]
    puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
    puts $vmdout_file "rotate stop"
    puts $vmdout_file "mol modstyle 0 0 CPK 1.800000 0.300000 8.000000 6.000000"
    puts $vmdout_file "mol modcolor 0 0 SegName"
    puts $vmdout_file "imd connect $HOSTNAME $port"
    puts $vmdout_file "imd transfer 1"
    puts $vmdout_file "imd keep 1"
     close $vmdout_file
    if { $start == 0 } {
	puts "Start VMD in the same directory on the machine you with :"
	puts "vmd -e vmd_start.script &"
	imd listen $wait
    } else {
	exec vmd -e vmd_start.script &
    }
}

