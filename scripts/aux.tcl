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
    if { "$cfg" != "-" } { blockfile $f write configs }

    # Close file
    flush $f; close $f
}


proc checkpoint_set { destination { cnt "all" } { tclvar "all" } { ia "all" } } {
    if { [string compare [lindex [split $destination "."] end] "gz"]==0 } {
	set f [open "|gzip -c - >$destination" w]
	set chk [open "[join [lrange [split $destination .] 0 [expr [llength [split $destination .]]-3]] .].chk" a]
    } else { 
	set f [open "$destination" "w"]
	set chk [open "[join [lrange [split $destination .] 0 [expr [llength [split $destination .]]-2]] .].chk" "a"]
    }
    blockfile $f write variable
    if { "$tclvar" != "-" } { foreach j $tclvar { blockfile $f write tclvariable $j } }
    if { "$ia" != "-" } { blockfile $f write interactions }
    blockfile $f write particles "id pos type q v f"
    blockfile $f write bonds
    blockfile $f write random
    flush $f
    if { "$cnt" == "all" } { blockfile $f write configs } else { blockfile $f write configs $cnt }
    flush $f; close $f
    puts $chk "$destination"; flush $chk; close $chk
}


proc checkpoint_read { origin } {
    if { [file exists "$origin.chk"] } { set chk [open "$origin.chk" "r"] 
    } elseif { [file exists "$origin"] } { set chk [open "$origin" "r"] 
    } else { puts "ERROR: Could not find checkpoint-list $origin!\nAborting..."; exit }
    while { [eof $chk]==0 } { if { [gets $chk source] > 0 } {
	if { [string compare [lindex [split $source "."] end] "gz"]==0 } { set f [open "|gzip -c - >$source" r]
	} else { set f [open "$source" "r"] }
	while { [blockfile $f read auto] != "eof" } {}
	close $f
    } }
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
# stopParticles
# -------------
# 
# Sets the velocities and forces of all particles currently
# stored in Espresso to zero.
#
#############################################################

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
    puts $vmdout_file "mol modstyle 0 0 CPK 1.000000 0.300000 8.000000 6.000000"
    puts $vmdout_file "mol modcolor 0 0 SegName"
    puts $vmdout_file "imd connect $HOSTNAME $port"
     close $vmdout_file
    if { $start == 0 } {
	puts "Start VMD in the same directory on the machine you with :"
	puts "vmd -e vmd_start.script &"
	imd listen $wait
    } else {
	exec vmd -e vmd_start.script &
    }
}