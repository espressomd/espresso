#############################################################
#                                                           #
# aux.tcl                                                   #
# =======                                                   #
#                                                           #
# Several additional auxiliary functions for tcl_md.        #
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
# - a list of 'tcl_md' parameters to be saved (out of node_grid|box_l|niatypes|time_step|skin|gamma|bjerrum|...
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

proc polyBlockWrite { destination {write_param {node_grid box_l niatypes time_step skin gamma bjerrum p3m_alpha p3m_r_cut p3m_mesh p3m_cao p3m_epsilon p3m_mesh_offset max_num_cells periodicity}} {write_part "id pos type q v f"} } {
    
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
    close $f
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


proc analysisInit { stat stat_out N_P MPC simtime } {
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
	puts $stat_out "\n$tmp_stat"; flush $stat_out
	puts "Done."
    }
}

proc analysis { stat stat_out N_P MPC simtime } {
    puts -nonewline $stat_out "$simtime "
    foreach i $stat {
	set tmp_stat [analyze $i]
	puts -nonewline ", $i=$tmp_stat"; flush stdout
	puts -nonewline $stat_out "$tmp_stat "
    }
    puts $stat_out " "; flush $stat_out
}


proc stopParticles { } {
    puts -nonewline "        Setting all particles' velocities and forces to zero... "; flush stdout
    set old_i 0; set old_e [setmd npart]; set old [expr 10*$old_e]
    for { set i 0} { $i < $old_e } {incr i} {
	if { [expr $old_i*100]>=$old } { set old [expr $old + 10*$old_e]; puts -nonewline ".,."; flush stdout }; incr old_i
	part $i v 0 0 0
	part $i f 0 0 0
    }
    puts ".,. Done (stopped $old_e particles and as many forces)."
}